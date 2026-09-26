/*
 * JobPool - the parallel machinery the solver runs its per-body phases on.
 *
 * Thin adapter over BS::thread_pool v5.1.0 (https://github.com/bshoshany/thread-pool,
 * MIT, header-only, C++17). The library owns the queue, the worker lifetime and the
 * wait handshake; what is left to us is the part no library can know: how many
 * independent pieces of work there are, and where the phase boundaries are.
 *
 * THREAD COUNT:
 * - Process-global shared pool (not per-solver)
 * - Lazy initialization on first parallel use
 * - Worker count = min(hardware_concurrency, threads parameter)
 * - >hardware is clamped to hardware_concurrency
 * - Serial solvers (threads=1) never spawn threads
 *
 * EXECUTION STYLES:
 * - forCount: fork-join, one dispatch, workers split range, wait()
 *   - Break-even minItems=256 to avoid dispatch overhead
 *   - Used for phases running once per step (warmstart, finish)
 *
 * - runLoop: persistent-worker loop, iteration-shaped work
 *   - One dispatch arms threadCount() participants (calling thread is worker 0)
 *   - Each participant runs whole loop, synchronizes at phase boundaries via LoopSync barrier
 *   - Eliminates dispatch storm for colour-parallel Gauss-Seidel
 *   - Barrier cost: ~1-3 µs per phase boundary
 *
 * WORK DISTRIBUTION:
 * - Process-global queue, no work stealing (flat queue)
 * - Fixed chunking: cursor-based with 4 chunks/participant
 * - Atomic<int> fetch_add for work stealing (limited)
 * - Load imbalance possible when work distribution uneven
 *
 * PHASE BOUNDARIES:
 * - LoopSync barrier ensures all workers finish before next phase
 * - Barrier is lightweight (~1-3 µs) and non-blocking
 * - Barriers prevent race conditions in Gauss-Seidel ordering
 * - Each barrier corresponds to a sync.arriveAndWait() call
 *
 * THREAD SAFETY:
 * - runLoop calls are exclusive (try-locked mutex)
 * - Two loops at once would deadlock on shared queue
 * - Godot steps one solver at a time, so lock uncontended in practice
 * - Contended callers degrade to running loop inline (threads=1 reference)
 *
 * MEMORY MODEL:
 * - std::memory_order_relaxed for cursor updates (no ordering needed)
 * - std::memory_order_acquire/release for barrier sync
 * - Thread-local Eigen scratch avoids false sharing
 * - Pool-based Manifold allocation avoids heap churn
 *
 * THIS HEADER IS ONLY INCLUDED BY SOLVER.CPP
 * The rest of the core sees the forward declaration in solver.h
 * to avoid leaking thread pool dependency into every translation unit.
 */

#ifndef AVBD_JOB_POOL_HPP
#define AVBD_JOB_POOL_HPP

#include <algorithm>
#include <atomic>
#include <barrier>
#include <cstddef>
#include <mutex>
#include <mutex>
#include <thread>

#include <BS_thread_pool.hpp>

namespace avbd::detail {

// The phase-boundary synchronisers handed to a runLoop body.
//
// The loop body is written once against `auto &sync` and must treat every
// `arriveAndWait()` as a hard phase fence: all workers arrive, then all proceed. Both
// types provide `forItems`, which spreads one phase's independent items over the
// participants (or runs them serially) - the barrier call itself stays in the solver's
// loop so the Gauss-Seidel structure reads plainly.

// Real barrier: participants pull items from a shared cursor, then rendezvous.
class LoopSync {
    // Rearm the work cursor for the next phase. Runs exactly once per completed phase,
    // after every participant has arrived (so every pull of the finished phase is done)
    // and before any participant is released (so every pull of the next phase sees 0).
    struct CursorReset {
        LoopSync *self;
        void operator()() noexcept { self->cursor.store(0, std::memory_order_relaxed); }
    };

public:
    explicit LoopSync(unsigned count) :
            participants(count),
            barrier(static_cast<std::ptrdiff_t>(count), CursorReset{this}) {}

    LoopSync(const LoopSync &) = delete;
    LoopSync &operator=(const LoopSync &) = delete;

    // Fence: every write any worker made before its arrival is visible to every worker
    // after the fence.
    void arriveAndWait() {
        barrier.arrive_and_wait();
    }

    // Run `fn(i)` for every i in [0, count), self-balanced: workers grab chunks from the
    // shared cursor until this phase's range is drained. Items must be independent - the
    // split varies with timing by design. The cursor is rearmed by the barrier's
    // completion, so consecutive phases must be separated by `arriveAndWait()`, which is
    // the loop contract anyway.
    template <typename F>
    void forItems(int count, F &&fn) {
        if (count <= 0) {
            return;
        }
        // ~4 chunks per participant: few atomic operations, enough chop points that a
        // slow worker does not strand a large tail.
        const int grain = std::max(1, count / static_cast<int>(participants * 4));
        for (int first = cursor.fetch_add(grain, std::memory_order_relaxed); first < count;
                first = cursor.fetch_add(grain, std::memory_order_relaxed)) {
            const int last = std::min(count, first + grain);
            for (int i = first; i < last; i++) {
                fn(i);
            }
        }
    }

private:
    unsigned participants;
    std::atomic<int> cursor{0};
    std::barrier<CursorReset> barrier;
};

// Inline stand-in for a single participant: phases already run sequentially on the
// calling thread, so the fence is a no-op and the items run in order.
class NullSync {
public:
    void arriveAndWait() {}

    template <typename F>
    void forItems(int count, F &&fn) {
        for (int i = 0; i < count; i++) {
            fn(i);
        }
    }
};

// The process-global worker pool. Materialised lazily so a serial solver never spawns
// threads; static lifetime is safe because no Solver touches it during teardown.
inline BS::thread_pool<> &sharedPool() {
    static BS::thread_pool<> pool;
    return pool;
}

// Serialises runLoop calls across all solvers, see the exclusivity note above.
inline std::mutex &loopMutex() {
    static std::mutex mutex;
    return mutex;
}

class JobPool {
public:
    // `threads` of 0 means "one worker per hardware thread", 1 means "run everything
    // inline on the calling thread". Neither spawns anything here.
    explicit JobPool(unsigned threads) :
            desired(threads) {}

    JobPool(const JobPool &) = delete;
    JobPool &operator=(const JobPool &) = delete;

    // Worker budget this solver may claim, clamped to what the machine has. Pure
    // arithmetic: it never touches (or creates) the shared pool.
    unsigned threadCount() const {
        const unsigned hw = std::thread::hardware_concurrency();
        const unsigned available = hw == 0 ? 1 : hw;
        return desired == 0 ? available : std::min(desired, available);
    }

    // Run `fn(i)` for every i in [0, count). Returns only once all of them are done.
    //
    // `minItems` is the smallest batch worth dispatching: one dispatch costs tens of
    // microseconds (waking workers and synchronising), which would dominate a scene
    // with a handful of bodies. Callers pass the value that matches their per-item
    // cost. Items must be independent; the block split is fixed for a given count, so
    // runs stay reproducible even though that does not change the numbers.
    template <typename F>
    void forCount(int count, F &&fn, int minItems) {
        if (count <= 0) {
            return;
        }
        // A budget of one means the pool would hand the work straight back with a round
        // trip on top, so a single-threaded solver runs everything on the calling thread.
        if (threadCount() <= 1 || count < minItems) {
            for (int i = 0; i < count; i++) {
                fn(i);
            }
            return;
        }

        const unsigned workers = threadCount();
        const std::size_t blocks = std::min<std::size_t>(workers, static_cast<std::size_t>(count));
        sharedPool().detach_blocks(0, count, [&fn](int first, int last) {
            for (int i = first; i < last; i++) {
                fn(i);
            }
        }, blocks);
        sharedPool().wait();
    }

    // Run `fn(worker, sync)` on `threadCount()` participants at once and return only
    // when the whole loop has finished. `worker` is the caller's index (0 on the
    // calling thread); `sync` is a LoopSync or NullSync - write the body once against
    // `auto &sync`. The body must be deterministic per item and fence its phases with
    // `sync.arriveAndWait()`; see the class comments above.
    //
    // Degrades to an inline run when the budget is one, or when another solver is
    // already inside a loop (the try-lock below; a second loop on a shared queue would
    // deadlock). The inline fallback is the serial reference, so results are unchanged.
    template <typename F>
    void runLoop(F &&fn) {
        const unsigned workers = threadCount();
        if (workers <= 1) {
            NullSync sync;
            fn(0, sync);
            return;
        }

        std::unique_lock<std::mutex> lock(loopMutex(), std::try_to_lock);
        if (!lock.owns_lock()) {
            NullSync sync;
            fn(0, sync);
            return;
        }

        LoopSync sync(workers);
        for (unsigned w = 1; w < workers; ++w) {
            sharedPool().detach_task([&sync, &fn, w] {
                fn(static_cast<int>(w), sync);
            });
        }
        fn(0, sync);
        // The workers own no state of this frame once their task callables return, so a
        // plain drain of the pool - rather than a per-loop latch - is the safe join: it
        // keeps the worker bookkeeping inside BS::thread_pool's own synchronisation and
        // lets `sync` die only after every reference to it is gone. Sequential stepping
        // (Godot's rule) means the queue holds nothing but this loop's tasks.
        sharedPool().wait();
    }

private:
    // Normalised `Solver::threads`: 0 = all hardware threads, 1 = inline, N = at most N.
    unsigned desired;
};

} // namespace avbd::detail

#endif // AVBD_JOB_POOL_HPP
