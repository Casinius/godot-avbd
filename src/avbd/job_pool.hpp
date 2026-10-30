/*
 * JobPool - the parallel machinery the solver runs its phases on.
 *
 * Thin adapter over BS::thread_pool v5.1.0 (https://github.com/bshoshany/thread-pool,
 * MIT, header-only, C++17). The library owns the worker queue and lifetime; what is left
 * to us is the part no library can know: how many independent pieces of work there are,
 * which of them may overlap, and where the phase boundaries are.
 *
 * TASK CLASSES (the solver sequences them directly):
 * - Independent work: per-item functions with no internal synchronization. Any worker may
 *   run any item, and several independent tasks may be in flight at once. Joined by a
 *   per-dispatch latch, never by a global queue drain.
 * - Region work: a body that runs on every participant and fences its own phases with
 *   LoopSync barriers (the Gauss-Seidel primal). Regions hold a budget of participant
 *   slots for their lifetime, so two solvers can run their regions concurrently without
 *   either silently falling back to serial. Regions must not nest.
 * - Serial work: runs on the calling thread; used for mutations (manifold linking,
 *   deletion, colouring) that must not race with anything else.
 *
 * WHY THE GLOBAL MUTEX IS GONE:
 * The previous design serialized all solvers with a try-locked mutex around the persistent
 * loop and joined with `sharedPool().wait()`. That `wait()` drains the whole queue, so it
 * cannot coexist with region participants blocked on a barrier: it would wait for them
 * forever. Every join here is a per-dispatch `std::latch`, which also lets independent
 * tasks from another solver run while a region is parked at a barrier.
 *
 * BUDGET:
 * - One process-wide worker pool, lazily created on first parallel use; serial solvers
 *   (threads=1) never create it.
 * - Each JobPool's `threadCount()` is its requested budget, clamped to the machine.
 * - A region needs `threadCount()-1` workers (the calling thread is participant 0) and
 *   reserves them from a counting semaphore. Two regions at full width therefore run one
 *   after the other (by budget, not by silent serialization); setting `threads=4` on two
 *   solvers on an 8-thread machine lets both run at once.
 *
 * MEMORY MODEL:
 * - `std::latch` release/acquire covers the dispatch join; `std::barrier` covers every
 *   region phase boundary. Work items are disjoint by contract, so no floating-point
 *   reduction happens across threads and the result is independent of the split.
 */

#ifndef AVBD_JOB_POOL_HPP
#define AVBD_JOB_POOL_HPP

#include <algorithm>
#include <atomic>
#include <barrier>
#include <condition_variable>
#include <cstddef>
#include <latch>
#include <memory>
#include <mutex>
#include <thread>

#include <BS_thread_pool.hpp>

namespace avbd::detail {

// The largest region width we can reserve. Hardware thread counts above this are clamped;
// 256 covers every desktop/console target this solver has run on.
inline constexpr std::ptrdiff_t kMaxRegionSlots = 256;

// The phase-boundary synchronisers handed to a region body.
//
// The loop body is written once against `auto &sync` and must treat every
// `arriveAndWait()` as a hard phase fence: all participants arrive, then all proceed. Both
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

// Count the dispatch parent down on every exit path (including an escaping exception);
// the parent then waits on the latch instead of a global queue drain.
struct LatchCountDown {
    std::latch &latch;
    ~LatchCountDown() { latch.count_down(); }
};

// All-or-nothing worker reservation for regions. A region asks for its whole helper
// width at once, before dispatching anything: holding a partial reservation while waiting
// for the rest is exactly the deadlock two concurrent regions can create, and a plain
// counting semaphore does not offer a bulk acquire. `notify_all` with a size predicate is
// not FIFO, which is fine here: regions are short and at most a handful exist. If a large
// region could be starved by a stream of small ones, add a ticket queue; the solver does
// not produce that pattern.
class SlotBudget {
public:
    explicit SlotBudget(unsigned capacity) :
            available(capacity) {}

    // The common case (slots are free) is a single atomic CAS, no lock. Only a region that
    // has to wait touches the mutex/condition variable, and the predicate is re-checked
    // under that mutex so a release in the check-to-wait window cannot be lost.
    void acquire(unsigned count) {
        if (count == 0) {
            return;
        }
        unsigned current = available.load(std::memory_order_relaxed);
        while (current >= count) {
            if (available.compare_exchange_weak(current, current - count,
                        std::memory_order_acq_rel, std::memory_order_relaxed)) {
                return;
            }
        }

        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this, count] {
            return available.load(std::memory_order_relaxed) >= count;
        });
        available.fetch_sub(count, std::memory_order_relaxed);
    }

    void release(unsigned count) {
        if (count == 0) {
            return;
        }
        {
            // The increment stays under the mutex so a waiter's predicate check and its
            // transition to sleep cannot interleave with it.
            std::lock_guard<std::mutex> lock(mutex);
            available.fetch_add(count, std::memory_order_relaxed);
        }
        condition.notify_all();
    }

private:
    std::mutex mutex;
    std::condition_variable condition;
    std::atomic<unsigned> available;
};

// The process-wide pool and the region-width budget. Materialised lazily so a serial
// solver never spawns threads; static lifetime is safe because no Solver touches it
// during teardown.
class Scheduler {
public:
    static Scheduler &shared() {
        static Scheduler scheduler;
        return scheduler;
    }

    Scheduler(const Scheduler &) = delete;
    Scheduler &operator=(const Scheduler &) = delete;

    [[nodiscard]] BS::thread_pool<> &pool() { return workers; }

    // Reserve `count` region workers. Blocking: a region that cannot get its full width
    // waits for another region to finish instead of silently running serially. A region
    // never holds slots while waiting for another region, so this cannot deadlock.
    void acquire(unsigned count) { slots.acquire(count); }

    void release(unsigned count) { slots.release(count); }

private:
    Scheduler() :
            workers(),
            slots(capacityOf(workers)) {}

    static unsigned capacityOf(const BS::thread_pool<> &pool) {
        return static_cast<unsigned>(std::min<std::ptrdiff_t>(
                kMaxRegionSlots, static_cast<std::ptrdiff_t>(pool.get_thread_count())));
    }

    BS::thread_pool<> workers;
    SlotBudget slots;
};

class JobPool {
public:
    // `threads` of 0 means "one worker per hardware thread", 1 means "run everything
    // inline on the calling thread". Neither spawns anything here.
    explicit JobPool(unsigned threads) :
            desired(threads) {}

    JobPool(const JobPool &) = delete;
    JobPool &operator=(const JobPool &) = delete;

    // Worker budget this solver may claim, clamped to what the machine has. Pure
    // arithmetic: it never touches (or creates) the shared scheduler.
    unsigned threadCount() const {
        const unsigned hw = std::thread::hardware_concurrency();
        const unsigned available = hw == 0 ? 1 : hw;
        return desired == 0 ? available : std::min(desired, available);
    }

    // Number of contiguous chunks `forBlocks` would split `count` items into: one per
    // participant, clamped to the item count (0 for an empty range).
    [[nodiscard]] int blockCount(int count) const {
        if (count <= 0)
            return 0;
        return static_cast<int>(std::min<unsigned>(threadCount(), static_cast<unsigned>(count)));
    }

    // Run `fn(i)` for every i in [0, count). Independent items: the split only affects which
    // worker computes an item, never the result. The calling thread runs one chunk and then
    // joins the others on a latch.
    //
    // `minItems` is the smallest batch worth dispatching: one dispatch costs tens of
    // microseconds (waking workers and synchronising), which would dominate a scene
    // with a handful of bodies.
    template <typename F>
    void forCount(int count, F &&fn, int minItems) {
        if (count <= 0) {
            return;
        }
        if (threadCount() <= 1 || count < minItems) {
            for (int i = 0; i < count; i++) {
                fn(i);
            }
            return;
        }
        // Target ~4 items per chunk: enough chop points for an uneven per-body cost, few
        // enough that the latch and task bookkeeping stay in the noise. Clamp to
        // [workers, 4*workers] so a batch near the break-even size does not turn into one
        // task per item while a large batch still splits evenly.
        const int workers = static_cast<int>(threadCount());
        const int blocks = std::min(std::clamp(count / 4, workers, workers * 4), count);
        runChunks(count, blocks, [&fn](int, int first, int last) {
            for (int i = first; i < last; i++) {
                fn(i);
            }
        });
    }

    // Run `fn(chunk, first, last)` over contiguous chunks of [0, count). Unlike `forCount`
    // this hands the chunk index to the callback, which is what lets a phase give each
    // participant its own output buffer (the broad phase's pair emission). Chunk boundaries
    // depend on the thread count, so the per-chunk order is only fixed for a given
    // `threads` value: only use this for work whose final order is imposed afterwards
    // (the pair list is integer-sorted), never for reduction order.
    template <typename F>
    void forBlocks(int count, F &&fn) {
        const int blocks = blockCount(count);
        if (blocks <= 0) {
            return;
        }
        runChunks(count, blocks, [&fn](int chunk, int first, int last) {
            fn(chunk, first, last);
        });
    }

    // Run a region body on `threadCount()` participants and return when every participant
    // has left the region. `body(worker, sync)` must fence its phases with
    // `sync.arriveAndWait()`. Regions reserve workers for their whole lifetime, so two
    // solvers can run concurrently (up to the machine's width) without the old global
    // mutex, and never silently degrade to serial.
    template <typename F>
    void runRegion(F &&body) {
        const unsigned participants = threadCount();
        if (participants <= 1) {
            NullSync sync;
            body(0, sync);
            return;
        }

        Scheduler &scheduler = Scheduler::shared();
        const unsigned helpers = participants - 1;
        scheduler.acquire(helpers);

        LoopSync sync(participants);
        std::latch done(static_cast<std::ptrdiff_t>(helpers));
        for (unsigned w = 1; w < participants; ++w) {
            scheduler.pool().detach_task([&body, &sync, &done, w] {
                LatchCountDown guard{done};
                body(static_cast<int>(w), sync);
            });
        }
        body(0, sync);
        done.wait();
        scheduler.release(helpers);
    }

    // The small-scene escape hatch: run the region body on the calling thread with the
    // serial synchroniser. Kept separate from `runRegion` so the flow can decide per step
    // whether the work justifies waking workers at all.
    template <typename F>
    void runRegionInline(F &&body) {
        NullSync sync;
        body(0, sync);
    }

private:
    template <typename ChunkFn>
    void runChunks(int count, int blocks, ChunkFn &&fn) {
        if (blocks <= 1) {
            fn(0, 0, count);
            return;
        }
        Scheduler &scheduler = Scheduler::shared();
        std::latch done(static_cast<std::ptrdiff_t>(blocks - 1));
        const auto firstOf = [count, blocks](int chunk) {
            return static_cast<int>((static_cast<long long>(count) * chunk) / blocks);
        };
        for (int b = 1; b < blocks; ++b) {
            scheduler.pool().detach_task([&fn, &done, b, firstOf] {
                LatchCountDown guard{done};
                fn(b, firstOf(b), firstOf(b + 1));
            });
        }
        fn(0, 0, firstOf(1));
        done.wait();
    }

    // Normalised `Solver::threads`: 0 = all hardware threads, 1 = inline, N = at most N.
    unsigned desired;
};

} // namespace avbd::detail

#endif // AVBD_JOB_POOL_HPP
