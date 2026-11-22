// JobPool: thin adapter over BS::thread_pool v5.1.0; the solver owns work splitting and phase fences.
// Regions reserve threadCount()-1 workers (never nest); joins are per-dispatch latches, never a global drain (which would deadlock barrier-blocked regions).

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

// Phase-boundary sync for region bodies: treat every `arriveAndWait()` as a hard fence (all arrive,
// then all proceed); `forItems` spreads one phase's independent items, barrier stays in the solver loop.

// Real barrier: participants pull items from a shared cursor, then rendezvous.
class LoopSync {
    // Rearm the cursor once per completed phase: after all arrivals (finished pulls done) and
    // before any release (next-phase pulls see 0).
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

    // Run `fn(i)` for i in [0, count): workers grab chunks from the shared cursor, so the split
    // is timing-dependent and items must be independent; the barrier rearms it between phases.
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

// All-or-nothing bulk helper reservation: a partial hold while waiting is the deadlock two regions
// can create (no bulk acquire on a plain semaphore); non-FIFO notify_all is fine for short regions.
class SlotBudget {
public:
    explicit SlotBudget(unsigned capacity) :
            available(capacity) {}

    // Fast path is one atomic CAS, no lock; waiters re-check the predicate under the mutex so a
    // release in the check-to-wait window cannot be lost.
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

// Process-wide pool + region-width budget, materialised lazily so serial solvers spawn nothing;
// static lifetime is safe (no Solver touches it at teardown).
class Scheduler {
public:
    static Scheduler &shared() {
        static Scheduler scheduler;
        return scheduler;
    }

    Scheduler(const Scheduler &) = delete;
    Scheduler &operator=(const Scheduler &) = delete;

    [[nodiscard]] BS::thread_pool<> &pool() { return workers; }

    // Reserve `count` region workers, blocking (never silently serial) until the full width is
    // free; a region never holds slots while waiting, so this cannot deadlock.
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

    // Run `fn(i)` for i in [0, count): independent items, so no cross-thread FP reduction (result
    // is split-independent); `minItems` skips dispatch when one costs tens of microseconds.
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
        // Target ~4 items/chunk (enough chop points, low bookkeeping), clamped to [workers,
        // 4*workers] so small batches aren't one task/item and large ones still split evenly.
        const int workers = static_cast<int>(threadCount());
        const int blocks = std::min(std::clamp(count / 4, workers, workers * 4), count);
        runChunks(count, blocks, [&fn](int, int first, int last) {
            for (int i = first; i < last; i++) {
                fn(i);
            }
        });
    }

    // Run `fn(chunk, first, last)`; chunk index gives each participant its own output buffer.
    // Order scales with `threads` (fixed per value): use only when sorted later (pairs int-sorted), never reductions.
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

    // Run `body(worker, sync)` on `threadCount()` participants; fence phases with `arriveAndWait()`.
    // Workers held for the region's lifetime; solvers run concurrently, never silently serial.
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

    // Small-scene escape hatch: run the body on the calling thread with the serial sync, kept
    // separate so the flow can decide per step whether waking workers is worthwhile.
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
