/*
 * JobPool - the parallel-for the solver uses for its per-body phases.
 *
 * Thin adapter over BS::thread_pool (https://github.com/bshoshany/thread-pool, MIT,
 * header-only). We use the library rather than hand-rolling a pool: it owns the queue,
 * the worker lifetime and the wait handshake, and `detach_blocks` + `wait` is exactly the
 * fork-join the solver needs.
 *
 * What is left to us is the part no library can know: how many independent pieces of work
 * there are. That is what `forCount` is for.
 *
 * Semantics that the solver relies on:
 *   - Work items are independent, so the result never depends on how the range was split.
 *     That is what makes the simulation bit-identical for every thread count (verified in
 *     the test suite by comparing digests).
 *   - Small batches run inline on the calling thread: a dispatch costs ~10-20 us, which
 *     would dominate a scene with a handful of bodies.
 *   - One pool per Solver. Jobs never migrate between simulations, so `wait()` can be a
 *     plain "all my work is done" - no cross-solver coupling, and two worlds run
 *     independently instead of serialising on a shared pool.
 *
 * This header is only included by solver.cpp: the rest of the core sees the forward
 * declaration in solver.h, so the dependency does not leak into every translation unit.
 */

#ifndef AVBD_JOB_POOL_HPP
#define AVBD_JOB_POOL_HPP

#include <algorithm>
#include <cstddef>

#include <BS_thread_pool.hpp>

namespace avbd::detail {

class JobPool {
public:
    // `threads` of 0 means "one worker per hardware thread".
    explicit JobPool(unsigned threads) :
            pool(threads) {}

    JobPool(const JobPool &) = delete;
    JobPool &operator=(const JobPool &) = delete;

    unsigned threadCount() const {
        return static_cast<unsigned>(pool.get_thread_count());
    }

    // Run `fn(i)` for every i in [0, count). Returns only once all of them are done.
    //
    // `minItems` is the smallest batch worth dispatching. One dispatch costs tens of
    // microseconds (waking workers and synchronising), so it only pays off above a number
    // that depends on how much work one item is: a broad-phase pair test is ~40 ns while a
    // body's primal update is ~2 us, and one break-even count cannot serve both. Callers
    // pass the value that matches their per-item cost.
    template <typename F>
    void forCount(int count, F &&fn, int minItems) {
        if (count <= 0) {
            return;
        }
        // One worker means the pool would hand the work straight back with a round trip on
        // top, so a single-threaded solver runs everything on the calling thread.
        if (pool.get_thread_count() <= 1 || count < minItems) {
            for (int i = 0; i < count; i++) {
                fn(i);
            }
            return;
        }

        // One block per worker, and the split of a given count never changes, so runs are
        // reproducible. (The items are independent, so even an unbalanced split would give
        // the same numbers - this just keeps the timing consistent too.)
        const std::size_t blocks = std::min<std::size_t>(pool.get_thread_count(), static_cast<std::size_t>(count));
        pool.detach_blocks(0, count, [&fn](int first, int last) {
            for (int i = first; i < last; i++) {
                fn(i);
            }
        }, blocks);
        pool.wait();
    }

private:
    // `tp::none`: no task priority, no pause support - the solver needs neither, and the
    // defaults keep the dispatch path as short as possible. The wait_deadlock_checks
    // instrumentation is also off, which matters because it makes `wait()` slower.
    BS::thread_pool<> pool;
};

} // namespace avbd::detail

#endif // AVBD_JOB_POOL_HPP
