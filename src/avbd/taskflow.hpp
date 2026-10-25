/*
 * TaskFlow - a minimal sequence of classified work for a physics step.
 *
 * The solver already knows which of its phases are independent (per-body / per-force
 * loops), which need every participant at a barrier (the colour Gauss-Seidel region), and
 * which mutate shared structure and must run alone (manifold linking, deletion, colouring).
 * This header makes that classification explicit instead of hiding it behind bare pool
 * calls, and gives the execution facility one place to sequence it.
 *
 * The flow is deliberately a linear sequence, not a general DAG: the physics step is
 * linear, and the value of the classification is in how each stage is allowed to use the
 * workers, not in reordering. If a stage grows a dependency on a later one, add the DAG
 * then; a linear flow keeps the determinism argument trivial.
 *
 *   Serial      - runs on the calling thread. Mutates state nothing else may touch.
 *   Independent - `fn(i)` for i in [0, count()). No internal synchronization; items are
 *                 disjoint, so any worker may take any item and the result does not depend
 *                 on the split. Joined by a per-dispatch latch.
 *   Region      - `body(worker, sync)` on every participant, with `sync.arriveAndWait()`
 *                 phase fences. Owns a participant budget for its lifetime. A small step
 *                 can run the same body inline with NullSync instead (`inlineWhen`).
 *
 * The flow is built before it runs (that is what lets stages be cheap closures), so a
 * later stage reads values produced by an earlier one through a captured reference.
 */

#ifndef AVBD_TASKFLOW_HPP
#define AVBD_TASKFLOW_HPP

#include <functional>
#include <utility>
#include <vector>

#include "avbd/job_pool.hpp"

namespace avbd::detail {

enum class TaskKind : unsigned char
{
    Serial,
    Independent,
    Region,
};

class TaskFlow {
public:
    explicit TaskFlow(JobPool &p_pool) :
            pool(p_pool) {
        entries.reserve(8);
    }

    TaskFlow(const TaskFlow &) = delete;
    TaskFlow &operator=(const TaskFlow &) = delete;

    // Drop every stage, keeping the vector's capacity: a step rebuilds the same handful of
    // stages every frame.
    void clear() { entries.clear(); }

    // A single-threaded stage: the calling thread runs it, and nothing else overlaps.
    void serial(std::function<void()> fn) {
        Entry entry;
        entry.kind = TaskKind::Serial;
        entry.run = [fn = std::move(fn)](bool) { fn(); };
        entries.push_back(std::move(entry));
    }

    // `count()` is evaluated when the flow runs, not when it is built: the update set and
    // the force list change earlier in the same step.
    void independent(std::function<int()> count, int minItems, std::function<void(int)> fn) {
        Entry entry;
        entry.kind = TaskKind::Independent;
        entry.minItems = minItems;
        entry.count = std::move(count);
        entry.item = std::move(fn);
        entries.push_back(std::move(entry));
    }

    // A barrier-synchronized region. `inlineWhen` (optional) lets a small step run the same
    // body on the calling thread with NullSync, skipping the dispatch and the barrier.
    template <typename F>
    void region(F &&fn, std::function<bool()> inlineWhen = {}) {
        Entry entry;
        entry.kind = TaskKind::Region;
        entry.inlineWhen = std::move(inlineWhen);
        entry.run = [&pool = pool, body = std::forward<F>(fn)](bool inlineNow) mutable {
            if (inlineNow) {
                pool.runRegionInline(body);
            } else {
                pool.runRegion(body);
            }
        };
        entries.push_back(std::move(entry));
    }

    // Execute every stage in order. A stage only starts after the previous one has fully
    // completed, which is the same fence the solver relied on when the phases were called
    // directly; the difference is that the joins are per-dispatch and a region owns its
    // budget instead of a process-global mutex.
    void run() {
        for (Entry &entry : entries) {
            switch (entry.kind) {
                case TaskKind::Serial:
                    entry.run(false);
                    break;
                case TaskKind::Independent:
                    pool.forCount(entry.count(), entry.item, entry.minItems);
                    break;
                case TaskKind::Region: {
                    const bool inlineNow = pool.threadCount() <= 1
                            || (entry.inlineWhen && entry.inlineWhen());
                    entry.run(inlineNow);
                    break;
                }
            }
        }
    }

private:
    struct Entry {
        TaskKind kind = TaskKind::Serial;
        int minItems = 0;
        std::function<int()> count;       // Independent: item count, evaluated at run time
        std::function<void(int)> item;    // Independent: one item
        std::function<bool()> inlineWhen; // Region: run inline instead of dispatching
        std::function<void(bool)> run;    // Serial/Region: bool = run inline
    };

    JobPool &pool;
    std::vector<Entry> entries;
};

} // namespace avbd::detail

#endif // AVBD_TASKFLOW_HPP
