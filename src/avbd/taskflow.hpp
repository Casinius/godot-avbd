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
 * Stages are a tagged variant of plain function pointers and a `void *ctx`, not
 * `std::function`: building a step allocates nothing, and there is no type-erased thunk
 * per stage. The context is an object owned by the caller (a stack `SubstepCtx`), so a
 * stage never outlives the values it reads.
 */

#ifndef AVBD_TASKFLOW_HPP
#define AVBD_TASKFLOW_HPP

#include <variant>
#include <vector>

#include "avbd/job_pool.hpp"

namespace avbd::detail {

struct SerialStage {
    void (*run)(void *ctx) noexcept;
    void *ctx;
};

struct IndependentStage {
    int (*count)(void *ctx) noexcept;
    void (*item)(void *ctx, int i) noexcept;
    int minItems;
    void *ctx;
};

struct RegionStage {
    void (*runInline)(void *ctx) noexcept;   // body(0, NullSync)
    void (*runParallel)(void *ctx) noexcept; // body(w, LoopSync) on every participant
    bool (*inlineWhen)(void *ctx) noexcept;  // may be null
    void *ctx;
};

using Stage = std::variant<SerialStage, IndependentStage, RegionStage>;

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
    void serial(void (*fn)(void *ctx) noexcept, void *ctx) {
        entries.push_back(SerialStage{fn, ctx});
    }

    // `count(ctx)` is evaluated when the flow runs, not when it is built: the update set
    // and the force list change earlier in the same step.
    void independent(int (*count)(void *ctx) noexcept, void (*item)(void *ctx, int i) noexcept,
            int minItems, void *ctx) {
        entries.push_back(IndependentStage{count, item, minItems, ctx});
    }

    // A barrier-synchronized region. `inlineWhen` (optional) lets a small step run the same
    // body on the calling thread with NullSync, skipping the dispatch and the barrier. The
    // two bodies are separate thunks because the sync type is part of the callee signature.
    void region(void (*runInline)(void *ctx) noexcept, void (*runParallel)(void *ctx) noexcept,
            bool (*inlineWhen)(void *ctx) noexcept, void *ctx) {
        entries.push_back(RegionStage{runInline, runParallel, inlineWhen, ctx});
    }

    // Execute every stage in order. A stage only starts after the previous one has fully
    // completed, which is the same fence the solver relied on when the phases were called
    // directly; the difference is that the joins are per-dispatch and a region owns its
    // budget instead of a process-global mutex.
    void run() {
        for (Stage &stage : entries) {
            if (auto *serial = std::get_if<SerialStage>(&stage)) {
                serial->run(serial->ctx);
            } else if (auto *independent = std::get_if<IndependentStage>(&stage)) {
                pool.forCount(independent->count(independent->ctx),
                        [item = independent->item, ctx = independent->ctx](int i) { item(ctx, i); },
                        independent->minItems);
            } else {
                auto *region = std::get_if<RegionStage>(&stage);
                const bool inlineNow = pool.threadCount() <= 1
                        || (region->inlineWhen && region->inlineWhen(region->ctx));
                if (inlineNow) {
                    region->runInline(region->ctx);
                } else {
                    region->runParallel(region->ctx);
                }
            }
        }
    }

private:
    JobPool &pool;
    std::vector<Stage> entries;
};

} // namespace avbd::detail

#endif // AVBD_TASKFLOW_HPP
