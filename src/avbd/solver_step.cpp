/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#include <climits>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "avbd/solver.h"

#include "avbd/job_pool.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "radixsort.hpp"

namespace avbd {

void Solver::clear()
{
    while (forces)
        delete forces;

    manifolds.clear();

    while (bodies)
        delete bodies;

    forceAdjacency.clear();
}

// One tick: broad phase / force warmstart / update-set selection once on the calling thread;
// runSubstep does warmstart, colour-parallel iteration and finish (stepInterleaved repeats it).
void Solver::step()
{
    ensurePool();
    collisionDt = dt;

    // Phase timings (diagnostics only, see Solver::PhaseTimes). Four clock reads per step are
    // noise next to a step; the numbers are what get_broadPhase_time() and friends report.
    using Clock = std::chrono::steady_clock;
    const auto msSince = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const Clock::time_point t0 = Clock::now();

    broadPhase();
    const Clock::time_point t1 = Clock::now();
    const int forceCount = warmstartForces();
    selectUpdateBodies();
    projectContacts(false);
    const Clock::time_point t2 = Clock::now();

    int targetIterations = iterations;
    roundsThisStep = std::max(1, targetIterations);
    runSubstep(forceCount, targetIterations, true);
    const Clock::time_point t3 = Clock::now();

    projectContacts(positionProjection);
    rebaseContacts();
    const Clock::time_point t4 = Clock::now();

    lastStepTimes = PhaseTimes{msSince(t0, t1), lastBroadIndexMs, lastBroadTreeMs, lastBroadEmitMs,
            msSince(t1, t2), msSince(t2, t3), msSince(t3, t4)};
}

// One substep: body warmstart, colour-parallel Gauss-Seidel region, BDF1 finish. Idle detection
// runs only on the tick's final substep, so substeps cannot sleep bodies faster than ticks.
void Solver::runSubstep(int forceCount, int targetIterations, bool idleDetection)
{
    const int warmstartCount = static_cast<int>(warmstartOrder.size());
    pool->forCount(warmstartCount, [this](int i) { warmstartBody(i); }, kMinBodiesPerDispatch);

    // Small scenes run the region inline: waking workers for a few dozen bodies costs
    // more than the work itself.
    const bool inlineNow = static_cast<int>(updateOrder.size()) < kMinBodiesPerDispatch
            && forceCount + manifolds.count() < kMinForcesPerDispatch;
    float prevResidual = 0.0f;
    std::atomic<bool> converged{false};
    if (inlineNow) {
        pool->runRegionInline([this, forceCount, targetIterations, &prevResidual, &converged](
                                      int workerIndex, auto &sync) {
            iterateBody(workerIndex, sync, targetIterations, forceCount, prevResidual, converged);
        });
    } else {
        pool->runRegion([this, forceCount, targetIterations, &prevResidual, &converged](
                                int workerIndex, auto &sync) {
            iterateBody(workerIndex, sync, targetIterations, forceCount, prevResidual, converged);
        });
    }

    pool->forCount(warmstartCount, [this, idleDetection](int i) {
        finishBody(i, idleDetection);
    }, kMinBodiesPerDispatch);
}

// Interleaved substeps: collision/forces/update-set once per tick, a few rebased rounds per
// substep; the contact set is not refreshed inside the tick. `dt` is already the substep dt.
void Solver::stepInterleaved(int substeps, int roundsPerSubstep)
{
    ensurePool();
    if (substeps <= 1) {
        step();
        return;
    }

#ifdef AVBD_DIAG
    diag = Diag{}; // warmstartForces fills reaped/soft/new below; reset first
#endif
    collisionDt = dt * substeps;
    using Clock = std::chrono::steady_clock;
    const auto msSince = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const Clock::time_point t0 = Clock::now();
    broadPhase();
    const Clock::time_point t1 = Clock::now();
    const int forceCount = warmstartForces();
    selectUpdateBodies();
    const Clock::time_point t2 = Clock::now();

    // `roundsPerSubstep` is an exact override, not a floor (0 = auto split of tick's iterations).
    // Stability knob: at 12 rounds/tick a welded chain diverges at >=4 substeps (tools/repro_joint.cpp).
    const int rounds = roundsPerSubstep > 0
            ? std::max(1, roundsPerSubstep)
            : std::max(2, (iterations + substeps - 1) / substeps);
#ifdef AVBD_DIAG
    diag.actualRounds = rounds;
#endif
    roundsThisStep = rounds;

    for (int sub = 0; sub < substeps; ++sub) {
        if (sub > 0) {
            rebaseContacts();
        }
        // Contacts refresh only once per tick, so each substep re-brakes against its own start-pose
        // gap over the time left in the tick.
        projectContacts(false, dt * (substeps - sub));
        runSubstep(forceCount, rounds, sub + 1 == substeps);
    }
    const Clock::time_point t3 = Clock::now();
    // Once per tick, not per substep: the correction is positional and velocities are
    // already final, so a second pass inside the same tick would only re-fight the first.
    projectContacts(positionProjection, collisionDt);
    // Rebase to the final pose so the host's penetration diagnostic reflects the tick's actual
    // geometry, not the pose at contact init; the next tick's initialize() replaces it anyway.
    rebaseContacts();
    lastStepTimes = PhaseTimes{msSince(t0, t1), lastBroadIndexMs, lastBroadTreeMs, lastBroadEmitMs,
            msSince(t1, t2), msSince(t2, t3), msSince(t3, Clock::now())};
}

// Rebase every manifold's constraint error at the current pose. Parallel over slots; a
// manifold only writes its own C0, so this is independent work.
void Solver::rebaseContacts()
{
    const int manifoldCount = manifolds.count();
    pool->forCount(manifoldCount, [this](int i) {
        manifolds.rebase(manifolds.live[i], maxPenetration);
    }, kMinForcesPerDispatch);
}

void Solver::ensurePool()
{
    // Pass `threads` straight through: 0 = one worker per hardware thread, 1 = single worker;
    // clamping 0 up to 1 here would silently make the default a serial run.
    const int requested = std::max(0, threads);
    if (pool != nullptr && poolThreads == requested)
        return;

    pool = std::make_unique<detail::JobPool>(static_cast<unsigned>(requested));
    poolThreads = requested;
}

int Solver::threadCount() const
{
    if (pool == nullptr)
        return 1;
    return static_cast<int>(pool->threadCount());
}



// Per iteration: colour-by-colour GS (barrier carries order), dual, block-Jacobi overflow.
// NullSync runs the same body inline in list order - the threads=1 digest-test reference.
template <typename Sync>
void Solver::iterateBody(int workerIndex, Sync &sync, int targetIterations, int forceCount,
        float &r_prevResidual, std::atomic<bool> &r_converged)
{
    // Newton rounds are a newtonRatio slice of the budget; later rounds skip primal passes and
    // decay the dual penalty. The convergence early exit is probed after each Newton round below.
    const int newtonRounds = std::clamp(static_cast<int>(std::ceil(
            newtonRatio * static_cast<float>(targetIterations))), 1, targetIterations);
    for (int it = 0; it < targetIterations; it++)
    {
        const bool newton = it < newtonRounds;
        // One pass per colour: assemble, solve, apply. Bodies of a colour share no force, so the pass
        // is race-free without a snapshot (manifolds stamp self-contained, other forces stamp directly).

        for (int c = 0; newton && c < colourCount; c++)
        {
            const int first = colourBegins[c];
            const int last = colourBegins[c + 1];
            if (last == first)
                continue; // empty colour: every worker skips the barrier alike
            sync.forItems(last - first, [this, first, it](int k) {
                Rigid *body = updateOrder[colourOrder[first + k]];
                Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
                        body->positionAng, body->inertialAng);
                for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
                {
                    if (ref.isManifold())
                        manifolds.updatePrimal(ref.slot(), *this, body, alpha, block);
                    else
                        ref.asForce()->updatePrimal(body, alpha, block);
                }
                block.apply(*body);
            });
            sync.arriveAndWait();
        }

        // Overflow (> kMaxColours cliques): 3-phase block Jacobi on the post-colour snapshot. Phase A
        // caches shared per-contact assembly; no pose is written until the combine pass.
        const int overflowCount = newton ? static_cast<int>(overflowOrder.size()) : 0;
        if (overflowCount > 0)
        {
            sync.forItems(manifolds.count(), [this](int i) {
                const uint32_t slot = manifolds.live[i];
                const bool aIn = manifolds.bodyA[slot]->inUpdateSet &&
                        bodyColour[manifolds.bodyA[slot]->updateIndex] == kOverflowColour;
                const bool bIn = manifolds.bodyB[slot]->inUpdateSet &&
                        bodyColour[manifolds.bodyB[slot]->updateIndex] == kOverflowColour;
                if (!aIn && !bIn)
                    return; // both sides solved in a colour, or frozen
                manifolds.preCompute(slot, *this, alpha);
            });
            sync.arriveAndWait();

            // Phase B: each bucket body assembles inertia + incident constraints (manifolds from the
            // phase-A cache) and solves for its delta.
            sync.forItems(overflowCount, [this](int k) {
                const int i = overflowOrder[k];
                Rigid *body = updateOrder[i];
                Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
                        body->positionAng, body->inertialAng);
                for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
                {
                    if (ref.isManifold())
                        manifolds.stampCached(ref.slot(), body, block);
                    else
                        ref.asForce()->updatePrimal(body, alpha, block);
                }
                block.solveDelta(bodyDeltaLin[i], bodyDeltaAng[i]);
            });
            sync.arriveAndWait();

            // Combine: apply each body's delta scaled by jacobiRelaxation; omega < 1 damps the overshoot a
            // snapshot solve can have where coloured Gauss-Seidel does not.
            sync.forItems(overflowCount, [this](int k) {
                const int i = overflowOrder[k];
                Rigid *body = updateOrder[i];
                body->positionLin = body->positionLin + bodyDeltaLin[i] * jacobiRelaxation;
                body->positionAng = body->positionAng + bodyDeltaAng[i] * jacobiRelaxation;
            });
            sync.arriveAndWait();
        }

        // Dual update is independent: each force reads frozen poses and writes only its own dual state.
        // Frozen constraints are skipped - every input they read is parked.
        const int manifoldCount = manifolds.count();
        const float decay = newton ? 1.0f : stiffnessDecay;
        sync.forItems(manifoldCount + forceCount,
                [this, manifoldCount, decay](int k) {
            if (k < manifoldCount)
            {
                const uint32_t slot = manifolds.live[k];
                if (!manifolds.frozen(slot))
                {
                    if (decay < 1.0f)
                        manifolds.decayPenalty(slot, decay);
                    manifolds.updateDual(slot, *this, alpha);
                }
                return;
            }
            forceOrder[k - manifoldCount]->updateDual(alpha);
        });
        sync.arriveAndWait();


        // Worker 0 alone probes the residual after a Newton round. Converged workers still join
        // every barrier (drain, don't return) - unequal participation deadlocks.
        if (workerIndex == 0 && convergenceThreshold > 0.0f && it + 1 < targetIterations)
        {
            float total = 0.0f;
            int count = 0;
            for (uint32_t slot : manifolds.live)
                manifolds.accumulatePenetration(slot, *this, total, count);
            const float residual = count > 0 ? total / count : 0.0f;
            r_converged.store(r_prevResidual - residual < convergenceThreshold,
                    std::memory_order_relaxed);
            r_prevResidual = residual;
        }
        const bool stop = r_converged.load(std::memory_order_relaxed);
        if (stop)
            break; // every participant breaks at the same iteration - barriers stay paired
    }
}

} // namespace avbd
