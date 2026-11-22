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

// One tick. The broad phase / force warmstart / update-set selection run once on the
// calling thread (they parallelise their own stages), then `runSubstep` runs the per-body
// warmstart, the colour-parallel iteration region and the velocity finish. Interleaved
// substepping calls `runSubstep` more than once per tick.
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

    lastStepTimes = PhaseTimes{msSince(t0, t1), lastBroadIndexMs, msSince(t1, t2), msSince(t2, t3),
            msSince(t3, t4)};
}

// One substep: body warmstart (independent), the colour-parallel Gauss-Seidel region,
// then the BDF1 velocity finish (independent). Idle detection only runs on the tick's
// final substep, so substeps cannot sleep bodies faster than ticks.
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

// Interleaved substepping: collision/forces/update-set once per tick, then a small number
// of iterations per substep with the contact error rebased at each substep start. This is
// the cheaper alternative to running the full iteration set once per substep; it trades
// some constraint freshness (the contact set is not refreshed inside the tick) for a
// much lower round count. `dt` is already the substep dt.
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

    // Rounds per substep. An explicit `roundsPerSubstep` is an override: the caller asked
    // for exactly that many primal rounds in each substep, and honouring it is what makes
    // the setting mean something (it used to be a floor, so asking for 2 while `iterations`
    // implied 3 silently ran 3 and the knob could not be lowered below the split). 0 asks
    // for the automatic split of the tick's `iterations` across the substeps.
    //
    // This number is a *stability* parameter, not just a cost one: at a fixed 12 rounds per
    // tick a four-link welded chain settles with 1 or 2 substeps and diverges with 4 or more
    // (measured, tools/repro_joint.cpp: |vz| = 0.07 / 0.06 / 0.31 / 140 / 413 / 1310 m/s for
    // 1x12 / 2x6 / 3x4 / 4x3 / 6x2 / 12x1). Fewer substeps means more rounds inside one
    // linearisation, and the augmented-Lagrangian dual converges inside that linearisation.
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
        // The contact set is only refreshed once per tick, so every substep re-brakes
        // approach speeds against the gap measured at its own start pose, over the time
        // still left in the tick.
        projectContacts(false, dt * (substeps - sub));
        runSubstep(forceCount, rounds, sub + 1 == substeps);
    }
    const Clock::time_point t3 = Clock::now();
    // Once per tick, not per substep: the correction is positional and velocities are
    // already final, so a second pass inside the same tick would only re-fight the first.
    projectContacts(positionProjection, collisionDt);
    // The tick's stored error is the one the solver saw when its contacts were initialized;
    // the bodies have moved since. Rebase it at the final pose so what the host reads back
    // (the penetration diagnostic) is the geometry the pair actually ended the tick with,
    // not a one-tick-old pose. The next tick's initialize() replaces it from its own
    // collision pass either way; this only feeds the diagnostics and the host.
    rebaseContacts();
    lastStepTimes = PhaseTimes{msSince(t0, t1), lastBroadIndexMs, msSince(t1, t2),
            msSince(t2, t3), msSince(t3, Clock::now())};
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

// Inverse inertia is diagonal in the body frame (`Rigid::moment` holds the inertia itself,
// the same convention GenericJoint::penaltyLimit reads), so both halves of a projection come
// from one rotation into that frame: the quadratic form t.(I^-1 t) is the effective angular
// compliance along the contact's torque axis, and I^-1 t is the rotation vector that axis's
// impulse produces.

void Solver::ensurePool()
{
    // Pass the value straight through: 0 means "one worker per hardware thread" and 1 means
    // a single worker, which both the pool and `threads` document the same way. (Clamping 0
    // up to 1 here would silently turn the default into a serial run.)
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



template <typename Sync>
void Solver::iterateBody(int workerIndex, Sync &sync, int targetIterations, int forceCount,
        float &r_prevResidual, std::atomic<bool> &r_converged)
{
    // Convergence early exit: measure the residual at the end of a round; when it
    // stops improving (delta below the threshold or rising), the remaining rounds would
    // burn cycles on float noise. 0 disables the check.
    const int newtonRounds = std::clamp(static_cast<int>(std::ceil(
            newtonRatio * static_cast<float>(targetIterations))), 1, targetIterations);
    for (int it = 0; it < targetIterations; it++)
    {
        const bool newton = it < newtonRounds;
        // Colour-parallel Gauss-Seidel. One pass per colour: assemble
        // (manifolds stamp self-contained against the live poses, other
        // forces stamp directly), solve, apply. Bodies of one colour
        // share no force, so the pass is race-free without a snapshot.

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

        // Overflow bucket (> kMaxColours cliques): block Jacobi against
        // the post-colour snapshot, in three phases. Phase A precomputes
        // each manifold's shared per-contact assembly (world arms,
        // angular jacobians, clamped forces) once - work the two body
        // stamps used to duplicate. Nobody writes a pose until the
        // combine pass, so every phase reads the same frozen snapshot.
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

            // Phase B: each bucket body assembles its full block system
            // (own inertia plus every incident constraint - manifolds
            // stamped from the phase-A cache, other forces stamping
            // directly) and solves for its delta.
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

            // Combine: each bucket body applies its own delta, scaled by
            // jacobiRelaxation. A snapshot solve can overshoot where the
            // coloured Gauss-Seidel updates do not; omega < 1 trades
            // per-iteration progress for stability.
            sync.forItems(overflowCount, [this](int k) {
                const int i = overflowOrder[k];
                Rigid *body = updateOrder[i];
                body->positionLin = body->positionLin + bodyDeltaLin[i] * jacobiRelaxation;
                body->positionAng = body->positionAng + bodyDeltaAng[i] * jacobiRelaxation;
            });
            sync.arriveAndWait();
        }

        // Dual update: one item per force, and each force only reads the (frozen)
        // body poses and writes its own dual state, so the whole pass is independent.
        // Frozen constraints have nothing to update: every input they read is parked.
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


        // Convergence probe after a Newton round: a full primal pass is expensive, so
        // skip further rounds when the residual stops shrinking. Only the calling
        // thread (worker 0) evaluates the reduction, and the verdict is stored in a
        // shared atomic that EVERY participant reads BEFORE the barrier - each worker
        // that sees `r_converged` skips its own remaining work but still participates
        // in every barrier (draining instead of returning), so the barrier
        // participation count stays equal and no worker deadlocks.
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

// One body's BDF1 velocity after the final iteration. Per-body and independent, so
// `runSubstep` dispatches it as an independent batch; it covers static bodies too, so
// their bookkeeping matches, though only movable ones get a new velocity.

} // namespace avbd
