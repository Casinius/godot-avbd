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

#include <cmath>
#include <algorithm>
#include <numeric>
#include <span>

#include "avbd/solver.h"

#include "avbd/job_pool.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "avbd/spatial_hash.hpp"

namespace avbd {

Solver::Solver()
{
}

Solver::~Solver()
{
    clear();
    // Release the contact-block cache this solver's churn filled. Blocks outlive their
    // solver on purpose (a destroyed contact is recycled by the next one), but when the
    // whole simulation shuts down the cache must not leak.
    Manifold::drainPool();
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local)
{
    return pick(origin, dir, local, 0xFFFFFFFFu);
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local, uint32_t p_mask)
{
    const float epsilon = 1.0e-6f;
    float bestT = INFINITY;
    Rigid *bestBody = 0;
    float3 bestLocal = {0, 0, 0};

    // Exact per-body test: ray against the body's local box. The broad phase is
    // a spatial hash rebuilt per step and not efficient to walk per ray, so
    // picking scans the body list; the exact test decides.
    const auto testBody = [&](Rigid *body) {
        if ((body->collisionLayer & p_mask) == 0)
            return;
        quat invRot = qconj(body->positionAng);
        float3 o = qrot(invRot, origin - body->positionLin);
        float3 d = qrot(invRot, dir);
        float3 half = body->size * 0.5f;

        float tEnter = 0.0f;
        float tExit = INFINITY;
        bool hit = true;

        for (int i = 0; i < 3; ++i)
        {
            if (std::fabs(d[i]) < epsilon)
            {
                if (o[i] < -half[i] || o[i] > half[i])
                {
                    hit = false;
                    break;
                }
                continue;
            }

            float invD = 1.0f / d[i];
            float t0 = (-half[i] - o[i]) * invD;
            float t1 = (half[i] - o[i]) * invD;
            if (t0 > t1)
            {
                float tmp = t0;
                t0 = t1;
                t1 = tmp;
            }

            tEnter = std::max(tEnter, t0);
            tExit = std::min(tExit, t1);
            if (tEnter > tExit)
            {
                hit = false;
                break;
            }
        }

        if (!hit)
            return;

        float tHit = tEnter >= 0.0f ? tEnter : tExit;
        if (tHit < 0.0f)
            return;

        if (tHit < bestT)
        {
            bestT = tHit;
            bestBody = body;
            bestLocal = o + d * tHit;
        }
    };

    // Unlike the upstream demo (which picks bodies for dragging) this also tests
    // static bodies, so rays hit the ground and ramps.
    for (Rigid *body = bodies; body != 0; body = body->next)
        testBody(body);

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
}

void Solver::collectAABB(float3 min, float3 max, uint32_t p_mask, std::vector<Rigid *> &out) const
{
    // The broad phase is a spatial hash rebuilt per step, so shape/point queries
    // walk the body list: each body's bound is centre plus its tight world
    // half-extent (the same formula the hash's spheres derive from).
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if ((body->collisionLayer & p_mask) == 0)
            continue;
        const float3 wx = aabbHalfExtent(*body);
        const float3 bmin = body->positionLin - wx;
        const float3 bmax = body->positionLin + wx;
        if (!(bmin.x > max.x || bmin.y > max.y || bmin.z > max.z ||
                min.x > bmax.x || min.y > bmax.y || min.z > bmax.z))
            out.push_back(body);
    }
}

void Solver::clear()
{
    while (forces)
        delete forces;

    while (bodies)
        delete bodies;
}

void Solver::step()
{
    ensurePool();

    broadPhase();
    const int forceCount = warmstartForces();
    selectUpdateBodies();
    warmstartBodies();
    solveIterations(forceCount);
    finishVelocities();
}

// A broad-phase candidate test: tight world AABBs plus Godot's layer/mask pair
// rule, minus pairs already joined by a force (IgnoreCollision, joints). The AABB
// is the tight box of the actual shape, unlike a bounding-sphere test whose radius
// is the shape's diagonal and therefore admits many far-apart candidates.
static bool pairOverlaps(Rigid *bodyA, const float3 &aHalf, Rigid *bodyB, const float3 &bHalf)
{
    const float3 aMin = bodyA->positionLin - aHalf;
    const float3 aMax = bodyA->positionLin + aHalf;
    const float3 bMin = bodyB->positionLin - bHalf;
    const float3 bMax = bodyB->positionLin + bHalf;
    const bool boundsOverlap = !(bMin.x > aMax.x || bMin.y > aMax.y || bMin.z > aMax.z ||
            aMin.x > bMax.x || aMin.y > bMax.y || aMin.z > bMax.z);
    // Godot semantics: the pair collides when either side's layer is in the other's
    // mask (defaults 1/1 keep every pair, as before these fields existed).
    const bool layersAllow = ((bodyA->collisionLayer & bodyB->collisionMask) != 0) ||
            ((bodyB->collisionLayer & bodyA->collisionMask) != 0);
    return boundsOverlap && layersAllow && !bodyA->constrainedTo(bodyB);
}

// Contact detection: hierarchical spatial hash (big-block screen + layered
// hash), rebuilt every step.
//
// The old x-axis sweep emitted 28k candidate pairs on a 40-layer resting wall,
// and the Morton LBVH that replaced it partitioned in 3D at the cost of a tree
// build and a recursive self-query. The hash bucketing also partitions in 3D,
// but with none of either: bodies go into dyadic levels by bounding-sphere
// size, each level is a sorted key table, and candidate pairs come from a
// 27-cell walk. The candidate set is a superset of the overlapping-sphere set,
// exactly as the tree's AABB overlap set was.
//
// Determinism: pair collection is integer-sorted and de-duplicated, and the
// final (min-index, other) sort fixes the manifold creation order - the same
// contract as the sweep and tree paths before it.
void Solver::broadPhase()
{
    bodiesInOrder.clear();
    bodiesInOrder.insert(bodiesInOrder.end(), next_range(bodies).begin(), next_range(bodies).end());

    const int count = static_cast<int>(bodiesInOrder.size());
    if (count < 2)
        return;

    // Cache each body's tight AABB half extents: the exact candidate filter runs on
    // every pair, so recomputing the rotation matrix there would dominate the phase.
    bodyAabbHalf.resize(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i)
        bodyAabbHalf[i] = aabbHalfExtent(*bodiesInOrder[i]);

    broadHash.buildAndCollectPairs(bodiesInOrder, 0.0f, sweepPairs);

    std::sort(sweepPairs.begin(), sweepPairs.end());
    sweepPairs.erase(std::unique(sweepPairs.begin(), sweepPairs.end()), sweepPairs.end());
    for (const auto &[i, j] : sweepPairs)
    {
        Rigid *bodyA = bodiesInOrder[i];
        Rigid *bodyB = bodiesInOrder[j];
        if (pairOverlaps(bodyA, bodyAabbHalf[i], bodyB, bodyAabbHalf[j]))
            new Manifold(this, bodyA, bodyB);
    }
}

// A constraint is frozen when neither of its bodies can move this step. Joints and springs
// stay out of this path: their initialize() does more than refresh contacts (initial
// error, spring setup), and the core suite pins their per-scene force counts. Manifolds
// between parked bodies (static or already asleep) are the pure case - the colouring
// already skips such bodies in the primal pass - and freezing their bookkeeping is the
// sleep win: a settled scene stops paying SAT refresh + dual updates for all its contacts.
static bool constraintFrozen(const Force *f)
{
    if (!f->isManifold)
        return false; // joints / springs / soft constraints keep their exact behaviour
    const auto awake = [](const Rigid *b) { return b != nullptr && b->mass > 0.0f && !b->sleeping; };
    return !awake(f->bodyA) && !awake(f->bodyB);
}

// Bring every force up to date for this step, and drop the ones that have gone inactive.
//
// Initialising is independent work - a force only reads the body poses and writes its own
// warm-started state - so it runs in parallel. Deleting is not: an inactive force (an empty
// manifold, a broken joint) unlinks itself from the solver and from both bodies' lists, so
// that mutation happens in one thread, in list order, once the parallel pass has finished.
//
// Returns the number of forces left, and leaves `forceOrder` indexing exactly those.
int Solver::warmstartForces()
{
    int forceCount = collectForces();
    forceActive.resize(forceCount);
    pool->forCount(forceCount, [this](int i) {
        Force *f = forceOrder[i];
        // Frozen constraints keep their state verbatim; nothing they read can move.
        forceActive[i] = constraintFrozen(f) || f->initialize() ? 1 : 0;
    }, kMinForcesPerDispatch);

    for (int i = 0; i < forceCount; i++)
    {
        if (!forceActive[i])
            delete forceOrder[i]; // unlinks itself from the solver and both bodies
    }

    forceCount = collectForces(); // the sweep above invalidated the index
    return forceCount;
}

// Compute every body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so it runs in parallel. It covers static bodies as well: their
// poses do not move, but the constraint code reads their cached `initialLin`/`initialAng`.
void Solver::warmstartBodies()
{
    const int count = static_cast<int>(warmstartOrder.size());
    pool->forCount(count, [this](int i) {
        Rigid *body = warmstartOrder[i];

        // Axis locks (Godot BodyAxis): a locked axis cannot carry velocity into the step,
        // and gravity does not act along a locked linear axis. Constraint coupling can
        // still displace a locked axis within a step; documented approximation.
        if (body->axisLockLinear)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockLinear & (1 << k))
                    body->velocityLin[k] = 0.0f;
            }
        }
        if (body->axisLockAngular)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockAngular & (1 << k))
                    body->velocityAng[k] = 0.0f;
            }
        }

        // Gravity runs along -Z only, so the lock on bit 2 (z) gates it entirely.
        const float g = (body->axisLockLinear & 0x4) ? 0.0f : body->gravity;

        body->inertialLin = body->positionLin + body->velocityLin * dt;
        if (body->mass > 0)
            body->inertialLin += float3{0, 0, g} * (dt * dt);
        body->inertialAng = body->positionAng + body->velocityAng * dt;

        // Adaptive warmstart (See original VBD paper)
        const float3 accel = (body->velocityLin - body->prevVelocityLin) / dt;
        const float accelExt = accel.z * std::copysign(1.0f, g);
        float accelWeight = std::clamp(accelExt / std::fabs(g), 0.0f, 1.0f);
        if (!std::isfinite(accelWeight))
            accelWeight = 0.0f;

        // Save the position at the start of the step (x-) and warm-start to the tentative
        // new position (See original VBD paper)
        body->initialLin = body->positionLin;
        body->initialAng = body->positionAng;
        if (body->mass > 0)
        {
            body->positionLin = body->positionLin + body->velocityLin * dt + float3{0, 0, g} * (accelWeight * dt * dt);
            body->positionAng = body->positionAng + body->velocityAng * dt;
        }
    }, kMinBodiesPerDispatch);
}

// The solver's main loop: `iterations` rounds of a primal pass followed by a dual pass.
void Solver::solveIterations(int forceCount)
{
    int targetIterations = iterations;

    if (adaptiveIterations)
    {
        // Estimate total penetration error from current manifolds
        float totalPenetration = 0.0f;
        int penetrationCount = 0;
        for (Force *f = forces; f != nullptr; f = f->next)
            f->accumulatePenetration(totalPenetration, penetrationCount);

        // Adaptive iteration count: more iterations needed if penetration is large
        if (penetrationCount > 0)
        {
            float avgPenetration = totalPenetration / penetrationCount;
            targetIterations = static_cast<int>(std::round(iterations * avgPenetration / maxPenetrationError));
        }
        targetIterations = std::max(2, targetIterations);  // Minimum 2 iterations
        targetIterations = std::min(targetIterations, adaptiveMaxIterations); // Impact-moment cost bound; see the member comment.
    }

    // One dispatch covers every round: the barrier fences inside `iterate` carry the
    // Gauss-Seidel ordering, so the loop body below runs the whole set, not one round.
    // newtonRatio splits the set: the leading rounds run primal + dual, the tail runs
    // dual-only relaxation (optionally with penalty decay via stiffnessDecay).
    const int newtonRounds = std::clamp(
            static_cast<int>(std::ceil(newtonRatio * static_cast<float>(targetIterations))), 0, targetIterations);
    iterate(targetIterations, forceCount, newtonRounds);
}

// One worker-loop body for the whole iteration set: `iterations` rounds of the
// colour-parallel Gauss-Seidel primal followed by a dual pass, with a barrier
// fence at every phase boundary.
//
// The primal runs colour by colour: every body of colour c assembles its full
// block system against the poses the earlier colours already applied, solves it
// and applies the update immediately - Gauss-Seidel across colours, parallel
// and race-free within one (a colour is an independent set of the constraint
// graph by construction). Manifolds take the self-contained updatePrimal path,
// which recomputes the shared assembly against the live pose; each manifold is
// seen by at most one body per colour, so nothing races.
//
// The colouring is capped at kMaxColours so the barrier count per round is
// bounded; bodies the greedy colouring could not fit run one block-Jacobi round
// after the colours (precompute, solve against the post-colour snapshot,
// combine), using jacobiRelaxation/Chebyshev as the combine controls.
//
// Both sync types run the identical body: NullSync (small scenes, or another solver
// holding the loop) executes it inline in list order, which is exactly the threads=1
// reference the digest tests compare against.
void Solver::iterate(int targetIterations, int forceCount, int newtonRounds)
{
    // Convergence early exit: measure the residual at the end of a Newton round; when it
    // stops improving (delta below the threshold or rising), the remaining rounds would
    // burn cycles on float noise. Measured only at the Newton boundary - the impulse tail
    // converges slowly by design, so its delta is checked against its own, much looser
    // scale. 0 disables the check.
    float prevResidual = 0.0f;
    int roundsRun = 0;
    std::atomic<bool> converged{false};
    const auto loop = [this, targetIterations, forceCount, newtonRounds, &prevResidual, &roundsRun, &converged](int workerIndex, auto &sync) {
        // Per-worker Chebyshev schedule state for the overflow bucket's combine:
        // the recurrence is deterministic, so every worker derives the same
        // weights without sharing anything.
        int chebRound = 0;
        float chebGammaPrev = 1.0f;
        for (int it = 0; it < targetIterations; it++)
        {
            const bool newton = it < newtonRounds;
            if (newton)
            {
                // Colour-parallel Gauss-Seidel. One pass per colour: assemble
                // (manifolds stamp self-contained against the live poses, other
                // forces stamp directly), solve, apply. Bodies of one colour
                // share no force, so the pass is race-free without a snapshot.
                for (int c = 0; c < colourCount; c++)
                {
                    const int first = colourBegins[c];
                    const int last = colourBegins[c + 1];
                    if (last == first)
                        continue; // empty colour: every worker skips the barrier alike
                    sync.forItems(last - first, [this, first, it](int k) {
                        Rigid *body = updateOrder[colourOrder[first + k]];
                        Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
                                body->positionAng, body->inertialAng);
                        for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
                            f->updatePrimal(body, alpha, block);
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
                const int overflowCount = static_cast<int>(overflowOrder.size());
                if (overflowCount > 0)
                {
                    sync.forItems(forceCount, [this](int i) {
                        Force *f = forceOrder[i];
                        if (!f->isManifold)
                            return;
                        Manifold *manifold = static_cast<Manifold *>(f);
                        const bool aIn = manifold->bodyA->inUpdateSet &&
                                bodyColour[manifold->bodyA->updateIndex] == kOverflowColour;
                        const bool bIn = manifold->bodyB->inUpdateSet &&
                                bodyColour[manifold->bodyB->updateIndex] == kOverflowColour;
                        if (!aIn && !bIn)
                            return; // both sides solved in a colour, or frozen
                        manifold->preCompute(alpha);
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
                        for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
                        {
                            if (f->isManifold)
                                static_cast<Manifold *>(f)->stampCached(body, block);
                            else
                                f->updatePrimal(body, alpha, block);
                        }
                        block.solveDelta(bodyDeltaLin[i], bodyDeltaAng[i]);
                    });
                    sync.arriveAndWait();

                    // Combine: each bucket body applies its own delta. With
                    // chebyshevRho > 0 the plain Jacobi estimate x + delta is
                    // convex-combined with the previous iterate through the
                    // Chebyshev weights (semi-iteration); otherwise the constant
                    // jacobiRelaxation scales it. The angular delta gets the same
                    // gamma but no prev-mixing: positions are a vector space,
                    // rotations are not, and the first-order quat integration
                    // stays exact as written.
                    float relaxGamma = jacobiRelaxation;
                    if (chebyshevRho > 0.0f)
                    {
                        const float rho2 = chebyshevRho * chebyshevRho;
                        if (chebRound == 0)
                            relaxGamma = 1.0f;
                        else if (chebRound == 1)
                            relaxGamma = 1.0f / (1.0f - 0.5f * rho2);
                        else
                            relaxGamma = 1.0f / (1.0f - 0.25f * rho2 * chebGammaPrev);
                        chebGammaPrev = relaxGamma;
                    }
                    chebRound++;
                    sync.forItems(overflowCount, [this, relaxGamma](int k) {
                        const int i = overflowOrder[k];
                        Rigid *body = updateOrder[i];
                        if (chebyshevRho <= 0.0f || !body->contactOnly)
                        {
                            // Plain relaxed Jacobi - the constant, or the safe fallback for
                            // bodies coupled by joints (extrapolation diverges those).
                            body->positionLin = body->positionLin + bodyDeltaLin[i] * jacobiRelaxation;
                            body->positionAng = body->positionAng + bodyDeltaAng[i] * jacobiRelaxation;
                            return;
                        }
                        // x_{k+1} = gamma * (x_k + delta_k) + (1 - gamma) * x_{k-1}
                        const float3 xLin = body->positionLin;
                        body->positionLin = relaxGamma * (xLin + bodyDeltaLin[i]) + (1.0f - relaxGamma) * bodyPrevLin[i];
                        bodyPrevLin[i] = xLin;
                        body->positionAng = body->positionAng + bodyDeltaAng[i] * relaxGamma;
                    });
                    sync.arriveAndWait();
                }
                else
                {
                    chebRound++;
                }
            }

            // Dual update: one item per force, and each force only reads the (frozen)
            // body poses and writes its own dual state, so the whole pass is independent.
            // Frozen constraints have nothing to update: every input they read is parked.
            // In the impulse tail the contact penalties decay first (stiffnessDecay < 1):
            // relaxation converges poorly against a stiff ramp, softening lets the cheap
            // rounds finish the job the Newton rounds started - fused into the same pass,
            // same force, same order, so the decay lands before the dual reads the
            // penalty without paying an extra pass and barrier per tail round.
            const float decay = stiffnessDecay;
            sync.forItems(forceCount, [this, newton, decay](int i) {
                Force *f = forceOrder[i];
                if (!newton && decay < 1.0f && f->isManifold) {
                    Manifold *m = static_cast<Manifold *>(f);
                    for (int ci = 0, n = m->numContacts; ci < n; ++ci)
                        m->contacts[ci].penalty = min(m->contacts[ci].penalty * decay, float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});
                }
                if (!constraintFrozen(f))
                    f->updateDual(alpha);
            });
            sync.arriveAndWait();
            ++roundsRun;

            // Convergence probe after a Newton round: a full primal pass is expensive, so
            // skip further rounds when the residual stops shrinking. Only the calling
            // thread (worker 0) evaluates the reduction, and the verdict is stored in a
            // shared atomic that EVERY participant reads BEFORE the barrier - each worker
            // that sees `converged` skips its own remaining work but still participates
            // in every barrier (draining instead of returning), so the barrier
            // participation count stays equal and no worker deadlocks.
            if (workerIndex == 0 && convergenceThreshold > 0.0f && it + 1 < targetIterations)
            {
                float total = 0.0f;
                int count = 0;
                for (Force *f = forces; f != nullptr; f = f->next)
                    f->accumulatePenetration(total, count);
                const float residual = count > 0 ? total / count : 0.0f;
                converged.store(prevResidual - residual < convergenceThreshold,
                        std::memory_order_relaxed);
                prevResidual = residual;
            }
            const bool stop = converged.load(std::memory_order_relaxed);
            if (stop)
                break; // every participant breaks at the same iteration - barriers stay paired
        }
    };

    // Below both break-even sizes every phase would run inline anyway, so skip the
    // dispatch entirely. Above them the loop dispatches once for the whole set.
    const bool tiny = static_cast<int>(updateOrder.size()) < kMinBodiesPerDispatch
            && forceCount < kMinForcesPerDispatch;
    if (tiny)
    {
        detail::NullSync sync;
        loop(0, sync);
    }
    else
    {
        pool->runLoop(loop);
    }
}

// BDF1 velocities after the final iteration. Per-body and independent; it covers static
// bodies too, so their bookkeeping matches, though only movable ones get a new velocity.
void Solver::finishVelocities()
{
    const int count = static_cast<int>(warmstartOrder.size());
    pool->forCount(count, [this](int i) {
        Rigid *body = warmstartOrder[i];
        body->prevVelocityLin = body->velocityLin;
        // Don't update velocity if sleeping (Godot 4.7 API)
        if (body->mass > 0 && !body->sleeping)
        {
            body->velocityLin = (body->positionLin - body->initialLin) / dt;
            body->velocityAng = (body->positionAng - body->initialAng) / dt;
        }
        // A locked axis reports no velocity even if constraint coupling nudged it: the
        // lock's promise is "this axis does not move" to whoever reads state back.
        if (body->axisLockLinear)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockLinear & (1 << k))
                    body->velocityLin[k] = 0.0f;
            }
        }
        if (body->axisLockAngular)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockAngular & (1 << k))
                    body->velocityAng[k] = 0.0f;
            }
        }
        // Idle detection: a movable, non-sleeping body whose velocities stayed under both
        // thresholds for sleepFrames consecutive steps falls asleep (velocity zeroed,
        // Godot semantics). Only touches this body's fields, so it is safe inside the
        // parallel dispatch and deterministic regardless of thread count.
        // Automatic sleep only for bodies whose host allows it (sleep_mode == SLEEP).
        // A body with can_sleep=false used to enter idle sleep here anyway; the next
        // selection then simulated it again (sleep_mode NEVER), but finishVelocities
        // skipped its velocity update because `sleeping` was still set - the body
        // moved with a permanently zero velocity, and contacts could pump it upward
        // indefinitely. That was the pyramid_dbl runaway.
        if (sleepFrames > 0 && body->mass > 0 && !body->sleeping && body->sleep_mode == 1)
        {
            if (length2(body->velocityLin) < sleepLinearThreshold * sleepLinearThreshold &&
                    length2(body->velocityAng) < sleepAngularThreshold * sleepAngularThreshold)
                ++body->stillFrames;
            else
                body->stillFrames = 0;
            if (body->stillFrames >= sleepFrames)
            {
                body->sleeping = true;
                body->velocityLin = float3{0, 0, 0};
                body->velocityAng = float3{0, 0, 0};
            }
        }
    }, kMinBodiesPerDispatch);
}

// Island-wide sleep bookkeeping. Movable, sleep-enabled bodies are grouped into islands
// through force connectivity (union-find over the flattened body list; static bodies and
// sleep_mode == NEVER bodies join nothing and never propagate). An island where every
// member cleared idle detection this step (stillFrames >= sleepFrames or already
// sleeping) is put to sleep as a unit: Godot semantics, and it removes the
// partially-assembled-chain edge case where the bottom box qualifies a step before the
// boxes above it. Deterministic: the union-find walks the force list in list order and
// the minimum-body-index island representative is independent of visit order.
void Solver::islandSleep()
{
    // Only bodies that can sleep join an island. With the default sleepFrames = 0 this
    // whole pass is skipped by the caller, so the core suite (which never sets
    // sleep_mode) is unaffected.
    islandBodies.clear();
    for (Rigid *body = bodies; body != 0; body = body->next)
        if (body->mass > 0 && body->sleep_mode == 1)
            islandBodies.push_back(body);

    const int count = static_cast<int>(islandBodies.size());
    if (count < 2)
        return;

    std::vector<int> parent(count);
    for (int i = 0; i < count; ++i)
        parent[i] = i;
    auto find = [&parent](int x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };

    // Slot lookup: map body pointer -> island index once (sorted pointers + binary
    // search), then each force is two O(log n) lookups. The map is rebuilt each call:
    // bodies enter and leave the island set as sleep_mode and mass change.
    islandForceSlots.clear();
    std::vector<Rigid *> sorted = islandBodies;
    std::sort(sorted.begin(), sorted.end());
    for (Force *f = forces; f != 0; f = f->next)
    {
        const auto slotOf = [&sorted](Rigid *b) -> int {
            const auto it = std::lower_bound(sorted.begin(), sorted.end(), b);
            return (it != sorted.end() && *it == b) ? static_cast<int>(it - sorted.begin()) : -1;
        };
        const int ia = slotOf(f->bodyA);
        const int ib = slotOf(f->bodyB);
        if (ia >= 0 && ib >= 0)
            islandForceSlots.emplace_back(ia, ib);
    }

    for (const auto &[ia, ib] : islandForceSlots)
        parent[find(ia)] = find(ib);

    // Compose islands: representative -> aggregate state.
    // islandAwake[r] = true iff any member is awake; islandStill[r] = min stillFrames.
    std::vector<char> islandAwake(count, 0);
    std::vector<int> islandStill(count, INT_MAX);
    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (!islandBodies[i]->sleeping)
            islandAwake[r] = 1;
        islandStill[r] = std::min(islandStill[r], islandBodies[i]->stillFrames);
    }

    // Sleep whole islands whose every member has cleared idle detection; leave mixed or
    // awake islands exactly as they are (the per-body wake-up pass above already ran).
    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (islandAwake[r] || islandStill[r] < sleepFrames)
            continue;
        Rigid *body = islandBodies[i];
        if (!body->sleeping)
        {
            body->sleeping = true;
            body->velocityLin = float3{0, 0, 0};
            body->velocityAng = float3{0, 0, 0};
        }
    }
}

// The primal half of one iteration for a single body (Eqs. 4-6): assemble the 6x6 system
// from every force acting on the body, then solve it and move the body.
//
// Solve the assembled system and move the body by the result (Eq. 4).
void Block::apply(Rigid &body) const noexcept
{
    float3 dxLin, dxAng;
    solveDelta(dxLin, dxAng);
    body.positionLin = body.positionLin + dxLin;
    body.positionAng = body.positionAng + dxAng;
}

// Colour-parallel Gauss-Seidel: bodies sharing a force never solve at the same
// instant, so each colour is an independent set - one barrier per colour carries
// the Gauss-Seidel ordering, and the result is deterministic for a fixed
// colouring no matter the thread count.
//
// The colouring is capped at kMaxColours: a body whose neighbours already used
// every colour goes to the overflow bucket, which the primal runs as one
// block-Jacobi round (every delta solved against the same post-colour snapshot,
// then combined). The bucket costs three barriers when non-empty and is empty in
// every ordinary scene; it exists so a pathological clique cannot make the
// barrier count unbounded.
//
// The wake-up and island-sleep logic is unchanged from the colouring era: a sleeping body
// sharing a force with a moving movable body re-enters the update set, and settled
// islands sleep together.
void Solver::selectUpdateBodies()
{
    warmstartOrder.clear();
    updateOrder.clear();

    // Wake-up pass: a sleeping body (sleep_mode enabled) sharing a force with a moving
    // movable body is woken before selection, so it re-enters this step's update set.
    // A neighbour that is itself below the idle thresholds is not a disturbance: stacks
    // must be allowed to fall asleep together, one body at a time. Serial and read-only
    // over the force adjacency apart from the examined body itself, so the outcome is
    // fixed regardless of thread count.
    if (sleepFrames > 0)
    {
        const float lin2 = sleepLinearThreshold * sleepLinearThreshold;
        const float ang2 = sleepAngularThreshold * sleepAngularThreshold;
        for (Rigid *body = bodies; body != 0; body = body->next)
        {
            if (!body->sleeping || body->sleep_mode != 1)
                continue;
            for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
            {
                Rigid *other = (f->bodyA == body) ? f->bodyB : f->bodyA;
                if (other != 0 && other->mass > 0 && !other->sleeping &&
                        (length2(other->velocityLin) > lin2 ||
                                length2(other->velocityAng) > ang2))
                {
                    body->sleeping = false;
                    body->stillFrames = 0;
                    break;
                }
            }
        }
        islandSleep();
    }

    // Every body is warmed up; only movable ones are updated. Static bodies are excluded
    // from the update set: they are never moved, so nothing solves for them.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        warmstartOrder.push_back(body);

        // Fresh per-step state; admitted bodies flip inUpdateSet below.
        body->inUpdateSet = false;
        body->updateIndex = -1;
        body->contactOnly = true;

        // Sleep mode support (Godot 4.7 API)
        // SLEEP_MODE_NEVER: always update
        // SLEEP_MODE_SLEEP: update unless idle detection has put it to sleep
        // SLEEP_MODE_START_IN_SLEEP: start in sleep
        if (body->mass > 0)
        {
            // Check if body should be updated based on sleep mode
            bool should_update = true;
            if (body->sleep_mode == 1) // SLEEP_MODE_SLEEP
            {
                // Idle detection (finishVelocities) is the only way into sleep here;
                // a sleeping body with contacts stays out of the update set until the
                // wake-up pass above re-admits it.
                if (body->sleeping)
                    should_update = false;
            }
            else if (body->sleep_mode == 2) // SLEEP_MODE_START_IN_SLEEP
            {
                // Start in sleep, wake up if velocity or forces are present
                if (!body->sleeping && 
                    length2(body->velocityLin) == 0.0f &&
                    body->forces == nullptr)
                {
                    should_update = false;
                }
                else
                {
                    body->sleeping = false; // Wake up
                }
            }
            
            if (should_update)
            {
                body->updateIndex = static_cast<int>(updateOrder.size());
                updateOrder.push_back(body);
                body->inUpdateSet = true;
            }
        }
        else
        {
            // Static body: excluded from the update set, so nothing solves for it.
        }
    }

    // Per-body delta buffers for the Jacobi primal's split-solve / combine pair, plus
    // the previous-iterate buffer the Chebyshev combination reads.
    bodyDeltaLin.resize(updateOrder.size());
    bodyDeltaAng.resize(updateOrder.size());
    bodyPrevLin.resize(updateOrder.size());

    // Chebyshev eligibility: a body carrying any non-manifold force (joint, spring, world
    // anchor) extrapolates unsafely - its gamma falls back to the constant relaxation in
    // the combine pass. Walk the force lists once, after warmstartForces() has pruned.
    for (Force *f = forces; f != 0; f = f->next)
    {
        if (f->isManifold)
            continue;
        if (f->bodyA != 0)
            f->bodyA->contactOnly = false;
        if (f->bodyB != 0)
            f->bodyB->contactOnly = false;
    }

    colourGraph();
}

// Greedy colouring of the update set's constraint graph, in updateOrder order:
// each body takes the smallest colour none of its already-coloured neighbours
// holds, capped at kMaxColours. The order is the body list order, so the colour
// assignment is a pure function of the scene - deterministic regardless of
// thread count. Neighbours are found through the bodies' own force lists; a
// body not yet coloured (later in the order) constrains nothing, and overflow
// neighbours never block a colour because the overflow bucket runs Jacobi
// against a snapshot and is safe to overlap with any colour.
void Solver::colourGraph()
{
    const int n = static_cast<int>(updateOrder.size());
    bodyColour.assign(n, kOverflowColour);
    colourOrder.clear();
    overflowOrder.clear();

    uint8_t used[kMaxColours];
    for (int i = 0; i < n; i++)
    {
        Rigid *body = updateOrder[i];
        for (int c = 0; c < kMaxColours; c++)
            used[c] = 0;

        for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
        {
            Rigid *other = (f->bodyA == body) ? f->bodyB : f->bodyA;
            if (other == 0 || other->updateIndex < 0 || other->updateIndex >= i)
                continue; // not in the update set, or coloured later in the pass
            const uint8_t c = bodyColour[other->updateIndex];
            if (c < kMaxColours)
                used[c] = 1;
        }

        uint8_t colour = kOverflowColour;
        for (int c = 0; c < kMaxColours; c++)
            if (!used[c])
            {
                colour = static_cast<uint8_t>(c);
                break;
            }
        bodyColour[i] = colour;
    }

    // Counting sort into per-colour runs (stable in updateOrder order).
    colourBegins.assign(kMaxColours + 1, 0);
    for (int i = 0; i < n; i++)
        if (bodyColour[i] < kMaxColours)
            colourBegins[bodyColour[i] + 1]++;
        else
            overflowOrder.push_back(i);
    for (int c = 0; c < kMaxColours; c++)
        colourBegins[c + 1] += colourBegins[c];
    colourOrder.resize(n - static_cast<int>(overflowOrder.size()));
    std::array<int, kMaxColours> cursor{};
    for (int c = 0; c < kMaxColours; c++)
        cursor[c] = colourBegins[c];
    for (int i = 0; i < n; i++)
        if (bodyColour[i] < kMaxColours)
            colourOrder[cursor[bodyColour[i]]++] = i;

    colourCount = 0;
    for (int c = 0; c < kMaxColours; c++)
        if (colourBegins[c + 1] > colourBegins[c])
            colourCount = c + 1;

    // The overflow bucket's Jacobi round needs the delta buffers.
    bodyDeltaLin.resize(updateOrder.size());
    bodyDeltaAng.resize(updateOrder.size());
    bodyPrevLin.resize(updateOrder.size());
}

// Fill `forceOrder` with the solver's forces, so the phases can address them by index.
int Solver::collectForces()
{
    forceOrder.clear();
    forceOrder.insert(forceOrder.end(), next_range(forces).begin(), next_range(forces).end());
    return static_cast<int>(forceOrder.size());
}

void Solver::ensurePool()
{
    if (pool != nullptr)
        return;

    // Pass the value straight through: 0 means "one worker per hardware thread" and 1 means
    // a single worker, which both the pool and `threads` document the same way. (Clamping 0
    // up to 1 here would silently turn the default into a serial run.)
    pool = std::make_unique<detail::JobPool>(static_cast<unsigned>(std::max(0, threads)));
}

int Solver::threadCount() const
{
    if (pool == nullptr)
        return 1;
    return static_cast<int>(pool->threadCount());
}

} // namespace avbd
