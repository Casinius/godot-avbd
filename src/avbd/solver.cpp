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
#include "avbd/bvh/bvh.hpp"
#include "avbd/bvh/node_storage.hpp"

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

    // Exact per-body test: ray against the body's local box, as before. Used at BVH
    // leaves and by the no-tree fallback below, so both paths agree bit for bit.
    const auto testBody = [&](Rigid *body) {
        if ((body->collisionLayer & p_mask) == 0)
            return;
        quat invRot = body->positionAng.conjugate();
        float3 o = invRot * (origin - body->positionLin);
        float3 d = invRot * dir;
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

    if (bvhNodes.size() == 0)
    {
        // No broad phase yet (no step has run): walk the list, as this always did.
        // Unlike the upstream demo (which picks bodies for dragging) this also tests
        // static bodies, so rays hit the ground and ramps.
        for (Rigid *body = bodies; body != 0; body = body->next)
            testBody(body);
    }
    else
    {
        // Descend the last built broad-phase tree. Node bounds are the tight
        // AABB-of-OBB leaves aggregated upward, so a ray/AABB interval test prunes
        // whole subtrees; pruning by the current best hit keeps the query output-
        // sensitive. The exact per-body test above decides at the leaves.
        const auto rayOverNode = [&](int slot, float &r_tEnter) {
            const float3 &bmin = bvhNodes.boundsMin()[slot];
            const float3 &bmax = bvhNodes.boundsMax()[slot];
            float tEnter = 0.0f;
            float tExit = bestT; // prune against the closest hit so far
            for (int i = 0; i < 3; ++i)
            {
                if (std::fabs(dir[i]) < epsilon)
                {
                    if (origin[i] < bmin[i] || origin[i] > bmax[i])
                        return false;
                    continue;
                }
                const float invD = 1.0f / dir[i];
                float t0 = (bmin[i] - origin[i]) * invD;
                float t1 = (bmax[i] - origin[i]) * invD;
                if (t0 > t1)
                {
                    const float tmp = t0;
                    t0 = t1;
                    t1 = tmp;
                }
                tEnter = std::max(tEnter, t0);
                tExit = std::min(tExit, t1);
                if (tEnter > tExit)
                    return false;
            }
            r_tEnter = tEnter;
            return true;
        };

        // Complete-tree depth bounds the descent (one pending node per level); 64
        // covers n up to 2^31 with room to spare.
        FixedStack<int, 64> stack;
        stack.push(0);
        while (!stack.empty())
        {
            const int slot = stack.pop();
            if (!bvhNodes.valid(slot))
                continue;
            float tEnter = 0.0f;
            if (!rayOverNode(slot, tEnter))
                continue;
            if (2 * slot + 1 < bvhNodes.size())
            {
                // Internal node (its children exist in the implicit layout).
                stack.push(2 * slot + 1);
                stack.push(2 * slot + 2);
                continue;
            }
            const int bodyIndex = bvhNodes.bodyIndex()[slot];
            if (bodyIndex < 0)
                continue;
            testBody(bodiesInOrder[bodyIndex]);
        }
    }

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
}

void Solver::collectAABB(float3 min, float3 max, uint32_t p_mask, std::vector<Rigid *> &out) const
{
    if (bvhNodes.size() == 0)
    {
        // No tree yet: the list walk answers with the same overlap rule.
        for (Rigid *body = bodies; body != 0; body = body->next)
        {
            if ((body->collisionLayer & p_mask) == 0)
                continue;
            // The body's own bound: centre plus its tight world half-extent (the same
            // formula the broad phase builds leaves with, evaluated inline).
            float3 bmin, bmax;
            float3 half;
            switch (body->shape)
            {
                case ShapeType::Sphere:
                    half = float3{body->radius, body->radius, body->radius};
                    break;
                case ShapeType::Cylinder:
                    half = float3{body->size.x(), body->size.z() * 0.5f, body->size.x()};
                    break;
                case ShapeType::Box:
                default:
                    half = body->size * 0.5f;
                    break;
            }
            const float3x3 rot(body->positionAng);
            const float3 wx{
                    std::abs(rot(0, 0)) * half.x() + std::abs(rot(0, 1)) * half.y() + std::abs(rot(0, 2)) * half.z(),
                    std::abs(rot(1, 0)) * half.x() + std::abs(rot(1, 1)) * half.y() + std::abs(rot(1, 2)) * half.z(),
                    std::abs(rot(2, 0)) * half.x() + std::abs(rot(2, 1)) * half.y() + std::abs(rot(2, 2)) * half.z()};
            bmin = body->positionLin - wx;
            bmax = body->positionLin + wx;
            if (!(bmin.x() > max.x() || bmin.y() > max.y() || bmin.z() > max.z() ||
                        min.x() > bmax.x() || min.y() > bmax.y() || min.z() > bmax.z()))
                out.push_back(body);
        }
        return;
    }

    const auto overlaps = [&](int slot) {
        const float3 &bmin = bvhNodes.boundsMin()[slot];
        const float3 &bmax = bvhNodes.boundsMax()[slot];
        return !(bmin.x() > max.x() || bmin.y() > max.y() || bmin.z() > max.z() ||
                min.x() > bmax.x() || min.y() > bmax.y() || min.z() > bmax.z());
    };

    FixedStack<int, 64> stack;
    stack.push(0);
    while (!stack.empty())
    {
        const int slot = stack.pop();
        if (!bvhNodes.valid(slot) || !overlaps(slot))
            continue;
        if (2 * slot + 1 < bvhNodes.size())
        {
            stack.push(2 * slot + 1);
            stack.push(2 * slot + 2);
            continue;
        }
        const int bodyIndex = bvhNodes.bodyIndex()[slot];
        if (bodyIndex < 0)
            continue;
        Rigid *body = bodiesInOrder[bodyIndex];
        if ((body->collisionLayer & p_mask) != 0)
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

// A broad-phase candidate test: bounding spheres plus Godot's layer/mask pair rule,
// minus pairs already joined by a force (IgnoreCollision, joints).
static bool pairOverlaps(Rigid *bodyA, Rigid *bodyB)
{
    const float3 dp = bodyA->positionLin - bodyB->positionLin;
    const float reach = bodyA->radius + bodyB->radius;
    // Godot semantics: the pair collides when either side's layer is in the other's
    // mask (defaults 1/1 keep every pair, as before these fields existed).
    const bool layersAllow = ((bodyA->collisionLayer & bodyB->collisionMask) != 0) ||
            ((bodyB->collisionLayer & bodyA->collisionMask) != 0);
    return dp.squaredNorm() <= reach * reach && layersAllow && !bodyA->constrainedTo(bodyB);
}

// Contact detection: Morton LBVH (implicit heap layout), rebuilt every step.
//
// The old x-axis sweep emitted 28k candidate pairs on a 40-layer resting wall (x cannot
// discriminate a z-pile); the Morton tree partitions in 3D so the candidate list tracks
// the real contact surface. Bodies are ordered by 64-bit Morton code - a pure function
// of the scene bounds and positions - and the code-sorted order is mapped into the
// implicit heap layout of a complete binary tree: node i has children 2i+1 / 2i+2, no
// left/right/parent arrays, best-possible cache behaviour for a pointer-free BVH.
//
// Determinism: leaf order is the Morton-sorted body order, pair collection is leaf
// order, and the final pair list is sorted ascending (min-index, other) - the same
// manifold creation order contract as the sweep path.
void Solver::broadPhase()
{
    bodiesInOrder.clear();
    bodiesInOrder.insert(bodiesInOrder.end(), next_range(bodies).begin(), next_range(bodies).end());

    const int count = static_cast<int>(bodiesInOrder.size());
    if (count < 2)
        return;

    std::vector<int> indices(count);
    for (int i = 0; i < count; ++i)
        indices[i] = i;

    std::vector<const Rigid *> bodiesPtr;
    bodiesPtr.reserve(count);
    for (Rigid *b : bodiesInOrder)
        bodiesPtr.push_back(b);

    bvh::builder::Builder builder(bvhNodes, indices, bodiesPtr);
    const int leafCount = builder.buildLBVH();
    if (leafCount < 2)
        return;
    const int padded = 1;
    int p = 1;
    while (p < leafCount)
        p <<= 1;
    const int firstLeaf = p - 1;
    (void)padded;

    // Leaf slots [firstLeaf, firstLeaf + leafCount) in Morton-sorted order. The pair
    // search is one top-down self-query: collideSubtrees walks every leaf pair between
    // two subtrees whose aggregated bounds overlap (descending the internal side), and
    // selfTraverse splits each internal node into its two children and recurses - each
    // internal-node overlap test runs once per ancestor path, where the per-leaf
    // from-root descent re-tested the top of the tree once per leaf. Empty subtrees
    // (the right tail when n is not a power of two) are marked invalid by the builder
    // and pruned before their uninitialized bounds are read. The body-index minmax and
    // the final sort keep the (min-index, other) manifold order.
    sweepPairs.clear();
    const auto boundsOverlap = [&](int sa, int sb) {
        const float3 &aMin = bvhNodes.boundsMin()[sa];
        const float3 &aMax = bvhNodes.boundsMax()[sa];
        const float3 &bMin = bvhNodes.boundsMin()[sb];
        const float3 &bMax = bvhNodes.boundsMax()[sb];
        // Same touching-counts-as-overlapping semantics as the old pairwise sweep.
        return !(bMin.x() > aMax.x() || bMin.y() > aMax.y() || bMin.z() > aMax.z() ||
                aMin.x() > bMax.x() || aMin.y() > bMax.y() || aMin.z() > bMax.z());
    };
    const auto collideSubtrees = [&](auto &&self, int sa, int sb) -> void {
        if (!bvhNodes.valid(sa) || !bvhNodes.valid(sb) || !boundsOverlap(sa, sb))
            return;
        const bool aLeaf = sa >= firstLeaf;
        const bool bLeaf = sb >= firstLeaf;
        if (aLeaf && bLeaf)
        {
            sweepPairs.emplace_back(std::minmax(bvhNodes.bodyIndex()[sa], bvhNodes.bodyIndex()[sb]));
            return;
        }
        if (!aLeaf)
        {
            self(self, 2 * sa + 1, sb);
            self(self, 2 * sa + 2, sb);
        }
        else
        {
            self(self, sa, 2 * sb + 1);
            self(self, sa, 2 * sb + 2);
        }
    };
    const auto selfTraverse = [&](auto &&self, int slot) -> void {
        if (slot >= firstLeaf)
            return;
        const int l = 2 * slot + 1;
        const int r = 2 * slot + 2;
        collideSubtrees(collideSubtrees, l, r);
        self(self, l);
        self(self, r);
    };
    selfTraverse(selfTraverse, 0);

    std::sort(sweepPairs.begin(), sweepPairs.end());
    sweepPairs.erase(std::unique(sweepPairs.begin(), sweepPairs.end()), sweepPairs.end());
    for (const auto &[i, j] : sweepPairs)
    {
        Rigid *bodyA = bodiesInOrder[i];
        Rigid *bodyB = bodiesInOrder[j];
        if (pairOverlaps(bodyA, bodyB))
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
    if (f->contactPointCount() <= 0)
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
        const float accelExt = accel.z() * std::copysign(1.0f, g);
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

// One worker-loop body for the whole iteration set: `iterations` rounds of the block
// Jacobi primal (one full block solve per body against the iteration-start pose, then a
// combine pass) followed by a dual pass, with a barrier fence at every phase boundary.
//
// The shape matters more than the phases: one dispatch arms the workers once per step and
// each boundary costs a ~1-3 us barrier. The Jacobi primal needs two fences per round
// regardless of how entangled the constraint graph is - the old colour-parallel
// Gauss-Seidel needed one per colour, which is why dense stacking scenes starved it (and
// why the incremental colouring's correctness was so fragile).
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
        // Per-worker Chebyshev schedule state: the recurrence is deterministic, so every
        // worker derives the same weights without sharing anything.
        int chebRound = 0;
        float chebGammaPrev = 1.0f;
        for (int it = 0; it < targetIterations; it++)
        {
            const bool newton = it < newtonRounds;
            if (newton)
            {
                // Block projected Jacobi in three phases. Phase A precomputes each
                // manifold's shared per-contact assembly (world arms, angular jacobians,
                // clamped forces) once from the round's pose snapshot - work the two
                // body stamps used to duplicate. Nobody writes a pose until the combine
                // pass, so every phase reads the same frozen snapshot: race-free by
                // construction, no colours, and no thread-count dependence.
                sync.forItems(forceCount, [this](int i) {
                    Force *f = forceOrder[i];
                    if (f->contactPointCount() <= 0)
                        return;
                    Manifold *manifold = static_cast<Manifold *>(f);
                    if (!manifold->bodyA->inUpdateSet && !manifold->bodyB->inUpdateSet)
                        return; // frozen: neither side moves this round
                    manifold->preCompute(alpha);
                });
                sync.arriveAndWait();

                // Phase B: each body assembles its full block system (own inertia plus
                // every incident constraint - manifolds stamped from the phase-A cache,
                // other forces stamping directly) and solves for its delta.
                sync.forItems(static_cast<int>(updateOrder.size()), [this](int i) {
                    Rigid *body = updateOrder[i];
                    Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
                            body->positionAng, body->inertialAng);
                    for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
                    {
                        if (f->contactPointCount() > 0)
                            static_cast<Manifold *>(f)->stampCached(body, block);
                        else
                            f->updatePrimal(body, alpha, block);
                    }
                    block.solveDelta(bodyDeltaLin[i], bodyDeltaAng[i]);
                });
                sync.arriveAndWait();

                // Combine: each body applies its own delta. With chebyshevRho > 0 the
                // plain Jacobi estimate x + delta is convex-combined with the previous
                // iterate through the Chebyshev weights (semi-iteration); otherwise the
                // constant jacobiRelaxation scales it. The angular delta gets the same
                // gamma but no prev-mixing: positions are a vector space, rotations are
                // not, and the first-order quat integration stays exact as written. A
                // body with no constraints solves to exactly its inertial pose - free
                // flight falls out of the same path.
                float gamma = jacobiRelaxation;
                if (chebyshevRho > 0.0f)
                {
                    const float rho2 = chebyshevRho * chebyshevRho;
                    if (chebRound == 0)
                        gamma = 1.0f;
                    else if (chebRound == 1)
                        gamma = 1.0f / (1.0f - 0.5f * rho2);
                    else
                        gamma = 1.0f / (1.0f - 0.25f * rho2 * chebGammaPrev);
                    chebGammaPrev = gamma;
                }
                chebRound++;
                sync.forItems(static_cast<int>(updateOrder.size()), [this, gamma](int i) {
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
                    body->positionLin = gamma * (xLin + bodyDeltaLin[i]) + (1.0f - gamma) * bodyPrevLin[i];
                    bodyPrevLin[i] = xLin;
                    body->positionAng = body->positionAng + bodyDeltaAng[i] * gamma;
                });
                sync.arriveAndWait();
            }

            // Dual update: one item per force, and each force only reads the (frozen)
            // body poses and writes its own dual state, so the whole pass is independent.
            // Frozen constraints have nothing to update: every input they read is parked.
            // In the impulse tail, decay the contact penalties first (stiffnessDecay < 1):
            // relaxation converges poorly against a stiff ramp, softening lets the cheap
            // rounds finish the job the Newton rounds started.
            if (!newton && stiffnessDecay < 1.0f)
            {
                const float decay = stiffnessDecay;
                const float cap = PENALTY_MAX;
                sync.forItems(forceCount, [this, decay, cap](int i) {
                    Force *f = forceOrder[i];
                    // contactPointCount() > 0 iff the force is a Manifold; other forces
                    // keep their own penalty ramp untouched.
                    if (f->contactPointCount() > 0) {
                        Manifold *m = static_cast<Manifold *>(f);
                        for (int ci = 0, n = m->numContacts; ci < n; ++ci)
                            m->contacts[ci].penalty = (m->contacts[ci].penalty * decay).cwiseMin(float3{cap, cap, cap});
                    }
                });
                sync.arriveAndWait();
            }
            sync.forItems(forceCount, [this](int i) {
                if (!constraintFrozen(forceOrder[i]))
                    forceOrder[i]->updateDual(alpha);
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
        if (sleepFrames > 0 && body->mass > 0 && !body->sleeping)
        {
            if (body->velocityLin.squaredNorm() < sleepLinearThreshold * sleepLinearThreshold &&
                    body->velocityAng.squaredNorm() < sleepAngularThreshold * sleepAngularThreshold)
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

// Block Jacobi needs no colours and no mass splitting: the primal is parallel over
// bodies, each solving its full block system (own inertia plus every incident
// constraint, stamped against the iteration-start pose snapshot) and holding the delta
// for the combine pass. Rigid-body constraint graphs have no shared-vertex structure to
// split around - each manifold couples exactly two bodies - so the body-level split is
// the natural one: race-free by construction, one 6x6 solve per body per round.
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
                        (other->velocityLin.squaredNorm() > lin2 ||
                                other->velocityAng.squaredNorm() > ang2))
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
                    body->velocityLin.squaredNorm() == 0.0f && 
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
        if (f->contactPointCount() > 0)
            continue;
        if (f->bodyA != 0)
            f->bodyA->contactOnly = false;
        if (f->bodyB != 0)
            f->bodyB->contactOnly = false;
    }
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
