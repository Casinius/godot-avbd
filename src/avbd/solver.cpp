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

    // Ray-cast against each OBB by transforming the ray into body local space.
    // Unlike the upstream demo (which uses this to pick bodies for dragging) this
    // also tests static bodies, so rays hit the ground and ramps.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if ((body->collisionLayer & p_mask) == 0)
            continue;
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
            continue;

        float tHit = tEnter >= 0.0f ? tEnter : tExit;
        if (tHit < 0.0f)
            continue;

        if (tHit < bestT)
        {
            bestT = tHit;
            bestBody = body;
            bestLocal = o + d * tHit;
        }
    }

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
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

    // Broad phase is parallelized via LoopSync barrier
    if (this->threads > 1)
    {
        // Parallel mode: use LoopSync barrier
        avbd::detail::LoopSync sync(this->threads);
        broadPhase(sync);
    }
    else
    {
        // Serial mode: use NullSync (no barrier, no overhead)
        avbd::detail::NullSync sync;
        broadPhase(sync);
    }

    const int forceCount = warmstartForces();
    colourGraph();
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
template <typename SyncType>
void Solver::broadPhase(SyncType &sync)
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

    // Leaf slots [firstLeaf, firstLeaf + leafCount) in Morton-sorted order; collect them
    // and pair each against all later slots (AABB overlap prune). The final sort keeps
    // the ascending (min-index, other) manifold order.
    sweepPairs.clear();
    // Parallel AABB overlap pruning using LoopSync barrier.
    // Each worker processes a chunk of leaf pairs; barrier ensures completion before sort.
    const int totalPairs = leafCount * (leafCount - 1) / 2;
    sync.forItems(totalPairs, [this, firstLeaf, leafCount](int i) {
        // Map 1D index to 2D (a, b) with a < b
        int a, b;
        int offset = 0;
        int a_temp = 0;
        while (offset + (leafCount - a_temp - 1) <= i)
        {
            offset += leafCount - a_temp - 1;
            a_temp++;
        }
        a = a_temp;
        b = a + 1 + (i - offset);

        const int na = firstLeaf + a;
        const int nb = firstLeaf + b;
        const float3 aMin = bvhNodes.boundsMin()[na];
        const float3 aMax = bvhNodes.boundsMax()[na];
        const float3 bMin = bvhNodes.boundsMin()[nb];
        const float3 bMax = bvhNodes.boundsMax()[nb];

        // AABB overlap prune (6 comparisons, early exit on fail)
        if (bMin.x() > aMax.x() || bMin.y() > aMax.y() || bMin.z() > aMax.z() ||
            bMax.x() < aMin.x() || bMax.y() < aMin.y() || bMax.z() < aMin.z())
            return; // No overlap, skip

        const int iIdx = bvhNodes.bodyIndex()[na];
        const int jIdx = bvhNodes.bodyIndex()[nb];
        sweepPairs.emplace_back(std::minmax(iIdx, jIdx));
    });
    sync.arriveAndWait();

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
        targetIterations = std::min(targetIterations, 16); // Cap at 16 iterations (Jolt/Godot use 8, but their impulse solver is 10x cheaper per round; AVBD's dense blocks need the room - the convergence early-exit is the intended way to skip cheap rounds)
    }

    // One dispatch covers every round: the barrier fences inside `iterate` carry the
    // Gauss-Seidel ordering, so the loop body below runs the whole set, not one round.
    // newtonRatio splits the set: the leading rounds run primal + dual, the tail runs
    // dual-only relaxation (optionally with penalty decay via stiffnessDecay).
    const int newtonRounds = std::clamp(
            static_cast<int>(std::ceil(newtonRatio * static_cast<float>(targetIterations))), 0, targetIterations);
    iterate(targetIterations, forceCount, newtonRounds);
}

// One worker-loop body for the whole iteration set: `iterations` rounds of a primal pass
// per colour followed by a dual pass, with a barrier fence at every phase boundary.
//
// The shape matters more than the phases: the old code re-entered the pool once per
// colour per iteration (~45 us a dispatch), which at 8 colours x 10 iterations burned
// ~4 ms of a 16 ms frame on synchronisation alone. Here one dispatch arms the workers
// once per step and each boundary costs a ~1-3 us barrier instead.
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
        for (int it = 0; it < targetIterations; it++)
        {
            const bool newton = it < newtonRounds;
            if (newton)
            {
                // Primal pass per colour: colour groups are independent (no cross-colour writes
                // within same iteration). Each colour's updates read previous colour's dual state,
                // which is published via the barrier before this colour starts. Barrier at iteration
                // boundary ensures Gauss-Seidel ordering across Newton rounds. Parallel per-colour
                // is safe and fast: same pattern as colour-parallel Gauss-Seidel.
                for (int colour = 0; colour < colours; colour++)
                {
                    const int first = colourStart[colour];
                    const int last = colourStart[colour + 1];
                    sync.forItems(last - first, [this, first](int i) {
                        updatePrimal(updateOrder[first + i]);
                    });
                    sync.arriveAndWait();
                }
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
            // Dual pass: each force reads frozen body poses and writes its own dual state.
            // Forces are independent (no cross-force writes in updateDual), so the whole
            // pass is parallelizable. The barrier ensures all writes are visible before
            // the next phase (primal pass reads updated duals).
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
// Only ever called for movable bodies: `updateOrder` holds exactly those, which is also
// why a static body needs no guard here.
// Solve the assembled system and move the body by the result (Eq. 4).
void Block::apply(Rigid &body) const noexcept
{
    float3 dxLin, dxAng;
    solve(lhsLin, lhsAng, lhsCross, -rhsLin, -rhsAng, dxLin, dxAng);
    body.positionLin = body.positionLin + dxLin;
    body.positionAng = body.positionAng + dxAng;
}

void Solver::updatePrimal(Rigid *body)
{
    Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
            body->positionAng, body->inertialAng);

    for (Force *force = body->forces; force != 0; force = (force->bodyA == body) ? force->nextA : force->nextB)
        force->updatePrimal(body, alpha, block);

    block.apply(*body);
}

// Greedy graph colouring of the constraint graph. A body's colour must differ from every
// body it shares a force with, which is precisely the condition for two bodies to be
// updatable at the same time.
//
// Incremental: a body in this step's update set keeps its previous colour unless a
// neighbour now holds it, in which case it re-runs the first-fit search. Contacts move
// slowly between steps, so the search is rarely paid. Bodies outside the update set
// (static or sleeping) carry `colour == -1`, keeping the invariant that `colour >= 0`
// marks exactly the bodies in `updateOrder` - a stale colour is never trusted.
void Solver::colourGraph()
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
    // from the colouring: they are never moved, so they take no part in the update order
    // and cannot conflict with anything.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        warmstartOrder.push_back(body);

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
            }
            else
            {
                // Not in this step's update set: drop any stale colour so the invariant
                // `colour >= 0` <=> "in updateOrder" holds.
                body->colour = -1;
            }
        }
        else
        {
            // Static body: excluded from the update set, so it must not contribute a
            // stale colour.
            body->colour = -1;
        }
    }

    std::vector<int> taken;
    int maxColour = -1;
    for (Rigid *body : updateOrder)
    {
        // Collect the neighbours' current colours.
        taken.clear();
        for (Force *force = body->forces; force != 0; force = (force->bodyA == body) ? force->nextA : force->nextB)
        {
            Rigid *other = (force->bodyA == body) ? force->bodyB : force->bodyA;
            if (other != 0 && other->colour >= 0)
                taken.push_back(other->colour);
        }

        // Incremental recolouring: keep the previous colour when no neighbour took it.
        // A kept colour is conflict-free: neighbours processed earlier in this loop hold
        // their final colours, and a later neighbour that would collide with a kept colour
        // sees it in `taken` and refits. Only bodies whose neighbourhood actually changed
        // pay for the first-fit search.
        int colour = body->colour;
        if (colour < 0 || std::ranges::find(taken, colour) != taken.end())
        {
            colour = 0;
            while (std::ranges::find(taken, colour) != taken.end())
                colour++;
        }

        body->colour = colour;
        if (colour > maxColour)
            maxColour = colour;
    }

    colours = maxColour + 1;

    // Group the bodies by colour, so each iteration walks contiguous, independent runs.
    std::vector<int> perColour;
    perColour.reserve(colours);
    perColour.assign(colours, 0);
    for (const Rigid *body : updateOrder)
        perColour[body->colour]++;

    colourStart.assign(colours + 1, 0);
    std::exclusive_scan(perColour.begin(), perColour.end(), colourStart.begin(), 0);
    colourStart[colours] = static_cast<int>(updateOrder.size());

    widest = perColour.empty() ? 0 : *std::ranges::max_element(perColour);

    std::vector<int> cursor;
    cursor.reserve(colours + 1);
    cursor = colourStart; // cursor[c] is the next free slot for colour c
    std::vector<Rigid *> grouped;
    grouped.reserve(updateOrder.size());
    for (Rigid *body : updateOrder)
        grouped.push_back(body);
    updateOrder = std::move(grouped);
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
