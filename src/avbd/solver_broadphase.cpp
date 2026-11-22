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
#include "avbd/solver_internal.hpp"

namespace avbd {

static bool pairOverlaps(const ForceAdjacency &adjacency, const ManifoldStore &manifolds,
                         Rigid *bodyA, const float3 &aHalf, const float3 &aCenter,
                         Rigid *bodyB, const float3 &bHalf, const float3 &bCenter,
                         bool verbose = false)
{
    const float3 aMin = aCenter - aHalf;
    const float3 aMax = aCenter + aHalf;
    const float3 bMin = bCenter - bHalf;
    const float3 bMax = bCenter + bHalf;
    const bool boundsOverlap = !(bMin.x > aMax.x || bMin.y > aMax.y || bMin.z > aMax.z ||
            aMin.x > bMax.x || aMin.y > bMax.y || aMin.z > bMax.z);
    // Godot semantics: the pair collides when either side's layer is in the other's
    // mask (defaults 1/1 keep every pair, as before these fields existed).
    const bool layersAllow = ((bodyA->collisionLayer & bodyB->collisionMask) != 0) ||
            ((bodyB->collisionLayer & bodyA->collisionMask) != 0);
    const bool constrained = pairConstrained(adjacency, manifolds, bodyA, bodyB);
    // Two static bodies can never need a contact: neither can move, and the solver has no
    // update for either. Without this the pair still got a manifold slot in broadPhase,
    // and `ManifoldStore::frozen()` then reported it alive every step (both sides are
    // immovable), so `warmstartForces` kept a zero-contact manifold forever - one store
    // slot and eight contact records per overlapping pair of static bodies, on a terrain
    // made of static boxes for instance. A pair that both *sleep* is a different thing: its
    // manifold carries the state that has to survive until it wakes, so only mass 0 counts.
    const bool bothStatic = bodyA->mass <= 0.0f && bodyB->mass <= 0.0f;
    if (verbose)
    {
        int heldBy = -1;
        if (constrained)
            for (const ForceRef &ref : adjacency.of(bodyA->broadIndex))
                if (otherBody(ref, bodyA, manifolds) == bodyB && ref.isManifold())
                    heldBy = static_cast<int>(ref.slot());
        std::printf("      pair %d/%d bounds=%d layers=%d constrained=%d heldSlot=%d\n"
                "        aC=(%.3f,%.3f,%.3f) aH=(%.3f,%.3f,%.3f)"
                " bC=(%.3f,%.3f,%.3f) bH=(%.3f,%.3f,%.3f)\n",
                bodyA->broadIndex, bodyB->broadIndex, (int)boundsOverlap, (int)layersAllow,
                (int)constrained, heldBy,
                aCenter.x, aCenter.y, aCenter.z, aHalf.x, aHalf.y, aHalf.z,
                bCenter.x, bCenter.y, bCenter.z, bHalf.x, bHalf.y, bHalf.z);
    }
    return boundsOverlap && layersAllow && !constrained && !bothStatic;
}

// A sleeping body whose swept AABB touches a moving movable body wakes up. AABBs are
// already grown by the speculative margin, so a touch here means "these two meet within the
// tick"; no real contact is needed. Serial by construction (called from broadPhase's pair
// loop): the fields are plain bools/ints and two pairs sharing a body used to write them
// from two threads.
static void wakeSleeperPair(Rigid *a, Rigid *b) noexcept
{
    if (a->sleeping && !b->sleeping && b->mass > 0.0f)
    {
        a->sleeping = false;
        a->stillFrames = 0;
    }
    else if (b->sleeping && !a->sleeping && a->mass > 0.0f)
    {
        b->sleeping = false;
        b->stillFrames = 0;
    }
}

// Contact detection: the Morton LBVH broad phase, rebuilt or refit every step.
//
// It indexes the swept world AABBs and emits packed pair candidates. The surrounding
// pipeline - pair sort, layer/constraint filter, manifold creation - is shared.
//
// Determinism: pair collection is integer-sorted and de-duplicated, and the final
// (min-index, other) sort fixes the manifold creation order.
void Solver::broadPhase()
{
    bodiesInOrder.clear();
    bodiesInOrder.insert(bodiesInOrder.end(), next_range(bodies).begin(), next_range(bodies).end());

    const int count = static_cast<int>(bodiesInOrder.size());
    for (int i = 0; i < count; i++)
        bodiesInOrder[i]->broadIndex = i;

    // Pair existence and AABB wake-up read the constraint graph, so it has to be current
    // before the candidate filter runs. Built again after warm-starting, when the live
    // manifold set has changed.
    buildForceAdjacency();

    if (count < 2)
        return;

    // Swept bounds: the union of the body's AABB at the step start and at the step end.
    // For linear motion that union is exact - centred at p + disp/2 with half extent
    // h + |disp|/2 - and the rotation term conservatively grows it by the arc a point at
    // the bounding radius sweeps. The speculative margin is the last growth, which is
    // what lets a pair whose surfaces are about to meet enter the narrow phase.
    bodyAabbHalf.resize(static_cast<size_t>(count));
    bodyAabbCenter.resize(static_cast<size_t>(count));
    pool->forCount(count, [this](int i) {
        Rigid *body = bodiesInOrder[i];
        const float3 tight = aabbHalfExtent(*body);

        // Same prediction warmstartBody applies to the pose; the adaptive warm-start's
        // extra gravity fraction is smaller than the full g*dt^2 term used here.
        float3 disp = body->velocityLin * collisionDt;
        if (body->mass > 0)
            disp += float3{0, 0, body->gravity * collisionDt * collisionDt};

        const float rotSweep = body->radius * length(body->velocityAng) * collisionDt;
        const float grow = rotSweep + speculativeMargin;

        bodyAabbCenter[i] = body->positionLin + disp * 0.5f;
        bodyAabbHalf[i] = tight
                + float3{std::fabs(disp.x), std::fabs(disp.y), std::fabs(disp.z)} * 0.5f
                + float3{grow, grow, grow};
    }, kMinBodiesPerCheapDispatch);

    {
        const std::chrono::steady_clock::time_point tIndex = std::chrono::steady_clock::now();
        broadLbvh.buildAndCollectPairs(bodiesInOrder, sweepPairs, bodyAabbHalf,
                bodyAabbCenter, *pool);
        lastBroadIndexMs = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - tIndex).count();
    }

    // Pairs are packed (min<<32 | max), so integer order is lexicographic order over the
    // old std::pair and the sort takes the radix path instead of std::sort's comparator.
    radix_std::sort_with_buffer(sweepPairs.begin(), sweepPairs.end(), pairSortScratch);
    sweepPairs.erase(std::unique(sweepPairs.begin(), sweepPairs.end()), sweepPairs.end());

    auto kMinPairsPerDispatch = kMinForcesPerDispatch;
    const size_t pairCount = sweepPairs.size();
    broadPairValid.resize(pairCount);
    pool->forCount(pairCount, [this](int k) {
        const int i = pairFirst(sweepPairs[k]);
        const int j = pairSecond(sweepPairs[k]);
        const bool traced = traceMatch(bodiesInOrder[i], bodiesInOrder[j]);
        auto overlaps = pairOverlaps(forceAdjacency, manifolds, bodiesInOrder[i], bodyAabbHalf[i],
                                    bodyAabbCenter[i], bodiesInOrder[j], bodyAabbHalf[j],
                                    bodyAabbCenter[j], traced);
        broadPairValid[k] = overlaps;
    }, kMinPairsPerDispatch);

    // 串行创建：manifold 顺序决定后续所有阶段的遍历顺序，必须单线程
    for (size_t k = 0; k < pairCount; k++) {
        if (!broadPairValid[k]) continue;
        const int i = pairFirst(sweepPairs[k]);

        // 睡眠 + awake 的 AABB 重叠就唤醒。不需要真实接触——AABB 已经含 margin，
        // 重叠意味着"下几个 substep 内会碰"。
        //
        // 这一步必须在串行阶段：`sleeping`/`stillFrames` 是普通字段，两个候选 pair
        // 在并行 filter 里同时写同一个 body 是数据竞争（C++ 里就是 UB，TSAN 会报）。
        // 放在这里还有一个好处：顺序由排好序的 pair 表决定，唤醒结果与线程数无关。
        wakeSleeperPair(bodiesInOrder[pairFirst(sweepPairs[k])],
                bodiesInOrder[pairSecond(sweepPairs[k])]);
        const int j = pairSecond(sweepPairs[k]);
        // A dormant manifold (contact grace) let this pair back into the candidate set
        // (pairConstrained skips dormant ones so overlap can wake sleepers); it will
        // revive in its own initialize(), so creating a second manifold here would double
        // the pair's constraints.
        bool exists = false;
        for (const ForceRef &ref : forceAdjacency.of(i))
        {
            if (!ref.isManifold())
                continue;
            if (otherBody(ref, bodiesInOrder[i], manifolds) == bodiesInOrder[j])
            {
                exists = true;
                break;
            }
        }
        if (exists)
        {
            if (traceMatch(bodiesInOrder[i], bodiesInOrder[j]))
                std::printf("TRACE cand %d/%d valid=1 exists=1 (slot held)\n", i, j);
            continue;
        }
        if (traceMatch(bodiesInOrder[i], bodiesInOrder[j]))
            std::printf("TRACE cand %d/%d valid=1 acquire\n", i, j);
        (void)manifolds.acquire(bodiesInOrder[i], bodiesInOrder[j]);
    }

}

// Bring every force up to date for this step, and drop the ones that have gone inactive.
//
// Initialising is independent work - a force only reads the body poses and writes its own
// warm-started state - so it runs in parallel. Deleting is not: an inactive force (a
// broken joint) unlinks itself from the solver and from both bodies' lists, so that
// mutation happens in one thread, in list order, once the parallel pass has finished.
// Manifolds live in the store, so "deleting" one is removing its slot from the live list
// and returning it to the free list.
//
// Returns the number of non-manifold forces left, and leaves `forceOrder` indexing those.

} // namespace avbd
