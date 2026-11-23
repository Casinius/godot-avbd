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
    // Two static (mass 0) bodies need no contact: without this they hold a frozen zero-contact
    // manifold forever. Sleeping pairs differ - their manifold must survive to wake, so only mass 0.
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

// Swept AABBs already include the speculative margin, so a touch wakes a sleeper within the tick
// without real contact. Serial by construction: parallel pairs sharing a body raced on bools.
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

// Morton LBVH over swept AABBs, rebuilt/refit each step; pair sort, filter and manifold creation
// are shared. Determinism: integer sort + dedup and (min-index, other) order fix creation order.
void Solver::broadPhase()
{
    bodiesInOrder.clear();
    bodiesInOrder.insert(bodiesInOrder.end(), next_range(bodies).begin(), next_range(bodies).end());

    const int count = static_cast<int>(bodiesInOrder.size());
    for (int i = 0; i < count; i++)
        bodiesInOrder[i]->broadIndex = i;

    // Pair filter and AABB wake-up read the constraint graph; rebuild it again after warm-starting
    // changes the live manifold set.
    buildForceAdjacency();

    if (count < 2)
        return;

    // Swept bounds = start/end AABB union: exact for linear motion (centre p + disp/2, extent
    // h + |disp|/2), grown for rotation then speculativeMargin so near-misses reach narrow phase.
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
        lastBroadTreeMs = broadLbvh.lastTreeMs;
        lastBroadEmitMs = broadLbvh.lastEmitMs;
    }

    // Pairs are packed (min<<32 | max), so integer order is lexicographic order over the
    // old std::pair and the sort takes the radix path instead of std::sort's comparator.
    radix_std::sort_with_buffer(sweepPairs.begin(), sweepPairs.end(), pairSortScratch);
    sweepPairs.erase(std::unique(sweepPairs.begin(), sweepPairs.end()), sweepPairs.end());

    auto kMinPairsPerDispatch = kMinBroadPairsPerDispatch;
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

        // 睡眠体与 awake 可动体的 swept AABB（已含 margin）重叠即唤醒，无需真实接触。
        // 必须在此串行执行：并行写 sleeping/stillFrames 是 UB；顺序由排序后的 pair 表决定，与线程数无关。
        wakeSleeperPair(bodiesInOrder[pairFirst(sweepPairs[k])],
                bodiesInOrder[pairSecond(sweepPairs[k])]);
        const int j = pairSecond(sweepPairs[k]);
        // A dormant contact-grace manifold is in the candidate set so overlap can wake sleepers;
        // it revives in its own initialize() - creating another here would double constraints.
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

} // namespace avbd
