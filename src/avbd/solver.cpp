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
#include <limits>
#include <algorithm>
#include <numeric>
#include <span>

#include "avbd/solver.h"

#include "avbd/job_pool.hpp"
#include "avbd/taskflow.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "avbd/bvh/bvh.hpp"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

namespace {

// Times one phase into a profile slot: construction starts the clock, `workDone` closes the item
// loop's span and `syncDone` the barrier's. A null slot - every worker but 0, and every build
// without AVBD_PROFILE - makes every call a no-op, which is what keeps the timers off the hot path.
class PhaseTimer
{
public:
    explicit PhaseTimer(Solver::PhaseProfile *p_slot) noexcept :
            slot(p_slot)
    {
#ifdef AVBD_PROFILE
        if (slot != nullptr)
            start = std::chrono::steady_clock::now();
#endif
    }

    void workDone() noexcept { stamp(true); }
    void syncDone() noexcept { stamp(false); }

private:
    void stamp(bool work) noexcept
    {
#ifdef AVBD_PROFILE
        if (slot == nullptr)
            return;
        const auto now = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(now - start).count();
        (work ? slot->work : slot->sync) += seconds;
        start = now;
#else
        (void)work;
#endif
    }

    Solver::PhaseProfile *slot;
#ifdef AVBD_PROFILE
    std::chrono::steady_clock::time_point start;
#endif
};

} // namespace

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

// Instrumentation (see Solver's "Instrumentation" block). Everything here folds away outside an
// AVBD_PROFILE / AVBD_DIAG build, which is why the phases can call it unconditionally.
Solver::PhaseProfile *Solver::phaseSlot(int id, int workerIndex) noexcept
{
#ifdef AVBD_PROFILE
    if (workerIndex == 0 && id >= 0 && id < kPhaseCount)
        return &profileData.phases[id];
#else
    (void)id;
    (void)workerIndex;
#endif
    return nullptr;
}

void Solver::StepProfile::reset() noexcept
{
    for (PhaseProfile &phase : phases)
        phase = PhaseProfile{};
    wall = 0.0;
    broad = 0.0;
    pairQuery = 0.0;
    bvhBuild = 0.0;
    pairSearch = 0.0;
    pairTests = 0.0;
    pairCandidates = 0.0;
    constrainedRebuilds = 0.0;
    movedBodies = 0.0;
    totalBodies = 0.0;
    incrementalSteps = 0.0;
    collide = 0.0;
    manifoldCreate = 0.0;
    forceInit = 0.0;
    selectStage = 0.0;
    step = 0.0;
    rounds = 0.0;
    samples = 0;
    samplesPerStep = 0;
}

#ifdef AVBD_DIAG
void Solver::Diag::reset() noexcept
{
    maxJointLambda.store(0.0f, std::memory_order_relaxed);
    maxJointC0.store(0.0f, std::memory_order_relaxed);
    chargeContact.store(0.0f, std::memory_order_relaxed);
    chargeJoint.store(0.0f, std::memory_order_relaxed);
    askedContact.store(0.0f, std::memory_order_relaxed);
    askedJoint.store(0.0f, std::memory_order_relaxed);
    maxBodyMove.store(0.0f, std::memory_order_relaxed);
}
#endif

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
            const float3 &bmin = bvhNodes.boundsMin(slot);
            const float3 &bmax = bvhNodes.boundsMax(slot);
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
            const int bodyIndex = bvhNodes.bodyIndex(slot);
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
        const float3 &bmin = bvhNodes.boundsMin(slot);
        const float3 &bmax = bvhNodes.boundsMax(slot);
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
        const int bodyIndex = bvhNodes.bodyIndex(slot);
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

// One tick, as a classified task flow. The broad phase / force warmstart / update-set
// selection run once on the calling thread, then `runSubstep` runs the per-body
// warmstart, the block-Jacobi iteration region and the velocity finish.
void Solver::step()
{
    ensurePool();

#ifdef AVBD_DIAG
    diag.reset();
#endif
#ifdef AVBD_PROFILE
    const auto stepStart = std::chrono::steady_clock::now();
#endif
    computeSleepFrames();
    broadPhase();
#ifdef AVBD_PROFILE
    const auto afterBroad = std::chrono::steady_clock::now();
    profileData.broad += std::chrono::duration<double>(afterBroad - stepStart).count();
#endif
    const int forceCount = warmstartForces();
#ifdef AVBD_PROFILE
    const auto afterForces = std::chrono::steady_clock::now();
    profileData.forceInit += std::chrono::duration<double>(afterForces - afterBroad).count();
#endif
    selectUpdateBodies();
#ifdef AVBD_PROFILE
    profileData.selectStage += std::chrono::duration<double>(std::chrono::steady_clock::now() - afterForces).count();
#endif

    int targetIterations = iterations;
    int newtonRounds = iterations;
    computeIterations(targetIterations, newtonRounds);
    runSubstep(forceCount, targetIterations, newtonRounds);
#ifdef AVBD_PROFILE
    profileData.step += std::chrono::duration<double>(std::chrono::steady_clock::now() - stepStart).count();
    profileData.samplesPerStep++;
#endif
}

// A broad-phase candidate test: bounding spheres plus Godot's layer/mask pair rule,
// minus pairs already joined by a force (IgnoreCollision, joints).
// The cheap exact filter: bounding spheres and Godot's layer/mask rule. Whether the pair is already
// constrained is a separate question, answered once per step for all candidates together (see
// Solver::constrainedPairKeys) - walking a body's force list per candidate is what made a settled pile's
// broad phase expensive.
static bool pairOverlaps(const Rigid *bodyA, const Rigid *bodyB)
{
    const float3 dp = bodyA->positionLin - bodyB->positionLin;
    const float reach = bodyA->radius + bodyB->radius;
    // Godot semantics: the pair collides when either side's layer is in the other's
    // mask (defaults 1/1 keep every pair, as before these fields existed).
    const bool layersAllow = ((bodyA->collisionLayer & bodyB->collisionMask) != 0) ||
            ((bodyB->collisionLayer & bodyA->collisionMask) != 0);
    return dp.squaredNorm() <= reach * reach && layersAllow;
}

// LSD radix sort for the broad phase's pair keys. A key is (smaller body id in the high 32 bits, larger
// in the low 32), so the bytes that can be nonzero are the low few of each half: ids are dense, which
// makes those passes a handful of linear passes over a few thousand keys. The comparison sort this
// replaces spent ~200 us a step on the constrained set alone. Stable, so the result is deterministic.
static void radixSortPairKeys(std::vector<uint64_t> &r_keys, std::vector<uint64_t> &scratch,
        uint32_t maxId) noexcept
{
    const size_t count = r_keys.size();
    if (count < 2)
        return;

    uint32_t significantBytes = 4;
    if (maxId < (1u << 8))
        significantBytes = 1;
    else if (maxId < (1u << 16))
        significantBytes = 2;
    else if (maxId < (1u << 24))
        significantBytes = 3;

    scratch.resize(count);
    std::vector<uint64_t> *source = &r_keys;
    std::vector<uint64_t> *target = &scratch;
    for (uint32_t half = 0; half < 2; ++half)
    {
        for (uint32_t byte = 0; byte < significantBytes; ++byte)
        {
            const uint32_t shift = half * 32 + byte * 8;
            size_t bucketCount[256] = {};
            for (size_t i = 0; i < count; ++i)
                ++bucketCount[(*source)[i] >> shift & 0xffu];
            size_t offset[256];
            size_t running = 0;
            for (size_t b = 0; b < 256; ++b)
            {
                offset[b] = running;
                running += bucketCount[b];
            }
            for (size_t i = 0; i < count; ++i)
                (*target)[offset[(*source)[i] >> shift & 0xffu]++] = (*source)[i];
            std::swap(source, target);
        }
    }
    if (source != &r_keys)
        r_keys.swap(scratch);
}

// Stable body ids for the broad phase's incremental pair set (see Rigid::id). Handed out on creation
// and recycled on destruction, so the id space stays dense even in a scene that streams bodies.
int32_t Solver::acquireBodyId()
{
    int32_t id;
    if (!freeBodyIds.empty())
    {
        id = freeBodyIds.back();
        freeBodyIds.pop_back();
    }
    else
    {
        id = nextBodyId++;
        bodiesById.resize(static_cast<size_t>(nextBodyId), nullptr);
        // An invalid sentinel box, so a fresh body's first comparison says "moved" whatever the epsilon
        // is (and a NaN never reads as "unmoved").
        constexpr float inf = std::numeric_limits<float>::infinity();
        broadAabbMin.resize(static_cast<size_t>(nextBodyId), float3{inf, inf, inf});
        broadAabbMax.resize(static_cast<size_t>(nextBodyId), float3{-inf, -inf, -inf});
        broadSweepMin.resize(static_cast<size_t>(nextBodyId));
        broadSweepMax.resize(static_cast<size_t>(nextBodyId));
        bodyMoved.resize(static_cast<size_t>(nextBodyId), 0);
    }
    bodiesById[static_cast<size_t>(id)] = nullptr; // the caller links it once it is constructed
    // A recycled id must not inherit the previous body's broad-phase reference, or the drift test would
    // measure the new body's box against a stranger's and could call an unchanged one "moved" - or, worse,
    // a moved one "settled" and never query it. The invalid sentinel forces the first comparison to say
    // "moved".
    constexpr float inf = std::numeric_limits<float>::infinity();
    broadAabbMin[static_cast<size_t>(id)] = float3{inf, inf, inf};
    broadAabbMax[static_cast<size_t>(id)] = float3{-inf, -inf, -inf};
    bodyMoved[static_cast<size_t>(id)] = 1;
    return id;
}

void Solver::releaseBodyId(int32_t id) noexcept
{
    if (id < 0 || static_cast<size_t>(id) >= bodiesById.size())
        return;
    bodiesById[static_cast<size_t>(id)] = nullptr;
    freeBodyIds.push_back(id);
}

// Tight world-space AABB of a body's shape at its current pose. The broad phase needs it for the tree,
// for the motion test and for the region a moved body has to re-query, and the builder takes it from
// here rather than walking the shapes again.
static inline void bodyAabb(Rigid &body, float3 &r_min, float3 &r_max) noexcept
{
    float3 half;
    switch (body.shape)
    {
        case ShapeType::Sphere:
            half = float3{body.radius, body.radius, body.radius};
            break;
        case ShapeType::Cylinder:
            half = float3{body.size.x(), body.size.z() * 0.5f, body.size.x()};
            break;
        case ShapeType::Box:
        default:
            half = body.size * 0.5f;
            break;
    }
    // The rotation is the body's cached one: the broad phase refreshes it from the pose just before
    // this call, and the narrow phase reads the same value later in the step.
    const float3x3 &rot = body.rot;
    const float3 wx{
            std::abs(rot(0, 0)) * half.x() + std::abs(rot(0, 1)) * half.y() + std::abs(rot(0, 2)) * half.z(),
            std::abs(rot(1, 0)) * half.x() + std::abs(rot(1, 1)) * half.y() + std::abs(rot(1, 2)) * half.z(),
            std::abs(rot(2, 0)) * half.x() + std::abs(rot(2, 1)) * half.y() + std::abs(rot(2, 2)) * half.z()};
    r_min = body.positionLin - wx;
    r_max = body.positionLin + wx;
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

#ifdef AVBD_PROFILE
    const auto pairStart = std::chrono::steady_clock::now();
#endif
    broadIndices.resize(static_cast<size_t>(count));
    broadBodies.resize(static_cast<size_t>(count));
    broadStepMin.resize(static_cast<size_t>(count));
    broadStepMax.resize(static_cast<size_t>(count));
    // Boxes, and which bodies moved since the last broad phase. A body's box is compared bit for bit:
    // an unchanged box means every overlap this body has is one it had last step, which is what lets the
    // previous candidate set stand for it. `broadMin/broadMax` is the box as of the last step; the union
    // of the two is the region a moved body has to re-query, so nothing it swept through is missed.
    int movedCount = 0;
    for (int i = 0; i < count; ++i)
    {
        Rigid *body = bodiesInOrder[static_cast<size_t>(i)];
        broadIndices[static_cast<size_t>(i)] = i;
        broadBodies[static_cast<size_t>(i)] = body;

        // Rotation first: the box needs it, the narrow phase reads it later in this same step, and a
        // body's pose is final from here until the warm start predicts the next one.
        body->rot = float3x3(body->positionAng);
        float3 mn, mx;
        bodyAabb(*body, mn, mx);
        broadStepMin[static_cast<size_t>(i)] = mn;
        broadStepMax[static_cast<size_t>(i)] = mx;
        const size_t id = static_cast<size_t>(body->id);
        bodiesById[id] = body;
        // Movement, measured against the box the last *query* ran from. Keeping that reference until
        // the drift exceeds the epsilon is what bounds what a missed pair can be: the box the carried-over
        // pairs were found with is never more than epsilon away from the live one, so a pair that closes
        // within the epsilon is looked at again within one epsilon of travel. Updating the reference every
        // step instead would let a slow body drift arbitrarily far without ever being re-queried.
        const float3 refMin = broadAabbMin[id];
        const float3 refMax = broadAabbMax[id];
        const float drift = std::max((mn - refMin).cwiseAbs().maxCoeff(), (mx - refMax).cwiseAbs().maxCoeff());
        const bool moved = !(drift <= params.collision.incrementalPairEpsilon);
        if (moved)
        {
            ++movedCount;
            bodyMoved[id] = 1;
            // The query covers where the body is and where its previous box was, so nothing it swept
            // through since the last look can be missed.
            broadSweepMin[id] = float3{std::min(mn.x(), refMin.x()), std::min(mn.y(), refMin.y()),
                    std::min(mn.z(), refMin.z())};
            broadSweepMax[id] = float3{std::max(mx.x(), refMax.x()), std::max(mx.y(), refMax.y()),
                    std::max(mx.z(), refMax.z())};
            broadAabbMin[id] = mn;
            broadAabbMax[id] = mx;
        }
        else
        {
            bodyMoved[id] = 0;
        }
    }

    bvh::builder::Builder builder(bvhNodes, broadIndices, broadBodies, broadStepMin, broadStepMax,
            broadCoded);
    const int leafCount = builder.buildLBVH();
#ifdef AVBD_PROFILE
    const auto afterBuild = std::chrono::steady_clock::now();
    profileData.bvhBuild += std::chrono::duration<double>(afterBuild - pairStart).count();
#endif
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
#ifdef AVBD_PROFILE
    long long pairTests = 0;
#endif
    // The nodes are packed (bounds and liveness in one record), so a test is one line per node rather
    // than three parallel-array reads, and the empty-subtree case - the right tail when the leaf count
    // is not a power of two - is answered by the same six compares instead of a fourth array.
    const std::span<const bvh::nodes::Node> nodes = bvhNodes.nodes();
    const auto collideSubtrees = [&](auto &&self, int sa, int sb) -> void {
        const bvh::nodes::Node &na = nodes[sa];
        const bvh::nodes::Node &nb = nodes[sb];
#ifdef AVBD_PROFILE
        ++pairTests;
#endif
        if (!na.valid() || !nb.valid())
            return;
        // Touching counts as overlapping, as it did in the pairwise sweep this replaced.
        if (nb.min.x() > na.max.x() || nb.min.y() > na.max.y() || nb.min.z() > na.max.z() ||
                na.min.x() > nb.max.x() || na.min.y() > nb.max.y() || na.min.z() > nb.max.z())
            return;
        const bool aLeaf = sa >= firstLeaf;
        const bool bLeaf = sb >= firstLeaf;
        if (aLeaf && bLeaf)
        {
            sweepPairs.emplace_back(std::minmax(bvhNodes.bodyIndex(sa), bvhNodes.bodyIndex(sb)));
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
    // Self-collision: test the two children of every internal node once. Each internal-node pair is
    // tested exactly once, which is what makes this the cheap shape - the alternative (walk each leaf up
    // to the root and collide it with each ancestor's sibling) re-tests the same subtree once per
    // ancestor level and measured 2x the node tests for the same pair set on a dense pile.
    const auto selfTraverse = [&](auto &&self, int slot) -> void {
        if (slot >= firstLeaf)
            return;
        const int l = 2 * slot + 1;
        const int r = 2 * slot + 2;
        collideSubtrees(collideSubtrees, l, r);
        self(self, l);
        self(self, r);
    };
    // The pair set, incremental when most of the scene stood still.
    //
    // A pair's boxes can only start overlapping when one of the pair moved, so a pair between two
    // unchanged bodies that was not a candidate last step cannot have become one - those candidates
    // simply carry over. Every pair involving a body that moved is re-derived by querying that body's
    // swept box against the tree, which covers both where it was and where it is. Above a quarter of the
    // scene moving, the per-body queries stop being the cheaper side and the whole self-collision runs
    // instead, which is also what the first step does (there is nothing to carry over yet).
    broadPairKeys.clear();
    newPairKeys.clear();
    const float movedFraction =
            count > 0 ? static_cast<float>(movedCount) / static_cast<float>(count) : 1.0f;
    const bool incremental = params.collision.incrementalPairs && !broadPairCache.empty() &&
            movedFraction <= params.collision.incrementalPairFraction;
    broadPairKeysCarried.clear();
    if (incremental)
    {
        for (uint64_t key : broadPairCache)
        {
            int32_t ia, ib;
            unpackPairIds(key, ia, ib);
            const size_t ida = static_cast<size_t>(ia);
            const size_t idb = static_cast<size_t>(ib);
            if (ida >= bodiesById.size() || idb >= bodiesById.size())
                continue;
            if (bodiesById[ida] == nullptr || bodiesById[idb] == nullptr)
                continue; // one of them has been deleted
            if (bodyMoved[ida] != 0 || bodyMoved[idb] != 0)
                continue; // re-derived below
            broadPairKeysCarried.push_back(key);
        }
#ifdef AVBD_PROFILE
        ++profileData.incrementalSteps;
#endif
        // One query per moved body: descend the tree with its swept box and emit every leaf whose box
        // overlaps. Leaf slots are the only place a body index is worth reading, so this reads one per
        // emitted pair rather than one per test.
        const auto queryBody = [&](auto &&self, int slot, const float3 &qmin, const float3 &qmax,
                                           int32_t queryId) -> void {
            const bvh::nodes::Node &node = nodes[slot];
#ifdef AVBD_PROFILE
            ++pairTests;
#endif
            if (!node.valid())
                return;
            if (qmin.x() > node.max.x() || qmin.y() > node.max.y() || qmin.z() > node.max.z() ||
                    node.min.x() > qmax.x() || node.min.y() > qmax.y() || node.min.z() > qmax.z())
                return;
            if (slot >= firstLeaf)
            {
                Rigid *other = bodiesInOrder[static_cast<size_t>(bvhNodes.bodyIndex(slot))];
                newPairKeys.push_back(packPairIds(queryId, other->id));
                return;
            }
            self(self, 2 * slot + 1, qmin, qmax, queryId);
            self(self, 2 * slot + 2, qmin, qmax, queryId);
        };
        for (Rigid *body : bodiesInOrder)
        {
            if (bodyMoved[static_cast<size_t>(body->id)] == 0)
                continue;
            queryBody(queryBody, 0, broadSweepMin[static_cast<size_t>(body->id)],
                    broadSweepMax[static_cast<size_t>(body->id)], body->id);
        }
    }
    else
    {
        selfTraverse(selfTraverse, 0);
        newPairKeys.reserve(sweepPairs.size());
        for (const auto &[i, j] : sweepPairs)
            newPairKeys.push_back(packPairIds(bodiesInOrder[static_cast<size_t>(i)]->id,
                    bodiesInOrder[static_cast<size_t>(j)]->id));
    }

#ifdef AVBD_PROFILE
    const auto afterTraverse = std::chrono::steady_clock::now();
    profileData.pairSearch += std::chrono::duration<double>(afterTraverse - afterBuild).count();
#endif

    // Only the newly derived keys need sorting: the carried-over ones are the previous step's list, which
    // was already in order and was only filtered. Sorting the whole set meant a comparison sort of several
    // thousand 64-bit keys every step (~170 us on a settled pile) where the new keys are a few hundred and
    // the two lists merge in one linear pass. The merge deduplicates as it goes, which is what `unique`
    // was for.
    radixSortPairKeys(newPairKeys, pairSortScratch, static_cast<uint32_t>(bodiesById.size()));
    broadPairKeys.swap(broadPairKeysCarried); // broadPairKeys now holds the carried list, in order
    broadPairMerged.clear();
    broadPairMerged.reserve(broadPairKeys.size() + newPairKeys.size());
    {
        size_t carried = 0;
        size_t fresh = 0;
        while (carried < broadPairKeys.size() && fresh < newPairKeys.size())
        {
            const uint64_t a = broadPairKeys[carried];
            const uint64_t b = newPairKeys[fresh];
            if (a < b)
            {
                broadPairMerged.push_back(a);
                ++carried;
            }
            else
            {
                if (a == b)
                    ++carried; // the same pair from both lists: emit once
                broadPairMerged.push_back(b);
                ++fresh;
            }
        }
        for (; carried < broadPairKeys.size(); ++carried)
            broadPairMerged.push_back(broadPairKeys[carried]);
        for (; fresh < newPairKeys.size(); ++fresh)
            broadPairMerged.push_back(newPairKeys[fresh]);
    }
    broadPairKeys.swap(broadPairMerged);
    broadPairCache = broadPairKeys;
#ifdef AVBD_PROFILE
    profileData.pairTests += static_cast<double>(pairTests);
    profileData.pairCandidates += static_cast<double>(broadPairKeys.size());
    profileData.movedBodies += static_cast<double>(movedCount);
    profileData.totalBodies += static_cast<double>(count);
    const auto afterPairs = std::chrono::steady_clock::now();
    profileData.pairQuery += std::chrono::duration<double>(afterPairs - pairStart).count();
#endif
    // Pairs that already carry a constraint (a manifold, a joint, an ignore pair), in one sorted list:
    // the candidates come sorted by the same key, so the creation loop is a merge and each candidate costs
    // a comparison instead of a walk over a body's force list.
    // Rebuilt only when a force or body has been created or destroyed since the last time. On a settled
    // scene that is nothing, so the set costs a comparison per step instead of a pass over every force
    // plus a sort of the result - which was the largest single item left in the broad phase.
    if (constrainedPairGeneration != structureGeneration)
    {
        constrainedPairKeys.clear();
        for (const Force *f = forces; f != nullptr; f = f->next)
        {
            if (f->bodyA == nullptr || f->bodyB == nullptr)
                continue;
            constrainedPairKeys.push_back(packPairIds(f->bodyA->id, f->bodyB->id));
        }
        radixSortPairKeys(constrainedPairKeys, constrainedSortScratch,
                static_cast<uint32_t>(bodiesById.size()));
        constrainedPairGeneration = structureGeneration;
#ifdef AVBD_PROFILE
        ++profileData.constrainedRebuilds;
#endif
    }

    size_t constrainedBase = 0;
    for (uint64_t key : broadPairKeys)
    {
        // Both lists are ascending, and the candidate list is strictly increasing (it was deduplicated),
        // so the cursor only moves forward.
        while (constrainedBase < constrainedPairKeys.size() && constrainedPairKeys[constrainedBase] < key)
            ++constrainedBase;
        if (constrainedBase < constrainedPairKeys.size() && constrainedPairKeys[constrainedBase] == key)
            continue; // already constrained

        int32_t ia, ib;
        unpackPairIds(key, ia, ib);
        Rigid *bodyA = bodiesById[static_cast<size_t>(ia)];
        Rigid *bodyB = bodiesById[static_cast<size_t>(ib)];
        if (bodyA == nullptr || bodyB == nullptr || bodyA == bodyB)
            continue;
        if (pairOverlaps(bodyA, bodyB))
            new Manifold(this, bodyA, bodyB);
    }
#ifdef AVBD_PROFILE
    profileData.manifoldCreate +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - afterPairs).count();
#endif
}

// A constraint is frozen when neither of its bodies can move this step. Joints and springs
// stay out of this path: their initialize() does more than refresh contacts (initial
// error, spring setup), and the core suite pins their per-scene force counts. Manifolds
// between parked bodies (static or already asleep) are the pure case - the colouring
// already skips such bodies in the primal pass - and freezing their bookkeeping is the
// sleep win: a settled scene stops paying SAT refresh + dual updates for all its contacts.
static bool constraintFrozen(const Force *f)
{
    if (f->kind != ForceKind::Contact)
        return false; // joints / springs / soft constraints keep their exact behaviour
    // A manifold that holds no contacts is not frozen but pending: a fresh one (the broad
    // phase news these every step) or one a previous step emptied. initialize() is what tells
    // the two apart, and it deletes the manifold that comes back with nothing.
    const Manifold *manifold = f->asManifold();
    if (manifold->numContacts <= 0)
        return false;
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

// Refresh a body's derived pose state (Rigid::rot / dqLin / dqAng). Every pose-writing phase calls
// this - the warmstart, the combine pass and the projection - so a constraint reading them never has
// to convert a quaternion or subtract a pose itself, and a new phase cannot forget to.
static inline void refreshDerivedPose(Rigid *body)
{
    body->rot = float3x3(body->positionAng);
    body->dqLin = body->positionLin - body->initialLin;
    body->dqAng = body->positionAng - body->initialAng;
}

// Compute one body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so it runs through the flow's Independent stage. It covers
// static bodies as well: their poses do not move, but the constraint code reads their
// cached `initialLin`/`initialAng`.
void Solver::warmstartBody(int i)
{
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
    refreshDerivedPose(body);

    // The projection's running total starts the step empty. It is summed over every pass of the
    // step rather than reset per pass, because what `finishBody` has to take back out of the BDF1
    // velocity is the whole distance the layer moved this body, not the last pass's share.
    if (body->updateSlot >= 0)
    {
        scratch.projLin[body->updateSlot] = float3{0, 0, 0};
        scratch.projAng[body->updateSlot] = float3{0, 0, 0};
    }
}

// The idle window the scene declared in seconds, in steps of the clock this tick runs on. Rounded
// rather than floored so that a window shorter than one step still means "one step", and a
// tick_rate change re-derives it on the next tick instead of sleeping on a stale count.
void Solver::computeSleepFrames()
{
    const float seconds = params.sleep.timeToSleep;
    sleepFrames = seconds > 0.0f && dt > 0.0f
            ? std::max(1, static_cast<int>(std::lround(seconds / dt)))
            : 0;
}

// Decide this step's round count once, on the caller, after warm-starting and before the
// region starts; every participant then reads the same values.
//
// `iterations` is the round budget, full stop. It used to be multiplied by the scene's mean
// penetration over a tolerance and then clamped to a second limit, which meant every scene resting
// at millimetre scale - i.e. every scene with a stack in it - ran the clamp and ignored the number
// the host asked for. A budget that silently becomes a constant is not a budget, and it also hid
// how much the round count was really costing. The estimate is gone.
void Solver::computeIterations(int &r_targetIterations, int &r_newtonRounds)
{
    const SolverParams::Primal &primal = params.primal;

    r_targetIterations = std::max(iterations, 0);
    // newtonRatio splits the set: the leading rounds run primal + dual, the tail runs
    // dual-only relaxation (optionally with penalty decay via stiffnessDecay).
    r_newtonRounds = std::clamp(
            static_cast<int>(std::ceil(primal.newtonRatio * static_cast<float>(r_targetIterations))),
            0, r_targetIterations);
}

// One step: body warmstart (independent), the block-Jacobi iteration region, then the BDF1
// velocity finish (independent). The stages are built every step because the update set and
// force count change earlier in the tick, and joined stage by stage, which is the same
// fence the phases relied on when they were called directly.
void Solver::runSubstep(int forceCount, int targetIterations, int newtonRounds)
{
#ifdef AVBD_PROFILE
    const auto profileStart = std::chrono::steady_clock::now();
#endif

    // One stage: the whole substep inside a single region.
    //
    // It used to be four dispatched stages - body warmstart, the iteration set, the closing
    // projection (two more) and the velocity finish - and at eight threads the dispatches themselves
    // costed ~340 us per substep regardless of scene size, which is 75% of a 137-body step: waking
    // the pool, joining on a latch, reserving the region's slots. Every one of those stages is
    // per-body work with no cross-item dependency, so they are phases of one region instead, fenced
    // by the same LoopSync the rounds already use.
    std::atomic<bool> converged{false};
    flow->clear();
    flow->region(
            [this, forceCount, targetIterations, newtonRounds, &converged](int workerIndex, auto &sync) {
                substepBody(workerIndex, sync, forceCount, targetIterations, newtonRounds, converged);
            },
            [this, forceCount] {
                // Small scene: the region's barrier and participant setup cost more than the phases they
                // fence, so the identical body runs inline on the calling thread.
                return static_cast<int>(updateOrder.size()) < params.primal.parallelMinBodies
                        && forceCount < kMinForcesPerDispatch;
            });
    flow->run();

#ifdef AVBD_PROFILE
    profileData.wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - profileStart).count();
    profileData.samples++;
#endif
}

// One substep, on the region's participants: the body warmstart, the iteration set, the closing
// projection and the BDF1 velocity finish. The order is the one the dispatched stages had, so the
// inline path (NullSync) is still exactly the serial reference the digest tests compare against.
template <typename Sync>
void Solver::substepBody(int workerIndex, Sync &sync, int forceCount, int targetIterations,
        int newtonRounds, std::atomic<bool> &converged)
{
    PhaseTimer warmTimer(phaseSlot(kPhaseWarmstart, workerIndex));
    sync.forItems(static_cast<int>(warmstartOrder.size()), [this](int i) { warmstartBody(i); });
    warmTimer.workDone();
    sync.arriveAndWait();
    warmTimer.syncDone();

    // The iteration set. Its convergence probe carries a residual across rounds, which is why it
    // lives here rather than in the round loop.
    float prevResidual = 0.0f;
    iterateBody(workerIndex, sync, targetIterations, forceCount, newtonRounds, prevResidual, converged);

    // The closing projection. It runs whether or not the iteration set finished, finished early on
    // convergence or never projected at all, and it is the pass that makes the promise of the layer
    // true of the pose actually handed to the host: what `finishBody` then measures the velocity
    // from is feasible.
    if (params.projection.rounds > 0)
        projectPass(sync, forceCount, workerIndex);

    PhaseTimer finishTimer(phaseSlot(kPhaseFinish, workerIndex));
    sync.forItems(static_cast<int>(warmstartOrder.size()), [this](int i) { finishBody(i); });
    finishTimer.workDone();
    sync.arriveAndWait();
    finishTimer.syncDone();
}

// One worker-loop body for the whole iteration set: `iterations` rounds of the block Jacobi primal
// (one full block solve per body against the iteration-start pose, applied immediately for a body
// whose forces are all contacts) followed by a dual pass, with a barrier fence at every phase
// boundary.
//
// The shape matters more than the phases: one dispatch arms the workers once per step and each
// boundary costs a ~1-3 us barrier. A round is two fences - one after the per-body assembly, one
// after the dual - regardless of how entangled the constraint graph is; the old colour-parallel
// Gauss-Seidel needed one per colour, which is why dense stacking scenes starved it (and why the
// incremental colouring's correctness was so fragile). Three shapes keep it at two:
//   * the assembly pass (per-contact world arms and clamped forces) runs only for the first round:
//     every later round reads what the previous round's dual emitted, and the dual measures exactly
//     the pose the next primal solves against;
//   * a body with no joint on it applies its delta in the same item that solved it, because a joint
//     is the only force that reads its partner's pose - and a joint makes both of its bodies
//     non-contactOnly, so those deltas are the only ones that need a fence between solve and apply;
//   * the combine pass therefore exists only for the bodies a joint or spring couples, and not at
//     all in a scene that has none.
//
// Both sync types run the identical body: NullSync (small scenes, or another solver
// holding the region) executes it inline in list order, which is exactly the threads=1
// reference the digest tests compare against.
template <typename Sync>
void Solver::iterateBody(int workerIndex, Sync &sync, int targetIterations, int forceCount,
        int newtonRounds, float &prevResidual, std::atomic<bool> &converged)
{
    // Per-worker Chebyshev schedule state: the recurrence is deterministic, so every
    // worker derives the same weights without sharing anything.
    int chebRound = 0;
    float chebGammaPrev = 1.0f;
    int roundsRun = 0;
    // Read once: every participant must agree on which rounds project, or the barriers inside
    // the pass stop being paired.
    const int projectionRounds = params.projection.rounds;
    const bool projectNewtonOnly = params.projection.trailingNewtonOnly;
    // Whether the next round has to assemble its manifolds from scratch: yes for the first round
    // (nothing has been measured yet) and after a projection pass, which moves poses in between the
    // dual that filled the cache and the stamps that read it.
    bool needAssemble = true;

    for (int it = 0; it < targetIterations; it++)
    {
        roundsRun++;
        const bool newton = it < newtonRounds;
        if (newton)
        {
            if (needAssemble)
            {
                // Block projected Jacobi, assembly pass: each manifold evaluates the shared
                // per-contact state (world arms, clamped forces) once from the round's pose
                // snapshot, work the two body stamps used to duplicate. Only the first round of a
                // step needs it - every later round reads what the previous round's dual emitted,
                // because the dual measures exactly the pose the next primal solves against. That
                // is one phase and one barrier per round that the old layout spent here.
                PhaseTimer preTimer(phaseSlot(kPhasePrecompute, workerIndex));
                sync.forItems(forceCount, [this](int i) {
                    Force *f = forceOrder[i];
                    if (f->kind != ForceKind::Contact)
                        return;
                    Manifold *manifold = f->asManifold();
                    if (manifold->bodyA->updateSlot < 0 && manifold->bodyB->updateSlot < 0)
                        return; // frozen: neither side moves this round
                    manifold->preCompute(params.dual.alpha);
                });
                preTimer.workDone();
                sync.arriveAndWait();
                preTimer.syncDone();
                needAssemble = false;
            }

            // The combine schedule, read once per round: loop-invariant inside both phases below.
            const float relax = params.primal.jacobiRelaxation;
            const float jointRelax = params.primal.jointRelaxation;
            const float rho = params.primal.chebyshevRho;
            float gamma = relax;
            if (rho > 0.0f)
            {
                const float rho2 = rho * rho;
                if (chebRound == 0)
                    gamma = 1.0f;
                else if (chebRound == 1)
                    gamma = 1.0f / (1.0f - 0.5f * rho2);
                else
                    gamma = 1.0f / (1.0f - 0.25f * rho2 * chebGammaPrev);
                chebGammaPrev = gamma;
            }
            chebRound++;

            // Phase B: each body assembles its full block system (own inertia plus every incident
            // constraint - manifolds stamped from the cache, other forces stamping directly), solves
            // it, and - when everything on it is a contact - applies the delta right here. A joint is
            // the only force that reads its partner's pose, and a joint makes both of its bodies
            // non-contactOnly, so a contact-only body's move is invisible to every other item in this
            // phase: the result still does not depend on how the items were split.
            PhaseTimer stampTimer(phaseSlot(kPhaseStamp, workerIndex));
            const auto stampBody = [this, gamma, rho, relax, jointRelax](int i) {
                        Rigid *body = updateOrder[i];
                        Block block(body->mass, body->moment, dt * dt, body->positionLin,
                                body->inertialLin, body->positionAng, body->inertialAng);
                        for (Force *f = body->forces; f != 0;
                                f = (f->bodyA == body) ? f->nextA : f->nextB)
                        {
                            if (f->kind == ForceKind::Contact)
                            {
                                Manifold *manifold = f->asManifold();
                                // Safety net for a manifold the previous round's dual skipped
                                // (frozen then, mobile now). A manifold the dual skipped is one no
                                // body stamps, so this cannot be reached from two threads at once.
                                if (!manifold->primalReady)
                                    manifold->preCompute(params.dual.alpha);
                                manifold->stampCached(body, block);
                            }
                            else
                            {
                                f->updatePrimal(body, params.dual.alpha, block);
                            }
                        }
                        float3 dLin;
                        float3 dAng;
                        block.solveDelta(dLin, dAng);
                        if (body->contactOnly)
                        {
                            // Applied on the spot: nothing else in this phase reads this body's pose,
                            // and the delta never has to reach memory.
                            applyDelta(i, body, dLin, dAng, gamma, rho, relax, jointRelax);
                        }
                        else
                        {
                            scratch.deltaLin[i] = dLin;
                            scratch.deltaAng[i] = dAng;
                        }
            };
            // Split across the participants only when there is enough work to pay for the fence; below
            // the threshold one worker walks the whole update set. The barrier is outside either way, so
            // every participant still arrives at the same points.
            const int bodyCount = static_cast<int>(updateOrder.size());
            if (bodyCount >= params.primal.parallelStampMinBodies)
            {
                sync.forItems(bodyCount, stampBody);
            }
            else if (workerIndex == 0)
            {
                for (int i = 0; i < bodyCount; ++i)
                    stampBody(i);
            }
            stampTimer.workDone();
            sync.arriveAndWait();
            stampTimer.syncDone();

            // The combine pass: only the bodies a joint or spring couples. Their stamps read each
            // other's poses, so their deltas have to land together after everything is assembled.
            // With no such body in the scene the pass does not exist, and a round is two barriers.
            if (hasJointedBodies)
            {
                PhaseTimer combineTimer(phaseSlot(kPhaseCombine, workerIndex));
                sync.forItems(static_cast<int>(updateOrder.size()),
                        [this, gamma, rho, relax, jointRelax](int i) {
                            Rigid *body = updateOrder[i];
                            if (body->contactOnly)
                                return; // applied in the assembly phase
                            applyDelta(i, body, scratch.deltaLin[i], scratch.deltaAng[i], gamma, rho,
                                    relax, jointRelax);
                        });
                combineTimer.workDone();
                sync.arriveAndWait();
                combineTimer.syncDone();
            }
        }

        // Dual update: one item per force, and each force only reads the (frozen) body poses and
        // writes its own dual state, so the whole pass is independent. Frozen constraints have
        // nothing to update: every input they read is parked. The relaxation tail's penalty decay
        // is folded into the same item - it touches the state the dual is about to read and nothing
        // else, so it does not need a phase (and a barrier) of its own.
        const bool decay = !newton && params.primal.stiffnessDecay < 1.0f;
        const float decayFactor = params.primal.stiffnessDecay;
        PhaseTimer dualTimer(phaseSlot(kPhaseDual, workerIndex));
        sync.forItems(forceCount, [this, newton, decay, decayFactor](int i) {
            Force *f = forceOrder[i];
            if (constraintFrozen(f))
                return;
            // A dual-only round relaxes a multiplier against a pose nothing is going to change. For a
            // contact row that is dual ascent on a clamped convex constraint, which is what the tail is
            // for. For a bilateral row it is pure windup: C cannot move, so lambda grows by K*C every
            // round and the next step's primal inherits a force no displacement ever paid for. Joints
            // therefore keep the multiplier the last Newton round gave them.
            if (!newton && f->kind != ForceKind::Contact)
                return;
            if (decay)
            {
                Manifold *manifold = f->asManifold();
                for (int ci = 0, n = manifold->numContacts; ci < n; ++ci)
                    manifold->contacts[ci].penalty =
                            (manifold->contacts[ci].penalty * decayFactor).cwiseMin(float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});
            }
            f->updateDual(params.dual.alpha);
        });
        dualTimer.workDone();
        sync.arriveAndWait();
        dualTimer.syncDone();

        // The trailing rounds of the set project. Leading rounds leave it to the Lagrangian,
        // which reaches the same target faster where it can - a projection has no force behind
        // it, so a stack held up by projections would separate the moment one contact breaks.
        // The last rounds are where the projection earns its keep: those are the ones whose
        // residual the set runs out of budget to pay.
        //
        // `trailingNewtonOnly` moves that window onto the last Newton rounds instead, which keeps
        // the layer out of the dual-only tail where stiffnessDecay is softening the very penalty
        // it measures against. Either way the window is a suffix of a fixed range, so every
        // participant agrees on it and the barriers stay paired.
        if (projectionRounds > 0)
        {
            const int lastProjected = projectNewtonOnly ? newtonRounds : targetIterations;
            if (it < lastProjected && it >= lastProjected - projectionRounds)
            {
                projectPass(sync, forceCount, workerIndex);
                // The pass moved poses after the dual filled this round's cache, so the next round
                // has to measure again instead of stamping what is now last round's news.
                needAssemble = true;
            }
        }

        // Convergence probe after a Newton round: a full primal pass is expensive, so
        // skip further rounds when the residual stops shrinking. Only the calling
        // thread (worker 0) evaluates the reduction, and the verdict is stored in a
        // shared atomic that EVERY participant reads BEFORE the barrier - each worker
        // that sees `converged` skips its own remaining work but still participates
        // in every barrier (draining instead of returning), so the barrier
        // participation count stays equal and no worker deadlocks.
        if (workerIndex == 0 && params.primal.convergenceThreshold > 0.0f && it + 1 < targetIterations)
        {
            float total = 0.0f;
            int count = 0;
            for (Force *f = forces; f != nullptr; f = f->next)
                f->accumulatePenetration(total, count);
            const float residual = count > 0 ? total / count : 0.0f;
            converged.store(prevResidual - residual < params.primal.convergenceThreshold,
                    std::memory_order_relaxed);
            prevResidual = residual;
        }
        const bool stop = converged.load(std::memory_order_relaxed);
        if (stop)
            break; // every participant breaks at the same iteration - barriers stay paired
    }

#ifdef AVBD_PROFILE
    // Diagnostics only, and only on the thread that owns the profile block.
    if (workerIndex == 0)
        profileData.rounds += roundsRun;
#else
    (void)roundsRun;
#endif
}

// Apply one body's solved delta with the combine pass's schedule. The relaxation constant depends on
// what the body's rows are (see SolverParams::Primal::jointRelaxation), and the derived pose state
// every constraint reads is refreshed here, which is the only place a pose change can be observed.
void Solver::applyDelta(int i, Rigid *body, float3 dLin, float3 dAng, float gamma, float rho,
        float relax, float jointRelax) noexcept
{
    if (rho <= 0.0f || !body->contactOnly)
    {
        const float step = body->contactOnly ? relax : jointRelax;
        body->positionLin = body->positionLin + dLin * step;
        body->positionAng = body->positionAng + dAng * step;
    }
    else
    {
        // Chebyshev semi-iteration: x_{k+1} = gamma (x_k + delta) + (1 - gamma) x_{k-1}. The angular
        // half gets the same gamma but no prev-mixing: positions are a vector space, rotations are
        // not, and the first-order quat integration stays exact as written.
        const float3 xLin = body->positionLin;
        body->positionLin = gamma * (xLin + dLin) + (1.0f - gamma) * scratch.prevLin[i];
        scratch.prevLin[i] = xLin;
        body->positionAng = body->positionAng + dAng * gamma;
    }
    refreshDerivedPose(body);
}

// The position projection's measure half: one constraint, reading the poses every other pass has
// left and writing only its own step-scoped state. Frozen constraints are skipped for the same
// reason the dual pass skips them - nothing they read can have moved.
void Solver::projectForce(int i)
{
    Force *f = forceOrder[i];
    if (constraintFrozen(f))
        return;
    f->projectMeasure(params.projection.scale, params.projection.maxCorrection);
}

// The projection's apply half: one body, summing the share of every constraint it is in and
// moving itself. Exactly the shape of the primal's phase B, which is what makes it a Jacobi step:
// no body is written by two items, so the result does not depend on how the work was split.
void Solver::projectBody(int i)
{
    Rigid *body = updateOrder[i];
    float3 deltaLin = float3{0, 0, 0};
    float3 deltaAng = float3{0, 0, 0};
    for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
        f->projectShare(body, deltaLin, deltaAng);

    // One ceiling per body per pass, on the sum: a body pinned by eight deep contacts is pulled
    // out of the worst of them, not teleported out of all of them at once. The angular half is
    // capped by the same number read as radians - the rotation is applied through the same
    // first-order integration the primal's delta uses, so its magnitude is an angle.
    const float ceiling = params.projection.maxCorrection;
    const float moved = deltaLin.norm();
#ifdef AVBD_DIAG
    // The magnitude *requested*, before the per-body ceiling below. Paired with the applied move it
    // is what makes the old accounting bug visible: a charge computed for the requested move while
    // the ceiling delivered a fraction of it.
    atomicMaxFloat(diag.maxBodyMove, moved);
#endif
    if (moved > ceiling && moved > 0.0f)
        deltaLin = deltaLin * (ceiling / moved);
    const float turned = deltaAng.norm();
    if (turned > ceiling && turned > 0.0f)
        deltaAng = deltaAng * (ceiling / turned);
    if (moved == 0.0f && turned == 0.0f)
        return;

    body->positionLin = body->positionLin + deltaLin;
    body->positionAng = body->positionAng + deltaAng;
    refreshDerivedPose(body);
    scratch.projLin[i] = scratch.projLin[i] + deltaLin;
    scratch.projAng[i] = scratch.projAng[i] + deltaAng;
}

// One projection: measure over the forces, fence, share and apply over the update set. The fence
// is the whole design - without it the share pass would read a constraint whose error was
// computed against a pose another thread had already corrected, and the layer would inherit the
// colour-parallel Gauss-Seidel's thread-count dependence, which is why that was removed.
//
// Both halves run through the same `Sync` the iteration rounds use, so a small scene or a solver
// that does not hold the region executes the identical body inline, in list order - the threads=1
// reference the digest tests compare against.
template <typename Sync>
void Solver::projectPass(Sync &sync, int forceCount, int workerIndex)
{
    PhaseTimer timer(phaseSlot(kPhaseProject, workerIndex));
    sync.forItems(forceCount, [this](int i) { projectForce(i); });
    timer.workDone();
    sync.arriveAndWait();
    timer.syncDone();
    sync.forItems(static_cast<int>(updateOrder.size()), [this](int i) { projectBody(i); });
    timer.workDone();
    sync.arriveAndWait();
    timer.syncDone();

    // Re-anchor the stabilisation reference onto the pose this pass delivered, when the layer is
    // configured to work that way. The two halves above move the bodies; this is what stops the
    // next round's Lagrangian from pulling them straight back to the pose the step started from.
    // The closing pass of a step needs no such stage: nothing reads a reference after it, and the
    // next step re-measures every reference at the pose it inherits.
    if (params.projection.rebase)
    {
        sync.forItems(forceCount, [this](int i) { projectRebase(i); });
        timer.workDone();
        sync.arriveAndWait();
        timer.syncDone();
    }
}

// One constraint's re-anchoring: it re-measures its own error at the pose the pass left behind and
// adopts that as the pose the Lagrangian will not fight. Nothing else is touched - the multiplier
// keeps carrying whatever load it was carrying, and the penalty keeps its value - because the pose
// is the state the projection changed, so the pose is the state the reference has to follow.
void Solver::projectRebase(int i)
{
    Force *f = forceOrder[i];
    if (constraintFrozen(f))
        return;
    f->rebaseStabilisation();
}

// The rest window's follow-up rate, and the credit a body in the hysteresis band loses per step.
// One tenth per step is a window of roughly ten steps, which is what makes the travel gate read as
// "where has this body been over the last sixth of a second" at any tick rate.
constexpr float kRestWindowDecay = 0.1f;
constexpr int kHysteresisLeak = 2;

// BDF1 velocities after the final iteration. Per-body and independent; it covers static
// bodies too, so their bookkeeping matches, though only movable ones get a new velocity.
void Solver::finishBody(int i)
{
    Rigid *body = warmstartOrder[i];
    body->prevVelocityLin = body->velocityLin;
    // Don't update velocity if sleeping (Godot 4.7 API)
    if (body->mass > 0 && !body->sleeping)
    {
        body->velocityLin = (body->positionLin - body->initialLin) / dt;
        body->velocityAng = (body->positionAng - body->initialAng) / dt;
        // BDF1 cannot tell a projection's move from a body's own: both are a change of pose
        // inside the step. projection_velocity_transfer says which of them the scene keeps. At 1
        // (the default) the correction is indistinguishable from motion, which is standard PBD; at
        // 0 the layer is purely positional - the pose is legal and the velocity is what it would
        // have been without it - and that mismatch is what measured badly on the jointed scenes.
        if (params.projection.rounds > 0 && params.projection.velocityTransfer < 1.0f
                && body->updateSlot >= 0)
        {
            const float giveBack = (1.0f - params.projection.velocityTransfer) / dt;
            body->velocityLin = body->velocityLin - scratch.projLin[body->updateSlot] * giveBack;
            body->velocityAng = body->velocityAng - scratch.projAng[body->updateSlot] * giveBack;
        }
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
    // Idle detection - a travel-only rest test, a loud band around it, and a leaky counter; the
    // counter is the point.
    //
    // The quantitative fact the whole of C was written against: a BDF1 velocity is a displacement
    // divided by the tick, so the same 0.4 mm of solver jitter reads 0.024 m/s at 60 Hz and
    // 0.096 m/s at 240 Hz - Godot's 0.1 m/s threshold, reached by pure noise. Judged on velocity
    // alone, and with a budget counted in steps (the old 30, which is 0.5 s at 60 Hz but 0.125 s at
    // 240), a faster clock made sleeping both harder to earn and shorter-lived. So the window is
    // seconds (computeSleepFrames), and the travel gate below is in metres: it asks "has this body
    // moved more than 5 mm over the last few ticks", which is the same question at any tick rate.
    //
    // Measured only for sleep_mode == SLEEP bodies. This used to run for every movable body, which
    // marked a can_sleep == false one `sleeping` as well - and the velocity update above is gated
    // on `!sleeping`, while the update-set gate keeps solving for it. The result was a body that
    // went on moving but stopped accelerating, forever, since nothing but the island pass clears
    // that flag for a body that never sleeps. That pass is the only writer of `sleeping` now.
    if (sleepFrames > 0 && body->mass > 0 && body->sleep_mode == 1 && !body->sleeping)
    {
        // The reference pose chases the body at one tenth per step, which makes the difference
        // below "travelled recently" rather than "not where it was parked": a box resting 5 mm
        // inside its collision margin is not travelling, because the reference sits there with it.
        // The angular half uses the solver's own rotation-delta convention (maths.h): subtract is
        // 2 sin(theta/2) about the axis, add integrates a rotation vector and renormalises.
        body->sleepRefLin = body->sleepRefLin + (body->positionLin - body->sleepRefLin) * kRestWindowDecay;
        body->sleepRefAng = body->sleepRefAng + (body->positionAng - body->sleepRefAng) * kRestWindowDecay;

        const SolverParams::Sleep &sp = params.sleep;
        const float lin2 = body->velocityLin.squaredNorm();
        const float ang2 = body->velocityAng.squaredNorm();
        const float moved2 = (body->positionLin - body->sleepRefLin).squaredNorm();
        const float turned2 = (body->positionAng - body->sleepRefAng).squaredNorm();
        const float ratio2 = sp.hysteresisRatio * sp.hysteresisRatio;
        // Rest is the two travel gates and nothing else. Deliberately not the velocity gates: the
        // jitter a settled stack produces is a fixed distance per STEP, not per second, so dividing
        // it by a shorter tick inflates it into a real speed. Measured on this scene - 60, 240 and
        // 480 Hz with the velocity gates in the quiet test - the stack fell asleep in 0.82 s at
        // 60 Hz and 2.2 s at both faster clocks, because at 480 Hz the same 0.4 mm of solver noise
        // is a 0.19 m/s body. The window says what the tick says nothing about: where has this body
        // actually been.
        const bool quiet = moved2 < sp.positionTolerance * sp.positionTolerance &&
                turned2 < sp.rotationTolerance * sp.rotationTolerance;
        // Out of the hysteresis band on any measure - rate or travel - and the credit is gone: this
        // is a body that is genuinely moving. Godot's thresholds still own the loudest half of the
        // decision, and the whole of the wake test, which is where a fast neighbour is the question
        // being asked anyway.
        const bool loud = lin2 > sp.linearThreshold * sp.linearThreshold * ratio2 ||
                ang2 > sp.angularThreshold * sp.angularThreshold * ratio2 ||
                moved2 > sp.positionTolerance * sp.positionTolerance * ratio2 ||
                turned2 > sp.rotationTolerance * sp.rotationTolerance * ratio2;

        if (quiet)
        {
            // Capped at the budget: the credit is a progress bar to the island pass's threshold, and
            // an island that will not commit should not be counting either.
            if (body->stillFrames < sleepFrames)
                ++body->stillFrames;
        }
        else if (loud)
            body->stillFrames = 0;
        else
            // In the band: leak. A stack whose top box twitches at one and a half times the
            // tolerance loses two steps of credit per step rather than its whole budget, so it
            // settles instead of oscillating between asleep and freshly awake.
            body->stillFrames = std::max(0, body->stillFrames - kHysteresisLeak);
    }
}

// Island-wide sleep bookkeeping, and the only thing that changes `sleeping` on the solver's own
// authority. Movable, sleep-enabled bodies are grouped into islands through force connectivity
// (union-find over the flattened body list; static bodies and sleep_mode == NEVER bodies join
// nothing and never propagate), and the island is then the unit of BOTH transitions.
//
// Both halves are in one pass because either one alone thrashes - that is the audit's finding.
// Sleeping boxes individually left a contact chain part static and part solved: an asleep body is
// out of the update set, so every constraint against it resolves as if it were a wall, and the
// neighbours above it are then solved against ground that was not there a step ago. Waking boxes
// individually, on one step of a neighbour's velocity, restarts the idle credit of every body the
// woken one touches - and at 240 Hz the noise floor of a settled stack is at Godot's velocity
// threshold, so that happened constantly: the transition itself was the disturbance.
//
// Deterministic: the union-find walks the force list in list order, and both aggregate flags fold
// over the whole island, so the outcome does not depend on which member is examined first.
void Solver::updateSleepIslands()
{
    // Only bodies that can sleep join an island. With the default timeToSleep = 0 this whole pass
    // is skipped by the caller, so the core suite (which never sets sleep_mode) is unaffected.
    islandBodies.clear();
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        // Every body is tagged here, so a non-member can never carry a stale slot from a
        // previous step: bodies enter and leave this set as sleep_mode and mass change.
        const bool member = body->mass > 0 && body->sleep_mode == 1;
        body->islandSlot = member ? static_cast<int32_t>(islandBodies.size()) : -1;
        if (member)
            islandBodies.push_back(body);
    }

    // A single box on the ground is an island of one, and it sleeps on its own budget: this pass is
    // what puts it there now that finishBody no longer does.
    const int count = static_cast<int>(islandBodies.size());
    if (count == 0)
        return;

    std::vector<int> parent(count);
    for (int i = 0; i < count; ++i)
        parent[i] = i;
    auto find = [&parent](int x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };

    // Each force is a union candidate. Both endpoints are read straight from the body, so
    // the walk is O(F) with no lookup structure; a force with a null endpoint (an anchor)
    // connects nothing. Unions run in force-list order, exactly as before.
    for (Force *f = forces; f != 0; f = f->next)
    {
        if (f->bodyA == nullptr || f->bodyB == nullptr)
            continue;
        const int ia = f->bodyA->islandSlot;
        const int ib = f->bodyB->islandSlot;
        if (ia >= 0 && ib >= 0)
            parent[find(ia)] = find(ib);
    }

    const SolverParams::Sleep &sp = params.sleep;
    const float ang2 = sp.angularThreshold * sp.angularThreshold;
    const float ratio2 = sp.hysteresisRatio * sp.hysteresisRatio;
    // "The neighbour actually moved", in the position domain: a body whose BDF1 velocity says
    // "moving" but whose pose did not change last step - a locked axis, or a projection whose
    // correction finishBody took back out of the velocity - is not a disturbance. Half the distance
    // one step of the sleep threshold implies, so that it means the same at any tick rate.
    const float movedMin = sp.linearThreshold * dt * 0.5f;
    const float movedMin2 = movedMin * movedMin;

    // Wake evidence, one step of it per sleeping body. The test is the RELATIVE motion along the
    // line joining the two, not the neighbour's bare speed: what matters to a sleeping box is
    // whether the gap to its neighbour is changing. A neighbour sliding tangentially past it does
    // not wake it, and the friction it applies is then solved against an immovable side - the same
    // approximation putting the island to sleep already makes. Rotation has no relative form to
    // take, so it stays a bare threshold.
    std::vector<char> islandWake(count, 0);
    for (int i = 0; i < count; ++i)
    {
        Rigid *body = islandBodies[i];
        if (!body->sleeping)
            continue;
        bool disturbance = false;
        for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
        {
            Rigid *other = (f->bodyA == body) ? f->bodyB : f->bodyA;
            if (other == nullptr || other->mass <= 0 || other->sleeping)
                continue;
            if ((other->positionLin - other->initialLin).squaredNorm() <= movedMin2)
                continue;
            const float3 d = other->positionLin - body->positionLin;
            const float dd = d.squaredNorm();
            if (dd < 1e-12f)
                continue; // coincident centres: there is no axis to measure along
            // The sleeper's own velocity is zero while it is asleep, so this is the neighbour's for
            // now; the subtraction is what keeps it honest if that stops being true.
            const float along = std::fabs((other->velocityLin - body->velocityLin).dot(d)) / std::sqrt(dd);
            const bool turning = other->velocityAng.squaredNorm() > ang2;
            if (along > sp.linearThreshold || turning)
            {
                disturbance = true;
                // Out of the hysteresis band is not noise, and an impact must not wait out the
                // debounce: the island would spend those steps solved against a wall.
                if (along > sp.linearThreshold * sp.hysteresisRatio ||
                        other->velocityAng.squaredNorm() > ang2 * ratio2)
                {
                    islandWake[find(i)] = 1;
                    break;
                }
            }
        }
        if (!disturbance)
            body->wakeFrames = 0;
        else if (++body->wakeFrames >= sp.wakeDebounce)
        {
            islandWake[find(i)] = 1;
            body->wakeFrames = 0;
        }
    }

    // A disturbed island wakes as a unit, and every member's rest window starts over: the credit
    // they accumulated was earned against a neighbourhood that has just moved.
    for (int i = 0; i < count; ++i)
    {
        if (islandWake[find(i)])
            islandBodies[i]->wake();
    }

    // Then the sleep half, measured after the wake-up so an island that has just been disturbed
    // cannot put itself back to sleep in the same step. An island sleeps once none of its members
    // is both awake and short of the budget.
    std::vector<char> islandWaiting(count, 0);
    for (int i = 0; i < count; ++i)
    {
        const Rigid *body = islandBodies[i];
        if (!body->sleeping && body->stillFrames < sleepFrames)
            islandWaiting[find(i)] = 1;
    }

    for (int i = 0; i < count; ++i)
    {
        Rigid *body = islandBodies[i];
        if (body->sleeping || islandWaiting[find(i)])
            continue;
        body->sleeping = true;
        body->velocityLin = float3{0, 0, 0};
        body->velocityAng = float3{0, 0, 0};
        body->wakeFrames = 0;
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
// Sleep state is settled before the set is built (Solver::updateSleepIslands): a sleeping body
// only re-enters it once its island's wake evidence says so, and an island that has earned its
// budget leaves it together.
void Solver::selectUpdateBodies()
{
    warmstartOrder.clear();
    updateOrder.clear();

    // Both sleep transitions are island-wide and both happen here, before the update set is chosen:
    // a body this pass wakes moves this step, and a body it puts to sleep does not. The evidence it
    // reads - the velocities and the displacements behind them - is all from the step that just
    // ran, which is why it has to be measured here rather than in finishBody, where the per-body
    // dispatch is parallel and a neighbour's numbers are whatever its worker left.
    // Serial and read-only apart from the examined bodies, so the outcome is fixed regardless of
    // thread count.
    if (sleepFrames > 0)
        updateSleepIslands();

    // Every body is warmed up; only movable ones are updated. Static bodies are excluded
    // from the update set: they are never moved, so nothing solves for them.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        warmstartOrder.push_back(body);

        // Fresh per-step state; admitted bodies get their update-set slot below.
        body->updateSlot = -1;
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
                // Idle detection (finishBody) is the only way into sleep here;
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
                body->updateSlot = static_cast<int32_t>(updateOrder.size());
                updateOrder.push_back(body);
            }
        }
        // Static bodies keep the -1 written above: they are never in the update set.
    }

    // The only sizing point for the primal's per-body columns: one call keeps every column
    // the same length as updateOrder, which is what Rigid::updateSlot indexes into.
    scratch.resize(updateOrder.size());

    // Chebyshev eligibility: a body carrying any non-manifold force (joint, spring, world
    // anchor) extrapolates unsafely - its gamma falls back to the constant relaxation in
    // the combine pass. Walk the force lists once, after warmstartForces() has pruned.
    for (Force *f = forces; f != 0; f = f->next)
    {
        if (f->kind == ForceKind::Contact)
            continue;
        if (f->bodyA != 0)
            f->bodyA->contactOnly = false;
        if (f->bodyB != 0)
            f->bodyB->contactOnly = false;
    }

    // Whether the combine pass has to exist this step. `contactOnly` is exactly the predicate the
    // assembly phase uses to decide if a body may apply its own delta, so both sides of that decision
    // read the same flag - a body it defers is a body the combine pass applies, always.
    hasJointedBodies = false;
    for (const Rigid *body : updateOrder)
        hasJointedBodies = hasJointedBodies || !body->contactOnly;

    // Load balance. The assembly phase hands bodies out in this list's order, and a body's cost is
    // dominated by how many constraints it stamps, so sorting dearest-first puts the long items at the
    // head where the dynamic cursor overlaps them instead of leaving one heavy body as a tail no
    // worker can help with. Items are independent, so the order cannot change the result - only which
    // participant takes which body - and the sort is on a stable key, so the order itself is fixed.
    if (updateOrder.size() > 1)
    {
        const size_t count = updateOrder.size();
        updateCost.resize(count);
        for (size_t i = 0; i < count; ++i)
        {
            Rigid *body = updateOrder[i];
            int incident = 0;
            for (Force *f = body->forces; f != 0; f = (f->bodyA == body) ? f->nextA : f->nextB)
                ++incident;
            updateCost[i] = incident;
        }
        updateOrderIndex.resize(count);
        for (size_t i = 0; i < count; ++i)
            updateOrderIndex[i] = static_cast<int>(i);
        std::stable_sort(updateOrderIndex.begin(), updateOrderIndex.end(), [this](int a, int b) {
            return updateCost[static_cast<size_t>(a)] > updateCost[static_cast<size_t>(b)];
        });
        updateOrderSorted.clear();
        for (int index : updateOrderIndex)
            updateOrderSorted.push_back(updateOrder[static_cast<size_t>(index)]);
        updateOrder.swap(updateOrderSorted);
        // Renumber the slots so Rigid::updateSlot keeps addressing the column the body now occupies;
        // every scratch array is indexed by it.
        for (size_t i = 0; i < count; ++i)
            updateOrder[i]->updateSlot = static_cast<int32_t>(i);
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
    flow = std::make_unique<detail::TaskFlow>(*pool);
}

int Solver::threadCount() const
{
    if (pool == nullptr)
        return 1;
    return static_cast<int>(pool->threadCount());
}

} // namespace avbd
