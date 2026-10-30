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
#include "avbd/taskflow.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "avbd/spatial_hash.hpp"
#include "radixsort.hpp"
namespace avbd {

Solver::Solver()
{
}

Solver::~Solver()
{
    clear();
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

    // Exact per-body test against the body's own shape. The broad phase is a spatial hash
    // rebuilt per step and not efficient to walk per ray, so picking scans the body list;
    // the exact test decides.
    //
    // The local half extents are per shape, and this used to be `size * 0.5f` for all of
    // them. For a sphere `size.x` is the *radius* (Rigid's constructor and aabbHalfExtent
    // both read it that way, and the server stores {radius, 0, 0}), so a sphere was tested
    // as a zero-thickness box of half-extent r/2 in y and z - a ray could not hit one at
    // all. A cylinder had its radius halved the same way.
    const auto localHalfExtent = [](const Rigid &body) noexcept {
        switch (body.shape)
        {
            case ShapeType::Sphere:
            {
                const float r = body.size.x;
                return float3{r, r, r};
            }
            case ShapeType::Cylinder:
                // Local +Z is the axis: radius spans X/Y, half the height spans Z.
                return float3{body.size.x, body.size.x, body.size.z * 0.5f};
            case ShapeType::Box:
            default:
                return body.size * 0.5f;
        }
    };

    const auto testBody = [&](Rigid *body) {
        if ((body->collisionLayer & p_mask) == 0)
            return;
        quat invRot = qconj(body->positionAng);
        float3 o = qrot(invRot, origin - body->positionLin);
        float3 d = qrot(invRot, dir);

        if (body->shape == ShapeType::Sphere)
        {
            // Exact ray/sphere. The box below would only be a proxy with corners that are
            // not on the surface.
            const float r = body->size.x;
            const float b = dot(o, d);
            const float c = length2(o) - r * r;
            const float disc = b * b - c;
            if (disc < 0.0f)
                return;
            const float sq = std::sqrt(disc);
            float t = -b - sq;
            if (t < 0.0f)
                t = -b + sq;
            if (t < 0.0f || t >= bestT)
                return;
            bestT = t;
            bestBody = body;
            bestLocal = o + d * t;
            return;
        }

        // Box, and the box that bounds a cylinder. `dir` is not normalised by the callers,
        // so `t` is in units of `dir` for every shape path.
        const float3 half = localHalfExtent(*body);

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

    manifolds.clear();

    while (bodies)
        delete bodies;

    forceAdjacency.clear();
}

// One tick, as a classified task flow. The broad phase / force warmstart / update-set
// selection run once on the calling thread (they already parallelise their own stages),
// then `runSubstep` runs the per-body warmstart, the colour-parallel iteration region and
// the velocity finish. `runSubstep` builds the flow because interleaved substepping calls
// it more than once per tick.
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

// The flow stages read one substep's arguments and scratch through this context. It is
// declared in solver.h and defined here so the storage and lifetime stay private to the
// implementation; runSubstep owns it on the stack and the flow never outlives it.
struct Solver::SubstepCtx {
    Solver *solver = nullptr;
    int forceCount = 0;
    int targetIterations = 0;
    bool idleDetection = false;
    float prevResidual = 0.0f;
    std::atomic<bool> converged{false};
};

int Solver::flowWarmstartCount(void *ctx) noexcept
{
    return static_cast<int>(static_cast<SubstepCtx *>(ctx)->solver->warmstartOrder.size());
}

void Solver::flowWarmstartItem(void *ctx, int i) noexcept
{
    static_cast<SubstepCtx *>(ctx)->solver->warmstartBody(i);
}

bool Solver::flowIterateInlineWhen(void *ctx) noexcept
{
    const SubstepCtx &c = *static_cast<SubstepCtx *>(ctx);
    return static_cast<int>(c.solver->updateOrder.size()) < kMinBodiesPerDispatch
            && c.forceCount + c.solver->manifolds.count() < kMinForcesPerDispatch;
}

void Solver::flowIterateInline(void *ctx) noexcept
{
    SubstepCtx &c = *static_cast<SubstepCtx *>(ctx);
    c.solver->pool->runRegionInline([&c](int workerIndex, auto &sync) {
        c.solver->iterateBody(workerIndex, sync, c.targetIterations, c.forceCount,
                c.prevResidual, c.converged);
    });
}

void Solver::flowIterateParallel(void *ctx) noexcept
{
    SubstepCtx &c = *static_cast<SubstepCtx *>(ctx);
    c.solver->pool->runRegion([&c](int workerIndex, auto &sync) {
        c.solver->iterateBody(workerIndex, sync, c.targetIterations, c.forceCount,
                c.prevResidual, c.converged);
    });
}

void Solver::flowFinishItem(void *ctx, int i) noexcept
{
    SubstepCtx &c = *static_cast<SubstepCtx *>(ctx);
    c.solver->finishBody(i, c.idleDetection);
}

// One substep: body warmstart (independent), the colour-parallel Gauss-Seidel region,
// then the BDF1 velocity finish (independent). Idle detection only runs on the tick's
// final substep, so substeps cannot sleep bodies faster than ticks.
void Solver::runSubstep(int forceCount, int targetIterations, bool idleDetection)
{
    SubstepCtx ctx;
    ctx.solver = this;
    ctx.forceCount = forceCount;
    ctx.targetIterations = targetIterations;
    ctx.idleDetection = idleDetection;

    flow->clear();
    flow->independent(&Solver::flowWarmstartCount, &Solver::flowWarmstartItem, kMinBodiesPerDispatch,
            &ctx);
    flow->region(&Solver::flowIterateInline, &Solver::flowIterateParallel,
            &Solver::flowIterateInlineWhen, &ctx);
    flow->independent(&Solver::flowWarmstartCount, &Solver::flowFinishItem, kMinBodiesPerDispatch,
            &ctx);

    flow->run();
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
namespace {

[[nodiscard]] float angularCompliance(const Rigid &body, const float3 &localAxis) noexcept
{
    const float3 &I = body.moment;
    float sum = 0.0f;
    for (int k = 0; k < 3; k++)
        if (I[k] > 0.0f)
            sum += localAxis[k] * localAxis[k] / I[k];
    return sum;
}

// World-space rotation vector produced by a unit impulse about `localAxis` (already in the
// body frame), returned in world space.
[[nodiscard]] float3 angularResponse(const Rigid &body, const float3 &localAxis) noexcept
{
    const float3 &I = body.moment;
    float3 out{0.0f, 0.0f, 0.0f};
    for (int k = 0; k < 3; k++)
        if (I[k] > 0.0f)
            out[k] = localAxis[k] / I[k];
    return qrot(body.positionAng, out);
}

// A locked axis does not move, so the projection drops that component. The effective mass
// was computed as if it did, which makes the correction slightly soft on locked bodies -
// the same documented approximation the rest of the solver makes for axis locks.
void applyLocked(Rigid &body, float3 &dLin, float3 &dAng) noexcept
{
    for (int k = 0; k < 3; k++)
    {
        if (body.axisLockLinear & (1 << k))
            dLin[k] = 0.0f;
        if (body.axisLockAngular & (1 << k))
            dAng[k] = 0.0f;
    }
}

} // namespace


void Solver::projectContacts(bool projectPositions, float brakeWindow)
{
    const int manifoldCount = manifolds.count();
    if (manifoldCount <= 0)
        return;

    const float slop = projectionSlop;
    const float beta = projectionBeta;
    const float cap = projectionMaxCorrection;
    const float window = brakeWindow > 0.0f ? brakeWindow : collisionDt;

    // Interception. A pair may only close the gap it has at this pose within the time left
    // before the contact set is refreshed again; otherwise a fast body travels a whole
    // substep into its neighbour and the solver spends the next round undoing that overlap,
    // which is what converts a collision into a launch.
    if (speculativeMargin > 0.0f)
    {
        for (uint32_t slot : manifolds.live)
        {
            Rigid *A = manifolds.bodyA[slot];
            Rigid *B = manifolds.bodyB[slot];
            const int first = ManifoldStore::base(slot);
            const int contactCount = manifolds.contactCount(slot);
            if (contactCount == 0)
                continue;

            const bool moveA = A->mass > 0.0f && A->inUpdateSet;
            const bool moveB = B->mass > 0.0f && B->inUpdateSet;
            if (!moveA && !moveB)
                continue;
            const float invMassA = moveA ? 1.0f / A->mass : 0.0f;
            const float invMassB = moveB ? 1.0f / B->mass : 0.0f;

            // Row 0 of the basis is the contact normal, pointing from B to A.
            const float3 n = row(manifolds.basis[slot], 0);

            // The whole patch is braked as one constraint, so a symmetric face cannot
            // acquire torque from the clamp.
            float3 patchRA{}, patchRB{};
            float velocityError = 0.0f;
            int braked = 0;

            for (int k = 0; k < contactCount; k++)
            {
                const float3 rA = qrot(A->positionAng, manifolds.contacts[first + k].rA);
                const float3 rB = qrot(B->positionAng, manifolds.contacts[first + k].rB);
                const float gap = dot(n, (A->positionLin + rA) - (B->positionLin + rB));
                if (!(gap > 0.0f))
                    continue;

                const float3 vAp = A->velocityLin + cross(A->velocityAng, rA);
                const float3 vBp = B->velocityLin + cross(B->velocityAng, rB);
                const float vn = dot(vAp - vBp, n);
                const float target = -gap / window;
                if (vn < target)
                {
                    patchRA += rA;
                    patchRB += rB;
                    velocityError += target - vn;
                    braked++;
                }
            }

            if (braked == 0)
                continue;

            const float invN = 1.0f / braked;
            const float3 tALocal = moveA ? qrot(qconj(A->positionAng), cross(patchRA * invN, n)) : float3{};
            const float3 tBLocal = moveB ? qrot(qconj(B->positionAng), cross(patchRB * invN, n)) : float3{};
            const float K = invMassA + invMassB
                    + (moveA ? angularCompliance(*A, tALocal) : 0.0f)
                    + (moveB ? angularCompliance(*B, tBLocal) : 0.0f);
            const float j = velocityError * invN / K;
            if (moveA)
            {
                float3 dvLin = n * (invMassA * j);
                float3 dvAng = angularResponse(*A, tALocal) * j;
                applyLocked(*A, dvLin, dvAng);
                A->velocityLin += dvLin;
                A->velocityAng += dvAng;
            }
            if (moveB)
            {
                float3 dvLin = n * (invMassB * j);
                float3 dvAng = angularResponse(*B, tBLocal) * j;
                applyLocked(*B, dvLin, dvAng);
                B->velocityLin -= dvLin;
                B->velocityAng -= dvAng;
            }
        }
    }

    if (!projectPositions || projectionIterations <= 0)
        return;

    projDeltaLin.resize(updateOrder.size());
    projDeltaAng.resize(updateOrder.size());
    projManifoldCount.resize(updateOrder.size());

    for (int iter = 0; iter < projectionIterations; iter++)
    {
        std::fill(projManifoldCount.begin(), projManifoldCount.end(), 0);
        std::fill(projDeltaLin.begin(), projDeltaLin.end(), float3{0, 0, 0});
        std::fill(projDeltaAng.begin(), projDeltaAng.end(), float3{0, 0, 0});

        for (uint32_t slot : manifolds.live)
        {
            Rigid *A = manifolds.bodyA[slot];
            Rigid *B = manifolds.bodyB[slot];
            const int first = ManifoldStore::base(slot);
            const int contactCount = manifolds.contactCount(slot);
            if (contactCount == 0)
                continue;

            const bool moveA = A->mass > 0.0f && A->inUpdateSet;
            const bool moveB = B->mass > 0.0f && B->inUpdateSet;
            if (!moveA && !moveB)
                continue;
            const float invMassA = moveA ? 1.0f / A->mass : 0.0f;
            const float invMassB = moveB ? 1.0f / B->mass : 0.0f;

            // Row 0 of the basis is the contact normal, pointing from B to A.
            const float3 n = row(manifolds.basis[slot], 0);

            // Average the active constraint before computing its mass so a symmetric face
            // cannot acquire torque.
            float3 patchRA{}, patchRB{};
            float patchDepth = 0.0f;
            int corrected = 0;

            for (int k = 0; k < contactCount; k++)
            {
                const float3 rA = qrot(A->positionAng, manifolds.contacts[first + k].rA);
                const float3 rB = qrot(B->positionAng, manifolds.contacts[first + k].rB);
                const float rawDepth = -dot(n, (A->positionLin + rA) - (B->positionLin + rB));

                const float depth = rawDepth - slop;
                if (depth > 0.0f)
                {
                    patchRA += rA;
                    patchRB += rB;
                    patchDepth += depth;
                    corrected++;
                }
            }

            if (corrected == 0)
                continue;

            const float invN = 1.0f / corrected;
            const float3 tALocal = moveA ? qrot(qconj(A->positionAng), cross(patchRA * invN, n)) : float3{};
            const float3 tBLocal = moveB ? qrot(qconj(B->positionAng), cross(patchRB * invN, n)) : float3{};
            const float K = invMassA + invMassB
                    + (moveA ? angularCompliance(*A, tALocal) : 0.0f)
                    + (moveB ? angularCompliance(*B, tBLocal) : 0.0f);
            const float dp = beta * patchDepth * invN / K;
            if (moveA)
            {
                const int i = A->updateIndex;
                projManifoldCount[i]++;
                projDeltaLin[i] += n * (invMassA * dp);
                projDeltaAng[i] += angularResponse(*A, tALocal) * dp;
            }
            if (moveB)
            {
                const int i = B->updateIndex;
                projManifoldCount[i]++;
                projDeltaLin[i] -= n * (invMassB * dp);
                projDeltaAng[i] -= angularResponse(*B, tBLocal) * dp;
            }
        }

        for (size_t i = 0; i < updateOrder.size(); i++)
        {
            Rigid *body = updateOrder[i];
            // Shared contacts must not multiply the correction at seams between colliders.
            const float invCount = 1.0f / std::max(projManifoldCount[i], 1);
            float3 dLin = projDeltaLin[i] * invCount;
            float3 dAng = projDeltaAng[i] * invCount;
            applyLocked(*body, dLin, dAng);
            if (cap > 0.0f)
            {
                const float moved = length(dLin);
                if (moved > cap)
                    dLin = dLin * (cap / moved);
                const float turned = length(dAng);
                if (turned > cap)
                    dAng = dAng * (cap / turned);
            }
            body->positionLin = body->positionLin + dLin;
            body->positionAng = body->positionAng + dAng;
        }
    }
}

// The body on the other side of an adjacency entry, or null when the entry does not
// touch `body`.
[[nodiscard]] static Rigid *otherBody(const ForceRef &ref, const Rigid *body,
        const ManifoldStore &manifolds) noexcept
{
    if (ref.isManifold())
    {
        const uint32_t slot = ref.slot();
        Rigid *a = manifolds.bodyA[slot];
        Rigid *b = manifolds.bodyB[slot];
        if (a == body)
            return b;
        if (b == body)
            return a;
        return nullptr;
    }
    const Force *f = ref.asForce();
    if (f->bodyA == body)
        return f->bodyB;
    if (f->bodyB == body)
        return f->bodyA;
    return nullptr;
}

// True when a force already joins the pair. Mirrors Rigid::constrainedTo: any joint,
// spring or ignore pair counts, and a live manifold counts, but a dormant one (contact
// grace) does not - the pair must still wake on AABB overlap and re-enter the narrow
// phase, or a sleeping body behind a dormant manifold is invisible to an incoming body.
[[nodiscard]] static bool pairConstrained(const ForceAdjacency &adjacency,
        const ManifoldStore &manifolds, Rigid *bodyA, Rigid *bodyB)
{
    for (const ForceRef &ref : adjacency.of(bodyA->broadIndex))
    {
        if (otherBody(ref, bodyA, manifolds) != bodyB)
            continue;
        if (ref.isManifold())
        {
            const uint32_t slot = ref.slot();
            if (manifolds.contactCount(slot) == 0 && manifolds.dormantTicks[slot] > 0)
                continue;
        }
        return true;
    }
    return false;
}

// A broad-phase candidate test: swept world AABBs plus Godot's layer/mask pair
// rule, minus pairs already joined by a force (IgnoreCollision, joints). The box
// starts from the tight bounds of the actual shape - unlike a bounding-sphere test
// whose radius is the shape's diagonal and therefore admits many far-apart
// candidates - and is grown to cover the body's motion for the step.
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

// Contact detection: hierarchical spatial hash (big-block screen + layered
// hash), rebuilt every step.
//
// The old x-axis sweep emitted 28k candidate pairs on a 40-layer resting wall,
// and the Morton LBVH that replaced it partitioned in 3D at the cost of a tree
// build and a recursive self-query. The hash bucketing also partitions in 3D,
// but with none of either: bodies go into dyadic levels by tight-AABB size,
// each level is a sorted key table, and candidate pairs come from a
// 27-cell walk. The candidate set is a superset of the tight-AABB overlap set,
// exactly as the tree's was.
//
// Determinism: pair collection is integer-sorted and de-duplicated, and the
// final (min-index, other) sort fixes the manifold creation order - the same
// contract as the sweep and tree paths before it.
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
        if (broadPhaseMode == 1)
            broadLbvh.buildAndCollectPairs(bodiesInOrder, sweepPairs, bodyAabbHalf,
                    bodyAabbCenter, *pool);
        else
            broadHash.buildAndCollectPairs(bodiesInOrder, 0.0f, sweepPairs, bodyAabbHalf,
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
int Solver::warmstartForces()
{
    const int forceCount = collectForces();
    forceActive.resize(forceCount);
    pool->forCount(forceCount, [this](int i) {
        forceActive[i] = forceOrder[i]->initialize() ? 1 : 0;
    }, kMinForcesPerDispatch);

    // One pass both reaps and compacts: `forceOrder` is list order, and deleting a force
    // only unlinks it from that list, so the survivors keep the exact order a second
    // collectForces() would produce.
    int write = 0;
    for (int i = 0; i < forceCount; i++)
    {
        Force *f = forceOrder[i];
        if (!forceActive[i])
        {
            delete f; // unlinks itself from the solver and both bodies
            continue;
        }
        forceOrder[write++] = f;
    }
    const int liveForceCount = write;

    const int manifoldCount = manifolds.count();
    manifoldActive.resize(static_cast<size_t>(manifoldCount));
    pool->forCount(manifoldCount, [this](int i) {
        const uint32_t slot = manifolds.live[i];
        // Frozen constraints keep their state verbatim; nothing they read can move.
        manifoldActive[i] = manifolds.frozen(slot) || manifolds.initialize(slot, *this) ? 1 : 0;
        if (traceMatch(manifolds.bodyA[slot], manifolds.bodyB[slot]))
            std::printf("TRACE init slot=%u active=%d contacts=%d frozen=%d\n", slot,
                    (int)manifoldActive[i], manifolds.contactCount(slot),
                    (int)manifolds.frozen(slot));
    }, kMinForcesPerDispatch);

    // Dormant manifolds (contact grace) keep their slot and hot state for a few ticks
    // instead of being reaped. initialize only enters dormancy while the grace setting is
    // positive, so with the default 0 this branch never fires.
    int liveWrite = 0;
    for (int i = 0; i < manifoldCount; i++)
    {
        const uint32_t slot = manifolds.live[i];
        if (!manifoldActive[i])
        {
            if (contactGraceTicks > 0 && manifolds.dormantTicks[slot] > 0)
            {
                manifolds.live[liveWrite++] = slot;
                continue;
            }
#ifdef AVBD_DIAG
            diag.reapedManifolds++;
#endif
            manifolds.release(slot);
            continue;
        }
        manifolds.live[liveWrite++] = slot;
    }
    manifolds.live.resize(static_cast<size_t>(liveWrite));

#ifdef AVBD_DIAG
    {
        int soft = 0;
        const float softLimit = 0.1f * betaLin;
        for (uint32_t slot : manifolds.live)
            for (int ci = 0, n = manifolds.contactCount(slot); ci < n; ci++)
                if (manifolds.contactPenalty[ManifoldStore::base(slot) + ci].x < softLimit)
                    soft++;
        diag.softContacts = soft;
    }
#endif

    // The live manifold set and force list are final for the step; the colouring, primal,
    // wake-up and island passes read this adjacency.
    buildForceAdjacency();
    return liveForceCount;
}

// Compute one body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so the step's flow runs it as an Independent stage. It covers
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
}

// One region body for the whole iteration set: `iterations` rounds of the
// colour-parallel Gauss-Seidel primal followed by a dual pass, with a barrier
// fence at every phase boundary.
//
// The primal runs colour by colour: every body of colour c assembles its full
// block system against the poses the earlier colours already applied, solves it
// and applies the update immediately - Gauss-Seidel across colours, parallel
// and race-free within one (a colour is an independent set of the constraint
// graph by construction). Manifolds take the self-contained updatePrimal path,
// which consumes the previous round's cache for the first touch and recomputes
// once for the second; each manifold is seen by at most one body per colour, so
// nothing races.
//
// The colouring is capped at kMaxColours so the barrier count per round is
// bounded; bodies the greedy colouring could not fit run one block-Jacobi round
// after the colours (precompute, solve against the post-colour snapshot,
// combine), using jacobiRelaxation as the combine control.
//
// NullSync (small scenes) runs the identical body inline in list order, which is
// exactly the threads=1 reference the digest tests compare against.
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

// One body's BDF1 velocity after the final iteration. Per-body and independent, so the
// step's flow runs it as an Independent stage; it covers static bodies too, so their
// bookkeeping matches, though only movable ones get a new velocity.
void Solver::finishBody(int i, bool idleDetection)
{
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
    // selection then simulated it again (sleep_mode NEVER), but the finish stage
    // skipped its velocity update because `sleeping` was still set - the body
    // moved with a permanently zero velocity, and contacts could pump it upward
    // indefinitely. That was the pyramid_dbl runaway.
    if (idleDetection && sleepFrames > 0 && body->mass > 0 && !body->sleeping && body->sleep_mode == 1)
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

    // Body -> island slot through the broad-phase index: one array write per island body,
    // then every constraint is two O(1) lookups (no pointer sort and binary search).
    bodyIsland.assign(bodiesInOrder.size(), -1);
    for (int i = 0; i < count; ++i)
        bodyIsland[static_cast<size_t>(islandBodies[i]->broadIndex)] = i;

    const auto join = [&](Rigid *a, Rigid *b) {
        if (a == nullptr || b == nullptr)
            return;
        if (a->broadIndex < 0 || b->broadIndex < 0)
            return;
        const int ia = bodyIsland[static_cast<size_t>(a->broadIndex)];
        const int ib = bodyIsland[static_cast<size_t>(b->broadIndex)];
        if (ia >= 0 && ib >= 0)
            parent[find(ia)] = find(ib);
    };
    for (uint32_t slot : manifolds.live)
        join(manifolds.bodyA[slot], manifolds.bodyB[slot]);
    for (Force *f = forces; f != 0; f = f->next)
        join(f->bodyA, f->bodyB);

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
            for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
            {
                Rigid *other = otherBody(ref, body, manifolds);
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
                    length2(body->velocityLin) == 0.0f &&
                    forceAdjacency.of(body->broadIndex).empty())
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

    // Per-body delta buffers for the Jacobi primal's split-solve / combine pair.
    bodyDeltaLin.resize(updateOrder.size());
    bodyDeltaAng.resize(updateOrder.size());

    // Overflow-bucket eligibility: a body carrying any non-manifold force (joint, spring,
    // world anchor) is routed into the overflow bucket by colourGraph - the contact map
    // is clamped, hard joints do not tolerate the snapshot solve. Every force in the
    // list is non-manifold; contacts do not count.
    for (Force *f = forces; f != 0; f = f->next)
    {
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

        for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
        {
            Rigid *other = otherBody(ref, body, manifolds);
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
    for (int i = 0; i < n; i++) {
        if (contactOnlyOverflow && updateOrder[i]->contactOnly) {
            bodyColour[i] = kOverflowColour;
            overflowOrder.push_back(i);
        }
        if (bodyColour[i] != kOverflowColour) {
            colourBegins[bodyColour[i] + 1]++;
        }
    }
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
}

// Fill `forceOrder` with the solver's non-manifold forces, so the phases can address them
// by index. Manifolds live in the store and are addressed by slot.
int Solver::collectForces()
{
    forceOrder.clear();
    forceOrder.insert(forceOrder.end(), next_range(forces).begin(), next_range(forces).end());
    return static_cast<int>(forceOrder.size());
}

void Solver::buildForceAdjacency()
{
    forceAdjacency.build(bodiesInOrder, forces, manifolds);
}

void Solver::ensurePool()
{
    // Pass the value straight through: 0 means "one worker per hardware thread" and 1 means
    // a single worker, which both the pool and `threads` document the same way. (Clamping 0
    // up to 1 here would silently turn the default into a serial run.)
    const int requested = std::max(0, threads);
    if (pool != nullptr && poolThreads == requested)
        return;

    pool = std::make_unique<detail::JobPool>(static_cast<unsigned>(requested));
    flow = std::make_unique<detail::TaskFlow>(*pool);
    poolThreads = requested;
}

int Solver::threadCount() const
{
    if (pool == nullptr)
        return 1;
    return static_cast<int>(pool->threadCount());
}

} // namespace avbd
