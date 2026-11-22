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

namespace {

// Effective angular compliance along `localAxis`: t.(I^-1 t), with I rotated into world space.
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


// Locked axes drop their correction component; effective mass still assumes mobility, so the
// correction is slightly soft - the solver's documented axis-lock approximation.
void applyLocked(Rigid &body, float3 &dLin, float3 &dAng) noexcept
{
    maskLockedLinear(body.axisLockLinear, dLin);
    maskLockedAngular(body.axisLockAngular, dAng);
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

    // Compact to movable pairs once: parked/sleeping pairs cannot change pose or velocity, so
    // testing them is a no-op; the brake and every projection pass share this small list.
    projectionActive.clear();
    for (uint32_t slot : manifolds.live)
    {
        if (manifolds.contactCount(slot) == 0)
            continue;
        const Rigid *A = manifolds.bodyA[slot];
        const Rigid *B = manifolds.bodyB[slot];
        if (!(A->mass > 0.0f && A->inUpdateSet) && !(B->mass > 0.0f && B->inUpdateSet))
            continue;
        projectionActive.push_back(slot);
    }
    if (projectionActive.empty())
        return;

    // Interception: close only the current gap within the time left before contacts refresh, or a
    // fast body tunnels a whole substep and the undo becomes a launch.
    if (speculativeMargin > 0.0f)
    {
        for (uint32_t slot : projectionActive)
        {
            Rigid *A = manifolds.bodyA[slot];
            Rigid *B = manifolds.bodyB[slot];
            const int first = ManifoldStore::base(slot);
            const int contactCount = manifolds.contactCount(slot);

            const bool moveA = A->mass > 0.0f && A->inUpdateSet;
            const bool moveB = B->mass > 0.0f && B->inUpdateSet;
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

    // Sub-slop residual corrections are below float noise, so resting/sleeping piles stop after the
    // first pass with nothing left to move.
    const float depthTolerance = std::max(slop * 0.01f, 1.0e-6f);

    for (int iter = 0; iter < projectionIterations; iter++)
    {
        std::fill(projManifoldCount.begin(), projManifoldCount.end(), 0);
        std::fill(projDeltaLin.begin(), projDeltaLin.end(), float3{0, 0, 0});
        std::fill(projDeltaAng.begin(), projDeltaAng.end(), float3{0, 0, 0});

        float maxDepth = 0.0f;
        for (uint32_t slot : projectionActive)
        {
            Rigid *A = manifolds.bodyA[slot];
            Rigid *B = manifolds.bodyB[slot];
            const int first = ManifoldStore::base(slot);
            const int contactCount = manifolds.contactCount(slot);

            const bool moveA = A->mass > 0.0f && A->inUpdateSet;
            const bool moveB = B->mass > 0.0f && B->inUpdateSet;
            const float invMassA = moveA ? 1.0f / A->mass : 0.0f;
            const float invMassB = moveB ? 1.0f / B->mass : 0.0f;

            // Row 0 of the basis is the contact normal, pointing from B to A.
            const float3 n = row(manifolds.basis[slot], 0);

            // Average the active constraint before computing its mass so a symmetric face
            // cannot acquire torque.
            float3 patchRA{}, patchRB{};
            float patchDepth = 0.0f;
            int corrected = 0;

            // Body rotations are loop-invariant across the contact patch, so one matrix
            // per body replaces a quaternion rotation per contact point.
            const float3x3 rotA = qmat(A->positionAng);
            const float3x3 rotB = qmat(B->positionAng);
            for (int k = 0; k < contactCount; k++)
            {
                const float3 rA = mul(rotA, manifolds.contacts[first + k].rA);
                const float3 rB = mul(rotB, manifolds.contacts[first + k].rB);
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
            maxDepth = std::max(maxDepth, patchDepth * invN);
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

        // The applied corrections were below the tolerance: another pass would only
        // re-measure float noise, so the remaining iterations are skipped.
        if (maxDepth <= depthTolerance)
            break;
    }
}

} // namespace avbd
