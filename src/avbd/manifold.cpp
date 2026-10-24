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
#include <array>
#include <atomic>
#include <mutex>
#include <vector>

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

Manifold::Manifold(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB)
    : Force(p_solver, p_bodyA, p_bodyB), numContacts(0)
{
    isManifold = true;
}

namespace {
// Manifold pool with simple mutex for correctness
static std::mutex poolMutex;
static std::vector<void *> poolFreeList;
} // namespace

void *Manifold::operator new(std::size_t count)
{
    if (count != sizeof(Manifold))
        return ::operator new(count);

    std::lock_guard<std::mutex> lock(poolMutex);
    if (!poolFreeList.empty())
    {
        void *ptr = poolFreeList.back();
        poolFreeList.pop_back();
        return ptr;
    }

    return ::operator new(sizeof(Manifold));
}

void Manifold::operator delete(void *ptr) noexcept
{
    if (ptr == nullptr)
        return;

    std::lock_guard<std::mutex> lock(poolMutex);
    poolFreeList.push_back(ptr);
}

void Manifold::drainPool() noexcept
{
    std::lock_guard<std::mutex> lock(poolMutex);
    for (void *ptr : poolFreeList)
        ::operator delete(ptr);
    poolFreeList.clear();
}

bool Manifold::initialize()
{
    // Compute friction
    friction = std::sqrt(bodyA->friction * bodyB->friction);

    // New contact set: cached jacobians describe the old geometry, drop them.
    jacobiansCached = false;

    // Compute new contacts
    std::array<Contact, 8> newContacts = {};
    int newNumContacts = collide(bodyA, bodyB, std::span<Contact>(newContacts), basis);

    // Merge old contact data with new contacts. prevErr[i] carries the matched
    // contact's error from the previous step so the adaptive-beta update below can
    // compare it against the fresh error.
    std::array<float, 8> prevErr{};
    prevErr.fill(-1.0f);
    for (int i = 0; i < newNumContacts; i++)
    {
        for (int j = 0; j < numContacts; j++)
        {
            if (newContacts[i].feature.key == contacts[j].feature.key)
            {
                float3 newRA = newContacts[i].rA;
                float3 newRB = newContacts[i].rB;
                newContacts[i] = contacts[j];
                prevErr[i] = contacts[j].C0.norm();
                newContacts[i].betaScale = contacts[j].betaScale;

                // If no static friction in last frame, use the new contact point locations
                if (!contacts[j].stick)
                {
                    newContacts[i].rA = newRA;
                    newContacts[i].rB = newRB;
                    newContacts[i].lambda = contacts[j].lambda;   // Keep lambda from prev frame
                    newContacts[i].penalty = contacts[j].penalty; // Keep penalty from prev frame
                }
                break;
            }
        }
    }

    // Copy new contacts to the manifold
    numContacts = newNumContacts;
    for (int i = 0; i < numContacts; i++)
        contacts[i] = newContacts[i];

    // Compute error at q- and update penalty and lambdas
    for (int i = 0; i < numContacts; i++)
    {
        // Error at q-
        float3 xA = bodyA->positionLin + bodyA->positionAng * contacts[i].rA;
        float3 xB = bodyB->positionLin + bodyB->positionAng * contacts[i].rB;
        contacts[i].C0 = basis * (xA - xB) + float3{COLLISION_MARGIN, 0, 0};

        // Clamp penetration to avoid excessive penalty forces
        float penetration = contacts[i].C0.norm();
        if (penetration > 0.02f) // Cap at 2cm penetration
        {
            contacts[i].C0 = contacts[i].C0.normalized() * 0.02f;
        }

        // Adaptive beta (the augmented-Lagrangian kappa update, per contact): a
        // sticking contact whose error did not shrink against last step doubles its
        // personal ramp multiplier; a well-converging one halves it back toward 1.
        // Unmatched (fresh) contacts start from the neutral multiplier.
        if (solver->adaptiveBeta && prevErr[i] > 0.0f)
        {
            const float err = contacts[i].C0.norm();
            float scale = contacts[i].betaScale;
            if (err > 0.5f * prevErr[i])
                scale = std::min(scale * 2.0f, ADAPTIVE_BETA_MAX);
            else if (err < 0.1f * prevErr[i])
                scale = std::max(scale * 0.5f, 1.0f);
            contacts[i].betaScale = scale;
        }
        else
        {
            contacts[i].betaScale = 1.0f;
        }

        // Initial penalty (Eq. 19): the original computed base * 10^(2*log10(pen)) —
        // algebraically pen^2 — and clamped to [PENALTY_MIN, PENALTY_MAX]. pen is capped
        // at 0.02 above, so the product is <= 4e-4 < PENALTY_MIN: the clamp always
        // saturates at the floor and the log10/pow chain was dead computation
        // (~300 us/step across a full manifold set). Penalty starts at PENALTY_MIN and is
        // ramped per iteration in updateDual (betaLin * |C|, capped by PENALTY_MAX).
        // The warmstart decay below is kept bit-for-bit.
        contacts[i].lambda = contacts[i].lambda * solver->alpha * solver->gamma;
        const float clampedLambdaPenalty = std::clamp(PENALTY_MIN * solver->gamma, PENALTY_MIN, PENALTY_MAX);
        contacts[i].penalty = float3{clampedLambdaPenalty, clampedLambdaPenalty, clampedLambdaPenalty};
    }

    return numContacts > 0;
}

// The primal runs in two phases so the two bodies of a manifold share one assembly:
// preCompute (parallel over forces, once per round) evaluates everything both stamps
// need - world arms, angular jacobians, and the clamped contact forces - from the
// round's pose snapshot; stampCached (parallel over bodies) stamps each body's side
// from that cache. Before the split, each of the two per-round updatePrimal calls
// recomputed the whole shared half and used only its own side of the stamps. Same
// inputs in the same order, so the cached values are bit-identical to what the stamps
// used to recompute.
void Manifold::preCompute(float alpha)
{
    float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    float3 dqBAng = bodyB->positionAng - bodyB->initialAng;

    // Same rotation-matrix hoist as updateDual: body quats are loop-invariant here.
    const float3x3 rotA(bodyA->positionAng);
    const float3x3 rotB(bodyB->positionAng);

    // Round-0 shared assembly. The angular jacobian is a function of the world-space
    // moment arms only (rows are rWorld x basis.row(i)); computed fresh because the
    // step's warmstart moved every pose, which always exceeds the drift threshold.
    // Rounds >= 1 get their jacobians from updateDual's drift-gated pass instead.
    const float3x3 jALin = basis;
    const float3x3 jBLin = -basis;
    for (int i = 0; i < numContacts; i++)
    {
        JacobianCache &cache = jacobianCache[i];

        const float3 rAWorld = rotA * contacts[i].rA;
        const float3 rBWorld = rotB * contacts[i].rB;

        // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
        // Row i of the angular jacobian is the moment arm crossed with row i of the linear one
        // (the original custom float3x3 was row-major; Eigen's .row() keeps the semantics).
        float3x3 jAAng, jBAng;
        for (int r = 0; r < 3; r++)
        {
            jAAng.row(r) = rAWorld.cross(jALin.row(r));
            jBAng.row(r) = rBWorld.cross(jBLin.row(r));
        }

        // Landed in the cache: stampCached reads these for this round.
        cache.rAWorld = rAWorld;
        cache.rBWorld = rBWorld;
        cache.jAAng = jAAng;
        cache.jBAng = jBAng;

        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = contacts[i].penalty.cwiseProduct(C) + contacts[i].lambda;

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction;
        float frictionScale = float2{F[1], F[2]}.norm();
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        contactForceCache[i] = F;
    }
}

void Manifold::stampCached(Rigid *body, Block &block)
{
    const bool fromIsA = body == bodyA;
    for (int i = 0; i < numContacts; i++)
    {
        const float3 &F = contactForceCache[i];
        // The penalty is a diagonal: every K product below is a row/column scaling,
        // not a generic 3x3 product (the diagonal zeros turned a third of the old
        // multiplies into bookkeeping).
        const auto kDiag = contacts[i].penalty.asDiagonal();
        const float3x3 jLin = fromIsA ? basis : -basis;
        const float3x3 &jAng = fromIsA ? jacobianCache[i].jAAng : jacobianCache[i].jBAng;

        // Stamp into LHS
        const float3x3 jLinT = jLin.transpose();
        const float3x3 jAngTk = jAng.transpose() * kDiag;

        block.lhsLin += (jLinT * kDiag) * jLin;
        block.lhsAng += jAngTk * jAng;
        block.lhsCross += jAngTk * jLin;

        // Stamp into RHS
        block.rhsLin += jLinT * F;
        block.rhsAng += jAng.transpose() * F;
    }
}

void Manifold::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // Self-contained fallback; the parallel primal drives the two phases directly.
    preCompute(alpha);
    stampCached(body, block);
}

void Manifold::updateDual(float alpha)
{
    float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    float3 dqBAng = bodyB->positionAng - bodyB->initialAng;

    // Both contact rotations share the two body quats, which never change inside this
    // loop: one rotation-matrix conversion (~20 ns) replaces two quat multiplies per
    // contact (~16 ns each at 8 contacts = ~256 ns/manifold).
    const float3x3 rotA(bodyA->positionAng);
    const float3x3 rotB(bodyB->positionAng);

    // Drift-gated jacobian reuse (the jacobianRebuildDistance semantics, now living in
    // the dual): when neither body's rotation moved further than the threshold since the
    // last capture, the cached jacobians are reused for C and the next round's stamp -
    // the arms and crosses are simply not computed this round.
    const bool cacheable = solver->jacobianRebuildDistance > 0.0f;
    const float rebuild2 = solver->jacobianRebuildDistance * solver->jacobianRebuildDistance;
    const float3 dA = rotA.col(0) - rotACaptured.col(0);
    const float3 dB = rotB.col(0) - rotBCaptured.col(0);
    const bool driftUsable = cacheable && jacobiansCached &&
            dA.squaredNorm() <= rebuild2 && dB.squaredNorm() <= rebuild2;

    const float3x3 jALin = basis;
    const float3x3 jBLin = -basis;
    for (int i = 0; i < numContacts; i++)
    {
        float3x3 jAAng, jBAng;
        if (driftUsable)
        {
            jAAng = jacobianCache[i].jAAng;
            jBAng = jacobianCache[i].jBAng;
        }
        else
        {
            const float3 rAWorld = rotA * contacts[i].rA;
            const float3 rBWorld = rotB * contacts[i].rB;

            // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
            // Row i of the angular jacobian is the moment arm crossed with row i of the linear one.
            for (int j = 0; j < 3; j++)
            {
                jAAng.row(j) = rAWorld.cross(jALin.row(j));
                jBAng.row(j) = rBWorld.cross(jBLin.row(j));
            }
            jacobianCache[i].jAAng = jAAng;
            jacobianCache[i].jBAng = jBAng;
            if (cacheable && i == numContacts - 1)
            {
                rotACaptured = rotA;
                rotBCaptured = rotB;
                jacobiansCached = true;
            }
        }

        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = contacts[i].penalty.cwiseProduct(C) + contacts[i].lambda;

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction;
        float frictionScale = float2{F[1], F[2]}.norm();
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        // Store updated force
        contacts[i].lambda = F;

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        const float beta = solver->betaLin * contacts[i].betaScale;
        if (F[0] < 0)
            contacts[i].penalty[0] = std::min(contacts[i].penalty[0] + beta * std::fabs(C[0]), PENALTY_MAX);
        if (frictionScale <= bounds)
        {
            contacts[i].penalty[1] = std::min(contacts[i].penalty[1] + beta * std::fabs(C[1]), PENALTY_MAX);
            contacts[i].penalty[2] = std::min(contacts[i].penalty[2] + beta * std::fabs(C[2]), PENALTY_MAX);
            contacts[i].stick = float2{C[1], C[2]}.norm() < STICK_THRESH;
        }

        // Leave the primal's cache ready for the next Newton round. The pose is frozen
        // between here and the next stampCached (the solve runs after the next
        // preCompute would), so the next preCompute would recompute exactly this
        // jacobian - and, with the ramped penalty and the just-stored lambda, exactly
        // this force. Doing it here removes the standalone preCompute from every round
        // after the first.
        jacobianCache[i].jAAng = jAAng;
        jacobianCache[i].jBAng = jBAng;
        float3 fNext = contacts[i].penalty.cwiseProduct(C) + contacts[i].lambda;
        fNext[0] = std::min(fNext[0], 0.0f);
        const float boundsNext = std::fabs(fNext[0]) * friction;
        const float frictionScaleNext = float2{fNext[1], fNext[2]}.norm();
        if (frictionScaleNext > boundsNext && frictionScaleNext > 0)
        {
            fNext[1] *= boundsNext / frictionScaleNext;
            fNext[2] *= boundsNext / frictionScaleNext;
        }
        contactForceCache[i] = fNext;
    }
}

// Accumulate this manifold's penetration error for the adaptive iteration estimate.
// Same loop, comparison and accumulation order as the dynamic_cast block it replaces
// in solveIterations, so the floating-point sequence is unchanged.
void Manifold::accumulatePenetration(float &r_total, int &r_count) const
{
    for (int i = 0; i < numContacts; i++)
    {
        const float penetration = contacts[i].C0.norm();
        if (penetration > 0.0f)
        {
            r_total += penetration;
            r_count++;
        }
    }
}

} // namespace avbd
