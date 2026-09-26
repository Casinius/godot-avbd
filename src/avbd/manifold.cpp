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
}

namespace {
// Lock-free stack using std::atomic pointer (MSPMC)
// Most manifolds are reused within same thread (thread-local reuse).
// Memory ordering: acquire for pop (read-then-use), release for push (publish-after-store).
static std::atomic<void *> freeList{nullptr};
} // namespace

void *Manifold::operator new(std::size_t count)
{
    if (count != sizeof(Manifold))
        return ::operator new(count);

    // Try to pop from free list (lock-free)
    void *ptr = freeList.exchange(nullptr, std::memory_order_acquire);

    if (ptr == nullptr)
    {
        // Free list empty, allocate from heap
        return ::operator new(sizeof(Manifold));
    }

    return ptr;
}

void Manifold::operator delete(void *ptr) noexcept
{
    if (ptr == nullptr)
        return;

    // Push to free list (lock-free)
    // Use release to publish to other threads
    freeList.store(ptr, std::memory_order_release);
}

void Manifold::drainPool() noexcept
{
    // Collect all free list pointers (single-threaded cleanup)
    void *ptr = freeList.exchange(nullptr, std::memory_order_acquire);
    while (ptr != nullptr)
    {
        void *next = freeList.load(std::memory_order_acquire);
        ::operator delete(ptr);
        ptr = next;
    }
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

    // Merge old contact data with new contacts
    for (int i = 0; i < newNumContacts; i++)
    {
        for (int j = 0; j < numContacts; j++)
        {
            if (newContacts[i].feature.key == contacts[j].feature.key)
            {
                float3 newRA = newContacts[i].rA;
                float3 newRB = newContacts[i].rB;
                newContacts[i] = contacts[j];

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

void Manifold::updatePrimal(Rigid *body, float alpha, Block &block)
{
    float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    float3 dqBAng = bodyB->positionAng - bodyB->initialAng;

    // Same rotation-matrix hoist as updateDual: body quats are loop-invariant here.
    const float3x3 rotA(bodyA->positionAng);
    const float3x3 rotB(bodyB->positionAng);

    // Incremental jacobian caching: the angular jacobian is a function of the world-space
    // moment arms only (rows are rWorld x basis.row(i)), and the arms move between
    // iterations solely because the body rotations do. The drift test is per BODY (the
    // two body rotations are loop-invariant): hoisted out of the contact loop so the
    // per-contact work in the cached path is a straight copy of the cached jacobians.
    const bool cacheable = solver->jacobianRebuildDistance > 0.0f;
    const float rebuild2 = solver->jacobianRebuildDistance * solver->jacobianRebuildDistance;
    const float3 dA = rotA.col(0) - rotACaptured.col(0);
    const float3 dB = rotB.col(0) - rotBCaptured.col(0);
    const bool cacheUsable = cacheable && jacobiansCached &&
            dA.squaredNorm() <= rebuild2 && dB.squaredNorm() <= rebuild2;

    for (int i = 0; i < numContacts; i++)
    {
        JacobianCache &cache = jacobianCache[i];

        float3 rAWorld, rBWorld;
        float3x3 jAAng, jBAng;
        if (cacheUsable)
        {
            rAWorld = cache.rAWorld;
            rBWorld = cache.rBWorld;
            jAAng = cache.jAAng;
            jBAng = cache.jBAng;
        }
        else
        {
            rAWorld = rotA * contacts[i].rA;
            rBWorld = rotB * contacts[i].rB;

            // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
            const float3x3 jALin = basis;
            const float3x3 jBLin = -basis;
            // Row i of the angular jacobian is the moment arm crossed with row i of the linear one
            // (the original custom float3x3 was row-major; Eigen's .row() keeps the semantics).
            for (int r = 0; r < 3; r++)
            {
                jAAng.row(r) = rAWorld.cross(jALin.row(r));
                jBAng.row(r) = rBWorld.cross(jBLin.row(r));
            }

            if (cacheable)
            {
                cache.rAWorld = rAWorld;
                cache.rBWorld = rBWorld;
                cache.jAAng = jAAng;
                cache.jBAng = jBAng;
                if (i == numContacts - 1) {
                    rotACaptured = rotA;
                    rotBCaptured = rotB;
                    jacobiansCached = true;
                }
            }
        }

        const float3x3 K = diagonal(contacts[i].penalty.x(), contacts[i].penalty.y(), contacts[i].penalty.z());
        const float3x3 jALin = basis;
        const float3x3 jBLin = -basis;
        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force
        float3 F = K * C + contacts[i].lambda;

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

        // Choose jacobian depending on input body
        float3x3 jLin = body == bodyA ? jALin : jBLin;
        float3x3 jAng = body == bodyA ? jAAng : jBAng;

        // Stamp into LHS
        float3x3 jLinT = jLin.transpose();
        float3x3 jAngT = jAng.transpose();
        float3x3 jAngTk = jAngT * K;

        block.lhsLin += jLinT * K * jLin;
        block.lhsAng += jAngTk * jAng;
        block.lhsCross += jAngTk * jLin;

        // Stamp into RHS
        block.rhsLin += jLinT * F;
        block.rhsAng += jAngT * F;
    }
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

    for (int i = 0; i < numContacts; i++)
    {
        float3 rAWorld = rotA * contacts[i].rA;
        float3 rBWorld = rotB * contacts[i].rB;

        // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
        float3x3 jALin = basis;
        float3x3 jBLin = -basis;
        // Row i of the angular jacobian is the moment arm crossed with row i of the linear one.
        float3x3 jAAng, jBAng;
        for (int j = 0; j < 3; j++)
        {
            jAAng.row(j) = rAWorld.cross(jALin.row(j));
            jBAng.row(j) = rBWorld.cross(jBLin.row(j));
        }

        float3x3 K = diagonal(contacts[i].penalty.x(), contacts[i].penalty.y(), contacts[i].penalty.z());
        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force
        float3 F = K * C + contacts[i].lambda;

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
        if (F[0] < 0)
            contacts[i].penalty[0] = std::min(contacts[i].penalty[0] + solver->betaLin * std::fabs(C[0]), PENALTY_MAX);
        if (frictionScale <= bounds)
        {
            contacts[i].penalty[1] = std::min(contacts[i].penalty[1] + solver->betaLin * std::fabs(C[1]), PENALTY_MAX);
            contacts[i].penalty[2] = std::min(contacts[i].penalty[2] + solver->betaLin * std::fabs(C[2]), PENALTY_MAX);
            contacts[i].stick = float2{C[1], C[2]}.norm() < STICK_THRESH;
        }
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
