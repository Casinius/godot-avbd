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

    // New contact set: cached moment arms describe the old geometry, drop them.
    armsCached = false;
    primalCacheReady = false;

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
                prevErr[i] = length(contacts[j].C0);
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
        float3 xA = bodyA->positionLin + qrot(bodyA->positionAng, contacts[i].rA);
        float3 xB = bodyB->positionLin + qrot(bodyB->positionAng, contacts[i].rB);
        contacts[i].C0 = mul(basis, xA - xB) + float3{COLLISION_MARGIN, 0, 0};

        // The reference implementation reports the full penetration depth here:
        // a deep impact then produces a proportionally large error and force, and
        // the primal pushes the body back out within the step. Capping the depth
        // (as this port once did, at 2 cm) removes the restoring gradient past the
        // cap - bodies settle anywhere below it and rebound off the stored spring
        // energy, which showed up as endless bouncing and sink-through on the
        // round-shape scenes. No clamp.

        // Adaptive beta (the augmented-Lagrangian kappa update, per contact): a
        // sticking contact whose error did not shrink against last step doubles its
        // personal ramp multiplier; a well-converging one halves it back toward 1.
        // Unmatched (fresh) contacts start from the neutral multiplier.
        if (solver->adaptiveBeta && prevErr[i] > 0.0f)
        {
            const float err = length(contacts[i].C0);
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
        // Warmstart the penalty too (Eq. 19). Resetting it to the floor every frame
        // threw away the ramp the previous step had already built, so a resting body
        // had to re-grow its support from PENALTY_MIN each frame and sank until the
        // dual happened to catch up. The reference implementation decays and clamps
        // the carried value instead; a fresh contact starts from the floor because
        // its Contact was value-initialised.
        contacts[i].penalty = clamp(contacts[i].penalty * solver->gamma, float3{PENALTY_MIN, PENALTY_MIN, PENALTY_MIN},
                float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});
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
    const float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    const float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    const float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    const float3 dqBAng = bodyB->positionAng - bodyB->initialAng;

    // Same rotation-matrix hoist as updateDual: body quats are loop-invariant here.
    const float3x3 rotA = qmat(bodyA->positionAng);
    const float3x3 rotB = qmat(bodyB->positionAng);

    // Constraint value in the basis frame. The linear jacobian rows are the basis rows
    // (±), and each angular row is rWorld x that basis row, so the four matrix-vector
    // products the old assembly built collapse into one bracket the basis then rotates:
    //   C = C0(1-alpha) + B * [ (dqA_lin - dqB_lin) + (dqA_ang x rA) - (dqB_ang x rB) ]
    // No angular jacobian matrices are needed: the basis is the only frame C lives in.
    for (int i = 0; i < numContacts; i++)
    {
        const float3 rAWorld = mul(rotA, contacts[i].rA);
        const float3 rBWorld = mul(rotB, contacts[i].rB);
        jacobianCache[i].rAWorld = rAWorld;
        jacobianCache[i].rBWorld = rBWorld;

        const float3 bracket = (dqALin - dqBLin) + cross(dqAAng, rAWorld) - cross(dqBAng, rBWorld);
        const float3 C = contacts[i].C0 * (1 - alpha) + mul(basis, bracket);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = cmul(contacts[i].penalty, C) + contacts[i].lambda;

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction;
        float frictionScale = length(float2{F[1], F[2]});
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        contactForceCache[i] = F;
    }
    primalCacheReady = true;
}

// void Manifold::stampCached(Rigid *body, Block &block)
// {
//     const bool fromIsA = body == bodyA;
//     const float side = fromIsA ? 1.0f : -1.0f;

//     // Both bodies share the linear jacobian up to sign, Jlin = ±basis, so the linear
//     // Hessian Jlin^T K Jlin = B^T K B is the same for A and B. K is diagonal, so that
//     // "basis sandwich" S is just a weighted sum of the three basis-row outer products.
//     // With R = skew(rWorld) the angular jacobian is Jang = ±B R^T, which gives the whole
//     // stamp in closed form:
//     //   lhsLin   = S
//     //   lhsAng   = Jang^T K Jang = R S R^T
//     //   lhsCross = Jang^T K Jlin = R S
//     //   rhsLin   = ±B^T F
//     //   rhsAng   = ±R (B^T F)
//     // This replaces the three generic 3x3 products (and their transposes) the old stamp
//     // built per contact, per body. The symmetric results are written mirrored so the
//     // pivot fallback in Block::solveDelta still sees a full matrix.
//     const float3 b0 = row(basis, 0);
//     const float3 b1 = row(basis, 1);
//     const float3 b2 = row(basis, 2);

//     for (int i = 0; i < numContacts; i++)
//     {
//         const float3 &F = contactForceCache[i];
//         const float3 k = contacts[i].penalty;
//         const float3 r = fromIsA ? jacobianCache[i].rAWorld : jacobianCache[i].rBWorld;

//         // S = B^T K B, symmetric. Six values (s01 = S[0][1], and so on).
//         const float s00 = k.x * b0.x * b0.x + k.y * b1.x * b1.x + k.z * b2.x * b2.x;
//         const float s11 = k.x * b0.y * b0.y + k.y * b1.y * b1.y + k.z * b2.y * b2.y;
//         const float s22 = k.x * b0.z * b0.z + k.y * b1.z * b1.z + k.z * b2.z * b2.z;
//         const float s01 = k.x * b0.x * b0.y + k.y * b1.x * b1.y + k.z * b2.x * b2.y;
//         const float s02 = k.x * b0.x * b0.z + k.y * b1.x * b1.z + k.z * b2.x * b2.z;
//         const float s12 = k.x * b0.y * b0.z + k.y * b1.y * b1.z + k.z * b2.y * b2.z;

//         // Column-major: column j of S is (S[0][j], S[1][j], S[2][j]).
//         const float3x3 S{float3{s00, s01, s02}, float3{s01, s11, s12}, float3{s02, s12, s22}};

//         // R S, one cross product per column of S.
//         const float3x3 RS{cross(r, S.x), cross(r, S.y), cross(r, S.z)};

//         // R S R^T, symmetric. Closed form in r and S (verified against the matrix
//         // product); h01 = H[0][1], h02 = H[0][2], h12 = H[1][2].
//         const float h00 = s11 * r.z * r.z - 2.0f * s12 * r.y * r.z + s22 * r.y * r.y;
//         const float h11 = s22 * r.x * r.x - 2.0f * s02 * r.x * r.z + s00 * r.z * r.z;
//         const float h22 = s00 * r.y * r.y - 2.0f * s01 * r.x * r.y + s11 * r.x * r.x;
//         const float h01 = s12 * r.x * r.z + s02 * r.y * r.z - s22 * r.x * r.y - s01 * r.z * r.z;
//         const float h02 = s12 * r.x * r.y + s01 * r.y * r.z - s11 * r.x * r.z - s02 * r.y * r.y;
//         const float h12 = s02 * r.x * r.y + s01 * r.x * r.z - s00 * r.y * r.z - s12 * r.x * r.x;

//         // B^T F and its moment about the contact arm: the only per-body part of the RHS.
//         const float3 fw = b0 * F.x + b1 * F.y + b2 * F.z;
//         const float3 fwAng = cross(r, fw);

//         block.lhsLin += S;
//         block.lhsAng += float3x3{float3{h00, h01, h02}, float3{h01, h11, h12}, float3{h02, h12, h22}};
//         block.lhsCross += RS;
//         block.rhsLin += fw * side;
//         block.rhsAng += fwAng * side;
//     }
// }


void Manifold::stampCached(Rigid *body, Block &block)
{
    const bool fromIsA = body == bodyA;
    const float side = fromIsA ? 1.0f : -1.0f;

    // Both bodies share the linear jacobian up to sign, Jlin = ±basis, so the linear
    // Hessian Jlin^T K Jlin = B^T K B is the same for A and B. K is diagonal, so that
    // "basis sandwich" S is just a weighted sum of the three basis-row outer products.
    // With R = skew(rWorld) the angular jacobian is Jang = ±B R^T, which gives the whole
    // stamp in closed form:
    //   lhsLin   = S
    //   lhsAng   = Jang^T K Jang = R S R^T
    //   lhsCross = Jang^T K Jlin = R S
    //   rhsLin   = ±B^T F
    //   rhsAng   = ±R (B^T F)
    // This replaces the three generic 3x3 products (and their transposes) the old stamp
    // built per contact, per body. The symmetric results are written mirrored so the
    // pivot fallback in Block::solveDelta still sees a full matrix.
    const float3 b0 = row(basis, 0);
    const float3 b1 = row(basis, 1);
    const float3 b2 = row(basis, 2);

    // Outer-product components of the basis rows: loop invariants for S.
    const float b0x2 = b0.x * b0.x, b0y2 = b0.y * b0.y, b0z2 = b0.z * b0.z;
    const float b0xy = b0.x * b0.y, b0xz = b0.x * b0.z, b0yz = b0.y * b0.z;
    const float b1x2 = b1.x * b1.x, b1y2 = b1.y * b1.y, b1z2 = b1.z * b1.z;
    const float b1xy = b1.x * b1.y, b1xz = b1.x * b1.z, b1yz = b1.y * b1.z;
    const float b2x2 = b2.x * b2.x, b2y2 = b2.y * b2.y, b2z2 = b2.z * b2.z;
    const float b2xy = b2.x * b2.y, b2xz = b2.x * b2.z, b2yz = b2.y * b2.z;

    // Accumulate in registers, then write back once.
    float3x3 accLin{};
    float3x3 accAng{};
    float3x3 accCross{};
    float3 accRhsLin{};
    float3 accRhsAng{};

    for (int i = 0; i < numContacts; i++)
    {
        const float3 &F = contactForceCache[i];
        const float3 k = contacts[i].penalty;
        const float3 r = fromIsA ? jacobianCache[i].rAWorld : jacobianCache[i].rBWorld;

        // S = B^T K B, symmetric. Six values (s01 = S[0][1], and so on).
        const float s00 = k.x * b0x2 + k.y * b1x2 + k.z * b2x2;
        const float s11 = k.x * b0y2 + k.y * b1y2 + k.z * b2y2;
        const float s22 = k.x * b0z2 + k.y * b1z2 + k.z * b2z2;
        const float s01 = k.x * b0xy + k.y * b1xy + k.z * b2xy;
        const float s02 = k.x * b0xz + k.y * b1xz + k.z * b2xz;
        const float s12 = k.x * b0yz + k.y * b1yz + k.z * b2yz;

        // Column-major: column j of S is (S[0][j], S[1][j], S[2][j]).
        const float3x3 S{float3{s00, s01, s02}, float3{s01, s11, s12}, float3{s02, s12, s22}};

        // R S, one cross product per column of S.
        const float3x3 RS{cross(r, S.x), cross(r, S.y), cross(r, S.z)};

        // Precompute second-order products of r.
        const float rx2 = r.x * r.x, ry2 = r.y * r.y, rz2 = r.z * r.z;
        const float rxy = r.x * r.y, rxz = r.x * r.z, ryz = r.y * r.z;

        // R S R^T, symmetric. Closed form in r and S (verified against the matrix
        // product); h01 = H[0][1], h02 = H[0][2], h12 = H[1][2].
        const float h00 = s11 * rz2 - 2.0f * s12 * ryz + s22 * ry2;
        const float h11 = s22 * rx2 - 2.0f * s02 * rxz + s00 * rz2;
        const float h22 = s00 * ry2 - 2.0f * s01 * rxy + s11 * rx2;
        const float h01 = s12 * rxz + s02 * ryz - s22 * rxy - s01 * rz2;
        const float h02 = s12 * rxy + s01 * ryz - s11 * rxz - s02 * ry2;
        const float h12 = s02 * rxy + s01 * rxz - s00 * ryz - s12 * rx2;

        // B^T F and its moment about the contact arm: the only per-body part of the RHS.
        const float3 fw = b0 * F.x + b1 * F.y + b2 * F.z;
        const float3 fwAng = cross(r, fw);

        accLin += S;
        accAng += float3x3{float3{h00, h01, h02}, float3{h01, h11, h12}, float3{h02, h12, h22}};
        accCross += RS;
        accRhsLin += fw;
        accRhsAng += fwAng;
    }

    block.lhsLin += accLin;
    block.lhsAng += accAng;
    block.lhsCross += accCross;
    block.rhsLin += accRhsLin * side;
    block.rhsAng += accRhsAng * side;
}

void Manifold::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // The first body to touch this manifold in a round consumes the force/jacobian cache
    // the previous round's dual left behind: no pose in this computation has moved since
    // then (the other body has not been solved yet either), so recomputing would produce
    // the same floats. The second body's colour runs after poses have moved, so it
    // recomputes - that is the round's one fresh evaluation.
    if (!primalCacheReady)
        preCompute(alpha);
    stampCached(body, block);
    primalCacheReady = false;
}

void Manifold::updateDual(float alpha)
{
    const float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    const float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    const float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    const float3 dqBAng = bodyB->positionAng - bodyB->initialAng;

    // Both contact rotations share the two body quats, which never change inside this
    // loop: one rotation-matrix conversion (~20 ns) replaces two quat multiplies per
    // contact (~16 ns each at 8 contacts = ~256 ns/manifold).
    const float3x3 rotA = qmat(bodyA->positionAng);
    const float3x3 rotB = qmat(bodyB->positionAng);

    // Drift-gated arm reuse (the jacobianRebuildDistance semantics, now living in the
    // dual): when neither body's rotation moved further than the threshold since the last
    // capture, the cached world moment arms are reused for C and for the next round's
    // stamp. Only the arms are cached now - the angular jacobian matrices are gone.
    const bool cacheable = solver->jacobianRebuildDistance > 0.0f;
    const float rebuild2 = solver->jacobianRebuildDistance * solver->jacobianRebuildDistance;
    const float3 dA = rotA[0] - rotACaptured[0];
    const float3 dB = rotB[0] - rotBCaptured[0];
    const bool driftUsable = cacheable && armsCached &&
            length2(dA) <= rebuild2 && length2(dB) <= rebuild2;

    for (int i = 0; i < numContacts; i++)
    {
        float3 rAWorld, rBWorld;
        if (driftUsable)
        {
            rAWorld = jacobianCache[i].rAWorld;
            rBWorld = jacobianCache[i].rBWorld;
        }
        else
        {
            rAWorld = mul(rotA, contacts[i].rA);
            rBWorld = mul(rotB, contacts[i].rB);
            jacobianCache[i].rAWorld = rAWorld;
            jacobianCache[i].rBWorld = rBWorld;
        }

        const float3 bracket = (dqALin - dqBLin) + cross(dqAAng, rAWorld) - cross(dqBAng, rBWorld);
        float3 C = contacts[i].C0 * (1 - alpha) + mul(basis, bracket);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = cmul(contacts[i].penalty, C) + contacts[i].lambda;

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction;
        float frictionScale = length(float2{F[1], F[2]});
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
            contacts[i].stick = length(float2{C[1], C[2]}) < STICK_THRESH;
        }

        // Leave the primal's cache ready for the next Newton round. The pose is frozen
        // between here and the next stampCached (the solve runs after the next
        // preCompute would), so the next preCompute would recompute exactly this bracket
        // - and, with the ramped penalty and the just-stored lambda, exactly this force.
        // Doing it here removes the standalone preCompute from every round after the first.
        float3 fNext = cmul(contacts[i].penalty, C) + contacts[i].lambda;
        fNext[0] = std::min(fNext[0], 0.0f);
        const float boundsNext = std::fabs(fNext[0]) * friction;
        const float frictionScaleNext = length(float2{fNext[1], fNext[2]});
        if (frictionScaleNext > boundsNext && frictionScaleNext > 0)
        {
            fNext[1] *= boundsNext / frictionScaleNext;
            fNext[2] *= boundsNext / frictionScaleNext;
        }
        contactForceCache[i] = fNext;
    }

    if (!driftUsable && cacheable)
    {
        rotACaptured = rotA;
        rotBCaptured = rotB;
        armsCached = true;
    }

    // Ready for the next round's first primal touch: nothing moves between here and it.
    primalCacheReady = true;
}

// Accumulate this manifold's penetration error for the adaptive iteration estimate.
// Same loop, comparison and accumulation order as the dynamic_cast block it replaces
// in solveIterations, so the floating-point sequence is unchanged.
void Manifold::accumulatePenetration(float &r_total, int &r_count) const
{
    for (int i = 0; i < numContacts; i++)
    {
        const float penetration = length(contacts[i].C0);
        if (penetration > 0.0f)
        {
            r_total += penetration;
            r_count++;
        }
    }
}

} // namespace avbd
