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
    : Force(p_solver, p_bodyA, p_bodyB, ForceKind::Contact), numContacts(0)
{
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
    // Local view of the knobs this function reads: warm-start decay and the adaptive-beta
    // ramp rule. Hoisting keeps the contact loop below free of the parameter path.
    const SolverParams::Dual &dual = solver->params.dual;

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
                // Penetration depth, not |C0|: C0's normal component carries the collision
                // margin, so a resting contact would otherwise look like a constant error.
                prevErr[i] = std::max(0.0f, COLLISION_MARGIN - contacts[j].C0.x());

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
        if (penetration > PENETRATION_MAX)
        {
            contacts[i].C0 = contacts[i].C0.normalized() * PENETRATION_MAX;
        }

        // Adaptive beta (the augmented-Lagrangian kappa update, per contact): a
        // sticking contact whose error did not shrink against last step doubles its
        // personal ramp multiplier; a well-converging one halves it back toward 1.
        // Unmatched (fresh) contacts start from the neutral multiplier. The residual is
        // penetration depth (the margin is not an error to converge), and a residual at
        // or below the floor is treated as solved so float-scale noise on a resting
        // contact cannot ratchet the multiplier to the cap.
        if (dual.adaptiveBeta && prevErr[i] >= 0.0f)
        {
            const float err = std::max(0.0f, COLLISION_MARGIN - contacts[i].C0.x());
            const float floor = dual.adaptiveBetaErrorFloor;
            float scale = contacts[i].betaScale;
            if (err <= floor)
                scale = std::max(scale * 0.5f, 1.0f);
            else if (prevErr[i] <= floor)
                scale = std::min(scale * 2.0f, ADAPTIVE_BETA_MAX);
            else if (err > 0.5f * prevErr[i])
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
        contacts[i].lambda = contacts[i].lambda * dual.alpha * dual.gamma;
        if (dual.adaptiveBeta)
        {
            // Adaptive beta is a per-contact ramp multiplier, so the penalty it produced
            // must survive into the next step: decay the carried value (Eq. 19) instead
            // of restarting from the floor, which threw the ramp away and made a resting
            // body re-grow its support every frame. Unmatched contacts are zero-initialised,
            // so this lands on PENALTY_MIN for a fresh contact.
            contacts[i].penalty = (contacts[i].penalty * dual.gamma)
                    .cwiseMax(float3{PENALTY_MIN, PENALTY_MIN, PENALTY_MIN})
                    .cwiseMin(float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});
        }
        else
        {
            const float clampedLambdaPenalty = std::clamp(PENALTY_MIN * dual.gamma, PENALTY_MIN, PENALTY_MAX);
            contacts[i].penalty = float3{clampedLambdaPenalty, clampedLambdaPenalty, clampedLambdaPenalty};
        }
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
    // Read once per round, not once per contact: the drift test below is per body, and the
    // value cannot change inside this loop.
    const float rebuildDistance = solver->params.primal.jacobianRebuildDistance;

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
    const bool cacheable = rebuildDistance > 0.0f;
    const float rebuild2 = rebuildDistance * rebuildDistance;
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

            // Always landed in the cache: stampCached reads these for this round even
            // when the cross-round drift feature is off (the capture bookkeeping below
            // is what jacobianRebuildDistance gates).
            cache.rAWorld = rAWorld;
            cache.rBWorld = rBWorld;
            cache.jAAng = jAAng;
            cache.jBAng = jBAng;
            if (cacheable && i == numContacts - 1) {
                rotACaptured = rotA;
                rotBCaptured = rotB;
                jacobiansCached = true;
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

        contactForceCache[i] = F;
    }
}

void Manifold::stampCached(Rigid *body, Block &block)
{
    const bool fromIsA = body == bodyA;
    for (int i = 0; i < numContacts; i++)
    {
        const float3 &F = contactForceCache[i];
        const float3x3 K = diagonal(contacts[i].penalty.x(), contacts[i].penalty.y(), contacts[i].penalty.z());
        const float3x3 jLin = fromIsA ? basis : -basis;
        const float3x3 &jAng = fromIsA ? jacobianCache[i].jAAng : jacobianCache[i].jBAng;

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

void Manifold::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // Self-contained fallback; the parallel primal drives the two phases directly.
    preCompute(alpha);
    stampCached(body, block);
}

// The projection's measure pass: for each contact row that is on the wrong side of its own
// equilibrium, the multiplier that would take it to zero, split across the two bodies by their
// compliance - and only that, written into this manifold's cache. Nothing here touches a body;
// the share pass reads these back and each body applies its own half.
//
// The error is recomputed here rather than taken from `contactForceCache`, and from the current
// poses rather than the round's snapshot: this pass runs after the dual update, and the whole
// point is to report what is left once every constraint has had its say.
//
// Unlike the primal and the dual it is measured *unstabilised* - C0 + J.dq instead of the
// C0.(1-alpha) + J.dq those two read. Both expressions share one zero, the speculative-margin
// plane a few millimetres inside touching where the solver wants a contact to rest anyway, so the
// two passes aim at the same surface; they differ only in how much of the distance to it one step
// is allowed to close. alpha is what makes a deep contact heal over several steps rather than
// one, and that delay is exactly what this layer exists to pay for: the guarantee a penalty method
// cannot give is a legal pose at the end of the step, and a projection that asked only for the
// percent the stabilisation leaves on the table could not give it.
void Manifold::projectMeasure(float scale, float maxCorrection)
{
    // The ceiling is on how far one *body* may be moved by one pass, and a body takes shares from
    // many contacts, so it is applied where those are summed - in the share pass.
    (void)maxCorrection;

    const Compliance cA = complianceOf(bodyA);
    const Compliance cB = complianceOf(bodyB);
    // Leaving the cache as it is here is safe, and it is worth saying why: a stale share is only
    // ever read by a body that can move, and a body that can move takes this early return away.
    if (cA.immobile() && cB.immobile())
        return;

    const float dtSq = solver->dt * solver->dt;
    const float3 dqALin = bodyA->positionLin - bodyA->initialLin;
    const float3 dqAAng = bodyA->positionAng - bodyA->initialAng;
    const float3 dqBLin = bodyB->positionLin - bodyB->initialLin;
    const float3 dqBAng = bodyB->positionAng - bodyB->initialAng;
    const float3x3 rotA(bodyA->positionAng);
    const float3x3 rotB(bodyB->positionAng);

    // A multiplier is a force once it is divided by the square of the step, because that is how
    // the primal reads one back: a body's response to a stamped force is dt^2 * F / m. Charging
    // the projected correction to the multiplier says "this much of the contact force has already
    // been delivered, by moving the pose instead of pushing through it", which is what stops the
    // next round's penalty ramp from re-buying the same correction.
    const bool feedback = solver->params.projection.dualFeedback && dtSq > 0.0f;

    // A contact patch is one constraint seen through up to eight points, and this is a Jacobi
    // shares one correction evenly instead - exact for a single contact, and for a patch it moves
    // the bodies by the mean of what its corners asked for. The disagreement between the corners
    // is left to the angular half, where it becomes the torque that levels the box.
    const float patch = 1.0f / static_cast<float>(numContacts > 0 ? numContacts : 1);

    for (int i = 0; i < numContacts; i++)
    {
        JacobianCache &cache = jacobianCache[i];
        cache.projLambda = float3{0, 0, 0};
        cache.projectedRows = 0;
        // Out of scope. Those two fields are the share pass's only window into this row, so
        // clearing them and stopping is what takes a contact out of the layer; skipping the clear
        // would re-apply whatever the previous pass measured.
        if (!solver->params.projection.contacts)
            continue;

        const float3 rAWorld = rotA * contacts[i].rA;
        const float3 rBWorld = rotB * contacts[i].rB;
        float3x3 jAAng, jBAng;
        for (int r = 0; r < 3; r++)
        {
            jAAng.row(r) = rAWorld.cross(basis.row(r));
            jBAng.row(r) = rBWorld.cross(-basis.row(r));
        }

        const float3 C = contacts[i].C0 + basis * dqALin - basis * dqBLin
                + jAAng * dqAAng + jBAng * dqBAng;

        // Normal first: the friction cone the tangential rows may use is sized by it.
        float3 mu = float3{0, 0, 0};
        uint8_t rows = 0;
        if (C[0] < 0.0f)
        {
            const float w = cA.row(basis.row(0), rAWorld) + cB.row(basis.row(0), rBWorld);
            if (w > 0.0f)
            {
                mu[0] = -C[0] * scale / w;
                rows = 1;
            }
        }

        // Tangential rows only carry static friction, and only up to the cone the normal
        // correction opened: a sliding contact is *supposed* to move along its tangents, and
        // projecting its error would be projecting the slide away.
        if (rows != 0 && contacts[i].stick)
        {
            const float cone = friction * std::fabs(mu[0]);
            for (int r = 1; r < 3; r++)
            {
                if (C[r] == 0.0f)
                    continue;
                const float w = cA.row(basis.row(r), rAWorld) + cB.row(basis.row(r), rBWorld);
                if (w > 0.0f)
                    mu[r] = -C[r] * scale / w;
            }
            const float slip = float2{mu[1], mu[2]}.norm();
            if (slip > cone && slip > 0.0f)
            {
                mu[1] *= cone / slip;
                mu[2] *= cone / slip;
            }
            if (mu[1] != 0.0f)
                rows |= 2;
            if (mu[2] != 0.0f)
                rows |= 4;
        }

        // Scaled after the cone clamp, so the ratio between the rows the cone allows is kept.
        mu = mu * patch;
        cache.projLambda = mu;
        cache.projectedRows = rows;

        if (feedback && rows != 0)
        {
            float3 &lambda = contacts[i].lambda;
            lambda = lambda - mu / dtSq;
            // The same two clamps updateDual applies to a force, for the same reasons: a contact
            // can only press, and friction cannot exceed what pressing provides.
            lambda[0] = std::min(lambda[0], 0.0f);
            const float bounds = std::fabs(lambda[0]) * friction;
            const float held = float2{lambda[1], lambda[2]}.norm();
            if (held > bounds && held > 0.0f)
            {
                lambda[1] *= bounds / held;
                lambda[2] *= bounds / held;
            }
        }
    }
}

void Manifold::projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const
{
    const bool isA = body == bodyA;
    if (!isA && body != bodyB)
        return;

    const Compliance c = complianceOf(isA ? bodyA : bodyB);
    if (c.immobile())
        return;

    const float3x3 rot(isA ? bodyA->positionAng : bodyB->positionAng);
    const float sign = isA ? 1.0f : -1.0f;
    for (int i = 0; i < numContacts; i++)
    {
        const uint8_t rows = jacobianCache[i].projectedRows;
        if (rows == 0)
            continue;
        const float3 &mu = jacobianCache[i].projLambda;
        const float3 arm = rot * (isA ? contacts[i].rA : contacts[i].rB);
        for (int r = 0; r < 3; r++)
        {
            if ((rows & (1 << r)) == 0)
                continue;
            // The row's Jacobian flips sign with the side it is stamped into, exactly as the
            // primal's does, so the two bodies are always moved towards the constraint, never
            // one towards and one away from the same error.
            const float3 n = basis.row(r);
            r_deltaLin = r_deltaLin + n * (sign * mu[r] * c.invMass);
            r_deltaAng = r_deltaAng + c.invInertiaTimes(arm.cross(n)) * (sign * mu[r]);
        }
    }
}

void Manifold::updateDual(float alpha)
{
    const float betaLinBase = solver->params.dual.betaLin;

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

        // A row the projection already moved to its target this round is not a residual the
        // penalty needs to grow against: ramping it would buy the correction a second time, and
        // the second bill is what turns a projection into an oscillator. The mask is consumed
        // here, so a row only skips the ramp for the one dual pass that follows the pass which
        // projected it. With the layer off the mask is always empty and this changes nothing.
        uint8_t projected = 0;
        if (solver->params.projection.dualFeedback)
        {
            projected = jacobianCache[i].projectedRows;
            jacobianCache[i].projectedRows = 0;
        }

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        const float beta = betaLinBase * contacts[i].betaScale;
        if (F[0] < 0 && (projected & 1) == 0)
            contacts[i].penalty[0] = std::min(contacts[i].penalty[0] + beta * std::fabs(C[0]), PENALTY_MAX);
        if (frictionScale <= bounds)
        {
            if ((projected & 2) == 0)
                contacts[i].penalty[1] = std::min(contacts[i].penalty[1] + beta * std::fabs(C[1]), PENALTY_MAX);
            if ((projected & 4) == 0)
                contacts[i].penalty[2] = std::min(contacts[i].penalty[2] + beta * std::fabs(C[2]), PENALTY_MAX);
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
