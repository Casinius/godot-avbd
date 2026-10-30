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

namespace {
// The symmetric linear Hessian S = B^T K B of one contact, as its six independent entries. The basis
// rows are the linear Jacobian rows (up to sign) and K is the penalty diagonal, so S is a weighted
// sum of the three basis-row outer products: six numbers, and the same six for both bodies of the
// manifold because the sign cancels in B^T K B. initialize() and updateDual() are the only places a
// penalty changes, so they are the only places this is recomputed.
inline Manifold::BasisHessian basisHessian(const float3x3 &basis, const float3 &k) noexcept
{
    const float3 b0 = basis.row(0);
    const float3 b1 = basis.row(1);
    const float3 b2 = basis.row(2);
    Manifold::BasisHessian s;
    s.s00 = k[0] * b0.x() * b0.x() + k[1] * b1.x() * b1.x() + k[2] * b2.x() * b2.x();
    s.s01 = k[0] * b0.x() * b0.y() + k[1] * b1.x() * b1.y() + k[2] * b2.x() * b2.y();
    s.s02 = k[0] * b0.x() * b0.z() + k[1] * b1.x() * b1.z() + k[2] * b2.x() * b2.z();
    s.s11 = k[0] * b0.y() * b0.y() + k[1] * b1.y() * b1.y() + k[2] * b2.y() * b2.y();
    s.s12 = k[0] * b0.y() * b0.z() + k[1] * b1.y() * b1.z() + k[2] * b2.y() * b2.z();
    s.s22 = k[0] * b0.z() * b0.z() + k[1] * b1.z() * b1.z() + k[2] * b2.z() * b2.z();
    return s;
}

// The constraint value of one contact in the basis frame, from the bracket the assembly builds:
//   C = C0 (1 - alpha) + B * [ (dqA_lin - dqB_lin) + dqA_ang x rA - dqB_ang x rB ]
// The linear Jacobian rows are the basis rows and each angular row is the world arm crossed with its
// basis row, so the four matrix-vector products the assembly used to evaluate collapse into one
// bracket the basis rotates - no angular Jacobian matrices at all. Shared by the primal and the dual
// so the two never disagree about which error they are looking at.
inline float3 constraintError(const float3x3 &basis, const float3 &c0, float alpha, const float3 &dqLin,
        const float3 &dqAAng, const float3 &dqBAng, const float3 &rAWorld, const float3 &rBWorld) noexcept
{
    const float3 bracket = dqLin + dqAAng.cross(rAWorld) - dqBAng.cross(rBWorld);
    return c0 * (1.0f - alpha) + basis * bracket;
}
} // namespace

bool Manifold::initialize()
{
    // Local view of the knobs this function reads: warm-start decay and the adaptive-beta
    // ramp rule. Hoisting keeps the contact loop below free of the parameter path.
    const SolverParams::Dual &dual = solver->params.dual;

    // Compute friction
    friction = std::sqrt(bodyA->friction * bodyB->friction);

    // New contact set: cached jacobians describe the old geometry, drop them.
    jacobiansCached = false;
    primalReady = false;

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

        // The penalty is final for this step's first round, so the shared linear Hessian can be
        // built here - the first round's stamps read it, and updateDual() refreshes it from then on.
        linHessian[i] = basisHessian(basis, contacts[i].penalty);
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
    // Read once per round, not once per contact: the drift test below is per body, and the value
    // cannot change inside this loop.
    const float rebuildDistance = solver->params.primal.jacobianRebuildDistance;

    // The two bodies' derived pose state, which the solver refreshes wherever a pose changes (see
    // Rigid::rot): two matrix loads and four vector loads, instead of a quaternion conversion and
    // four pose subtractions per contact.
    const float3 dqALin = bodyA->dqLin;
    const float3 dqAAng = bodyA->dqAng;
    const float3 dqBLin = bodyB->dqLin;
    const float3 dqBAng = bodyB->dqAng;
    const float3x3 rotA = bodyA->rot;
    const float3x3 rotB = bodyB->rot;
    const float3 dqLin = dqALin - dqBLin;

    // Incremental arm caching: the world arms move between rounds only because the body rotations
    // do, so later rounds reuse them while neither rotation has drifted further than the threshold
    // from the capture. The drift test is per BODY (the two rotations are loop-invariant here).
    const bool cacheable = rebuildDistance > 0.0f;
    const float rebuild2 = rebuildDistance * rebuildDistance;
    const float3 dA = rotA.col(0) - rotACaptured.col(0);
    const float3 dB = rotB.col(0) - rotBCaptured.col(0);
    const bool cacheUsable = cacheable && jacobiansCached &&
            dA.squaredNorm() <= rebuild2 && dB.squaredNorm() <= rebuild2;

    for (int i = 0; i < numContacts; i++)
    {
        JacobianCache &cache = jacobianCache[i];
        const float3 rAWorld = cacheUsable ? cache.rAWorld : rotA * contacts[i].rA;
        const float3 rBWorld = cacheUsable ? cache.rBWorld : rotB * contacts[i].rB;
        cache.rAWorld = rAWorld;
        cache.rBWorld = rBWorld;

        const float3 C = constraintError(basis, contacts[i].C0, alpha, dqLin, dqAAng, dqBAng, rAWorld, rBWorld);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        const float3 &k = contacts[i].penalty;
        float3 F{k[0] * C[0] + contacts[i].lambda[0], k[1] * C[1] + contacts[i].lambda[1],
                k[2] * C[2] + contacts[i].lambda[2]};

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        const float bounds = std::fabs(F[0]) * friction;
        const float frictionScale = float2{F[1], F[2]}.norm();
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        contactForceCache[i] = F;
    }

    if (cacheable && !cacheUsable)
    {
        rotACaptured = rotA;
        rotBCaptured = rotB;
        jacobiansCached = true;
    }

    // The stamps of the round this ran for can read the cache now. updateDual() keeps it current
    // from round to round, which is what removes this pass from every round after the first.
    primalReady = true;
}

void Manifold::stampCached(Rigid *body, Block &block)
{
    const bool fromIsA = body == bodyA;
    const float side = fromIsA ? 1.0f : -1.0f;
    const float3 b0 = basis.row(0);
    const float3 b1 = basis.row(1);
    const float3 b2 = basis.row(2);

    // Both bodies share the linear Jacobian up to sign (J_lin = +-B), so the linear Hessian
    // J_lin^T K J_lin = B^T K B = S is six numbers both of them can read (Manifold::BasisHessian),
    // and with R = skew(rWorld) the angular Jacobian is Jang = +-B R^T. The whole stamp then has a
    // closed form:
    //   lhsLin   = S
    //   lhsAng   = Jang^T K Jang   = R S R^T
    //   lhsCross = Jang^T K J_lin  = R S          (sign-invariant: both Jacobians flip together)
    //   rhsLin   = +-B^T F
    //   rhsAng   = +-R (B^T F)
    // which replaces three generic 3x3 triple products and their transposes per contact per body.
    // Accumulated in registers and written into the block once, so a manifold's eight contacts cost
    // one load and one store of the block instead of eight.
    float3x3 accLin = float3x3::Zero();
    float3x3 accAng = float3x3::Zero();
    float3x3 accCross = float3x3::Zero();
    float3 accRhsLin = float3::Zero();
    float3 accRhsAng = float3::Zero();

    for (int i = 0; i < numContacts; i++)
    {
        const float3 &F = contactForceCache[i];
        const float3 &r = fromIsA ? jacobianCache[i].rAWorld : jacobianCache[i].rBWorld;
        const BasisHessian &s = linHessian[i];

        float3x3 S;
        S << s.s00, s.s01, s.s02,
             s.s01, s.s11, s.s12,
             s.s02, s.s12, s.s22;

        // R S: one cross product per column of S.
        float3x3 RS;
        RS.col(0) = r.cross(S.col(0));
        RS.col(1) = r.cross(S.col(1));
        RS.col(2) = r.cross(S.col(2));

        // R S R^T, symmetric, expanded in r: 18 products instead of two 3x3 multiplies.
        const float rx2 = r.x() * r.x(), ry2 = r.y() * r.y(), rz2 = r.z() * r.z();
        const float rxy = r.x() * r.y(), rxz = r.x() * r.z(), ryz = r.y() * r.z();
        const float h00 = s.s11 * rz2 - 2.0f * s.s12 * ryz + s.s22 * ry2;
        const float h11 = s.s22 * rx2 - 2.0f * s.s02 * rxz + s.s00 * rz2;
        const float h22 = s.s00 * ry2 - 2.0f * s.s01 * rxy + s.s11 * rx2;
        const float h01 = s.s12 * rxz + s.s02 * ryz - s.s22 * rxy - s.s01 * rz2;
        const float h02 = s.s12 * rxy + s.s01 * ryz - s.s11 * rxz - s.s02 * ry2;
        const float h12 = s.s02 * rxy + s.s01 * rxz - s.s00 * ryz - s.s12 * rx2;
        float3x3 H;
        H << h00, h01, h02,
             h01, h11, h12,
             h02, h12, h22;

        // B^T F and its moment about the arm: the only per-body part of the right-hand side.
        const float3 fw = b0 * F[0] + b1 * F[1] + b2 * F[2];

        accLin += S;
        accCross += RS;
        accAng += H;
        accRhsLin += fw;
        accRhsAng += r.cross(fw);
    }

    block.lhsLin += accLin;
    block.lhsCross += accCross;
    block.lhsAng += accAng;
    block.rhsLin += accRhsLin * side;
    block.rhsAng += accRhsAng * side;
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
    const float3 dqALin = bodyA->dqLin;
    const float3 dqAAng = bodyA->dqAng;
    const float3 dqBLin = bodyB->dqLin;
    const float3 dqBAng = bodyB->dqAng;
    const float3x3 rotA = bodyA->rot;
    const float3x3 rotB = bodyB->rot;

    // A multiplier is a force once it is divided by the square of the step, because that is how
    // the primal reads one back: a body's response to a stamped force is dt^2 * F / m. Charging
    // the projected correction to the multiplier says "this much of the contact force has already
    // been delivered, by moving the pose instead of pushing through it", which is what stops the
    // next round's penalty ramp from re-buying the same correction.
    //
    // One ledger per kind, because the two kinds need opposite things. A contact is unilateral and
    // clamped, so charging its multiplier is how the lift a pass delivered is *held* across rounds
    // (the clamps keep the charge in the cone, so it cannot run away). A joint is bilateral and
    // unbounded, so a charge there is an oscillator - the reference has to move instead. Contacts
    // therefore keep `dualFeedback` regardless of this knob, and the joint classes read it as
    // "charge the multiplier or re-anchor instead".
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

        // Unstabilised - alpha = 0 - because this layer reports what is left today, not what the
        // step's stabilisation is willing to leave. Same bracket form as the primal's.
        const float3 C = constraintError(basis, contacts[i].C0, 0.0f, dqALin - dqBLin, dqAAng, dqBAng,
                rAWorld, rBWorld);

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

#ifdef AVBD_DIAG
        // What the pass asked to remove, and - when it is charged - what that came to as a force.
        atomicAddFloat(solver->diag.askedContact,
                std::fabs(mu[0]) + std::fabs(mu[1]) + std::fabs(mu[2]));
#endif

        if (feedback && rows != 0)
        {
#ifdef AVBD_DIAG
            atomicAddFloat(solver->diag.chargeContact,
                    (std::fabs(mu[0]) + std::fabs(mu[1]) + std::fabs(mu[2])) / dtSq);
#endif
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

    // The same per-body derived pose state preCompute reads (see Rigid::rot), refreshed by the
    // combine pass, which is what moved the bodies since that read.
    const float3 dqALin = bodyA->dqLin;
    const float3 dqAAng = bodyA->dqAng;
    const float3 dqBLin = bodyB->dqLin;
    const float3 dqBAng = bodyB->dqAng;
    const float3x3 rotA = bodyA->rot;
    const float3x3 rotB = bodyB->rot;
    const float3 dqLin = dqALin - dqBLin;

    for (int i = 0; i < numContacts; i++)
    {
        const float3 rAWorld = rotA * contacts[i].rA;
        const float3 rBWorld = rotB * contacts[i].rB;

        // The same error the primal solves against, at the pose this round leaves.
        const float3 C = constraintError(basis, contacts[i].C0, alpha, dqLin, dqAAng, dqBAng, rAWorld, rBWorld);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        const float3 &k = contacts[i].penalty;
        float3 F{k[0] * C[0] + contacts[i].lambda[0], k[1] * C[1] + contacts[i].lambda[1],
                k[2] * C[2] + contacts[i].lambda[2]};

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        const float bounds = std::fabs(F[0]) * friction;
        const float frictionScale = float2{F[1], F[2]}.norm();
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

        // The penalty just changed, so the shared linear Hessian is stale for the next round's
        // stamps. Refreshing it here, where the new penalty is already in hand, is what keeps
        // stampCached from having to build the sandwich itself.
        linHessian[i] = basisHessian(basis, contacts[i].penalty);

        // Hand the next round everything its stamps need. The pose is frozen between here and then
        // (only the combine pass moves a body, and that has already run), so the arms and the force
        // this pass measured are exactly what preCompute would have produced - which is why the
        // assembly pass does not have to run again. The force is the primal's own expression,
        // re-evaluated with the multiplier and the penalty this pass just stored.
        JacobianCache &cache = jacobianCache[i];
        cache.rAWorld = rAWorld;
        cache.rBWorld = rBWorld;
        const float3 &kNext = contacts[i].penalty;
        const float3 &lambdaNext = contacts[i].lambda;
        float3 next{kNext[0] * C[0] + lambdaNext[0], kNext[1] * C[1] + lambdaNext[1],
                kNext[2] * C[2] + lambdaNext[2]};
        next[0] = std::min(next[0], 0.0f);
        const float boundsNext = std::fabs(next[0]) * friction;
        const float scaleNext = float2{next[1], next[2]}.norm();
        if (scaleNext > boundsNext && scaleNext > 0)
        {
            next[1] *= boundsNext / scaleNext;
            next[2] *= boundsNext / scaleNext;
        }
        contactForceCache[i] = next;
    }

    primalReady = true;
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
