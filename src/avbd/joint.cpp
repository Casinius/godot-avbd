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

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

namespace {
// The three rows of the anchor error are the three world axes.
const float3 kAxis[3] = {float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)};
} // namespace

inline float3x3 geometricStiffnessBallSocket(int k, float3 v)
{
    float3x3 m = diagonal(-v[k], -v[k], -v[k]);
    m(0, k) += v[0];
    m(1, k) += v[1];
    m(2, k) += v[2];
    return m;
}

Joint::Joint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB,
        float p_stiffnessLin, float p_stiffnessAng, float p_fracture)
    : Force(p_solver, p_bodyA, p_bodyB, ForceKind::Joint), rA(p_rA), rB(p_rB), stiffnessLin(p_stiffnessLin),
    stiffnessAng(p_stiffnessAng), fracture(p_fracture), broken(false)
{
    // Dual variables and penalties start at zero and are built up over the first steps.
    penaltyLin = penaltyAng = float3{0, 0, 0};
    lambdaLin = lambdaAng = float3{0, 0, 0};

    // Scale the angular constraint by the size of the bodies it connects, so that its
    // stiffness is meaningful independently of how big they are.
    torqueArm = ((p_bodyA ? p_bodyA->size : float3{0, 0, 0}) + p_bodyB->size).squaredNorm();
}

bool Joint::initialize()
{
    // Store constraint function at beginnning of timestep C(x-)
    // Note: if bodyA is null, it is assumed that the joint connects a body to the world space position rA
    C0Lin = (bodyA ? bodyA->positionLin + bodyA->positionAng * rA : rA) - (bodyB->positionLin + bodyB->positionAng * rB);
    C0Ang = ((bodyA ? bodyA->positionAng : quat::Identity()) - bodyB->positionAng) * torqueArm;

    // Warmstart the dual variables and penalty parameters (Eq. 19)
    // Penalty is safely clamped to a minimum and maximum value
    lambdaLin = lambdaLin * solver->params.dual.alpha * solver->params.dual.gamma;
    lambdaAng = lambdaAng * solver->params.dual.alpha * solver->params.dual.gamma;
    penaltyLin = float3{
        std::clamp(penaltyLin.x() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyLin.y() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyLin.z() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX)
    };
    penaltyAng = float3{
        std::clamp(penaltyAng.x() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyAng.y() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyAng.z() * solver->params.dual.gamma, PENALTY_MIN, PENALTY_MAX)
    };

    // Clamp penalty to material stiffness
    penaltyLin = float3{
        std::min(penaltyLin.x(), stiffnessLin),
        std::min(penaltyLin.y(), stiffnessLin),
        std::min(penaltyLin.z(), stiffnessLin)
    };
    penaltyAng = float3{
        std::min(penaltyAng.x(), stiffnessAng),
        std::min(penaltyAng.y(), stiffnessAng),
        std::min(penaltyAng.z(), stiffnessAng)
    };

    return !broken;
}

void Joint::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // Linear constraint
    if (penaltyLin.squaredNorm() > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x(), penaltyLin.y(), penaltyLin.z());
        float3 C = (bodyA ? bodyA->positionLin + bodyA->positionAng * rA : rA) - (bodyB->positionLin + bodyB->positionAng * rB);
        
        // Stabilization
        if (std::isinf(stiffnessLin))
            C -= C0Lin * alpha;

        // Compute force
        float3 F = K * C + lambdaLin;

        // Choose jacobian depending on input body
        float3x3 negIdentity(-float3x3::Identity());
        float3x3 jLin = body == bodyA ? float3x3::Identity() : negIdentity;
        // Cross-product (skew-symmetric) matrix of the world-space moment arm.
        float3x3 jAng;
        if (body == bodyA)
        {
            const float3 r = -(bodyA->positionAng * rA);
            jAng <<    0, -r.z(),  r.y(),
                   r.z(),      0, -r.x(),
                  -r.y(),  r.x(),      0;
        }
        else
        {
            const float3 r = bodyB->positionAng * rB;
            jAng <<    0, -r.z(),  r.y(),
                   r.z(),      0, -r.x(),
                  -r.y(),  r.x(),      0;
        }

        // Stamp into LHS
        float3x3 jLinT = jLin.transpose();
        float3x3 jAngT = jAng.transpose();
        float3x3 jAngTk = jAngT * K;

        block.lhsLin += jLinT * K * jLin;
        block.lhsAng += jAngTk * jAng;
        block.lhsCross += jAngTk * jLin;

        // Diagonal approximation for higher order terms
        float3 r = body == bodyA ? bodyA->positionAng * rA : -(bodyB->positionAng * rB);
        float3x3 H = 
            geometricStiffnessBallSocket(0, r) * F[0] +
            geometricStiffnessBallSocket(1, r) * F[1] +
            geometricStiffnessBallSocket(2, r) * F[2];
        block.lhsAng += diagonal(H(0, 0), H(1, 1), H(2, 2));

        // Stamp into RHS
        block.rhsLin += jLinT * F;
        block.rhsAng += jAngT * F;
    }

    // Angular constraint
    if (penaltyAng.squaredNorm() > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x(), penaltyAng.y(), penaltyAng.z());
        float3 C = ((bodyA ? bodyA->positionAng : quat::Identity()) - bodyB->positionAng) * torqueArm;

        // Stabilization
        if (std::isinf(stiffnessAng))
            C -= C0Ang * alpha;

        // Compute force
        float3 F = K * C + lambdaAng;

        // Choose jacobian depending on input body
        float3x3 identity = float3x3::Identity();
        float3x3 negIdentity;
        negIdentity << 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f, 0.0f, -1.0f;
        float3x3 jAng = (body == bodyA ? identity : negIdentity) * torqueArm;

        // Stamp into LHS
        block.lhsAng += jAng.transpose() * K * jAng;

        // Stamp into RHS
        block.rhsAng += jAng.transpose() * F;
    }
}

// The projection's measure pass: the anchor error, row by row, as the multiplier that would take
// each row to zero. Bilateral, unlike a contact: an anchor is wrong in both directions.
//
// Only a hard anchor is projected. A finite stiffnessLin caps its own penalty (updateDual clamps
// to it), so the error left under that cap is the compliance the joint was built with, and taking
// it out by moving a body would turn a spring into a rod. Infinitely stiff is different: there the
// method only *approaches* the anchor, one percent of the remaining error per step because that is
// what the alpha stabilisation asks, and a positional layer is where the rest can be paid for.
//
// Only the linear half is projected. The angular lock of this joint is stamped through a
// Jacobian built on `torqueArm` - a body-size scale factor, with the identity for A and a
// degenerate matrix for B - which is a first-order stand-in rather than a rotation Jacobian, so
// a compliance-weighted split of it would move the bodies by an amount that is not the error.
// The Lagrangian keeps that row; see GenericJoint::projectMeasure for the same decision.
void Joint::projectMeasure(float scale, float maxCorrection)
{
    (void)maxCorrection; // a per-body ceiling, applied where the shares are summed

    // Cleared before every way out of this function: the share pass reads it for whichever body is
    // in the update set, and a left-over multiplier from an earlier step would be applied to a
    // constraint that has since decided it is soft.
    projLin = float3{0, 0, 0};

    // Out of scope, after the clear: the share pass reads projLin and nothing else, so a joint left
    // measuring nothing must leave it holding nothing.
    if (!solver->params.projection.joints)
        return;

    if (!std::isinf(stiffnessLin))
        return;

    const Compliance cA = complianceOf(bodyA);
    const Compliance cB = complianceOf(bodyB);
    if (cA.immobile() && cB.immobile())
        return;

    const float3 anchorA = bodyA ? bodyA->positionLin + bodyA->positionAng * rA : rA;
    const float3 anchorB = bodyB->positionLin + bodyB->positionAng * rB;
    // Measured at the current pose, against the constraint's own zero, unstabilised - see
    // Manifold::projectMeasure for the same choice and why.
    const float3 C = anchorA - anchorB;

    const float3 rAWorld = bodyA ? bodyA->positionAng * rA : float3{0, 0, 0};
    const float3 rBWorld = bodyB->positionAng * rB;
    const float dtSq = solver->dt * solver->dt;
    // Skipped when the layer re-anchors the reference instead: the two are alternative ledgers for
    // the same correction, and running both books the same displacement twice (see the knob).
    const bool feedback = solver->params.projection.dualFeedback &&
            !solver->params.projection.rebase && dtSq > 0.0f;

    float3 mu = float3{0, 0, 0};
    for (int r = 0; r < 3; r++)
    {
        if (C[r] == 0.0f)
            continue;
        const float w = cA.row(kAxis[r], rAWorld) + cB.row(kAxis[r], rBWorld);
        if (w > 0.0f)
            mu[r] = -C[r] * scale / w;
    }

    projLin = mu;

#ifdef AVBD_DIAG
    atomicAddFloat(solver->diag.askedJoint, std::fabs(mu[0]) + std::fabs(mu[1]) + std::fabs(mu[2]));
#endif

    if (feedback)
    {
#ifdef AVBD_DIAG
        atomicAddFloat(solver->diag.chargeJoint,
                (std::fabs(mu[0]) + std::fabs(mu[1]) + std::fabs(mu[2])) / dtSq);
#endif
        lambdaLin = lambdaLin - mu / dtSq;
    }
}

// Adopt the pose the pass delivered as the anchor's reference, leaving lambdaLin alone: the anchor
// error is stabilised against C0Lin (SolverParams::Projection::rebase), so moving C0Lin is what
// makes the next round accept the projected pose, while the multiplier keeps carrying the load the
// chain was already transferring through this joint.
void Joint::rebaseStabilisation()
{
    if (!solver->params.projection.rebase || projLin.isZero(0.0f))
        return;

    C0Lin = (bodyA ? bodyA->positionLin + bodyA->positionAng * rA : rA) -
            (bodyB->positionLin + bodyB->positionAng * rB);
}

void Joint::projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const
{
    if (projLin.isZero(0.0f))
        return;

    const bool isA = body == bodyA;
    if (!isA && body != bodyB)
        return;

    const Compliance c = complianceOf(isA ? bodyA : bodyB);
    if (c.immobile())
        return;

    const float sign = isA ? 1.0f : -1.0f;
    const float3 arm = isA ? (bodyA ? bodyA->positionAng * rA : float3{0, 0, 0})
                           : bodyB->positionAng * rB;
    for (int r = 0; r < 3; r++)
    {
        if (projLin[r] == 0.0f)
            continue;
        r_deltaLin = r_deltaLin + kAxis[r] * (sign * projLin[r] * c.invMass);
        r_deltaAng = r_deltaAng + c.invInertiaTimes(arm.cross(kAxis[r])) * (sign * projLin[r]);
    }
}

void Joint::breakNow()
{
    penaltyLin = { 0, 0, 0 };
    penaltyAng = { 0, 0, 0 };
    lambdaLin = { 0, 0, 0 };
    lambdaAng = { 0, 0, 0 };
    broken = true;
}

void Joint::updateDual(float alpha)
{
    // Linear constraint
    if (penaltyLin.squaredNorm() > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x(), penaltyLin.y(), penaltyLin.z());
        float3 C = (bodyA ? bodyA->positionLin + bodyA->positionAng * rA : rA) - (bodyB->positionLin + bodyB->positionAng * rB);

        if (std::isinf(stiffnessLin))
        {
            // Stabilization
            C -= C0Lin * alpha;

            // Compute force
            float3 F = K * C + lambdaLin;

            // Store updated force
            lambdaLin = F;
        }

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        penaltyLin = (penaltyLin + C.cwiseAbs() * solver->params.dual.betaLin).cwiseMin(float3{stiffnessLin, stiffnessLin, stiffnessLin}.cwiseMin(float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX}));

        // Angular
        penaltyAng = (penaltyAng + C.cwiseAbs() * solver->params.dual.betaAng).cwiseMin(float3{stiffnessAng, stiffnessAng, stiffnessAng}.cwiseMin(float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX}));
    }
    // Angular constraint
    if (penaltyAng.squaredNorm() > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x(), penaltyAng.y(), penaltyAng.z());
        float3 C = ((bodyA ? bodyA->positionAng : quat::Identity()) - bodyB->positionAng) * torqueArm;

        if (std::isinf(stiffnessAng))
        {
            // Stabilization
            C -= C0Ang * alpha;

            // Compute force
            float3 F = K * C + lambdaAng;

            // Store updated force
            lambdaAng = F;
        }

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        penaltyAng = (penaltyAng + C.cwiseAbs() * solver->params.dual.betaAng).cwiseMin(float3{stiffnessAng, stiffnessAng, stiffnessAng}.cwiseMin(float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX}));
    }

    // Fracture test
    if (lambdaAng.squaredNorm() > fracture * fracture)
    {
        penaltyLin = { 0, 0, 0 };
        penaltyAng = { 0, 0, 0 };
        lambdaLin = { 0, 0, 0 };
        lambdaAng = { 0, 0, 0 };
        broken = true;
    }

#ifdef AVBD_DIAG
    // The two numbers that say whether the Lagrangian and the projection are aiming at the same
    // pose: the force this joint is carrying, and the reference the stabilisation holds it to. A
    // lambda that grew while |C0| did not shrink is an injection, not a pose correction.
    atomicMaxFloat(solver->diag.maxJointC0, C0Lin.norm());
    atomicMaxFloat(solver->diag.maxJointLambda, lambdaLin.norm() + lambdaAng.norm());
#endif
}

} // namespace avbd
