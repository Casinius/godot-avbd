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

namespace avbd {

inline float3x3 geometricStiffnessBallSocket(int k, float3 v)
{
    float3x3 m = diagonal(-v[k], -v[k], -v[k]);
    m[k][0] += v[0];
    m[k][1] += v[1];
    m[k][2] += v[2];
    return m;
}

Joint::Joint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB,
        float p_stiffnessLin, float p_stiffnessAng, float p_fracture)
    : Force(p_solver, p_bodyA, p_bodyB), rA(p_rA), rB(p_rB), stiffnessLin(p_stiffnessLin),
    stiffnessAng(p_stiffnessAng), fracture(p_fracture), broken(false)
{
    // Dual variables and penalties start at zero and are built up over the first steps.
    penaltyLin = penaltyAng = float3{0, 0, 0};
    lambdaLin = lambdaAng = float3{0, 0, 0};

    // Scale the angular constraint by the size of the bodies it connects, so that its
    // stiffness is meaningful independently of how big they are.
    torqueArm = length2((p_bodyA ? p_bodyA->size : float3{0, 0, 0}) + p_bodyB->size);
}

bool Joint::initialize()
{
    // Store constraint function at beginnning of timestep C(x-)
    // Note: if bodyA is null, it is assumed that the joint connects a body to the world space position rA
    C0Lin = (bodyA ? bodyA->positionLin + qrot(bodyA->positionAng, rA) : rA) - (bodyB->positionLin + qrot(bodyB->positionAng, rB));
    C0Ang = ((bodyA ? bodyA->positionAng : quatIdentity) - bodyB->positionAng) * torqueArm;

    // Warmstart the dual variables and penalty parameters (Eq. 19)
    // Penalty is safely clamped to a minimum and maximum value
    lambdaLin = lambdaLin * solver->alpha * solver->gamma;
    lambdaAng = lambdaAng * solver->alpha * solver->gamma;
    penaltyLin = float3{
        std::clamp(penaltyLin.x * solver->gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyLin.y * solver->gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyLin.z * solver->gamma, PENALTY_MIN, PENALTY_MAX)
    };
    penaltyAng = float3{
        std::clamp(penaltyAng.x * solver->gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyAng.y * solver->gamma, PENALTY_MIN, PENALTY_MAX),
        std::clamp(penaltyAng.z * solver->gamma, PENALTY_MIN, PENALTY_MAX)
    };

    // Clamp penalty to material stiffness
    penaltyLin = float3{
        std::min(penaltyLin.x, stiffnessLin),
        std::min(penaltyLin.y, stiffnessLin),
        std::min(penaltyLin.z, stiffnessLin)
    };
    penaltyAng = float3{
        std::min(penaltyAng.x, stiffnessAng),
        std::min(penaltyAng.y, stiffnessAng),
        std::min(penaltyAng.z, stiffnessAng)
    };

    return !broken;
}

void Joint::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // Linear constraint
    if (length2(penaltyLin) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x, penaltyLin.y, penaltyLin.z);
        float3 C = (bodyA ? bodyA->positionLin + qrot(bodyA->positionAng, rA) : rA) - (bodyB->positionLin + qrot(bodyB->positionAng, rB));
        
        // Stabilization
        if (std::isinf(stiffnessLin))
            C -= C0Lin * alpha;

        // Compute force
        float3 F = mul(K, C) + lambdaLin;

        // Choose jacobian depending on input body
        const float3x3 negIdentity = -identity3x3;
        float3x3 jLin = body == bodyA ? identity3x3 : negIdentity;
        // Cross-product (skew-symmetric) matrix of the world-space moment arm: [r]x with
        // [r]x * v = r x v. Columns are the images of the unit axes.
        float3x3 jAng;
        if (body == bodyA)
        {
            const float3 r = -qrot(bodyA->positionAng, rA);
            jAng = float3x3{{0, r.z, -r.y}, {-r.z, 0, r.x}, {r.y, -r.x, 0}};
        }
        else
        {
            const float3 r = qrot(bodyB->positionAng, rB);
            jAng = float3x3{{0, r.z, -r.y}, {-r.z, 0, r.x}, {r.y, -r.x, 0}};
        }

        // Stamp into LHS
        const float3x3 jLinT = transpose(jLin);
        const float3x3 jAngT = transpose(jAng);
        const float3x3 jAngTk = mul(jAngT, K);

        block.lhsLin += mul(jLinT, K, jLin);
        block.lhsAng += mul(jAngTk, jAng);
        block.lhsCross += mul(jAngTk, jLin);

        // Diagonal approximation for higher order terms
        float3 r = body == bodyA ? qrot(bodyA->positionAng, rA) : -qrot(bodyB->positionAng, rB);
        float3x3 H = 
            geometricStiffnessBallSocket(0, r) * F[0] +
            geometricStiffnessBallSocket(1, r) * F[1] +
            geometricStiffnessBallSocket(2, r) * F[2];
        block.lhsAng += diagonal(H[0][0], H[1][1], H[2][2]);

        // Stamp into RHS
        block.rhsLin += mul(jLinT, F);
        block.rhsAng += mul(jAngT, F);
    }

    // Angular constraint
    if (length2(penaltyAng) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x, penaltyAng.y, penaltyAng.z);
        float3 C = ((bodyA ? bodyA->positionAng : quatIdentity) - bodyB->positionAng) * torqueArm;

        // Stabilization
        if (std::isinf(stiffnessAng))
            C -= C0Ang * alpha;

        // Compute force
        float3 F = mul(K, C) + lambdaAng;

        // Choose jacobian depending on input body. The reference implementation's
        // body-B half is diag(0, -1, -1) (the x row is carried by A's identity half);
        // changing it to diag(1, -1, -1) doubles the x correction and tears stiff
        // joint lattices apart, so keep the reference matrix.
        //
        // Known limitation: that split assumes A exists. With bodyA == nullptr (a joint to
        // the world) nothing carries the x row, so the x component of the angular constraint
        // produces no force while its multiplier keeps ramping. Use GenericJoint for
        // world-anchored angular constraints - it reads the world frame as body A's frame
        // and constrains all three axes. This class is kept for the core tests' joints.
        const float3x3 jAng = (body == bodyA ? identity3x3
                : float3x3{float3{0, 0, 0}, float3{0, -1, 0}, float3{0, 0, -1}}) * torqueArm;

        // Stamp into LHS
        block.lhsAng += mul(transpose(jAng), K, jAng);

        // Stamp into RHS
        block.rhsAng += mul(transpose(jAng), F);
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
    if (length2(penaltyLin) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x, penaltyLin.y, penaltyLin.z);
        float3 C = (bodyA ? bodyA->positionLin + qrot(bodyA->positionAng, rA) : rA) - (bodyB->positionLin + qrot(bodyB->positionAng, rB));

        if (std::isinf(stiffnessLin))
        {
            // Stabilization
            C -= C0Lin * alpha;

            // Compute force
            float3 F = mul(K, C) + lambdaLin;

            // Store updated force
            lambdaLin = F;
        }

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        penaltyLin = min(penaltyLin + abs(C) * solver->betaLin, min(float3{stiffnessLin, stiffnessLin, stiffnessLin}, float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX}));

        // The reference implementation also ramps penaltyAng here, with *this* block's C -
        // a length in metres fed to the angular beta, and then the angular block below ramps
        // it again with the real angular error, so a joint with both halves enabled grew its
        // angular penalty at up to twice the rate and used the wrong units for the first
        // half. Only the angular block may touch penaltyAng.
    }
    // Angular constraint
    if (length2(penaltyAng) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x, penaltyAng.y, penaltyAng.z);
        float3 C = ((bodyA ? bodyA->positionAng : quatIdentity) - bodyB->positionAng) * torqueArm;

        if (std::isinf(stiffnessAng))
        {
            // Stabilization
            C -= C0Ang * alpha;

            // Compute force
            float3 F = mul(K, C) + lambdaAng;

            // Store updated force
            lambdaAng = F;
        }

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        penaltyAng = min(penaltyAng + abs(C) * solver->betaAng, min(float3{stiffnessAng, stiffnessAng, stiffnessAng}, float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX}));
    }

    // Fracture test
    if (length2(lambdaAng) > fracture * fracture)
    {
        penaltyLin = { 0, 0, 0 };
        penaltyAng = { 0, 0, 0 };
        lambdaLin = { 0, 0, 0 };
        lambdaAng = { 0, 0, 0 };
        broken = true;
    }
}

} // namespace avbd
