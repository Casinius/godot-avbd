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

inline float3x3 geometricStiffnessBallSocket(int k, float3 v)
{
    float3x3 m = diagonal(-v[k], -v[k], -v[k]);

    m[0][k] += v[0];
    m[1][k] += v[1];
    m[2][k] += v[2];

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
    torqueArm = lengthSq((p_bodyA ? p_bodyA->size : float3{0, 0, 0}) + p_bodyB->size);
}

bool Joint::initialize()
{
    // Store constraint function at beginnning of timestep C(x-)
    // Note: if bodyA is null, it is assumed that the joint connects a body to the world space position rA
    C0Lin = (bodyA ? transform(bodyA->positionLin, bodyA->positionAng, rA) : rA) - transform(bodyB->positionLin, bodyB->positionAng, rB);
    C0Ang = ((bodyA ? bodyA->positionAng : quat{ 0, 0, 0, 1 }) - bodyB->positionAng) * torqueArm;

    // Warmstart the dual variables and penalty parameters (Eq. 19)
    // Penalty is safely clamped to a minimum and maximum value
    lambdaLin = lambdaLin * solver->alpha * solver->gamma;
    lambdaAng = lambdaAng * solver->alpha * solver->gamma;
    penaltyLin = clamp(penaltyLin * solver->gamma, PENALTY_MIN, PENALTY_MAX);
    penaltyAng = clamp(penaltyAng * solver->gamma, PENALTY_MIN, PENALTY_MAX);

    // Clamp penalty to material stiffness
    penaltyLin = min(penaltyLin, stiffnessLin);
    penaltyAng = min(penaltyAng, stiffnessAng);

    return !broken;
}

void Joint::updatePrimal(Rigid *body, float alpha, Block &block)
{
    // Linear constraint
    if (lengthSq(penaltyLin) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x, penaltyLin.y, penaltyLin.z);
        float3 C = (bodyA ? transform(bodyA->positionLin, bodyA->positionAng, rA) : rA) - transform(bodyB->positionLin, bodyB->positionAng, rB);
        
        // Stabilization
        if (std::isinf(stiffnessLin))
            C -= C0Lin * alpha;

        // Compute force
        float3 F = K * C + lambdaLin;

        // Choose jacobian depending on input body
        float3x3 jLin = body == bodyA ? float3x3{ 1, 0, 0, 0, 1, 0, 0, 0, 1 } : float3x3{ -1, 0, 0, 0, -1, 0, 0, 0, -1 };
        float3x3 jAng = body == bodyA ? skew(-rotate(bodyA->positionAng, rA)) : skew(rotate(bodyB->positionAng, rB));

        // Stamp into LHS
        float3x3 jLinT = transpose(jLin);
        float3x3 jAngT = transpose(jAng);
        float3x3 jAngTk = jAngT * K;

        block.lhsLin += jLinT * K * jLin;
        block.lhsAng += jAngTk * jAng;
        block.lhsCross += jAngTk * jLin;

        // Diagonal approximation for higher order terms
        float3 r = body == bodyA ? rotate(bodyA->positionAng, rA) : -rotate(bodyB->positionAng, rB);
        float3x3 H = 
            geometricStiffnessBallSocket(0, r) * F[0] +
            geometricStiffnessBallSocket(1, r) * F[1] +
            geometricStiffnessBallSocket(2, r) * F[2];
        block.lhsAng += diagonalize(H);

        // Stamp into RHS
        block.rhsLin += jLinT * F;
        block.rhsAng += jAngT * F;
    }

    // Angular constraint
    if (lengthSq(penaltyAng) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x, penaltyAng.y, penaltyAng.z);
        float3 C = ((bodyA ? bodyA->positionAng : quat{ 0, 0, 0, 1 }) - bodyB->positionAng) * torqueArm;

        // Stabilization
        if (std::isinf(stiffnessAng))
            C -= C0Ang * alpha;

        // Compute force
        float3 F = K * C + lambdaAng;

        // Choose jacobian depending on input body
        float3x3 jAng = (body == bodyA ? float3x3{ 1, 0, 0, 0, 1, 0, 0, 0, 1 } : float3x3{ -1, 0, 0, 0, -1, 0, 0, 0, -1 }) * torqueArm;

        // Stamp into LHS
        block.lhsAng += transpose(jAng) * K * jAng;

        // Stamp into RHS
        block.rhsAng += transpose(jAng) * F;
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
    if (lengthSq(penaltyLin) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyLin.x, penaltyLin.y, penaltyLin.z);
        float3 C = (bodyA ? transform(bodyA->positionLin, bodyA->positionAng, rA) : rA) - transform(bodyB->positionLin, bodyB->positionAng, rB);

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
        penaltyLin = min(penaltyLin + abs(C) * solver->betaLin, min(stiffnessLin, PENALTY_MAX));
    }

    // Angular constraint
    if (lengthSq(penaltyAng) > 0)
    {
        // Compute constraint and jacobians
        float3x3 K = diagonal(penaltyAng.x, penaltyAng.y, penaltyAng.z);
        float3 C = ((bodyA ? bodyA->positionAng : quat{ 0, 0, 0, 1 }) - bodyB->positionAng) * torqueArm;

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
        penaltyAng = min(penaltyAng + abs(C) * solver->betaAng, min(stiffnessAng, PENALTY_MAX));
    }

    // Fracture test
    if (lengthSq(lambdaAng) > fracture * fracture)
    {
        penaltyLin = { 0, 0, 0 };
        penaltyAng = { 0, 0, 0 };
        lambdaLin = { 0, 0, 0 };
        lambdaAng = { 0, 0, 0 };
        broken = true;
    }
}

} // namespace avbd
