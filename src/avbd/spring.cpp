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

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

Spring::Spring(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB,
        float p_stiffness, float p_rest)
    : Force(p_solver, p_bodyA, p_bodyB), rA(p_rA), rB(p_rB), rest(p_rest), stiffness(p_stiffness)
{
    // A negative rest length means "derive it from the bodies' current pose".
    if (rest < 0.0f)
    {
        const float3 pA = p_bodyA->positionLin + p_bodyA->positionAng * rA;
        const float3 pB = p_bodyB->positionLin + p_bodyB->positionAng * rB;
        rest = (pA - pB).norm();
    }
}

// A spring has no stabilised constraint term and no dual variables, so it ignores `alpha`
// and its `updateDual` does nothing: hence the unnamed parameters.
void Spring::updatePrimal(Rigid *body, float /*alpha*/, Block &block)
{
    float3 pA = bodyA->positionLin + bodyA->positionAng * rA;
    float3 pB = bodyB->positionLin + bodyB->positionAng * rB;
    float3 d = pA - pB;
    float dLen = d.norm();
    if (dLen <= 1.0e-6f)
        return;

    float3 n = d / dLen;
    float C = dLen - rest;
    float f = stiffness * C;

    float3 rWorld;
    float3 jLin;
    float3 jAng;
    if (body == bodyA)
    {
        rWorld = bodyA->positionAng * rA;
        jLin = n;
        jAng = rWorld.cross(n);
    }
    else
    {
        rWorld = bodyB->positionAng * rB;
        jLin = -n;
        jAng = -rWorld.cross(n);
    }

    float3 F = jLin * f;
    float3 Tau = jAng * f;
    float3x3 Kll = jLin * jLin.transpose() * stiffness;
    float3x3 Kla = jAng * jLin.transpose() * stiffness;
    float3x3 Kaa = jAng * jAng.transpose() * stiffness;

    block.lhsLin += Kll;
    block.lhsAng += Kaa;
    block.lhsCross += Kla;
    block.rhsLin += F;
    block.rhsAng += Tau;
}

void Spring::updateDual(float /*alpha*/)
{
}

} // namespace avbd
