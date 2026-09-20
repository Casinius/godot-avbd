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

Manifold::Manifold(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB)
    : Force(p_solver, p_bodyA, p_bodyB), numContacts(0)
{
}

bool Manifold::initialize()
{
    // Compute friction
    friction = std::sqrt(bodyA->friction * bodyB->friction);

    // Compute new contacts
    Contact newContacts[8] = {};
    int newNumContacts = collide(bodyA, bodyB, newContacts, basis);

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
        float3 xA = transform(bodyA->positionLin, bodyA->positionAng, contacts[i].rA);
        float3 xB = transform(bodyB->positionLin, bodyB->positionAng, contacts[i].rB);
        contacts[i].C0 = basis * (xA - xB) + float3{COLLISION_MARGIN, 0, 0};

        // Clamp penetration to avoid excessive penalty forces
        float penetration = length(contacts[i].C0);
        if (penetration > 0.02f) // Cap at 2cm penetration
        {
            contacts[i].C0 = normalize(contacts[i].C0) * 0.02f;
            penetration = 0.02f;
        }

        // Adaptive penalty based on initial penetration depth
        // Use log scale for large penetration variations, but avoid division by zero
        float basePenalty = 1.0f; // PENALTY_MIN
        float logPenetration = std::log10(std::max(penetration, 0.001f));
        float adaptivePenalty = basePenalty * std::pow(10.0f, logPenetration * 2.0f);
        float clampedPenalty = clamp(adaptivePenalty, PENALTY_MIN, PENALTY_MAX);
        contacts[i].penalty = float3{clampedPenalty, clampedPenalty, clampedPenalty};

        // Warmstart the dual variables and penalty parameters (Eq. 19)
        // Penalty is safely clamped to a minimum and maximum value
        contacts[i].lambda = contacts[i].lambda * solver->alpha * solver->gamma;
        float clampedLambdaPenalty = clamp(clampedPenalty * solver->gamma, PENALTY_MIN, PENALTY_MAX);
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

    for (int i = 0; i < numContacts; i++)
    {
        float3 rAWorld = rotate(bodyA->positionAng, contacts[i].rA);
        float3 rBWorld = rotate(bodyB->positionAng, contacts[i].rB);

        // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
        float3x3 jALin = basis;
        float3x3 jBLin = -basis;
        float3x3 jAAng = float3x3{cross(rAWorld, jALin[0]), cross(rAWorld, jALin[1]), cross(rAWorld, jALin[2])};
        float3x3 jBAng = float3x3{cross(rBWorld, jBLin[0]), cross(rBWorld, jBLin[1]), cross(rBWorld, jBLin[2])};

        float3x3 K = diagonal(contacts[i].penalty.x, contacts[i].penalty.y, contacts[i].penalty.z);
        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force
        float3 F = K * C + contacts[i].lambda;

        // Clamp normal force
        F[0] = min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction;
        float frictionScale = length(float2{F[1], F[2]});
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        // Choose jacobian depending on input body
        float3x3 jLin = body == bodyA ? jALin : jBLin;
        float3x3 jAng = body == bodyA ? jAAng : jBAng;

        // Stamp into LHS
        float3x3 jLinT = transpose(jLin);
        float3x3 jAngT = transpose(jAng);
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

    for (int i = 0; i < numContacts; i++)
    {
        float3 rAWorld = rotate(bodyA->positionAng, contacts[i].rA);
        float3 rBWorld = rotate(bodyB->positionAng, contacts[i].rB);

        // Compute the Taylor series approximation of the constraint function C(x) (Sec 4)
        float3x3 jALin = basis;
        float3x3 jBLin = -basis;
        float3x3 jAAng = float3x3{cross(rAWorld, jALin[0]), cross(rAWorld, jALin[1]), cross(rAWorld, jALin[2])};
        float3x3 jBAng = float3x3{cross(rBWorld, jBLin[0]), cross(rBWorld, jBLin[1]), cross(rBWorld, jBLin[2])};

        float3x3 K = diagonal(contacts[i].penalty.x, contacts[i].penalty.y, contacts[i].penalty.z);
        float3 C = contacts[i].C0 * (1 - alpha) + jALin * dqALin + jBLin * dqBLin + jAAng * dqAAng + jBAng * dqBAng;

        // Compute force
        float3 F = K * C + contacts[i].lambda;

        // Clamp normal force
        F[0] = min(F[0], 0.0f);

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
        if (F[0] < 0)
            contacts[i].penalty[0] = min(contacts[i].penalty[0] + solver->betaLin * std::fabs(C[0]), PENALTY_MAX);
        if (frictionScale <= bounds)
        {
            contacts[i].penalty[1] = min(contacts[i].penalty[1] + solver->betaLin * std::fabs(C[1]), PENALTY_MAX);
            contacts[i].penalty[2] = min(contacts[i].penalty[2] + solver->betaLin * std::fabs(C[2]), PENALTY_MAX);
            contacts[i].stick = length(float2{C[1], C[2]}) < STICK_THRESH;
        }
    }
}

} // namespace avbd
