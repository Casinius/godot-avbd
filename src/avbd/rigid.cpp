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

namespace avbd {

Rigid::Rigid(Solver *p_solver, float3 p_size, float p_density, float p_friction, float3 p_position, float3 p_velocity)
    : solver(p_solver), forces(0), next(0), positionLin(p_position), positionAng({ 0, 0, 0, 1 }),
    velocityLin(p_velocity), velocityAng({ 0, 0, 0 }), prevVelocityLin(p_velocity), size(p_size),
    friction(p_friction)
{
    // Add to linked list
    next = p_solver->bodies;
    p_solver->bodies = this;

    // Compute mass properties and bounding radius
    mass = p_size.x * p_size.y * p_size.z * p_density;
    moment = float3 {
        (p_size.y * p_size.y + p_size.z * p_size.z) / 12.0f * mass,
        (p_size.x * p_size.x + p_size.z * p_size.z) / 12.0f * mass,
        (p_size.x * p_size.x + p_size.y * p_size.y) / 12.0f * mass
    };
    radius = length(p_size * 0.5f);
}

Rigid::~Rigid()
{
    // Remove from linked list
    Rigid** p = &solver->bodies;
    while (*p != this)
        p = &(*p)->next;
    *p = next;
}

bool Rigid::constrainedTo(Rigid* other) const
{
    // Walk this body's own constraint list, following nextA/nextB. `Force::next` is the
    // solver-wide chain, so following it here would make every broad phase pair test
    // O(every force in the scene) instead of O(this body's constraints) - which, with the
    // pair test running once per candidate pair per step, costs more than the whole rest of
    // the solver put together.
    for (Force* f = forces; f != 0; f = (f->bodyA == this) ? f->nextA : f->nextB)
        if ((f->bodyA == this && f->bodyB == other) || (f->bodyA == other && f->bodyB == this))
            return true;
    return false;
}

} // namespace avbd
