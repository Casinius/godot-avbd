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

Force::Force(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB)
    : solver(p_solver), bodyA(p_bodyA), bodyB(p_bodyB), nextA(0), nextB(0)
{
    // Add to the solver's linked list
    next = p_solver->forces;
    p_solver->forces = this;

    // Add to body linked lists
    if (p_bodyA)
    {
        nextA = p_bodyA->forces;
        p_bodyA->forces = this;
    }
    if (p_bodyB)
    {
        nextB = p_bodyB->forces;
        p_bodyB->forces = this;
    }
}


Force::~Force()
{
    // Remove from solver linked list
    Force** p = &solver->forces;
    while (*p != this)
        p = &(*p)->next;
    *p = next;

    // Remove from body linked lists
    if (bodyA)
    {
        p = &bodyA->forces;
        while (*p != this)
            p = (*p)->bodyA == bodyA ? &(*p)->nextA : &(*p)->nextB;
        *p = nextA;
    }

    if (bodyB)
    {
        p = &bodyB->forces;
        while (*p != this)
            p = (*p)->bodyA == bodyB ? &(*p)->nextA : &(*p)->nextB;
        *p = nextB;
    }
}

} // namespace avbd
