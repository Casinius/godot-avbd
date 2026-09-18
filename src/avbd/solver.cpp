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

#include <algorithm>
#include <numeric>

#include "avbd/job_pool.hpp"
#include "avbd/bvh/bvh.hpp"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

Solver::Solver()
{
}

Solver::~Solver()
{
    clear();
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local)
{
    return pick(origin, dir, local, 0xFFFFFFFFu);
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local, uint32_t p_mask)
{
    const float epsilon = 1.0e-6f;
    float bestT = INFINITY;
    Rigid *bestBody = 0;
    float3 bestLocal = {0, 0, 0};

    // Ray-cast against each OBB by transforming the ray into body local space.
    // Unlike the upstream demo (which uses this to pick bodies for dragging) this
    // also tests static bodies, so rays hit the ground and ramps.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if ((body->collisionLayer & p_mask) == 0)
            continue;
        quat invRot = conjugate(body->positionAng);
        float3 o = rotate(invRot, origin - body->positionLin);
        float3 d = rotate(invRot, dir);
        float3 half = body->size * 0.5f;

        float tEnter = 0.0f;
        float tExit = INFINITY;
        bool hit = true;

        for (int i = 0; i < 3; ++i)
        {
            if (std::fabs(d[i]) < epsilon)
            {
                if (o[i] < -half[i] || o[i] > half[i])
                {
                    hit = false;
                    break;
                }
                continue;
            }

            float invD = 1.0f / d[i];
            float t0 = (-half[i] - o[i]) * invD;
            float t1 = (half[i] - o[i]) * invD;
            if (t0 > t1)
            {
                float tmp = t0;
                t0 = t1;
                t1 = tmp;
            }

            tEnter = max(tEnter, t0);
            tExit = min(tExit, t1);
            if (tEnter > tExit)
            {
                hit = false;
                break;
            }
        }

        if (!hit)
            continue;

        float tHit = tEnter >= 0.0f ? tEnter : tExit;
        if (tHit < 0.0f)
            continue;

        if (tHit < bestT)
        {
            bestT = tHit;
            bestBody = body;
            bestLocal = o + d * tHit;
        }
    }

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
}

void Solver::clear()
{
    while (forces)
        delete forces;

    while (bodies)
        delete bodies;
}

void Solver::step()
{
    ensurePool();

    broadPhase();
    const int forceCount = warmstartForces();
    colourGraph();
    warmstartBodies();
    solveIterations(forceCount);
    finishVelocities();
}

// Contact detection. The naive O(n^2) scan over bounding spheres, as in the reference
// implementation: enough for the low thousands of bodies, and the obvious place to add a
// spatial hash when it is not.
//
// Stays on the calling thread: it appends to the force list, so running it in parallel would
// make the manifold order - and with it the result - depend on scheduling.
void Solver::broadPhase()
{
    for (Rigid *bodyA = bodies; bodyA != 0; bodyA = bodyA->next)
    {
        for (Rigid *bodyB = bodyA->next; bodyB != 0; bodyB = bodyB->next)
        {
            const float3 dp = bodyA->positionLin - bodyB->positionLin;
            const float reach = bodyA->radius + bodyB->radius;
            // Godot semantics: the pair collides when either side's layer is in the other's
            // mask (defaults 1/1 keep every pair, as before these fields existed).
            const bool layersAllow = ((bodyA->collisionLayer & bodyB->collisionMask) != 0) ||
                    ((bodyB->collisionLayer & bodyA->collisionMask) != 0);
            if (dot(dp, dp) <= reach * reach && layersAllow && !bodyA->constrainedTo(bodyB))
                new Manifold(this, bodyA, bodyB);
        }
    }
}

// Bring every force up to date for this step, and drop the ones that have gone inactive.
//
// Initialising is independent work - a force only reads the body poses and writes its own
// warm-started state - so it runs in parallel. Deleting is not: an inactive force (an empty
// manifold, a broken joint) unlinks itself from the solver and from both bodies' lists, so
// that mutation happens in one thread, in list order, once the parallel pass has finished.
//
// Returns the number of forces left, and leaves `forceOrder` indexing exactly those.
int Solver::warmstartForces()
{
    int forceCount = collectForces();
    forceActive.resize(forceCount);
    pool->forCount(forceCount, [this](int i) {
        forceActive[i] = forceOrder[i]->initialize() ? 1 : 0;
    }, kMinForcesPerDispatch);

    for (int i = 0; i < forceCount; i++)
    {
        if (!forceActive[i])
            delete forceOrder[i]; // unlinks itself from the solver and both bodies
    }

    forceCount = collectForces(); // the sweep above invalidated the index
    return forceCount;
}

// Compute every body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so it runs in parallel. It covers static bodies as well: their
// poses do not move, but the constraint code reads their cached `initialLin`/`initialAng`.
void Solver::warmstartBodies()
{
    const int count = static_cast<int>(warmstartOrder.size());
    pool->forCount(count, [this](int i) {
        Rigid *body = warmstartOrder[i];

        // Axis locks (Godot BodyAxis): a locked axis cannot carry velocity into the step,
        // and gravity does not act along a locked linear axis. Constraint coupling can
        // still displace a locked axis within a step; documented approximation.
        if (body->axisLockLinear)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockLinear & (1 << k))
                    body->velocityLin[k] = 0.0f;
            }
        }
        if (body->axisLockAngular)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockAngular & (1 << k))
                    body->velocityAng[k] = 0.0f;
            }
        }

        // Gravity runs along -Z only, so the lock on bit 2 (z) gates it entirely.
        const float g = (body->axisLockLinear & 0x4) ? 0.0f : body->gravity;

        body->inertialLin = body->positionLin + body->velocityLin * dt;
        if (body->mass > 0)
            body->inertialLin += float3{0, 0, g} * (dt * dt);
        body->inertialAng = body->positionAng + body->velocityAng * dt;

        // Adaptive warmstart (See original VBD paper)
        const float3 accel = (body->velocityLin - body->prevVelocityLin) / dt;
        const float accelExt = accel.z * sign(g);
        float accelWeight = clamp(accelExt / abs(g), 0.0f, 1.0f);
        if (!std::isfinite(accelWeight))
            accelWeight = 0.0f;

        // Save the position at the start of the step (x-) and warm-start to the tentative
        // new position (See original VBD paper)
        body->initialLin = body->positionLin;
        body->initialAng = body->positionAng;
        if (body->mass > 0)
        {
            body->positionLin = body->positionLin + body->velocityLin * dt + float3{0, 0, g} * (accelWeight * dt * dt);
            body->positionAng = body->positionAng + body->velocityAng * dt;
        }
    }, kMinBodiesPerDispatch);
}

// The solver's main loop: `iterations` rounds of a primal pass followed by a dual pass.
void Solver::solveIterations(int forceCount)
{
    for (int it = 0; it < iterations; it++)
    {
        // Primal update, group by group. Every body inside a group is independent, so the
        // group is spread over the worker threads; between groups the updates stay ordered,
        // which is the Gauss-Seidel coupling the method relies on.
        for (int colour = 0; colour < colours; colour++)
        {
            const int first = colourStart[colour];
            const int last = colourStart[colour + 1];
            pool->forCount(last - first, [this, first](int i) {
                updatePrimal(updateOrder[first + i]);
            }, kMinBodiesPerDispatch);
        }

        // Dual update: one item per force, and each force only reads the (frozen) body poses
        // and writes its own dual state, so the whole pass is independent.
        pool->forCount(forceCount, [this](int i) {
            forceOrder[i]->updateDual(alpha);
        }, kMinForcesPerDispatch);
    }
}

// BDF1 velocities after the final iteration. Per-body and independent; it covers static
// bodies too, so their bookkeeping matches, though only movable ones get a new velocity.
void Solver::finishVelocities()
{
    const int count = static_cast<int>(warmstartOrder.size());
    pool->forCount(count, [this](int i) {
        Rigid *body = warmstartOrder[i];
        body->prevVelocityLin = body->velocityLin;
        // Don't update velocity if sleeping (Godot 4.7 API)
        if (body->mass > 0 && !body->sleeping)
        {
            body->velocityLin = (body->positionLin - body->initialLin) / dt;
            body->velocityAng = (body->positionAng - body->initialAng) / dt;
        }
        // A locked axis reports no velocity even if constraint coupling nudged it: the
        // lock's promise is "this axis does not move" to whoever reads state back.
        if (body->axisLockLinear)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockLinear & (1 << k))
                    body->velocityLin[k] = 0.0f;
            }
        }
        if (body->axisLockAngular)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (body->axisLockAngular & (1 << k))
                    body->velocityAng[k] = 0.0f;
            }
        }
    }, kMinBodiesPerDispatch);
}

// The primal half of one iteration for a single body (Eqs. 4-6): assemble the 6x6 system
// from every force acting on the body, then solve it and move the body.
//
// Only ever called for movable bodies: `updateOrder` holds exactly those, which is also
// why a static body needs no guard here.
// Solve the assembled system and move the body by the result (Eq. 4).
void Block::apply(Rigid &body) const noexcept
{
    float3 dxLin, dxAng;
    solve(lhsLin, lhsAng, lhsCross, -rhsLin, -rhsAng, dxLin, dxAng);
    body.positionLin = body.positionLin + dxLin;
    body.positionAng = body.positionAng + dxAng;
}

void Solver::updatePrimal(Rigid *body)
{
    Block block(body->mass, body->moment, dt * dt, body->positionLin, body->inertialLin,
            body->positionAng, body->inertialAng);

    for (Force *force = body->forces; force != 0; force = (force->bodyA == body) ? force->nextA : force->nextB)
        force->updatePrimal(body, alpha, block);

    block.apply(*body);
}

// Greedy graph colouring of the constraint graph. A body's colour must differ from every
// body it shares a force with, which is precisely the condition for two bodies to be
// updatable at the same time.
//
// Incremental: a body in this step's update set keeps its previous colour unless a
// neighbour now holds it, in which case it re-runs the first-fit search. Contacts move
// slowly between steps, so the search is rarely paid. Bodies outside the update set
// (static or sleeping) carry `colour == -1`, keeping the invariant that `colour >= 0`
// marks exactly the bodies in `updateOrder` - a stale colour is never trusted.
void Solver::colourGraph()
{
    warmstartOrder.clear();
    updateOrder.clear();

    // Every body is warmed up; only movable ones are updated. Static bodies are excluded
    // from the colouring: they are never moved, so they take no part in the update order
    // and cannot conflict with anything.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        warmstartOrder.push_back(body);
        
        // Sleep mode support (Godot 4.7 API)
        // SLEEP_MODE_NEVER: always update
        // SLEEP_MODE_SLEEP: update if velocity is zero and no forces
        // SLEEP_MODE_START_IN_SLEEP: start in sleep
        if (body->mass > 0)
        {
            // Check if body should be updated based on sleep mode
            bool should_update = true;
            if (body->sleep_mode == 1) // SLEEP_MODE_SLEEP
            {
                // Skip if sleeping and no velocity or forces
                if (body->sleeping && 
                    std::sqrt(body->velocityLin.x*body->velocityLin.x + body->velocityLin.y*body->velocityLin.y + body->velocityLin.z*body->velocityLin.z) == 0.0f && 
                    body->forces == nullptr)
                {
                    should_update = false;
                }
            }
            else if (body->sleep_mode == 2) // SLEEP_MODE_START_IN_SLEEP
            {
                // Start in sleep, wake up if velocity or forces are present
                if (!body->sleeping && 
                    std::sqrt(body->velocityLin.x*body->velocityLin.x + body->velocityLin.y*body->velocityLin.y + body->velocityLin.z*body->velocityLin.z) == 0.0f && 
                    body->forces == nullptr)
                {
                    should_update = false;
                }
                else
                {
                    body->sleeping = false; // Wake up
                }
            }
            
            if (should_update)
            {
                updateOrder.push_back(body);
            }
            else
            {
                // Not in this step's update set: drop any stale colour so the invariant
                // `colour >= 0` <=> "in updateOrder" holds.
                body->colour = -1;
            }
        }
        else
        {
            // Static body: excluded from the update set, so it must not contribute a
            // stale colour.
            body->colour = -1;
        }
    }

    std::vector<int> taken;
    int maxColour = -1;
    for (Rigid *body : updateOrder)
    {
        // Collect the neighbours' current colours.
        taken.clear();
        for (Force *force = body->forces; force != 0; force = (force->bodyA == body) ? force->nextA : force->nextB)
        {
            Rigid *other = (force->bodyA == body) ? force->bodyB : force->bodyA;
            if (other != 0 && other->colour >= 0)
                taken.push_back(other->colour);
        }

        // Incremental recolouring: keep the previous colour when no neighbour took it.
        // A kept colour is conflict-free: neighbours processed earlier in this loop hold
        // their final colours, and a later neighbour that would collide with a kept colour
        // sees it in `taken` and refits. Only bodies whose neighbourhood actually changed
        // pay for the first-fit search.
        int colour = body->colour;
        if (colour < 0 || std::find(taken.begin(), taken.end(), colour) != taken.end())
        {
            colour = 0;
            while (std::find(taken.begin(), taken.end(), colour) != taken.end())
                colour++;
        }

        body->colour = colour;
        if (colour > maxColour)
            maxColour = colour;
    }

    colours = maxColour + 1;

    // Group the bodies by colour, so each iteration walks contiguous, independent runs.
    std::vector<int> perColour(colours, 0);
    for (const Rigid *body : updateOrder)
        perColour[body->colour]++;

    colourStart.assign(colours + 1, 0);
    std::exclusive_scan(perColour.begin(), perColour.end(), colourStart.begin(), 0);
    colourStart[colours] = static_cast<int>(updateOrder.size());

    widest = perColour.empty() ? 0 : *std::max_element(perColour.begin(), perColour.end());

    std::vector<int> cursor = colourStart; // cursor[c] is the next free slot for colour c
    std::vector<Rigid *> grouped(updateOrder.size());
    for (Rigid *body : updateOrder)
        grouped[cursor[body->colour]++] = body;
    updateOrder = std::move(grouped);
}

// Fill `forceOrder` with the solver's forces, so the phases can address them by index.
int Solver::collectForces()
{
    forceOrder.clear();
    for (Force *force = forces; force != 0; force = force->next)
        forceOrder.push_back(force);
    return static_cast<int>(forceOrder.size());
}

void Solver::ensurePool()
{
    if (pool != nullptr)
        return;

    // Pass the value straight through: 0 means "one worker per hardware thread" and 1 means
    // a single worker, which both the pool and `threads` document the same way. (Clamping 0
    // up to 1 here would silently turn the default into a serial run.)
    pool = std::make_unique<detail::JobPool>(static_cast<unsigned>(std::max(0, threads)));
}

int Solver::threadCount() const
{
    if (pool == nullptr)
        return 1;
    return static_cast<int>(pool->threadCount());
}

} // namespace avbd
