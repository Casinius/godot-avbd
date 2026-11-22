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

#include <climits>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "avbd/solver.h"

#include "avbd/job_pool.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "radixsort.hpp"
#include "avbd/solver_internal.hpp"

namespace avbd {

// Warm every force for the step and drop inactive ones; returns the non-manifold force count.
int Solver::warmstartForces()
{
    const int forceCount = collectForces();
    forceActive.resize(forceCount);
    pool->forCount(forceCount, [this](int i) {
        forceActive[i] = forceOrder[i]->initialize() ? 1 : 0;
    }, kMinForcesPerDispatch);

    // Reap and compact in one pass: forceOrder is list order and deletion only unlinks, so
    // survivors keep the exact order a second collectForces() would produce.
    int write = 0;
    for (int i = 0; i < forceCount; i++)
    {
        Force *f = forceOrder[i];
        if (!forceActive[i])
        {
            delete f; // unlinks itself from the solver and both bodies
            continue;
        }
        forceOrder[write++] = f;
    }
    const int liveForceCount = write;

    const int manifoldCount = manifolds.count();
    manifoldActive.resize(static_cast<size_t>(manifoldCount));
    pool->forCount(manifoldCount, [this](int i) {
        const uint32_t slot = manifolds.live[i];
        // Frozen constraints keep their state verbatim; nothing they read can move.
        manifoldActive[i] = manifolds.frozen(slot) || manifolds.initialize(slot, *this) ? 1 : 0;
        if (traceMatch(manifolds.bodyA[slot], manifolds.bodyB[slot]))
            std::printf("TRACE init slot=%u active=%d contacts=%d frozen=%d\n", slot,
                    (int)manifoldActive[i], manifolds.contactCount(slot),
                    (int)manifolds.frozen(slot));
    }, kMinForcesPerDispatch);

    // Contact-grace dormancy keeps a manifold's slot and hot state; grace defaults to 0, so
    // this reaping exception never fires with default settings.
    int liveWrite = 0;
    for (int i = 0; i < manifoldCount; i++)
    {
        const uint32_t slot = manifolds.live[i];
        if (!manifoldActive[i])
        {
            if (contactGraceTicks > 0 && manifolds.dormantTicks[slot] > 0)
            {
                manifolds.live[liveWrite++] = slot;
                continue;
            }
#ifdef AVBD_DIAG
            diag.reapedManifolds++;
#endif
            manifolds.release(slot);
            continue;
        }
        manifolds.live[liveWrite++] = slot;
    }
    manifolds.live.resize(static_cast<size_t>(liveWrite));

#ifdef AVBD_DIAG
    {
        int soft = 0;
        const float softLimit = 0.1f * betaLin;
        for (uint32_t slot : manifolds.live)
            for (int ci = 0, n = manifolds.contactCount(slot); ci < n; ci++)
                if (manifolds.contactPenalty[ManifoldStore::base(slot) + ci].x < softLimit)
                    soft++;
        diag.softContacts = soft;
    }
#endif

    // The live manifold set and force list are final for the step; the colouring, primal,
    // wake-up and island passes read this adjacency.
    buildForceAdjacency();
    return liveForceCount;
}

// Inertial position and adaptive warm start (VBD paper Eq. 2), dispatched as an independent
// batch; static bodies are included because constraint code reads their cached initial pose.
void Solver::warmstartBody(int i)
{
    Rigid *body = warmstartOrder[i];

    // Godot BodyAxis locks zero velocity on locked axes; coupling can still displace them within
    // a step - the documented approximation.
    maskLockedLinear(body->axisLockLinear, body->velocityLin);
    maskLockedAngular(body->axisLockAngular, body->velocityAng);

    // Gravity runs along -Z only, so the lock on bit 2 (z) gates it entirely.
    const float g = (body->axisLockLinear & 0x4) ? 0.0f : body->gravity;

    body->inertialLin = body->positionLin + body->velocityLin * dt;
    if (body->mass > 0)
        body->inertialLin += float3{0, 0, g} * (dt * dt);
    body->inertialAng = body->positionAng + body->velocityAng * dt;

    // Adaptive warmstart (See original VBD paper)
    const float3 accel = (body->velocityLin - body->prevVelocityLin) / dt;
    const float accelExt = accel.z * std::copysign(1.0f, g);
    float accelWeight = std::clamp(accelExt / std::fabs(g), 0.0f, 1.0f);
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
}

// One body's BDF1 velocity after the final iteration. Independent per body, so it runs as a
// parallel batch; static bodies are included only for bookkeeping.
void Solver::finishBody(int i, bool idleDetection)
{
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
    maskLockedLinear(body->axisLockLinear, body->velocityLin);
    maskLockedAngular(body->axisLockAngular, body->velocityAng);
    // Idle sleep needs sleep_mode == SLEEP; it touches only per-body fields, so it is thread-count
    // independent. Dropping the mode gate caused the pyramid_dbl zero-velocity runaway.
    if (idleDetection && sleepFrames > 0 && body->mass > 0 && !body->sleeping && body->sleep_mode == 1)
    {
        if (length2(body->velocityLin) < sleepLinearThreshold * sleepLinearThreshold &&
                length2(body->velocityAng) < sleepAngularThreshold * sleepAngularThreshold)
            ++body->stillFrames;
        else
            body->stillFrames = 0;
        if (body->stillFrames >= sleepFrames)
        {
            body->sleeping = true;
            body->velocityLin = float3{0, 0, 0};
            body->velocityAng = float3{0, 0, 0};
        }
    }
}
// Force-connected sleep-enabled bodies sleep as one island (fixes the bottom-of-chain edge case);
// static/NEVER bodies join none. Union-find in list order + min-index rep make it order-independent.
void Solver::islandSleep()
{
    // Only sleep-enabled bodies join; sleepFrames = 0 makes the caller skip this pass entirely.
    islandBodies.clear();
    for (Rigid *body = bodies; body != 0; body = body->next)
        if (body->mass > 0 && body->sleep_mode == 1)
            islandBodies.push_back(body);

    const int count = static_cast<int>(islandBodies.size());
    if (count < 2)
        return;

    std::vector<int> parent(count);
    for (int i = 0; i < count; ++i)
        parent[i] = i;
    auto find = [&parent](int x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };

    // Body -> island slot through the broad-phase index: one array write per island body,
    // then every constraint is two O(1) lookups (no pointer sort and binary search).
    bodyIsland.assign(bodiesInOrder.size(), -1);
    for (int i = 0; i < count; ++i)
        bodyIsland[static_cast<size_t>(islandBodies[i]->broadIndex)] = i;

    const auto join = [&](Rigid *a, Rigid *b) {
        if (a == nullptr || b == nullptr)
            return;
        if (a->broadIndex < 0 || b->broadIndex < 0)
            return;
        const int ia = bodyIsland[static_cast<size_t>(a->broadIndex)];
        const int ib = bodyIsland[static_cast<size_t>(b->broadIndex)];
        if (ia >= 0 && ib >= 0)
            parent[find(ia)] = find(ib);
    };
    for (uint32_t slot : manifolds.live)
        join(manifolds.bodyA[slot], manifolds.bodyB[slot]);
    for (Force *f = forces; f != 0; f = f->next)
        join(f->bodyA, f->bodyB);

    // Compose islands: representative -> aggregate state.
    // islandAwake[r] = true iff any member is awake; islandStill[r] = min stillFrames.
    std::vector<char> islandAwake(count, 0);
    std::vector<int> islandStill(count, INT_MAX);
    
    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (!islandBodies[i]->sleeping)
            islandAwake[r] = 1;
        islandStill[r] = std::min(islandStill[r], islandBodies[i]->stillFrames);
    }

    // Sleep whole islands whose every member has cleared idle detection; leave mixed or
    // awake islands exactly as they are (the per-body wake-up pass above already ran).
    
    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (islandAwake[r] || islandStill[r] < sleepFrames)
            continue;
        Rigid *body = islandBodies[i];
        if (!body->sleeping)
        {
            body->sleeping = true;
            body->velocityLin = float3{0, 0, 0};
            body->velocityAng = float3{0, 0, 0};
        }
    }
}

// Solve the assembled 6x6 block and move the body by the result (Eq. 4); this function does not
// assemble - each force stamps its own block in the primal pass.
void Block::apply(Rigid &body) const noexcept
{
    float3 dxLin, dxAng;
    solveDelta(dxLin, dxAng);
    body.positionLin = body.positionLin + dxLin;
    body.positionAng = body.positionAng + dxAng;
}

// The step's movable, awake bodies, in deterministic list order; everything else is skipped.
void Solver::selectUpdateBodies()
{
    warmstartOrder.clear();
    updateOrder.clear();

    // A sleeping body sharing a force with a moving neighbour re-enters the update set; a
    // below-threshold neighbour is no disturbance (stacks sleep together). Serial: thread-count independent.
    if (sleepFrames > 0)
    {
        const float lin2 = sleepLinearThreshold * sleepLinearThreshold;
        const float ang2 = sleepAngularThreshold * sleepAngularThreshold;
        for (Rigid *body = bodies; body != 0; body = body->next)
        {
            if (!body->sleeping || body->sleep_mode != 1)
                continue;
            for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
            {
                Rigid *other = otherBody(ref, body, manifolds);
                if (other != 0 && other->mass > 0 && !other->sleeping &&
                        (length2(other->velocityLin) > lin2 ||
                                length2(other->velocityAng) > ang2))
                {
                    body->sleeping = false;
                    body->stillFrames = 0;
                    break;
                }
            }
        }
        islandSleep();
    }

    // Every body is warmed up; only movable ones are updated. Static bodies are excluded
    // from the update set: they are never moved, so nothing solves for them.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        warmstartOrder.push_back(body);

        // Fresh per-step state; admitted bodies flip inUpdateSet below.
        body->inUpdateSet = false;
        body->updateIndex = -1;
        body->contactOnly = true;

        // Godot 4.7 sleep modes: NEVER always updates, SLEEP skips when idle-sleeping,
        // START_IN_SLEEP starts asleep and wakes on motion or forces.
        if (body->mass > 0)
        {
            // Check if body should be updated based on sleep mode
            bool should_update = true;
            if (body->sleep_mode == 1) // SLEEP_MODE_SLEEP
            {
                // finishBody's idle detection is the only way in; a sleeping body stays out
                // until the wake-up pass above re-admits it.
                if (body->sleeping)
                    should_update = false;
            }
            else if (body->sleep_mode == 2) // SLEEP_MODE_START_IN_SLEEP
            {
                // Start in sleep, wake up if velocity or forces are present
                if (!body->sleeping && 
                    length2(body->velocityLin) == 0.0f &&
                    forceAdjacency.of(body->broadIndex).empty())
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
                body->updateIndex = static_cast<int>(updateOrder.size());
                updateOrder.push_back(body);
                body->inUpdateSet = true;
            }
        }
        else
        {
            // Static body: excluded from the update set, so nothing solves for it.
        }
    }

    // Per-body delta buffers for the Jacobi primal's split-solve / combine pair.
    bodyDeltaLin.resize(updateOrder.size());
    bodyDeltaAng.resize(updateOrder.size());

    // Any non-manifold force (joint, spring, anchor) routes a body to colourGraph's overflow
    // bucket: hard joints cannot tolerate the clamped contact snapshot solve.
    for (Force *f = forces; f != 0; f = f->next)
    {
        if (f->bodyA != 0)
            f->bodyA->contactOnly = false;
        if (f->bodyB != 0)
            f->bodyB->contactOnly = false;
    }

    colourGraph();
}

// Independent sets of the constraint graph, capped at kMaxColours; the rest go to the
// block-Jacobi overflow bucket. Colour order is the Gauss-Seidel order and is deterministic.
void Solver::colourGraph()
{
    const int n = static_cast<int>(updateOrder.size());
    bodyColour.assign(n, kOverflowColour);
    colourOrder.clear();
    overflowOrder.clear();

    uint8_t used[kMaxColours];
    for (int i = 0; i < n; i++)
    {
        Rigid *body = updateOrder[i];
        for (int c = 0; c < kMaxColours; c++)
            used[c] = 0;

        for (const ForceRef &ref : forceAdjacency.of(body->broadIndex))
        {
            Rigid *other = otherBody(ref, body, manifolds);
            if (other == 0 || other->updateIndex < 0 || other->updateIndex >= i)
                continue; // not in the update set, or coloured later in the pass
            const uint8_t c = bodyColour[other->updateIndex];
            if (c < kMaxColours)
                used[c] = 1;
        }

        uint8_t colour = kOverflowColour;
        for (int c = 0; c < kMaxColours; c++)
            if (!used[c])
            {
                colour = static_cast<uint8_t>(c);
                break;
            }
        bodyColour[i] = colour;
    }

    // Counting sort into per-colour runs (stable in updateOrder order).
    colourBegins.assign(kMaxColours + 1, 0);
    for (int i = 0; i < n; i++) {
        if (contactOnlyOverflow && updateOrder[i]->contactOnly) {
            bodyColour[i] = kOverflowColour;
            overflowOrder.push_back(i);
        }
        if (bodyColour[i] != kOverflowColour) {
            colourBegins[bodyColour[i] + 1]++;
        }
    }
    for (int c = 0; c < kMaxColours; c++)
        colourBegins[c + 1] += colourBegins[c];
    colourOrder.resize(n - static_cast<int>(overflowOrder.size()));
    std::array<int, kMaxColours> cursor{};
    for (int c = 0; c < kMaxColours; c++)
        cursor[c] = colourBegins[c];
    for (int i = 0; i < n; i++)
        if (bodyColour[i] < kMaxColours)
            colourOrder[cursor[bodyColour[i]]++] = i;

    colourCount = 0;
    for (int c = 0; c < kMaxColours; c++)
        if (colourBegins[c + 1] > colourBegins[c])
            colourCount = c + 1;
}

// Fill `forceOrder` with the solver's non-manifold forces, so the phases can address them
// by index. Manifolds live in the store and are addressed by slot.
int Solver::collectForces()
{
    forceOrder.clear();
    forceOrder.insert(forceOrder.end(), next_range(forces).begin(), next_range(forces).end());
    return static_cast<int>(forceOrder.size());
}

void Solver::buildForceAdjacency()
{
    forceAdjacency.build(bodiesInOrder, forces, manifolds);
}


} // namespace avbd
