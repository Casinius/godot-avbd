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

int Solver::warmstartForces()
{
    const int forceCount = collectForces();
    forceActive.resize(forceCount);
    pool->forCount(forceCount, [this](int i) {
        forceActive[i] = forceOrder[i]->initialize() ? 1 : 0;
    }, kMinForcesPerDispatch);

    // One pass both reaps and compacts: `forceOrder` is list order, and deleting a force
    // only unlinks it from that list, so the survivors keep the exact order a second
    // collectForces() would produce.
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

    // Dormant manifolds (contact grace) keep their slot and hot state for a few ticks
    // instead of being reaped. initialize only enters dormancy while the grace setting is
    // positive, so with the default 0 this branch never fires.
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

// Compute one body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so `runSubstep` dispatches it as an independent batch. It
// covers static bodies as well: their poses do not move, but the constraint code reads
// their cached `initialLin`/`initialAng`.
void Solver::warmstartBody(int i)
{
    Rigid *body = warmstartOrder[i];

    // Axis locks (Godot BodyAxis): a locked axis cannot carry velocity into the step,
    // and gravity does not act along a locked linear axis. Constraint coupling can
    // still displace a locked axis within a step; documented approximation.
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

// One region body for the whole iteration set: `iterations` rounds of the
// colour-parallel Gauss-Seidel primal followed by a dual pass, with a barrier
// fence at every phase boundary.
//
// The primal runs colour by colour: every body of colour c assembles its full
// block system against the poses the earlier colours already applied, solves it
// and applies the update immediately - Gauss-Seidel across colours, parallel
// and race-free within one (a colour is an independent set of the constraint
// graph by construction). Manifolds take the self-contained updatePrimal path,
// which consumes the previous round's cache for the first touch and recomputes
// once for the second; each manifold is seen by at most one body per colour, so
// nothing races.
//
// The colouring is capped at kMaxColours so the barrier count per round is
// bounded; bodies the greedy colouring could not fit run one block-Jacobi round
// after the colours (precompute, solve against the post-colour snapshot,
// combine), using jacobiRelaxation as the combine control.
//
// NullSync (small scenes) runs the identical body inline in list order, which is
// exactly the threads=1 reference the digest tests compare against.
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
    // Idle detection: a movable, non-sleeping body whose velocities stayed under both
    // thresholds for sleepFrames consecutive steps falls asleep (velocity zeroed,
    // Godot semantics). Only touches this body's fields, so it is safe inside the
    // parallel dispatch and deterministic regardless of thread count.
    // Automatic sleep only for bodies whose host allows it (sleep_mode == SLEEP).
    // A body with can_sleep=false used to enter idle sleep here anyway; the next
    // selection then simulated it again (sleep_mode NEVER), but the finish stage
    // skipped its velocity update because `sleeping` was still set - the body
    // moved with a permanently zero velocity, and contacts could pump it upward
    // indefinitely. That was the pyramid_dbl runaway.
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
// Island-wide sleep bookkeeping. Movable, sleep-enabled bodies are grouped into islands
// through force connectivity (union-find over the flattened body list; static bodies and
// sleep_mode == NEVER bodies join nothing and never propagate). An island where every
// member cleared idle detection this step (stillFrames >= sleepFrames or already
// sleeping) is put to sleep as a unit: Godot semantics, and it removes the
// partially-assembled-chain edge case where the bottom box qualifies a step before the
// boxes above it. Deterministic: the union-find walks the force list in list order and
// the minimum-body-index island representative is independent of visit order.
void Solver::islandSleep()
{
    // Only bodies that can sleep join an island. With the default sleepFrames = 0 this
    // whole pass is skipped by the caller, so the core suite (which never sets
    // sleep_mode) is unaffected.
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

// The primal half of one iteration for a single body (Eqs. 4-6): assemble the 6x6 system
// from every force acting on the body, then solve it and move the body.
//
// Solve the assembled system and move the body by the result (Eq. 4).
void Block::apply(Rigid &body) const noexcept
{
    float3 dxLin, dxAng;
    solveDelta(dxLin, dxAng);
    body.positionLin = body.positionLin + dxLin;
    body.positionAng = body.positionAng + dxAng;
}

// Colour-parallel Gauss-Seidel: bodies sharing a force never solve at the same
// instant, so each colour is an independent set - one barrier per colour carries
// the Gauss-Seidel ordering, and the result is deterministic for a fixed
// colouring no matter the thread count.
//
// The colouring is capped at kMaxColours: a body whose neighbours already used
// every colour goes to the overflow bucket, which the primal runs as one
// block-Jacobi round (every delta solved against the same post-colour snapshot,
// then combined). The bucket costs three barriers when non-empty and is empty in
// every ordinary scene; it exists so a pathological clique cannot make the
// barrier count unbounded.
//
// The wake-up and island-sleep logic is unchanged from the colouring era: a sleeping body
// sharing a force with a moving movable body re-enters the update set, and settled
// islands sleep together.
void Solver::selectUpdateBodies()
{
    warmstartOrder.clear();
    updateOrder.clear();

    // Wake-up pass: a sleeping body (sleep_mode enabled) sharing a force with a moving
    // movable body is woken before selection, so it re-enters this step's update set.
    // A neighbour that is itself below the idle thresholds is not a disturbance: stacks
    // must be allowed to fall asleep together, one body at a time. Serial and read-only
    // over the force adjacency apart from the examined body itself, so the outcome is
    // fixed regardless of thread count.
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

        // Sleep mode support (Godot 4.7 API)
        // SLEEP_MODE_NEVER: always update
        // SLEEP_MODE_SLEEP: update unless idle detection has put it to sleep
        // SLEEP_MODE_START_IN_SLEEP: start in sleep
        if (body->mass > 0)
        {
            // Check if body should be updated based on sleep mode
            bool should_update = true;
            if (body->sleep_mode == 1) // SLEEP_MODE_SLEEP
            {
                // Idle detection (finishBody) is the only way into sleep here;
                // a sleeping body with contacts stays out of the update set until the
                // wake-up pass above re-admits it.
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

    // Overflow-bucket eligibility: a body carrying any non-manifold force (joint, spring,
    // world anchor) is routed into the overflow bucket by colourGraph - the contact map
    // is clamped, hard joints do not tolerate the snapshot solve. Every force in the
    // list is non-manifold; contacts do not count.
    for (Force *f = forces; f != 0; f = f->next)
    {
        if (f->bodyA != 0)
            f->bodyA->contactOnly = false;
        if (f->bodyB != 0)
            f->bodyB->contactOnly = false;
    }

    colourGraph();
}

// Greedy colouring of the update set's constraint graph, in updateOrder order:
// each body takes the smallest colour none of its already-coloured neighbours
// holds, capped at kMaxColours. The order is the body list order, so the colour
// assignment is a pure function of the scene - deterministic regardless of
// thread count. Neighbours are found through the bodies' own force lists; a
// body not yet coloured (later in the order) constrains nothing, and overflow
// neighbours never block a colour because the overflow bucket runs Jacobi
// against a snapshot and is safe to overlap with any colour.
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
