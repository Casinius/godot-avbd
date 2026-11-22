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

#pragma once

#include "avbd/solver.h"

namespace avbd {

// Godot BodyAxis bits: one helper for every place a locked component must read back as
// motionless (the positional projection, the warm start and the final velocity). Keeping the
// bit tests here means the three sites cannot silently disagree about which axis is locked.
inline void maskLockedLinear(uint8_t locks, float3 &v) noexcept
{
    if (locks)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (locks & (1 << k))
                v[k] = 0.0f;
        }
    }
}

inline void maskLockedAngular(uint8_t locks, float3 &v) noexcept
{
    if (locks)
    {
        for (int k = 0; k < 3; ++k)
        {
            if (locks & (1 << k))
                v[k] = 0.0f;
        }
    }
}

// The body on the other side of an adjacency entry, or null when the entry does not
// touch `body`.
[[nodiscard]] static Rigid *otherBody(const ForceRef &ref, const Rigid *body,
        const ManifoldStore &manifolds) noexcept
{
    if (ref.isManifold())
    {
        const uint32_t slot = ref.slot();
        Rigid *a = manifolds.bodyA[slot];
        Rigid *b = manifolds.bodyB[slot];
        if (a == body)
            return b;
        if (b == body)
            return a;
        return nullptr;
    }
    const Force *f = ref.asForce();
    if (f->bodyA == body)
        return f->bodyB;
    if (f->bodyB == body)
        return f->bodyA;
    return nullptr;
}

// True when a force already joins the pair. Mirrors Rigid::constrainedTo: any joint,
// spring or ignore pair counts, and a live manifold counts, but a dormant one (contact
// grace) does not - the pair must still wake on AABB overlap and re-enter the narrow
// phase, or a sleeping body behind a dormant manifold is invisible to an incoming body.
[[nodiscard]] static bool pairConstrained(const ForceAdjacency &adjacency,
        const ManifoldStore &manifolds, Rigid *bodyA, Rigid *bodyB)
{
    for (const ForceRef &ref : adjacency.of(bodyA->broadIndex))
    {
        if (otherBody(ref, bodyA, manifolds) != bodyB)
            continue;
        if (ref.isManifold())
        {
            const uint32_t slot = ref.slot();
            if (manifolds.contactCount(slot) == 0 && manifolds.dormantTicks[slot] > 0)
                continue;
        }
        return true;
    }
    return false;
}

} // namespace avbd
