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

#include <cstdint>
#include <span>
#include <vector>

#include "avbd/types.hpp"

namespace avbd {

struct Force;

// Maximum contact points a manifold can carry. The collision paths emit at most eight.
inline constexpr int kMaxManifoldContacts = 8;

// Collision manifold types. Manifolds themselves live in the solver's ManifoldStore
// (structure of arrays, stable slots) rather than in per-pair objects, so this struct is
// a namespace for the shared records plus the solver-free pair query.
struct Manifold
{
    // Used to track contact features between frames
    union FeaturePair
    {
        struct
        {
            char inR;
            char outR;
            char inI;
            char outI;
        };

        int key;
    };

    // Contact point geometry. This is what collision detection writes; the hot solver
    // state lives in ManifoldStore's parallel arrays, so the per-contact stride the primal
    // and dual walk stays small.
    struct Contact
    {
        FeaturePair feature;
        float3 rA; // contact offset in A's local space (relative to center)
        float3 rB; // contact offset in B's local space (relative to center)
    };

    // The symmetric linear Hessian B^T K B, six independent entries: K is diagonal and
    // the basis is orthonormal, so every pair of bodies of this manifold shares it. It
    // changes only when the penalty changes, so initialize()/updateDual() refresh it once
    // per round and both body stamps read it.
    struct BasisHessian
    {
        float s00, s01, s02, s11, s12, s22;
    };

    // Shape-pair collision query at the bodies' current poses: the solver's narrow phase
    // entry, also used by space queries.
    static int collide(Rigid *bodyA, Rigid *bodyB, std::span<Contact> contacts, float3x3 &basis);
};

// Every collision manifold of one solver, as parallel arrays over stable slots.
//
// A slot is allocated from `freeList` (or appended) and stays put until the manifold is
// released; `live` lists the slots that currently own a manifold, in creation order, so
// every per-manifold phase is a compact walk over `live`. Contact state is indexed
// `slot * kMaxManifoldContacts + k`; the fixed stride wastes a partly unused tail on
// manifolds with few contacts, but keeps warm-start matching, parallel per-slot writes
// and the per-slot passes trivial.
struct ManifoldStore
{
    // Per-slot state.
    std::vector<Rigid *> bodyA, bodyB;
    std::vector<float3x3> basis; // normal in the first row (B -> A), tangents in rows 2/3
    std::vector<uint8_t> numContacts;
    std::vector<float> friction;
    std::vector<uint8_t> flags;
    std::vector<int32_t> dormantTicks;
    std::vector<int32_t> dormantSavedCount;
    // Drift-gated arm cache: world arms captured at the last full assembly, plus the body
    // rotations at capture time (see Solver::jacobianRebuildDistance).
    std::vector<float3x3> rotACaptured, rotBCaptured;

    // Per-contact state, stride kMaxManifoldContacts.
    std::vector<Manifold::Contact> contacts;
    std::vector<float3> contactC0; // constraint value at q-
    std::vector<float3> contactPenalty, contactLambda;
    std::vector<float3> contactForceCache;
    std::vector<float3> rAWorld, rBWorld; // cached world arms
    std::vector<uint8_t> contactStick;
    std::vector<float> contactBetaScale; // per-contact betaLin multiplier (adaptive beta)
    std::vector<Manifold::BasisHessian> linHessian;

    std::vector<uint32_t> live;     // slots owning a manifold, walk order
    std::vector<uint32_t> freeList; // recycled slots, not in `live`

    // flags
    static constexpr uint8_t kArmsCached = 1;
    static constexpr uint8_t kPrimalCacheReady = 2;

    [[nodiscard]] uint32_t acquire(Rigid *p_bodyA, Rigid *p_bodyB);
    // Return a slot to the free list. The caller removes it from `live`.
    void release(uint32_t slot) noexcept;
    void clear();

    [[nodiscard]] int count() const noexcept { return static_cast<int>(live.size()); }
    [[nodiscard]] static constexpr int base(uint32_t slot) noexcept
    {
        return static_cast<int>(slot) * kMaxManifoldContacts;
    }
    [[nodiscard]] int contactCount(uint32_t slot) const noexcept { return numContacts[slot]; }
    [[nodiscard]] bool dormant(uint32_t slot) const noexcept { return dormantTicks[slot] > 0; }

    // Neither body can move this step: the manifold's bookkeeping is parked verbatim.
    [[nodiscard]] bool frozen(uint32_t slot) const noexcept;
    // Scale every contact's penalty in this slot (stiffness decay for the dual-only tail).
    void decayPenalty(uint32_t slot, float decay) noexcept;

    // Refresh the contact set for the step and warm-start the hot state. Returns false
    // when the pair separated (and is not parked in the contact-grace window); the slot
    // stays allocated either way and the caller decides whether to release it.
    bool initialize(uint32_t slot, Solver &solver);
    // Recompute C0 at the current pose without re-running collision detection. Used between
    // interleaved substeps; invalidates the prepared primal force cache.
    void rebase(uint32_t slot, float maxPenetration) noexcept;
    // Two consumers share one cache. updatePrimal (the coloured Gauss-Seidel path) uses
    // the force/moment-arm cache the previous round's dual left behind: the first body to
    // touch the manifold this round consumes it for free, and the second colour - by then
    // poses have moved - recomputes once (preCompute) and stamps its side. The overflow
    // bucket's Jacobi path instead runs preCompute over slots, then stampCached over
    // bodies, so one assembly serves both sides against the same snapshot.
    void preCompute(uint32_t slot, const Solver &solver, float alpha) noexcept;
    void stampCached(uint32_t slot, Rigid *body, Block &block) const noexcept;
    void updatePrimal(uint32_t slot, const Solver &solver, Rigid *body, float alpha, Block &block);
    void updateDual(uint32_t slot, const Solver &solver, float alpha) noexcept;
    void accumulatePenetration(uint32_t slot, const Solver &solver, float &r_total,
            int &r_count) const noexcept;
};

// One incident constraint in the per-body adjacency: a non-manifold Force* (low bit 0,
// pointers are aligned) or a manifold slot tag (low bit 1).
struct ForceRef
{
    uintptr_t value = 0;

    [[nodiscard]] static ForceRef manifold(uint32_t slot) noexcept
    {
        return { (static_cast<uintptr_t>(slot) << 1) | 1u };
    }
    [[nodiscard]] static ForceRef force(Force *f) noexcept
    {
        return { reinterpret_cast<uintptr_t>(f) };
    }
    [[nodiscard]] bool isManifold() const noexcept { return (value & 1u) != 0; }
    [[nodiscard]] uint32_t slot() const noexcept { return static_cast<uint32_t>(value >> 1); }
    [[nodiscard]] Force *asForce() const noexcept
    {
        return reinterpret_cast<Force *>(value);
    }
};

// Per-body constraint adjacency in CSR form, keyed by Rigid::broadIndex. Rebuilt every
// step (before the broad phase for the pair filter, after warmstart for the solve): the
// colouring, the primal walk, the sleep wake-up and the island grouping read it instead
// of following the intrusive body lists. Manifold references carry a store slot.
struct ForceAdjacency
{
    std::vector<uint32_t> start; // bodyCount + 1 offsets
    std::vector<ForceRef> refs;

    void clear();
    // `bodies` is the step's flattened body list; every body in it has a valid
    // Rigid::broadIndex. Bodies outside it are skipped (they are not part of the step).
    void build(std::span<Rigid *const> bodies, const Force *forceHead, const ManifoldStore &manifolds);
    [[nodiscard]] std::span<const ForceRef> of(int bodyIndex) const noexcept
    {
        if (bodyIndex < 0 || static_cast<size_t>(bodyIndex) + 1 >= start.size())
            return {};
        return std::span<const ForceRef>(refs.data() + start[bodyIndex],
                start[bodyIndex + 1] - start[bodyIndex]);
    }

private:
    std::vector<uint32_t> cursor; // fill positions, reused across builds
};

// Discrete shape-pair collision query: up to 8 contact points between two `Shape`s,
// independent of any solver. `contacts` must hold at least 8 slots (`Manifold`'s
// capacity) - exported because space queries reuse the exact path `Manifold` walks.
// `basis` receives an orthonormal frame whose first row is the contact normal (A -> B).
// `speculativeMargin` lets separated shapes report contacts up to that distance apart
// (the constraint error is then positive); 0 keeps the query strictly overlap-only.
// Returns the contact count (0 = separated).
int collideShapes(const Shape &a, const Shape &b,
                  std::span<Manifold::Contact> contacts,
                  float3x3 &basis,
                  float speculativeMargin = 0.0f);

} // namespace avbd
