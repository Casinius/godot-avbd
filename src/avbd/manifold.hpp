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

// Manifolds live in the solver's ManifoldStore (SoA, stable slots), not per-pair objects; this struct
// just namespaces the shared records and solver-free pair query.
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

    // Contact geometry written by collision detection; hot solver state is in ManifoldStore's parallel
    // arrays, keeping the per-contact stride the primal/dual walk small.
    struct Contact
    {
        FeaturePair feature;
        float3 rA; // contact offset in A's local space (relative to center)
        float3 rB; // contact offset in B's local space (relative to center)
    };

    // Symmetric B^T K B (six entries); K diagonal and basis orthonormal, so both bodies' stamps share it.
    // It changes only with penalty, so initialize()/updateDual() refresh it once per round.
    struct BasisHessian
    {
        float s00, s01, s02, s11, s12, s22;
    };

    // Shape-pair collision query at the bodies' current poses: the solver's narrow phase
    // entry, also used by space queries.
    static int collide(Rigid *bodyA, Rigid *bodyB, std::span<Contact> contacts, float3x3 &basis);
};

// SoA over stable slots: a slot stays until released; `live` gives a compact creation-order walk.
// Contacts at slot*8+k - fixed stride wastes a partial tail but keeps warm-start matching and parallel writes trivial.
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

    // Refresh contacts and warm-start hot state. False if separated (and not in the grace window);
    // the slot stays allocated either way - the caller decides whether to release it.
    bool initialize(uint32_t slot, Solver &solver);
    // Recompute C0 at the current pose without re-running collision detection. Used between
    // interleaved substeps; invalidates the prepared primal force cache.
    void rebase(uint32_t slot, float maxPenetration) noexcept;
    // One cache, two consumers: colour GS reuses last round's dual arms (first body free, second recomputes
    // after poses moved); the Jacobi overflow preComputes all slots then stamps both bodies against one snapshot.
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

// CSR adjacency per broadIndex, rebuilt each step (pre-broad-phase for pair filter, post-warmstart for
// solve); colouring, primal, wake-up and islands read it. Manifold refs carry a store slot.
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

// Solver-free shape-pair query, <=8 contacts; `contacts` needs 8 slots (Manifold capacity, reused by space queries).
// `basis` row 0 = normal (A->B); speculativeMargin makes near-misses positive-error; 0 = overlap-only; returns count.
int collideShapes(const Shape &a, const Shape &b,
                  std::span<Manifold::Contact> contacts,
                  float3x3 &basis,
                  float speculativeMargin = 0.0f);

} // namespace avbd
