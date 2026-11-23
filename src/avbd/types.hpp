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

#include <array>
#include <cstdint>
#include <limits>
#include <span>

#include "avbd/maths.h"

namespace avbd {

// Minimum penalty parameter
inline constexpr float PENALTY_MIN = 1.0f;
// Maximum penalty parameter
inline constexpr float PENALTY_MAX = 10000000000.0f;
// Cap of the adaptive per-contact beta multiplier (see Solver::adaptiveBeta).
inline constexpr float ADAPTIVE_BETA_MAX = 64.0f;
// 1 mm skin stops contact flicker at the surface. Deliberately small: it adds to contact error,
// so a large value becomes a resting penetration floor.
inline constexpr float COLLISION_MARGIN_BASE = 0.001f;
// Base collision margin for compatibility
inline constexpr float COLLISION_MARGIN = COLLISION_MARGIN_BASE;
// Position threshold for sticking contacts (ie static friction)
inline constexpr float STICK_THRESH = 0.00001f;
// Frozen sticky anchors stay glued while sliding (STICK_THRESH never grows; C0 goes phantom), so this
// drift re-adopts them from the collision pass and resets force history; smaller, at noise scale, churns.
inline constexpr float STICK_DRIFT_HEAL = 0.002f;
// Separation that drops a warm-started dual: ~0.99/tick decay leaks impulse for hundreds of ticks
// (measured 10+ m/s pop from a still pile). Above skin+tolerance = real gap; below it, noise must not churn.
inline constexpr float CONTACT_RELEASE_GAP = 2.0e-3f;

// Penetration the solve sees per step (m). Uncapped, v = dp/dt launches deep impacts at tens of m/s
// and the pile never settles; cap bounds expulsion to 2 cm/step, residual depth resolves over later steps.
inline constexpr float CONTACT_MAX_PENETRATION = 0.02f;

struct Rigid;
struct Force;
struct Manifold;
struct Solver;

// Collision/inertia shape. `size`: Box = full xyz extents; Sphere = size.x radius (yz ignored);
// Cylinder = x/y radius, z height, axis local +Z (Godot's Y-up is rotated by the server; solver is always +Z).
enum class ShapeType : int
{
    Box = 0,
    Sphere = 1,
    Cylinder = 2,
};

// Shape-only collision query (no solver needed, for space queries). `half`/`radius`/`halfHeight` per type:
// Box = half extents; Sphere = radius only; Cylinder = radius + halfHeight, axis unit +Z in `center`'s space.
struct Shape
{
    ShapeType type = ShapeType::Box;
    float3 center;
    quat rotation;
    float3 half;          // box
    float radius = 0;     // sphere, cylinder
    float halfHeight = 0; // cylinder
    float3 axis = float3{0, 0, 0}; // cylinder axis, world space, unit length
};

// Holds all the state for a single rigid body that is needed by AVBD
struct Rigid
{
    Solver *solver;
    Rigid *next;
    float3 positionLin;
    quat positionAng;
    float3 initialLin;
    quat initialAng;
    float3 inertialLin;
    quat inertialAng;
    float3 velocityLin;
    float3 velocityAng;
    float3 prevVelocityLin;
    float3 size; // See ShapeType for what this means per shape
    ShapeType shape;
    float mass;
    float3 moment;
    float friction;
    float radius;
    // Per-body gravity: signed magnitude along -Z, defaulting to the solver's global value; a host may
    // stamp a different value per body (how per-body gravity modes are expressed).
    float gravity;
    bool sleeping = false;
    // Consecutive steps under both sleep thresholds; reaching sleepFrames sets `sleeping`. Only
    // meaningful for movable bodies with sleep_mode enabled.
    int stillFrames = 0;

    // Godot-style: a pair collides if either layer is in the other's mask. Defaults (1/1) collide all,
    // matching pre-filtering behaviour.
    uint32_t collisionLayer = 1;
    uint32_t collisionMask = 1;

    // Per-axis locks (bits 0..2 = x/y/z): no gravity, no carried-in velocity. Constraint coupling can still
    // move a locked axis within a step - a documented approximation.
    uint8_t axisLockLinear = 0;
    uint8_t axisLockAngular = 0;

    // Mass/inertia follow from shape, hence constructor args rather than setters: a body with wrong
    // inertia is a latent bug. The five-argument form is a box.
    Rigid(Solver *p_solver, float3 p_size, ShapeType p_shape, float p_density, float p_friction,
            float3 p_position, float3 p_velocity = float3{0, 0, 0});
    Rigid(Solver *p_solver, float3 p_size, float p_density, float p_friction, float3 p_position,
            float3 p_velocity = float3{0, 0, 0});
    ~Rigid();

    // A body is a node in the solver's intrusive list; copying it would corrupt the links.
    Rigid(const Rigid &) = delete;
    Rigid &operator=(const Rigid &) = delete;
    Rigid(Rigid &&) = delete;
    Rigid &operator=(Rigid &&) = delete;

    // Rewritten by selectUpdateBodies() each step: `inUpdateSet` marks movable, awake bodies the primal may
    // move; static/sleeping bodies are never solved.
    bool inUpdateSet = false;
    // Position of this body in Solver::updateOrder, or -1 when out of the update
    // set. The colouring reads its neighbours' colours through it.
    int updateIndex = -1;
    // Position in Solver::bodiesInOrder (rebuilt each broad phase); adjacency and pair checks address
    // bodies by it, avoiding pointer hashing.
    int broadIndex = -1;
    // All forces are contact manifolds. Only these enter colourGraph's Jacobi overflow: contacts are
    // clamped, but hard joints do not tolerate the snapshot solve.
    bool contactOnly = true;
    int sleep_mode = 0; // 0=NEVER, 1=SLEEP, 2=START_IN_SLEEP
};

// Tight world-space axis-aligned half extents of a body's collision shape. Used by
// the broad phase and by the space queries, so both agree on the same bound.
[[nodiscard]] float3 aabbHalfExtent(const Rigid &body) noexcept;

// One body's 6x6 Newton system: two 3x3 diagonal blocks + coupling + RHS. Constraints stamp their
// gradient/Hessian here; passing it as one object replaced six `updatePrimal` reference params.
struct Block
{
    float3x3 lhsLin;   // lower triangle used
    float3x3 lhsAng;   // lower triangle used
    float3x3 lhsCross; // coupling between the linear and angular halves
    float3 rhsLin;
    float3 rhsAng;

    // M/dt^2 with error as RHS; keep the division - a reciprocal multiply is not bit-identical.
    // Angular mass = world inertia R I R^T, not body-frame I: all Jacobians are world-frame; only isotropic OK.
    Block(float mass, float3 moment, float dtSq, float3 positionLin, float3 inertialLin,
          quat positionAng, quat inertialAng) noexcept {
        // The M/dt^2 diagonal is one value per matrix, then componentwise products.
        const float3 massOverDtSq = float3{mass, mass, mass} / dtSq;
        lhsLin = diagonal(massOverDtSq.x, massOverDtSq.y, massOverDtSq.z);
        lhsCross = float3x3{}; // zero
        rhsLin = cmul(massOverDtSq, positionLin - inertialLin);

        // I_w = R diag(moment) R^T, then M_ang = I_w / dt^2 and the right-hand side is the
        // same matrix acting on the world rotation vector since the step started.
        const float3x3 rot = qmat(positionAng);
        const float3x3 inertiaWorld = mul(rot,
                mul(diagonal(moment.x, moment.y, moment.z), transpose(rot)));
        lhsAng = inertiaWorld / dtSq;
        rhsAng = mul(lhsAng, positionAng - inertialAng);
    }

    // Solve without applying (Jacobi path: store delta, combine applies it). Assembles the 6x6 in place,
    // no by-value matrix copies on the hot path.
    void solveDelta(float3 &r_dxLin, float3 &r_dxAng) const noexcept {
        std::array<std::array<float, 6>, 6> a{};
        
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                a[r][c] = lhsLin[c][r];
                a[r + 3][c + 3] = lhsAng[c][r];
                a[r + 3][c] = lhsCross[c][r];
                a[r][c + 3] = lhsCross[r][c];
            }
        }
        const std::array<float, 6> b{-rhsLin[0], -rhsLin[1], -rhsLin[2], -rhsAng[0], -rhsAng[1], -rhsAng[2]};
        std::array<float, 6> x{};
        if (ldlt6_solve(a, b, x)) {
            r_dxLin = float3{x[0], x[1], x[2]};
            r_dxAng = float3{x[3], x[4], x[5]};
            return;
        }
        // Degenerate block: the pivoted dynamic retry (its own fixed attempt fails the
        // same way, so this lands straight in the dynamic path).
        solve(lhsLin, lhsAng, lhsCross, -rhsLin, -rhsAng, r_dxLin, r_dxAng);
    }

    // Solve for and apply the update (Eq. 4).
    void apply(Rigid &body) const noexcept;
};

} // namespace avbd
