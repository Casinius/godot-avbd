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
// Base collision margin: a 1 mm skin keeps contacts from flickering on/off at the
// exact surface. It is deliberately small - the margin is added to the contact
// error, so a large value is a resting penetration floor.
inline constexpr float COLLISION_MARGIN_BASE = 0.001f;
// Base collision margin for compatibility
inline constexpr float COLLISION_MARGIN = COLLISION_MARGIN_BASE;
// Position threshold for sticking contacts (ie static friction)
inline constexpr float STICK_THRESH = 0.00001f;
// A sticky contact's anchors are material points frozen when static friction engaged, and
// because the two anchors stay glued to each other while the pair slides, their separation
// (what STICK_THRESH measures) never grows - stick never clears on its own, even as the
// frozen anchors slide metres off the touching surfaces and turn C0 into phantom depth.
// The merge heals instead: when a sticky pair's anchor world positions drift this far from
// where the current collision pass puts the contact, the anchors are re-adopted and the
// force history resets. The threshold sits at the collision-margin scale: below it the
// "drift" is numerical noise on a genuinely static contact and re-anchoring would churn
// resting manifolds for nothing.
inline constexpr float STICK_DRIFT_HEAL = 0.002f;
// Geometric separation at which a contact's warm-started dual is dropped. The dual is a
// force, and the alpha*gamma warm-start decay (~0.99/tick) bleeds an impact impulse off
// over hundreds of ticks; while any of it is left the contact keeps applying it even
// though the pair has visibly parted - measured as a box popping out of an otherwise
// still pile at 10+ m/s seconds after the impact that loaded it. Anything above the
// collision-margin skin plus this tolerance is a real gap, and a real gap carries no
// contact force. Below it, the skin's own noise must not churn resting contacts.
inline constexpr float CONTACT_RELEASE_GAP = 2.0e-3f;

// Cap on the penetration the constraint solve sees in one step (metres). A position-level
// step turns displacement into velocity at finishBody (v = dp/dt), so removing a deep
// impact overlap in a single step launches the pair at tens of m/s; that body then
// crosses other boxes before their next collision pass and the pile never settles. The
// cap bounds the per-step expulsion to 2 cm, i.e. <= 2 cm/dt, and the residual depth is
// resolved over the following steps.
inline constexpr float CONTACT_MAX_PENETRATION = 0.02f;

struct Rigid;
struct Force;
struct Manifold;
struct Solver;

// The collision and inertia shape of a body.
//
// `size` means different things per shape, and this documents it once:
//   Box      - full extents in x, y and z
//   Sphere   - size.x is the radius; y and z are ignored
//   Cylinder - size.x and size.y are the radius, size.z is the full height, and the axis is the
//              body's local +Z. (Godot's CylinderShape3D is Y-up; the server rotates it into
//              sim space, see godot_convert.hpp - the solver itself is always +Z.)
enum class ShapeType : int
{
    Box = 0,
    Sphere = 1,
    Cylinder = 2,
};

// A shape collision query: the collision and inertia shape of one body, decoupled from
// the body itself so shape pairs can be tested without a solver (space queries do this).
//
// `half` / `radius` / `halfHeight` mean different things per type, documented once:
//   Box      - half is half the full extents in x, y and z
//   Sphere   - radius; half and halfHeight unused
//   Cylinder - radius and halfHeight; the axis is the shape's local +Z, unit length,
//              expressed in the same space as `center` (the solver's cylinder is Z-axis)
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
    Force *forces;
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
    // Per-body gravity: the signed magnitude applied along -Z. Defaults to the
    // solver's global gravity; a host may stamp a different value per body (this
    // is how per-body gravity modes are expressed).
    float gravity;
    bool sleeping = false;
    // Consecutive steps with both velocities under the solver's sleep
    // thresholds; reaching sleepFrames flips `sleeping` on. Only meaningful when the body
    // is movable and sleep_mode enables sleeping.
    int stillFrames = 0;

    // Godot-style collision filtering: a pair collides when either side's layer is in the
    // other's mask. Defaults keep every pair colliding, which is what the solver did
    // before these existed.
    uint32_t collisionLayer = 1;
    uint32_t collisionMask = 1;

    // Per-axis motion locks, bit 0..2 = solver x, y, z. A locked axis gets no gravity and
    // cannot carry velocity into a step, which is what keeps a locked body hanging in
    // place. Constraint coupling can still move it within a step; that residual is a
    // documented approximation.
    uint8_t axisLockLinear = 0;
    uint8_t axisLockAngular = 0;

    // Mass and inertia follow from the shape, so it is a constructor argument rather than
    // something set afterwards: a body that exists with the wrong inertia is a bug waiting to
    // happen. The five-argument form is a box.
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

    // Block-Jacobi update-set bookkeeping, rewritten by Solver::selectUpdateBodies() each
    // step: `inUpdateSet` marks the movable, awake bodies whose pose the primal phases may
    // change. Bodies outside the set (static, sleeping) are never solved or moved.
    bool inUpdateSet = false;
    // Position of this body in Solver::updateOrder, or -1 when out of the update
    // set. The colouring reads its neighbours' colours through it.
    int updateIndex = -1;
    // Position of this body in Solver::bodiesInOrder, rewritten by the broad phase every
    // step. The constraint adjacency and the pair-existence check address bodies through
    // it, so a rebuild never has to hash pointers.
    int broadIndex = -1;
    // True when every force on the body is a contact manifold. Only such bodies are
    // routed into the overflow bucket's block-Jacobi round by colourGraph: the contact
    // map is clamped, hard joints do not tolerate the snapshot solve.
    bool contactOnly = true;
    int sleep_mode = 0; // 0=NEVER, 1=SLEEP, 2=START_IN_SLEEP
};

// Tight world-space axis-aligned half extents of a body's collision shape. Used by
// the broad phase and by the space queries, so both agree on the same bound.
[[nodiscard]] float3 aabbHalfExtent(const Rigid &body) noexcept;

// The 6x6 Newton system a single body solves each iteration, split into its two 3x3
// diagonal blocks plus the coupling block, and the matching right-hand side.
//
// This is the output of the primal phase: every constraint the body takes part in stamps
// its gradient and Hessian into one of these, and then it is solved for the update. Passing
// it as one object replaced six separate reference parameters on `updatePrimal`.
struct Block
{
    float3x3 lhsLin;   // lower triangle used
    float3x3 lhsAng;   // lower triangle used
    float3x3 lhsCross; // coupling between the linear and angular halves
    float3 rhsLin;
    float3 rhsAng;

    // Mass term of the system: M / dt^2, with the body's current error as the right-hand
    // side. `dtSq` is dt*dt: the division is kept as a division on purpose, since
    // multiplying by a reciprocal would not be bit-identical.
    //
    // The angular half is the *world* inertia R I R^T, not the body-frame diagonal I. The
    // rest of the system is world-frame: the constraint angle Jacobians are `r x n` with r
    // in world space (ManifoldStore::stampCached builds R S R^T, Joint stamps `r x jLin`),
    // and the update `positionAng + dAng` is a world-frame rotation vector. Pairing that with
    // a body-frame diagonal inertia makes the effective angular mass wrong by whatever
    // rotation the body happens to have: only an isotropic inertia (a sphere, a cube) is
    // unaffected, and the demo's own 1 x 0.5 x 1 box is not (0.104 / 0.167 / 0.104).
    // projectContacts' angularCompliance/angularResponse already did this correctly.
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

    // Solve for the update without applying it (block-Jacobi path: each body stores its
    // delta and a later combine pass applies it). Assembles the 6x6 directly from the
    // block's members - no by-value matrix copies on the hot path.
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
