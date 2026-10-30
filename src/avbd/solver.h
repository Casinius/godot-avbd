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

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>
#include <array>
#include <span>

#ifdef AVBD_PROFILE
#include <chrono>
#endif

#include "avbd/maths.h"
#include "avbd/bvh/node_storage.hpp"

// Forward declarations
namespace bvh::nodes { class NodeStorage; }

namespace avbd {

// Defined in job_pool.hpp / taskflow.hpp, which are only included by solver.cpp.
namespace detail {
class JobPool;
class TaskFlow;
}

// Minimum penalty parameter
inline constexpr float PENALTY_MIN = 1.0f;
// Maximum penalty parameter
inline constexpr float PENALTY_MAX = 10000000000.0f;
// Cap of the adaptive per-contact betaScale multiplier (see SolverParams::Dual::adaptiveBeta).
inline constexpr float ADAPTIVE_BETA_MAX = 64.0f;
// Speculative contact margin: the normal constraint is satisfied at the margin surface,
// not at true zero separation, so a resting contact does not flicker in and out of the
// broad phase.
inline constexpr float COLLISION_MARGIN = 0.005f;
// Cap on the contact error a manifold will model: deeper than this, the residual is
// clamped (a shape that fell through its partner is not going to be rescued by one more
// unit of penalty).
inline constexpr float PENETRATION_MAX = 0.02f;
// Position threshold for sticking contacts (ie static friction)
inline constexpr float STICK_THRESH = 0.00001f;

struct Rigid;
struct Force;
struct Manifold;
struct Solver;

// Monotonic max / sum on a float held in an atomic, for the AVBD_DIAG counters. Best-effort
// reporting, not solver state: several workers can raise the same counter, and a lost race only
// means a slightly staler number.
inline void atomicMaxFloat(std::atomic<float> &slot, float p_value) noexcept
{
    float current = slot.load(std::memory_order_relaxed);
    while (p_value > current &&
            !slot.compare_exchange_weak(current, p_value, std::memory_order_relaxed)) {
    }
}

inline void atomicAddFloat(std::atomic<float> &slot, float p_value) noexcept
{
    slot.fetch_add(p_value, std::memory_order_relaxed);
}

// The collision and inertia shape of a body.
//
// `size` means different things per shape, and this documents it once:
//   Box      - full extents in x, y and z
//   Sphere   - size.x is the radius; y and z are ignored
//   Cylinder - size.x and size.z are the radius, size.y is the full height, and the axis is the
//              body's local +Y (matching Godot's CylinderMesh)
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
    float3 axis;          // cylinder axis, world space, unit length
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
    // Derived per-body state: the rotation as a matrix, and the offset from the pose this step
    // started at. Every constraint reads both, once per contact, and both follow from the pose - so
    // they are refreshed wherever the pose changes (Solver::warmstartBody, the combine pass, the
    // projection) instead of being recomputed by every constraint that needs them. Per body, not per
    // constraint: a manifold's eight contacts used to convert the same two quaternions eight times.
    // Solver state, not host state: never read across a step.
    float3x3 rot;
    float3 dqLin; // positionLin - initialLin
    float3 dqAng; // positionAng - initialAng, as a rotation vector
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

    // The one-byte flags, kept in one run so their padding is spent once instead of once per
    // field.
    bool sleeping = false;
    // True when every force on the body is a contact manifold. Only such bodies take the
    // Chebyshev extrapolation: the contact map is clamped and tolerates gamma > 1, hard
    // joints do not.
    bool contactOnly = true;
    // Per-axis motion locks, bit 0..2 = solver x, y, z. A locked axis gets no gravity and
    // cannot carry velocity into a step, which is what keeps a locked body hanging in
    // place. Constraint coupling can still move it within a step; that residual is a
    // documented approximation.
    uint8_t axisLockLinear = 0;
    uint8_t axisLockAngular = 0;
    int8_t sleep_mode = 0; // 0=NEVER, 1=SLEEP, 2=START_IN_SLEEP
    // Three bytes of this word are still free: the next flag lands here rather than in a new
    // cache line.

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

    // Sleep bookkeeping. All of it intrusive rather than a scratch column on purpose: a sleeping
    // body has no update slot (`updateSlot` is -1, it is out of the primal set), and the wake pass
    // reads this state for exactly those bodies; it also has to survive between steps, which the
    // per-step scratch does not.
    //
    // Idle credit: steps in which finishBody() found the body quiet, capped at the budget in
    // Solver::sleepFrames. Solver::updateSleepIslands() is the only thing that turns it into
    // `sleeping`, and it does so for a whole island at once.
    int32_t stillFrames = 0;
    // Consecutive steps the wake evidence has held for this sleeping body. One step's velocity
    // spike is a numerical nuisance, not a disturbance, so a body has to earn its wake-up over
    // SolverParams::Sleep::wakeDebounce steps.
    int32_t wakeFrames = 0;
    // The pose the rest test measures drift against: a one-tenth-per-step low-pass of the body's
    // own pose, so what it reports is roughly "how far has this body travelled over the last ten
    // steps". That is the position-domain, and therefore dt-independent, version of the idle test -
    // the velocity half of the test divides the same travel by the tick, which at 240 Hz reads four
    // times the speed it does at 60 Hz for a body that is moving exactly as slowly.
    float3 sleepRefLin;
    quat sleepRefAng;

    // Snap the rest window onto the current pose and drop the idle credit. A reference pose from
    // before a teleport, a wake-up or a rebuild says nothing about where the body is now, and
    // leaving it would either hold the body awake while the window catches up or, worse, hand it
    // back credit measured against a pose it no longer has.
    void resetRestWindow()
    {
        sleepRefLin = positionLin;
        sleepRefAng = positionAng;
        stillFrames = 0;
        wakeFrames = 0;
    }

    // Forced wake, the one that does not go through the island pass: a host that moves a body,
    // gives it a velocity or asks for it awake has broken the assumption the sleep state was built
    // on, and its rest window has to start over.
    void wake()
    {
        sleeping = false;
        resetRestWindow();
    }

    // Scratch-set membership, rewritten by Solver::selectUpdateBodies() and
    // Solver::updateSleepIslands() each step. These are handles into the solver's per-step arrays,
    // not state of the body: `updateSlot` is the index of this body in `Solver::updateOrder`
    // (and therefore into the columns of Solver::BodyScratch), and -1 marks a body the primal
    // phases may not move this step (static, sleeping). `islandSlot` is the same idea for the
    // union-find's island body list. A handle instead of a flag: the array index the phases
    // need is then free, and reading one says nothing about the other set.
    int32_t updateSlot = -1;
    int32_t islandSlot = -1;

    // Godot-style collision filtering: a pair collides when either side's layer is in the
    // other's mask. Defaults keep every pair colliding, which is what the solver did
    // before these existed.
    //
    // These two, the flags above and `host` below are host vocabulary: the primal never reads
    // them. They are read once each per step - pairOverlaps() for the layer/mask pair,
    // warmstartBody() for the axis locks, the update-set selection for sleep_mode - so they
    // live on the body they describe instead of costing a lookup per body per step. That is
    // the line this struct holds: solver state and host mapping stay here, per-pass working
    // sets go in Solver's scratch, and a new field has to say which side it belongs to.
    uint32_t collisionLayer = 1;
    uint32_t collisionMask = 1;

    // Host-side back-pointer: the solver never reads or writes it, it exists so a host resolves
    // a body - and, with it, its own record - from the solver object the broad phase hands back.
    // A body the host has not registered keeps nullptr. Forces get no such pointer: a host that
    // reaches a constraint has already reached at least one of its bodies, and the six pointers
    // plus the vtable fill this object's first 64 bytes exactly, so a byte here would cost a
    // whole word on every Manifold.
    void *host = nullptr;

    bool constrainedTo(Rigid *other) const;
};

// How much a body gives way along a constraint direction, for the position projection. The
// projection splits one correction between a constraint's two sides in proportion to these, so
// a side that cannot move this step (static, or sleeping and therefore out of the update set)
// takes none of it and the other side takes all of it.
//
// A body's inertia is diagonal (`Rigid::moment`), which is what makes the split arithmetic a
// componentwise product rather than a 3x3 solve.
struct Compliance
{
    float invMass = 0.0f;
    float3 invInertia = float3{0, 0, 0};

    [[nodiscard]] bool immobile() const noexcept { return invMass == 0.0f && invInertia.isZero(0.0f); }

    // The body's inverse inertia applied to a world-frame vector. Componentwise by definition:
    // `Rigid::moment` is diagonal, which is what makes the projection's split arithmetic a
    // multiply instead of a 3x3 solve. The inertia is stored in the body's own frame and the
    // solver's angular block treats the rotation vector as if it were in the world one (the
    // first-order approximation the primal stamps are built on), so the projection follows the
    // same convention rather than inventing a rotation of its own.
    [[nodiscard]] float3 invInertiaTimes(float3 v) const noexcept
    {
        return float3{invInertia.x() * v.x(), invInertia.y() * v.y(), invInertia.z() * v.z()};
    }

    // w = J M^-1 J^T for a row whose linear half is `n` and whose angular half is `arm x n`.
    // Invariant under flipping `n`, which is why one call serves both sides of a constraint.
    [[nodiscard]] float row(float3 n, float3 arm) const noexcept
    {
        const float3 jAng = arm.cross(n);
        return invMass + jAng.dot(invInertiaTimes(jAng));
    }
};

// Defined after Rigid, which is what it reads. `updateSlot >= 0` is the primal's own answer to
// "may this body be moved this step": the same test the iteration phases use, so a projection
// can never move a body the rest of the step has left frozen - that would break both the sleep
// bookkeeping and the thread-count independence the update set is there to protect.
inline Compliance complianceOf(const Rigid *body)
{
    Compliance c;
    if (body == nullptr || body->mass <= 0.0f || body->updateSlot < 0)
        return c;
    c.invMass = 1.0f / body->mass;
    // A zero or negative principal inertia means the body cannot rotate about that axis
    // (a cylinder's own axis, for instance); the projection leaves those alone.
    c.invInertia = float3{body->moment.x() > 0.0f ? 1.0f / body->moment.x() : 0.0f,
            body->moment.y() > 0.0f ? 1.0f / body->moment.y() : 0.0f,
            body->moment.z() > 0.0f ? 1.0f / body->moment.z() : 0.0f};
    return c;
}

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
    Block(float mass, float3 moment, float dtSq, float3 positionLin, float3 inertialLin,
          quat positionAng, quat inertialAng) noexcept :
            lhsLin(diagonal(mass, mass, mass) / dtSq),
            lhsAng(diagonal(moment.x(), moment.y(), moment.z()) / dtSq),
            lhsCross(float3x3::Zero()),
            rhsLin(diagonal(mass, mass, mass) / dtSq * (positionLin - inertialLin)),
            rhsAng(diagonal(moment.x(), moment.y(), moment.z()) / dtSq * (positionAng - inertialAng)) {}

    // Solve for the update without applying it (block-Jacobi path: each force stores its
    // partial delta and a later combine pass applies the sum over a body's constraints).
    void solveDelta(float3 &r_dxLin, float3 &r_dxAng) const noexcept {
        solve(lhsLin, lhsAng, lhsCross, -rhsLin, -rhsAng, r_dxLin, r_dxAng);
    }

    // Solve for and apply the update (Eq. 4).
    void apply(Rigid &body) const noexcept;
};

// Holds all user defined and derived constraint parameters, and provides a common interface for all forces.

// Which constraint a Force actually is, fixed at construction. The solver dispatches on this
// rather than calling contactPointCount() as a type test: that query answers "how many
// contacts does this manifold hold right now", and every one of its dispatch uses wanted the
// type. A fresh or emptied manifold is exactly where the two differ - it holds no contacts
// yet but is still a manifold - so the count stays a reporting query and the kind carries the
// dispatch.
enum class ForceKind : uint8_t
{
    Contact,      // Manifold
    Joint,        // Joint (revolute plus angle constraint)
    GenericJoint, // the 6-DOF joint
    Spring,
    Ignore,       // exists only so the broad phase skips a pair
};

struct Force
{
    Solver *solver;
    Rigid *bodyA;
    Rigid *bodyB;
    Force *nextA;
    Force *nextB;
    Force *next;
    const ForceKind kind;

    Force(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, ForceKind p_kind);
    virtual ~Force();

    // Same reason as Rigid: a force is a node in three intrusive lists at once.
    Force(const Force &) = delete;
    Force &operator=(const Force &) = delete;
    Force(Force &&) = delete;
    Force &operator=(Force &&) = delete;

    // Downcast, valid exactly when `kind == ForceKind::Contact`. Checking the kind is the
    // caller's job; this does not verify it, and there is nothing cheap to verify.
    [[nodiscard]] Manifold *asManifold();
    [[nodiscard]] const Manifold *asManifold() const;

    virtual bool initialize() = 0;
    virtual void updatePrimal(Rigid *body, float alpha, Block &block) = 0;
    virtual void updateDual(float alpha) = 0;

    // Adopt the current pose as this constraint's stabilisation reference
    // (SolverParams::Projection::rebase). Called after a projection pass has moved the bodies, and
    // only when that knob is on. A force with no stabilisation reference has nothing to adopt and
    // leaves the default, and a *contact* keeps the default on purpose: a unilateral row is held up
    // across rounds by the multiplier charge, and re-anchoring it as well tells the penalty that the
    // sunken pose is acceptable (measured on the 10-box stack: worst interface 0.00199 m with the
    // charge, 0.00845 m with both).
    virtual void rebaseStabilisation() {}

    // --- position projection (see SolverParams::Projection) -------------------
    //
    // The projection runs as two passes, and the split is the whole point: pass one is
    // dispatched over forces, so a constraint may write only its own step-scoped state; pass two
    // is dispatched over bodies, so a body may write only itself. That is what makes the layer
    // thread-count independent, and it is why the measured correction is stored here rather than
    // applied on the spot - applying it would be a second body writing this one's pose.
    //
    // Both default to doing nothing: a constraint that models a soft pull (a spring) or no pull
    // at all (an ignore pair) has no feasibility to guarantee.
    //
    // `scale` is the fraction of the error to remove, `maxCorrection` the per-body distance
    // ceiling in metres.
    virtual void projectMeasure(float scale, float maxCorrection)
    {
        (void)scale;
        (void)maxCorrection;
    }
    // This body's own share of what projectMeasure stored, summed into the two accumulators by
    // the caller. Zero for a body that is not one of this constraint's, and for one the pass may
    // not move.
    virtual void projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const
    {
        (void)body;
        (void)r_deltaLin;
        (void)r_deltaAng;
    }

    // How many contact points this force currently holds (0 for every non-contact force).
    // A reporting query - penetration statistics, the host's contact list. Dispatch uses
    // `kind`, which is also not a virtual call.
    virtual int contactPointCount() const { return 0; }

    // Accumulate penetration error (total) and positive-penetration contact count, in
    // contact order, for the adaptive iteration estimate. No-op for non-contact forces.
    virtual void accumulatePenetration(float &r_total, int &r_count) const { (void)r_total; (void)r_count; }

    // Whether this constraint has broken and stopped acting on its bodies. Only
    // breakable joints ever report true; every other force is permanent.
    virtual bool isBroken() const { return false; }
};

// Revolute joint + angle constraint between two rigid bodies, with optional fracture
struct Joint : Force
{
    float3 rA, rB;
    float3 C0Lin, C0Ang;
    float3 penaltyLin, penaltyAng;
    float3 lambdaLin, lambdaAng;
    float stiffnessLin, stiffnessAng, fracture;
    float torqueArm;
    bool broken;

    // Position projection state: the multiplier of each of the three anchor-error rows, in world
    // axes. Step-scoped for the same reason a manifold's cache is - it describes the correction
    // being applied right now. The share pass recomputes the arms and rows from the poses, which
    // have not moved between the two passes, so nothing else has to be stored.
    float3 projLin = float3{0, 0, 0};

    Joint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB, float p_stiffnessLin = INFINITY, float p_stiffnessAng = 0.0f, float p_fracture = INFINITY);

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;
    void projectMeasure(float scale, float maxCorrection) override;
    void projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const override;
    void rebaseStabilisation() override;

    // Break the joint immediately, discarding its dual state. Used when a host
    // application wants the fracture decision to survive being re-created (the
    // solver has no way to remember a broken joint once it is deleted).
    void breakNow();

    bool isBroken() const override { return broken; }
};

// Standard spring force
struct Spring : Force
{
    float3 rA, rB;
    float rest;
    float stiffness;

    Spring(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB, float p_stiffness, float p_rest = -1);

    bool initialize() override { return true; }
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;
};

// Force which has no physical effect, but is used to ignore collisions between two bodies
struct IgnoreCollision : Force
{
    // No state of its own, so it overrides nothing: this force exists purely so that two
    // bodies skip each other in the broad phase.
    IgnoreCollision(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB)
        : Force(p_solver, p_bodyA, p_bodyB, ForceKind::Ignore) {}

    bool initialize() override { return true; }
    void updatePrimal(Rigid *, float, Block &) override {}
    void updateDual(float) override {}
};

// How one degree of freedom of a GenericJoint is constrained.
enum class AxisMode : int
{
    Free = 0,    // unconstrained, apart from a spring if one is configured
    Locked = 1,  // held at the joint's zero
    Limited = 2, // free between the bounds, blocked outside them
};

// One degree of freedom of a GenericJoint: how it is constrained, plus the dual state the
// augmented Lagrangian keeps for it.
struct JointAxis
{
    AxisMode mode = AxisMode::Free;
    float lower = 0.0f; // used by Limited
    float upper = 0.0f;

    float springStiffness = 0.0f; // 0 = no spring
    float springEquilibrium = 0.0f;
    float springDamping = 0.0f;

    float lambda = 0.0f; // dual variable of the lock, or of the bound being held
    float penalty = 0.0f;
    float initial = 0.0f; // constraint value at the start of the step, for stabilisation
};

// A general 6-DOF joint: three linear and three angular degrees of freedom, each Free, Locked
// or Limited, each with an optional spring. The degrees of freedom are the axes of body A's
// local frame (world axes when bodyA is null).
//
// Angular degrees of freedom are measured from the pose the joint was created in. Before the
// locked axes are evaluated, the twist about every non-locked axis is removed: without that, a
// wheel spinning about its axle would drag the locked axes with it, because a rotation vector
// wraps every half turn and so cannot be read directly off the deviation. Removing one twist is
// exact; with several non-locked axes it is applied axis by axis, which is the usual approximation.
struct GenericJoint : Force
{
    float3 rA, rB; // anchors, in the local frame of each body
    // Orientation of B relative to A, in A's frame, at the moment the joint was created. Every
    // degree of freedom reads its displacement from this arrangement.
    quat rest;

    JointAxis linear[3];
    JointAxis angular[3];

    // Position projection state: the multiplier of each linear axis's row, along that axis in
    // world space. Step-scoped, like Joint::projLin and the manifold's cache.
    std::array<float, 3> projLinear = {};

    GenericJoint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB);

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;

    // Position projection of the three linear axes only; see projectMeasure for why the angular
    // degrees of freedom are left to the Lagrangian.
    void projectMeasure(float scale, float maxCorrection) override;
    void projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const override;
    void rebaseStabilisation() override;

    // Anchor of each body in world space.
    [[nodiscard]] float3 anchorA() const;
    [[nodiscard]] float3 anchorB() const;
    // Body A's orientation, or the world's if bodyA is null.
    [[nodiscard]] quat frameOrientation() const;


    // Current displacement of one degree of freedom, measured in the joint frame: metres for a
    // linear axis, radians for an angular one, both relative to the pose the joint was created
    // in. Reporting only - the solver never reads these back.
    [[nodiscard]] float linearValue(int p_axis) const;
    [[nodiscard]] float angularValue(int p_axis) const;
    // The whole angular measurement vector at once, cheaper than three separate calls.
    void angularValues(float p_out[3]) const;
    // The largest penalty this degree of freedom can use and still be solvable in the iterations
    // available. See the definition for why an unbounded penalty is a trap.
    [[nodiscard]] float penaltyLimit(const Rigid *p_body, float3 p_axis, bool p_angular) const;

private:
    // Constraint value of one degree of freedom for a given measurement, or false when it
    // applies no constraint this step (a Free axis, or a Limited one inside its range).
    [[nodiscard]] static bool axisConstraint(const JointAxis &p_axis, float p_value, float &r_constraint);
    // The force one degree of freedom applies this pass, and the constraint value behind it.
    // Shared by both passes so the primal and the dual never disagree about the error.
    [[nodiscard]] bool axisForce(const JointAxis &p_axis, float p_value, float p_alpha, float &r_constraint,
            float &r_force) const;
    // Stamp one degree of freedom into the body's system: its constraint, and its spring if it
    // has one.
    //
    // A constraint is stamped through *both* halves of the system: `p_jLin` is how it responds
    // to the body moving and `p_jAng` to the body turning. They are the two Jacobian rows of one
    // constraint, so a linear axis anchored off the centre of mass needs both - it moves the
    // anchor point as much by rotating the body as by translating it.
    void stampAxis(Block &block, const JointAxis &p_axis, float p_value, float p_alpha, float3 p_jLin,
            float3 p_jAng, float p_rate) const;
};

// Collision manifold between two rigid bodies, which contains up to eight frictional contact points
struct Manifold : Force
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

    // Contact point information for a single contact.
    //
    // This is the cross-step half of a manifold's state: warm-started dual variables that a
    // resting contact carries from one step to the next, matched into the new contact set by
    // feature id in initialize(). Anything that only means something inside the current step
    // belongs in JacobianCache below instead - the split is the rule, and the reason it is
    // worth a rule: Contact lives in a fixed eight-slot array, so a new field here grows every
    // pooled manifold block and the per-slot copy loop initialize() runs to warm the set up.
    // Before adding one, ask whether it survives the step; only then does it come here.
    struct Contact
    {
        FeaturePair feature;
        float3 rA; // contact offset in A's local space (relative to center)
        float3 rB; // contact offset in B's local space (relative to center)
        // The contact's linear constraint at the previous step's geometry, restamped by
        // initialize() for the current one: the stabilised form blends the live error against
        // it, and the adaptive-beta rule reads the value just before it is overwritten to
        // know how far this contact's error moved. Recomputed every step, but its predecessor
        // has to survive the boundary - which is what puts it here rather than in the cache.
        float3 C0;
        float3 penalty;
        float3 lambda;
        bool stick;
        // Personal multiplier on the betaLin penalty ramp, carried across steps with
        // the warmstart. Adaptive beta (SolverParams::Dual::adaptiveBeta) raises it for contacts
        // whose error stopped shrinking and lets it fall back for converging ones.
        float betaScale = 1.0f;
    };

    std::array<Contact, 8> contacts;
    float3x3 basis; // Normal in the first row (pointing from B to A), and tangents in the second and third rows
    int numContacts;
    float friction;

    // Per-contact cache of what the primal and the dual both need, refreshed once per round by
    // preCompute.
    //
    // Step-scoped state, the other half of the split documented on Contact: a field that only means
    // something inside the current iteration set goes here, not into the eight-slot persistent
    // array. The angular Jacobians used to live here as two 3x3 matrices per contact; they do not
    // any more, because the stamp is computed from the closed form (see stampCached) and the only
    // thing either half needs from the arms is the arms themselves.
    struct JacobianCache
    {
        float3 rAWorld;
        float3 rBWorld;
        // Position projection (SolverParams::Projection): the correction the measure pass worked
        // out for each row of this contact, as a multiplier in the contact basis, and which rows
        // it actually corrected. Step-scoped for the same reason the arms are - it describes the
        // correction being applied right now, and nothing about it survives to the next step. The
        // row mask is what the dual ramp reads to avoid ramping against an error already gone.
        float3 projLambda = float3{0, 0, 0};
        uint8_t projectedRows = 0;
    };
    std::array<JacobianCache, 8> jacobianCache;
    float3x3 rotACaptured = float3x3::Identity(); // body rotations at capture
    float3x3 rotBCaptured = float3x3::Identity();
    bool jacobiansCached = false;
    // The symmetric linear Hessian B^T K B of one contact, as its six independent entries. The basis
    // is orthonormal and the penalty diagonal, so the linear half of every stamp is a weighted sum of
    // the three basis-row outer products - the same six numbers for both bodies of the manifold, and
    // the same six for every round whose penalty did not change. initialize() and updateDual() are
    // the only places a penalty changes, so they are the only places this is recomputed.
    struct BasisHessian
    {
        float s00, s01, s02, s11, s12, s22;
    };
    std::array<BasisHessian, 8> linHessian;
    // Clamped contact forces for the round's stamps, plus the arms and the linear Hessians they need.
    // Written by preCompute for the first round of a step and by updateDual from then on - the dual
    // measures the same pose the next round's primal will, so it can hand it the forces directly and
    // the assembly pass disappears from every round after the first. `primalReady` says whether what
    // is in the cache belongs to the current pose; a projection pass moves poses after the dual that
    // filled it, so the round after one assembles again.
    std::array<float3, 8> contactForceCache;
    bool primalReady = false;

    Manifold(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB);

    // Class-specific pooled allocation: manifolds churn with contact formation and
    // separation (broad phase news them, warmstartForces deletes them), so recycled blocks
    // avoid the heap on the hot path. The constructor still fully runs on every new
    // expression, so a recycled block is indistinguishable from a fresh one.
    static void *operator new(std::size_t count);
    static void operator delete(void *ptr) noexcept;
    // Release every cached block (called when a solver dies; blocks outlive solvers by
    // design so a surviving solver's next contact pays no allocation).
    static void drainPool() noexcept;

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    // Two-phase primal: preCompute (parallel over forces) evaluates the shared per-
    // contact arms/jacobians/clamped forces once; stampCached (parallel over bodies)
    // stamps one body's side from that cache. updatePrimal is the self-contained
    // fallback that runs both.
    void preCompute(float alpha);
    void stampCached(Rigid *body, Block &block);
    void updateDual(float alpha) override;
    // Position projection. The measure pass reads the contact's error at the poses the rest of
    // the step has left, unstabilised, works out the multiplier that would take each violating row
    // to zero - shared between the two bodies in proportion to their compliance, the row's total
    // being the sum of theirs - and keeps it in the cache for the share pass. With dual feedback
    // it also charges that correction to the contact's multiplier, which is what stops the next
    // round re-fighting an error the projection has already paid for.
    void projectMeasure(float scale, float maxCorrection) override;
    void projectShare(const Rigid *body, float3 &r_deltaLin, float3 &r_deltaAng) const override;
    int contactPointCount() const override { return numContacts; }
    void accumulatePenetration(float &r_total, int &r_count) const override;

    static int collide(Rigid *bodyA, Rigid *bodyB, std::span<Contact> contacts, float3x3 &basis);
};

// Defined here rather than in Force, which is where they are declared: the downcast needs
// Manifold to be complete, and that happens further down this header.
inline Manifold *Force::asManifold() { return static_cast<Manifold *>(this); }
inline const Manifold *Force::asManifold() const { return static_cast<const Manifold *>(this); }

// Discrete shape-pair collision query: up to 8 contact points between two `Shape`s,
// independent of any solver. `contacts` must hold at least 8 slots (`Manifold`'s
// capacity) - exported because space queries reuse the exact path `Manifold` walks.
// `basis` receives an orthonormal frame whose first row is the contact normal (A -> B).
// Returns the contact count (0 = separated).
int collideShapes(const Shape &a, const Shape &b,
                  std::span<Manifold::Contact> contacts,
                  float3x3 &basis);

// Work per item differs by orders of magnitude between the phases below, so each one gets
// its own break-even batch size. Measured on a 12-core machine: one dispatch into the pool
// costs ~45 us, a body's primal update ~1 us, one manifold's contact solve ~85 ns, and a
// dual update ~200 ns. Break-even is therefore roughly 45 bodies, 500 contact solves and
// 250 dual updates - i.e. a scene has to be a few hundred bodies before splitting a phase
// across threads pays for the dispatch itself.
constexpr int kMinBodiesPerDispatch = 48;
constexpr int kMinForcesPerDispatch = 256;

// Tuning knobs, grouped by the phase that reads them, so a host can hand the whole set over
// in one assignment (Solver::setParams) instead of copying field by field, and so every knob
// has one name and one clamp site.
//
// Deliberately NOT in here: dt, gravity, iterations and threads. Those are per-tick inputs
// the phases read bare and often - dt alone appears in the inertial prediction, the Block
// construction and the BDF1 finish - and nesting them would add a hop to the hottest lines
// for no clarity.
struct SolverParams
{
    // The force law and the warm-start decay: read by every constraint, in every primal and
    // dual pass. `beta*` in particular is unit dependent, which is why these are exposed
    // rather than fixed.
    struct Dual
    {
        float alpha = 0.99f;       // Stabilization parameter, in (0, 1]
        float betaLin = 10000.0f;  // Penalty ramping for linear constraints
        float betaAng = 100.0f;    // Penalty ramping for angular constraints
        float gamma = 0.999f;      // Warmstarting decay, < 1

        // Adaptive beta: per-contact penalty-ramp multipliers carried across steps,
        // raised (x2, capped) for contacts whose error did not shrink against last step
        // and decayed (x0.5, floored at 1) for converging ones - the classic augmented
        // Lagrangian kappa update, per constraint instead of a global betaLin.
        bool adaptiveBeta = false;
        // Residual penetration depth (metres) at or below which a contact counts as solved
        // by the adaptive-beta update: the multiplier decays toward 1 instead of doubling,
        // so float-scale noise on a resting contact cannot ratchet it to the cap.
        float adaptiveBetaErrorFloor = 1e-3f;
    } dual;

    // How the iteration set is scheduled and relaxed: read by computeIterations and by the
    // primal/dual rounds, never by the constraint code.
    struct Primal
    {
        // Real-time phase split: the first ceil(newtonRatio * iterations) rounds run the full
        // Newton primal (dense 6x6 solve); the remaining rounds skip the primal entirely and
        // run dual-only relaxation on the existing jacobians. 1.0 = original AVBD behaviour.
        // Newton carries the stiff corrections; the dual-only tail is where the convergence
        // is cheapest to give up.
        float newtonRatio = 1.0f;
        // Per-iteration multiplicative decay applied to contact penalties during the
        // dual-only tail: a stiff penalty converges slowly under relaxation, softening it
        // lets the cheap rounds actually converge instead of spinning. 1.0 = no decay.
        float stiffnessDecay = 1.0f;

        // Upper bound of the adaptive iteration boost. The boost multiplies the round count
        // by the penetration ratio, so at impact moments a scene pays this many dense rounds
        // per substep - the frame spike the profiler blames. Lower it for smoother frames;
        // the healing those rounds did is recovered by substeps at high tick rates instead.
        //
        // Early-exit threshold on the mean penetration residual: when a Newton round
        // improves the average residual by less than this (or worsens it), the iteration set
        // ends early - remaining rounds would polish float noise at full price. 0 disables
        // the check (default: bit-identical iteration counts). Realistic value ~1e-5.
        float convergenceThreshold = 0.0f;

        // Combine-pass relaxation for the block-Jacobi primal. Every body solves its full
        // block system against the iteration-start pose snapshot (Jacobi, not Gauss-Seidel),
        // which for stiff densely-coupled scenes can overshoot where the sequential update
        // would not; scaling the applied delta by omega < 1 trades per-iteration progress
        // for stability. 1.0 = unrelaxed. The Chebyshev variant replaces this constant with
        // a per-iteration schedule. Suite sweep: 0.5-0.7 hold every rigid scene; 1.0 drifts
        // the 32-layer pyramid. Applies to the bodies whose forces are all contacts; see
        // jointRelaxation for the other half.
        float jacobiRelaxation = 0.7f;

        // Combine-pass relaxation for a body carrying any non-contact force (a joint, spring or
        // world anchor; Rigid::contactOnly is the test). A contact-only body's rows are unilateral
        // and clamped, and its Jacobi map tolerates the plain relaxation; a stiffness-coupled chain
        // does not. Block-Jacobi on a stiff chain with a hard end has a spectral radius above 1 at
        // omega = 0.7 - the error grows, the penalty ramps against the growing error, and the chain
        // is expelled (measured on the 20-link Rope: omega 0.7 -> |x| = 1618 m after 1000 steps,
        // 0.5 -> 20.0 m, where the chain's own length is 20). This is the damping that keeps the
        // chain's map inside the unit circle, and it is why the jointed scenes needed the projection
        // layer to look stable at all. It costs convergence per round, so it is only paid by the
        // bodies that need it.
        float jointRelaxation = 0.5f;

        // Chebyshev semi-iteration over the block-Jacobi primal (0 = off, the default). The
        // plain Jacobi estimate x + delta is convex-combined with the previous iterate
        // through the Chebyshev weights for a symmetric iteration matrix with spectral
        // radius <= rho: gamma_0 = 1, gamma_1 = 1/(1 - rho^2/2), then
        // gamma_{k+1} = 1/(1 - rho^2 * gamma_k / 4). Only pure-contact bodies take it
        // (Rigid::contactOnly): extrapolation diverges joint chains. Suite sweeps: rho=0.9
        // recovers deep-stack penetration (0.058 -> 0.027 m on the 10-box stack) but
        // diverges the 32-layer stack's interfaces; rho=0.7 is safe but the prev-mixing
        // smoothing gives back what the extrapolation gains. The textbook schedule assumes
        // a linear map; this one is clamped with ramping penalties, so it ships default-off
        // as a per-scene experiment rather than a win.
        float chebyshevRho = 0.0f;

        // Incremental jacobian rebuild distance (metres): after the first Newton round
        // assembles a contact's angular jacobian, later rounds reuse it while the contact's
        // world-space moment arms have drifted less than this from the captured values.
        // 0 disables caching (full assembly every round; default, bit-identical digests).
        // Realistic value 1e-4 (0.1 mm): in-iteration arm drift is millimetre-scale.
        float jacobianRebuildDistance = 0.0f;
    } primal;

    // The position projection layer (PBD-style), which sits on top of the augmented Lagrangian
    // and does what a penalty method structurally cannot: make the delivered pose feasible.
    // The Lagrangian converges towards zero error over rounds, its penalty is capped and its
    // residual is clamped, so whatever the round budget did not pay for stays as penetration.
    // Projection takes the error that is left and moves the bodies out of it directly, split
    // between the two sides by their own compliance so it stays a Jacobi step and therefore
    // thread-count independent.
    //
    // `rounds` is the whole switch: 0 runs none of the passes, and every digest in the suite is
    // the one measured before this layer existed.
    struct Projection
    {
        // Trailing rounds of the iteration set that end in a projection pass, and the closing
        // pass after the set. 0 = off, 8 = every round.
        int rounds = 0;
        // How much of the measured error one pass removes. 1.0 would be exact projection;
        // below 1 it is a damped step, which is what keeps a stack of projections from
        // oscillating when several constraints fight over the same body.
        float scale = 0.9f;
        // Ceiling, in metres, on the correction one pass applies to one body. Measured down from
        // the 0.01 first shipped: a ceiling-sized move is a teleport, and the jointed Bridge
        // absorbs it at 3 mm (|x| stays 20 at every round count) but throws itself apart at 10 mm
        // with only the trailing two rounds projected (|x| 1382).
        float maxCorrection = 0.003f;
        // Fraction of the projected displacement that shows up as velocity. 1 - standard PBD -
        // because 0 is measured to be the unstable half: subtracting the move back out of BDF1
        // leaves the pose saying one thing and its velocity another, and the jointed scenes diverge
        // at every ceiling tried (Bridge 1141 - 2152 at transfer 0, held at 20 at transfer 1).
        float velocityTransfer = 1.0f;
        // Charge the projected correction to the constraint's dual variable (the multiplier)
        // instead of leaving it as a purely geometric move, and let the penalty ramp skip the
        // rows the projection already resolved. On, the correction survives the next round
        // instead of being re-fought by the Lagrangian from scratch, which is what lets a
        // smaller iteration budget reach the same residual.
        // Which of the two ledgers the layer uses for the correction it delivers. Both are on the
        // same pose change, but only one of them can be booked per constraint kind:
        //
        // Contacts keep `dualFeedback`: the charge is what holds the lift across rounds (a resting
        // stack that is only re-anchored sinks back between passes), and the clamps make it safe -
        // a contact can only press and friction cannot exceed the cone, so the charge cannot
        // outgrow the load.
        //
        // Joints take `rebase`: their rows are bilateral and unbounded, so a charge there has no
        // clamp to hold it and no equilibrium to converge to. Re-anchoring the stabilisation
        // reference instead moves the pose the Lagrangian accepts, and leaves the multiplier
        // carrying the load it was already carrying. See `rebase` below for the full argument.
        bool dualFeedback = true;

        // Re-anchor each *bilateral* constraint's stabilisation reference onto the pose a pass just
        // delivered (Joint::C0Lin / JointAxis::initial), instead of charging the multiplier. On by
        // default: it is the only one of the two ledgers that is stable for a row with no clamp, and
        // while `rounds` is 0 - the shipped default - the knob cannot reach any scene at all.
        // Contacts are excluded (see Force::rebaseStabilisation): their charge is what holds the lift,
        // and re-anchoring them too costs a factor of four in resting penetration.
        //
        // This is the two-ledger problem in one sentence. The stabilised row is satisfied at
        // `C = alpha * C0`, so C0 - not lambda - is the state that declares which pose the
        // Lagrangian is willing to accept; lambda is only the force at that pose. A projection pass
        // moves the pose. Charging the move to lambda leaves C0 behind, and the next round's Newton
        // step pulls the body straight back towards `alpha * C0` with a force of order `K * C`,
        // while the charged lambda pushes with an unrelated force of its own - two injections per
        // pass that do not cancel. On a bilateral chain - the jointed Bridge, whose rows carry a
        // gravity load and never reach zero error, because the stabilisation keeps 1% of C0 per
        // step - that is an oscillator. Measured on the Bridge (1000 steps, threads=8,
        // newtonRatio=1, jacobiRelaxation=0.7, projection rounds=2, iterations=12): charging the
        // multiplier leaves |x|=1100 m and lambda in the thousands, moving the reference holds
        // |x|=20.0 m with lambda three times lower and no charge at all.
        bool rebase = true;

        // Run the projection passes on the trailing Newton rounds only, not on the dual-only tail
        // of the set. In the tail `primal.stiffnessDecay` is softening the contact penalty round by
        // round, so the surface the dual relaxes against is moving away from the surface the
        // projection measures - and joint penalties are not decayed at all, so the two halves of
        // the set drift apart. Off, the window is the last `rounds` rounds of the set.
        bool trailingNewtonOnly = false;

        // Which half of the constraint set the passes may touch. Separate switches, because the two
        // halves measured nothing alike (run_scene over the core scenes, 300 steps, threads=8, with
        // AVBD_PROJECTION_ROUNDS and the two scope knobs below):
        //   joints: the baseline loses its jointed chains - by step 300 the Bridge reaches
        //     |x|=617 and the Rope 1075 - and projection holds the Bridge at 20 with either scope
        //     and the Rope at 20 with both. Nothing else in the suite moves at all, because nothing
        //     else has joints. A hard anchor the Lagrangian leaves 2.7 cm open closes to zero
        //     (avbd_core_test, "projection closes a hard anchor").
        //   contacts: pays for the residual in scene motion. The stack's worst box-box interface
        //     goes 0.1145 -> 0.0020 m and the 512-box Pyramid XL improves too (314 -> 152), but the
        //     16-layer Pyramid, which the baseline holds at |x|=8.008 for 1000 steps, loses a box
        //     at every scale (0.25 - 1), ceiling (0.001 - 0.01), depth gate, round count and
        //     transfer tried, the most conservative configuration included (|x| 17 - 184). The
        //     mass-ratio stack (6.97 -> 250) and the breakable beam (5.24 -> 18.05) go the same way.
        // That asymmetry is why `rounds` stays 0 by default: the joint half has no measured cost,
        // the contact half trades residual for stability, and a project has to decide which of the
        // two it is buying.
        bool joints = true;
        bool contacts = true;
    } projection;

    // Idle detection: a movable body whose recent travel (see Rigid::sleepRefLin) stays under the
    // tolerances falls asleep - velocity zeroed, Godot semantics - once that has held for
    // `timeToSleep` seconds. 0 disables automatic sleeping, which is the default the core tests
    // measure against; the Godot server opens it at Godot's own space default.
    //
    // The velocity thresholds are deliberately NOT part of the rest test: solver jitter is a fixed
    // distance per step, so its BDF1 velocity is a distance divided by the tick and no velocity
    // gate can be tick-rate independent. They still own the loudest half of the idle decision
    // (hysteresisRatio) and the whole of the wake test in Solver::updateSleepIslands, which is
    // where a fast neighbour is the question being asked anyway.
    struct Sleep
    {
        float linearThreshold = 0.1f;       // m/s
        float angularThreshold = 0.139626f; // rad/s (Godot default: 8 deg/s)
        // Seconds of uninterrupted idle before a body sleeps. Seconds, not steps, because the
        // clock the scene runs on is a setting: at tick_rate 120 a step is half as long, and a
        // budget counted in steps would sleep everything twice as early. Solver::sleepFrames is
        // the same budget converted for the dt of the tick being run.
        float timeToSleep = 0.0f;
        // The position half of the rest test: how far the body may have travelled over the window
        // its reference pose defines (roughly ten steps) and still count as idle. Metres and
        // radians, so neither depends on the tick the scene runs on. Both sit well above the
        // collision margin's own 5 mm of resting overlap because that overlap is a static offset
        // the reference pose follows, not travel - only movement accumulates against it.
        float positionTolerance = 0.005f;
        float rotationTolerance = 0.02f;
        // Hysteresis. A body inside the tolerances earns a step of idle credit; a body beyond
        // `hysteresisRatio` times a tolerance loses all of it; in between it only leaks credit, so
        // a stack whose top box twitches at one and a half times the tolerance does not have to
        // restart its whole budget from nothing. The ratio is the width of that middle band, and
        // the reason it exists is the thing this section is about: the transition itself used to be
        // the disturbance.
        float hysteresisRatio = 3.0f;
        // Steps a wake condition has to hold for before the island is woken. One step of a
        // neighbour's velocity is the same noise the rest test rides out, and an island that wakes
        // on it restarts every budget inside it - which is what a settled stack did at high tick
        // rates, over and over. A neighbour moving `hysteresisRatio` times the threshold is not
        // noise, and wakes the island on the spot.
        int wakeDebounce = 2;
    } sleep;
};

// Core solver class which holds all the rigid bodies and forces, and has logic to step the simulation forward in time
struct Solver
{
    // Per-tick inputs, tuned for metre/kilogram/second scenes.
    float dt = 1.0f / 60.0f; // Timestep
    float gravity = -10.0f;  // Gravity, along -Z
    // The round budget: how many primal+dual rounds one step runs, split by
    // SolverParams::Primal::newtonRatio. This is the number the host asked for and the number the
    // solver spends - 16 because that is what the pinned scenes need at their tolerance (a box
    // resting on the slab lands at 1.6 cm of skin; below 16 the 10-box stack's worst interface breaks
    // the suite's 5 cm bound). It used to be multiplied by a penetration estimate and clamped to a
    // second limit, which is why it was effectively a constant: this is that constant, stated where a
    // host can see it. Cheaper settings are a measurement, not a default - with the projection layer
    // on, 12 rounds beat 16 without it by a wide margin on the same scenes (mean penetration 0.0033
    // vs 0.0200 m on the 10-box stack, 0.0004 vs 0.0058 m on the Pyramid), for a third less time.
    int iterations = 16;

    // Worker threads for the per-body phases. 0 = one per hardware thread, 1 = run
    // everything inline on the calling thread. The result does not depend on this value:
    // every body solves its own block system against the same pose snapshot, and applies
    // its own delta.
    int threads = 0;

    SolverParams params;

    // The host's entry point for the tuning knobs. Range fixing happens here, once, so the
    // phases can read their knobs bare and a bad project setting cannot reach the math.
    void setParams(SolverParams p)
    {
        p.primal.newtonRatio = std::clamp(p.primal.newtonRatio, 0.0f, 1.0f);
        p.primal.stiffnessDecay = std::clamp(p.primal.stiffnessDecay, 0.0f, 1.0f);
        p.primal.jacobiRelaxation = std::clamp(p.primal.jacobiRelaxation, 0.0f, 1.0f);
        p.primal.jointRelaxation = std::clamp(p.primal.jointRelaxation, 0.0f, 1.0f);
        p.primal.jacobianRebuildDistance = std::max(p.primal.jacobianRebuildDistance, 0.0f);
        p.primal.chebyshevRho = std::max(p.primal.chebyshevRho, 0.0f);
        p.projection.rounds = std::clamp(p.projection.rounds, 0, 8);
        p.projection.scale = std::clamp(p.projection.scale, 0.0f, 1.0f);
        p.projection.maxCorrection = std::max(p.projection.maxCorrection, 0.0f);
        p.projection.velocityTransfer = std::clamp(p.projection.velocityTransfer, 0.0f, 1.0f);
        p.sleep.linearThreshold = std::max(p.sleep.linearThreshold, 0.0f);
        p.sleep.angularThreshold = std::max(p.sleep.angularThreshold, 0.0f);
        p.sleep.timeToSleep = std::max(p.sleep.timeToSleep, 0.0f);
        p.sleep.positionTolerance = std::max(p.sleep.positionTolerance, 0.0f);
        p.sleep.rotationTolerance = std::max(p.sleep.rotationTolerance, 0.0f);
        // Below 1 the leaky band is empty and the hard-reset threshold sits inside the quiet
        // threshold, which would reset the credit of a body the same step it earned it.
        p.sleep.hysteresisRatio = std::max(p.sleep.hysteresisRatio, 1.0f);
        p.sleep.wakeDebounce = std::max(p.sleep.wakeDebounce, 1);
        params = p;
    }

    // --- Instrumentation ---------------------------------------------------------------
    // Both surfaces are compiled out of the shipped build. They are configure options rather than
    // local hacks because the defines have to match across the library, the test binary and the
    // GDExtension: `xmake f --profile=y --diag=y`.
    enum PhaseId : int
    {
        kPhaseWarmstart = 0, // body warmstart: inertial prediction + derived pose state
        kPhasePrecompute,    // manifold assembly over forces (first round of a step only)
        kPhaseStamp,         // per-body block assembly, solve, and apply for a contact-only body
        kPhaseCombine,       // delta application for the bodies a joint couples
        kPhaseDual,          // per-force dual update, and the tail's penalty decay
        kPhaseProject,       // projection measure + share
        kPhaseFinish,        // BDF1 velocity finish
        kPhaseCount,
    };

    struct PhaseProfile
    {
        double work = 0.0; // seconds inside the phase's item loop, worker 0
        double sync = 0.0; // seconds at the barrier that fences it, worker 0
    };

    struct StepProfile
    {
        // Totals over the window since the last reset, not the last substep's numbers: a single
        // step's barrier wait swings by an order of magnitude with what else the machine is doing,
        // so a mean over a few hundred is the only figure worth reading. Divide by `samples`.
        PhaseProfile phases[kPhaseCount];
        double wall = 0.0; // whole substep (the region), measured on the calling thread
        // The stages outside the iteration region, which its phase timers cannot see: broad phase,
        // force warm-start (mostly collision detection inside initialize), and update-set selection.
        // `step` is the whole tick, so these plus `wall` explain it.
        double broad = 0.0;
        // Broad phase, split: the candidate-pair search (Morton build + traversal + sort) and the
        // manifold objects the surviving pairs become. The second half is the constraint lifecycle,
        // not the tree, and it is charged here because that is where it runs.
        double pairQuery = 0.0;
        double bvhBuild = 0.0;
        double pairSearch = 0.0;
        double manifoldCreate = 0.0;
        double forceInit = 0.0;
        double selectStage = 0.0;
        double step = 0.0;
        double rounds = 0.0; // rounds actually run (Newton + relaxation), total
        long long samples = 0; // substeps accumulated
        long long samplesPerStep = 0; // ticks accumulated (the stage timings are per tick)
        void reset() noexcept;
    };

    // The slot worker 0 should time into, or null. Null for every other worker - one thread's view
    // of, and wait at, the phase is the measurement - and for every build without AVBD_PROFILE,
    // where the callers' clock calls fold away entirely.
    [[nodiscard]] PhaseProfile *phaseSlot(int id, int workerIndex) noexcept;
    [[nodiscard]] const StepProfile &stepProfile() const noexcept { return profileData; }
    // Start a fresh measurement window (the profile accumulates until this is called).
    void resetStepProfile() noexcept { profileData.reset(); }

#ifdef AVBD_DIAG
    // What the projection layer actually did this step, for the bridge investigation. The point is
    // to separate "the layer delivered pose" from "the layer injected force": a pose correction
    // shows up as `maxBodyMove` and shrinks `maxJointC0`; an injection shows up as `charge*` and
    // raises `maxJointLambda` without the load having changed.
    struct Diag
    {
        std::atomic<float> maxJointLambda{0.0f}; // largest |lambda| on any joint row
        std::atomic<float> maxJointC0{0.0f};     // largest |C0| any joint row stabilised against
        std::atomic<float> chargeContact{0.0f};  // sum |d(lambda)| the projection charged, contacts
        std::atomic<float> chargeJoint{0.0f};    // ... and joints
        std::atomic<float> askedContact{0.0f};   // sum |C| the measure pass asked to remove
        std::atomic<float> askedJoint{0.0f};
        std::atomic<float> maxBodyMove{0.0f}; // largest single-body move a share pass applied
        void reset() noexcept;
    } diag;
#endif

    // Heads of the two intrusive lists. Owned: `clear()` deletes through them.
    Rigid *bodies = nullptr;
    Force *forces = nullptr;

    Solver();
    ~Solver();

    // Bodies and forces are linked into the solver and into each other by raw pointers, and
    // a step deletes forces (see warmstartForces), so neither can be copied or moved: a copy
    // would share the list heads and double-delete them.
    Solver(const Solver &) = delete;
    Solver &operator=(const Solver &) = delete;
    Solver(Solver &&) = delete;
    Solver &operator=(Solver &&) = delete;

    // Closest body along a ray, and the hit point in that body's local space (null if none).
    [[nodiscard]] Rigid *pick(float3 origin, float3 dir, float3 &local);
    // Masked variant with Godot's pair semantics: a body is a candidate when
    // (body.layer & p_mask) != 0. The two-argument form above forwards with an
    // accept-everything mask. Both walk the last built broad-phase tree (ray/AABB
    // pruned); before any step has built one, they fall back to a body-list walk.
    [[nodiscard]] Rigid *pick(float3 origin, float3 dir, float3 &local, uint32_t p_mask);
    // Bodies whose last-built broad-phase bounds overlap the query AABB, filtered by
    // (layer & p_mask) != 0. Backs the direct space state's shape and point queries;
    // the same no-tree fallback applies.
    void collectAABB(float3 min, float3 max, uint32_t p_mask, std::vector<Rigid *> &out) const;
    // Delete every body and force.
    void clear();
    void step();

    // Decide which bodies the primal phases may move this step (wake-up, island sleep and
    // the sleep-mode gate), and size the per-body delta buffers.
    void selectUpdateBodies();

    // Number of worker threads actually in use (1 when running inline).
    [[nodiscard]] int threadCount() const;

private:
    // Every body in list order. The inertial/warmstart and BDF1 phases run over this,
    // because they also have to keep the cached state of static bodies current: contact
    // and joint code reads `initialLin`/`initialAng` of *both* bodies of a constraint, and
    // the ground is usually one of them.
    std::vector<Rigid *> warmstartOrder;

    // Broad-phase: the body list flattened to index-addressable order (rebuilt every
    // step), plus the BVH built over those bodies and the candidate pair list.
    std::vector<Rigid *> bodiesInOrder;
    bvh::nodes::NodeStorage bvhNodes;
    std::vector<std::pair<int, int>> sweepPairs;

    // Bodies that take part in the primal update (movable, awake ones), in list order.
    // The Jacobi primal solves one block system per body against the iteration-start
    // pose, holding the delta in the scratch columns below until the combine pass.
    std::vector<Rigid *> updateOrder;
    // Load-balance scratch for the update set (see selectUpdateBodies): the per-body constraint count
    // the sort reads, and the permutation it is applied through.
    std::vector<int> updateCost;
    std::vector<int> updateOrderIndex;
    std::vector<Rigid *> updateOrderSorted;

    // The per-step, per-body working set, addressed by Rigid::updateSlot. One struct rather
    // than loose members so the columns that travel together are declared and sized
    // together - one resize call, in selectUpdateBodies - and an index can never be valid in
    // one array and stale in the next.
    struct BodyScratch
    {
        std::vector<float3> deltaLin; // the block solve's linear half, per body
        std::vector<float3> deltaAng; // ... and its angular half
        // Previous linear iterate, the Chebyshev combination's x_{k-1}; only read when
        // chebyshevRho > 0. Angular iterates are not mixed: quats are not a vector space, so
        // the angular delta only takes the schedule's gamma.
        std::vector<float3> prevLin;
        // What the position projection has moved this body so far this step, summed over every
        // pass (SolverParams::Projection). `finishBody` gives this back out of the BDF1 velocity
        // only when projection_velocity_transfer asks for a purely positional layer; at the
        // default the move stays in the velocity, as standard PBD. Zeroed once per step in
        // warmstartBody, not per pass.
        std::vector<float3> projLin;
        std::vector<float3> projAng;

        void resize(size_t bodyCount)
        {
            deltaLin.resize(bodyCount);
            deltaAng.resize(bodyCount);
            prevLin.resize(bodyCount);
            projLin.resize(bodyCount);
            projAng.resize(bodyCount);
        }
    } scratch;

    // The force list rebuilt each step so the phases can address forces by index, plus the
    // per-force result of initialize() so that deletion can be done in one thread after the
    // parallel initialisation.
    std::vector<Force *> forceOrder;
    std::vector<uint8_t> forceActive;

    // Island scratch (rebuilt when sleep is active): the union-find's member list, kept as a
    // member so its capacity survives between steps. Rigid::islandSlot is the index side of
    // the deal, which is why the per-force index pairs no longer need storing.
    std::vector<Rigid *> islandBodies;

    // The idle window expressed in steps, derived from Sleep::timeToSleep and `dt` at the top of
    // every tick (computeSleepFrames). Idle detection counts steps because the state it compares
    // against, Rigid::stillFrames, is a step counter; the budget is seconds because the tick is a
    // setting. 0 means idle detection is off.
    int sleepFrames = 0;

    std::unique_ptr<detail::JobPool> pool;
    // Sequence of classified stages for one step, rebuilt each time `runSubstep` runs.
    std::unique_ptr<detail::TaskFlow> flow;

    // Adaptive iteration count parameters

    // Create the pool on first use so that a serial solver never spawns threads.
    void ensurePool();
    // Fill `forceOrder` from the solver's force list and return the count.
    int collectForces();

    // The phases of a step, in the order `step()` runs them.
    void broadPhase();
    int warmstartForces();
    // One body's inertial/warm-start update (independent across the list).
    void warmstartBody(int i);
    // Decide this step's round count once, before the region starts.
    void computeIterations(int &r_targetIterations, int &r_newtonRounds);
    // Convert the seconds-based idle window into steps for the dt of this tick.
    void computeSleepFrames();
    // The whole iteration set as one persistent-worker region: `targetIterations` rounds of the
    // block-Jacobi primal (one full block solve per body against the iteration-start pose, applied
    // on the spot when the body has no joint on it) plus dual passes, fenced by LoopSync barriers.
    // Two fences per round, and the combine pass only exists for the bodies a joint couples.
    // Both sync types run the identical body; NullSync is the inline/reference path.
    template <typename Sync>
    void iterateBody(int workerIndex, Sync &sync, int targetIterations, int forceCount,
            int newtonRounds, float &r_prevResidual, std::atomic<bool> &r_converged);
    // Apply one body's solved delta with the combine pass's schedule: the per-kind relaxation, or
    // the Chebyshev convex combination for a contact-only body, plus the derived pose state every
    // constraint reads. Shared by the assembly phase (a contact-only body applies its own delta
    // there) and the combine pass (the bodies a joint couples).
    void applyDelta(int i, Rigid *body, float3 dLin, float3 dAng, float gamma, float rho, float relax,
            float jointRelax) noexcept;
    // One body's BDF1 velocity update and idle detection (independent across the list).
    void finishBody(int i);
    // The position projection layer: two fenced passes, measure over forces then share and apply
    // over the update set. Runs as part of the iteration region (the trailing rounds) and again
    // as the closing stage of a step, so both call sites drive it through the same `Sync`.
    // Does nothing when projection_rounds is 0, and issues no barriers in that case.
    template <typename Sync>
    void projectPass(Sync &sync, int forceCount, int workerIndex);
    // The projection's per-item halves: one constraint's measure, one body's share and apply, and
    // - when SolverParams::Projection::rebase is on - one constraint's re-anchoring.
    void projectForce(int i);
    void projectBody(int i);
    void projectRebase(int i);
    // Per-phase wall clock, filled only under AVBD_PROFILE (see Instrumentation above).
    StepProfile profileData;
    // True when the update set holds a body with a non-contact force. Those are the bodies whose
    // stamp reads another body's pose (a joint measures both anchors), so they are the ones that
    // cannot apply their delta inside the assembly phase; when there are none, the combine pass does
    // not exist and a round is two barriers instead of three.
    bool hasJointedBodies = false;
    // One substep, as a flow: a single region holding every phase (see runSubstep).
    void runSubstep(int forceCount, int targetIterations, int newtonRounds);
    // The substep's phases, on the region's participants.
    template <typename Sync>
    void substepBody(int workerIndex, Sync &sync, int forceCount, int targetIterations,
            int newtonRounds, std::atomic<bool> &converged);

    // Island-wide sleep bookkeeping, and the only thing that changes `sleeping` on the solver's own
    // authority (it runs when sleeping is enabled): movable sleep-enabled bodies are grouped into
    // islands through force connectivity, an island whose every member has earned its idle credit
    // sleeps as a unit, and an island whose sleeper has held wake evidence for
    // Sleep::wakeDebounce steps wakes as a unit. Either half alone thrashes: boxes sleeping one at a
    // time leaves a chain that is part static and part solved, and waking one box of a settled
    // island restarts the credit of every body it touches.
    void updateSleepIslands();
};

} // namespace avbd
