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
#include <memory>
#include <vector>
#include <array>
#include <span>

#include "avbd/maths.h"
#include "avbd/spatial_hash.hpp"

namespace avbd {

// Defined in job_pool.hpp, which is only included by solver.cpp.
namespace detail {
class JobPool;
}

// Minimum penalty parameter
inline constexpr float PENALTY_MIN = 1.0f;
// Maximum penalty parameter
inline constexpr float PENALTY_MAX = 10000000000.0f;
// Cap of the adaptive per-contact betaScale multiplier (see Solver::adaptiveBeta).
inline constexpr float ADAPTIVE_BETA_MAX = 64.0f;
// Base collision margin: a 1 mm skin keeps contacts from flickering on/off at the
// exact surface. It is deliberately small - the margin is added to the contact
// error, so a large value is a resting penetration floor.
inline constexpr float COLLISION_MARGIN_BASE = 0.001f;
// Margin damping factor per iteration
inline constexpr float COLLISION_MARGIN_DAMPING = 0.5f;
// Cap on penetration margin
inline constexpr float COLLISION_MARGIN_MAX = 0.02f;
// Base collision margin for compatibility
inline constexpr float COLLISION_MARGIN = COLLISION_MARGIN_BASE;
// Position threshold for sticking contacts (ie static friction)
inline constexpr float STICK_THRESH = 0.00001f;

struct Rigid;
struct Force;
struct Manifold;
struct Solver;

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
    // Consecutive finishVelocities() steps with both velocities under the solver's sleep
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
    // True when every force on the body is a contact manifold. Only such bodies take the
    // Chebyshev extrapolation: the contact map is clamped and tolerates gamma > 1, hard
    // joints do not.
    bool contactOnly = true;
    int sleep_mode = 0; // 0=NEVER, 1=SLEEP, 2=START_IN_SLEEP

    bool constrainedTo(Rigid *other) const;
};

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
          quat positionAng, quat inertialAng) noexcept {
        // The M/dt^2 diagonal is one value per matrix, then componentwise products.
        const float3 massOverDtSq = float3{mass, mass, mass} / dtSq;
        const float3 momentOverDtSq = moment / dtSq;
        lhsLin = diagonal(massOverDtSq.x, massOverDtSq.y, massOverDtSq.z);
        lhsAng = diagonal(momentOverDtSq.x, momentOverDtSq.y, momentOverDtSq.z);
        lhsCross = float3x3{}; // zero
        rhsLin = cmul(massOverDtSq, positionLin - inertialLin);
        rhsAng = cmul(momentOverDtSq, positionAng - inertialAng);
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

// Holds all user defined and derived constraint parameters, and provides a common interface for all forces.
struct Force
{
    Solver *solver;
    Rigid *bodyA;
    Rigid *bodyB;
    Force *nextA;
    Force *nextB;
    Force *next;

    Force(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB);
    virtual ~Force();

    // Same reason as Rigid: a force is a node in three intrusive lists at once.
    Force(const Force &) = delete;
    Force &operator=(const Force &) = delete;
    Force(Force &&) = delete;
    Force &operator=(Force &&) = delete;


    virtual bool initialize() = 0;
    virtual void updatePrimal(Rigid *body, float alpha, Block &block) = 0;
    virtual void updateDual(float alpha) = 0;

    // Number of contact points this force contributes (0 for non-contact forces).
    // Used for statistics/reporting only; the solver never reads it.
    virtual int contactPointCount() const { return 0; }

    // Accumulate penetration error (total) and positive-penetration contact count, in
    // contact order, for the adaptive iteration estimate. No-op for non-contact forces.
    virtual void accumulatePenetration(float &r_total, int &r_count) const { (void)r_total; (void)r_count; }

    // Whether this constraint has broken and stopped acting on its bodies. Only
    // breakable joints ever report true; every other force is permanent.
    virtual bool isBroken() const { return false; }

    // True for collision manifolds. The primal's per-body dispatch used a virtual
    // contactPointCount() call per incident force per round to learn this; a plain
    // flag keeps the hot walk free of indirect calls (5% of serial time at 1653
    // bodies was the vtable tax).
    bool isManifold = false;

    // Host-side ownership tag. The solver never reads or writes it; it exists so a
    // host can map a force back to its own object, including one the solver removed.
    void *owner = nullptr;
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

    Joint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB, float p_stiffnessLin = INFINITY, float p_stiffnessAng = 0.0f, float p_fracture = INFINITY);

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;

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
    // No state of its own, so it inherits the base constructor and overrides nothing:
    // this force exists purely so that two bodies skip each other in the broad phase.
    using Force::Force;

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

    GenericJoint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB);

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;

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

    // Contact point information for a single contact
    struct Contact
    {
        FeaturePair feature;
        float3 rA; // contact offset in A's local space (relative to center)
        float3 rB; // contact offset in B's local space (relative to center)
        float3 C0;
        float3 penalty;
        float3 lambda;
        bool stick;
        // Personal multiplier on the betaLin penalty ramp, carried across steps with
        // the warmstart. Adaptive beta (Solver::adaptiveBeta) raises it for contacts
        // whose error stopped shrinking and lets it fall back for converging ones.
        float betaScale = 1.0f;
    };

    std::array<Contact, 8> contacts;
    float3x3 basis; // Normal in the first row (pointing from B to A), and tangents in the second and third rows
    int numContacts;
    float friction;

    // Incremental jacobian cache: world-space moment arms + angular jacobians captured
    // at the last full assembly, plus the body rotations at capture time. The drift test
    // is per BODY (rotations are loop-invariant across a manifold's contacts): rebuild
    // everything when either rotation moved more than jacobianRebuildDistance (unit-
    // column distance ~ 2 sin(angle/2)). Invalidated by initialize() (new contact set).
    struct JacobianCache
    {
        float3 rAWorld;
        float3 rBWorld;
        float3x3 jAAng;
        float3x3 jBAng;
    };
    std::array<JacobianCache, 8> jacobianCache;
    float3x3 rotACaptured = identity3x3; // body rotations at capture
    float3x3 rotBCaptured = identity3x3;
    bool jacobiansCached = false;
    // Clamped contact forces from the primal's shared precompute pass; the two body
    // stamps read these instead of each recomputing the same C/F from the snapshot.
    std::array<float3, 8> contactForceCache;

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
    int contactPointCount() const override { return numContacts; }
    void accumulatePenetration(float &r_total, int &r_count) const override;

    static int collide(Rigid *bodyA, Rigid *bodyB, std::span<Contact> contacts, float3x3 &basis);
};

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

// Core solver class which holds all the rigid bodies and forces, and has logic to step the simulation forward in time
struct Solver
{
    // Parameter defaults. Tuned for metre/kilogram/second scenes; `beta*` in particular is
    // unit dependent, which is why they are exposed rather than fixed.
    float dt = 1.0f / 60.0f; // Timestep
    float gravity = -10.0f;  // Gravity, along -Z
    int iterations = 10;     // Solver iterations per step

    float alpha = 0.99f;    // Stabilization parameter, in (0, 1]
    float betaLin = 10000.0f;  // Penalty ramping for linear constraints
    float betaAng = 100.0f;    // Penalty ramping for angular constraints
    float gamma = 0.999f;      // Warmstarting decay, < 1

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

    // Adaptive beta: per-contact penalty-ramp multipliers carried across steps,
    // raised (x2, capped) for contacts whose error did not shrink against last step
    // and decayed (x0.5, floored at 1) for converging ones - the classic augmented
    // Lagrangian kappa update, per constraint instead of a global betaLin.
    bool adaptiveBeta = false;

    // Upper bound of the adaptive iteration boost. The boost multiplies the round count
    // by the penetration ratio, so at impact moments a scene pays this many dense rounds
    // per substep - the frame spike the profiler blames. Lower it for smoother frames;
    // the healing those rounds did is recovered by substeps at high tick rates instead.
    int adaptiveMaxIterations = 16;

    // Early-exit threshold on the mean penetration residual: when a Newton round
    // improves the average residual by less than this (or worsens it), the iteration set
    // ends early - remaining rounds would polish float noise at full price. 0 disables
    // the check (default: bit-identical iteration counts). Realistic value ~1e-5.
    float convergenceThreshold = 0.0f;

    // Combine-pass relaxation for the overflow bucket of the coloured primal - the
    // bodies the greedy colouring could not fit inside kMaxColours, which run as one
    // block-Jacobi round against a pose snapshot. Such a snapshot solve can overshoot
    // where the coloured Gauss-Seidel updates do not; scaling the applied delta by
    // omega < 1 trades per-iteration progress for stability. 1.0 = unrelaxed. The
    // Chebyshev variant replaces this constant with a per-iteration schedule. The
    // bucket is empty in ordinary scenes, so these are tuning knobs for pathological
    // constraint cliques. Suite sweep: 0.5-0.7 hold every rigid bucket scene; 1.0
    // drifts a fully-cliqued stack.
    float jacobiRelaxation = 0.7f;

    // Chebyshev semi-iteration over the overflow bucket's block-Jacobi round
    // (0 = off, the default). The
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

    // Worker threads for the per-body phases. 0 = one per hardware thread, 1 = run
    // everything inline on the calling thread. The result does not depend on this value:
    // every body solves its own block system against the same pose snapshot, and applies
    // its own delta.
    int threads = 0;

    // Idle detection: a movable body whose BDF1 velocities stay under both thresholds for
    // this many consecutive steps falls asleep (velocity zeroed, Godot semantics).
    // 0 disables automatic sleeping; the core tests rely on the default for bit-identical
    // digests. Thresholds follow Godot's space-parameter defaults.
    float sleepLinearThreshold = 0.1f;       // m/s
    float sleepAngularThreshold = 0.139626f; // rad/s (Godot default: 8 deg/s)
    int sleepFrames = 0;                     // steps below thresholds before sleeping

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
    // the sleep-mode gate), size the per-body delta buffers, and colour the constraint
    // graph for the Gauss-Seidel primal.
    void selectUpdateBodies();
    // Greedy kMaxColours-bounded colouring of the update set's constraint graph;
    // bodies that do not fit join the block-Jacobi overflow bucket.
    void colourGraph();

    // Number of worker threads actually in use (1 when running inline).
    [[nodiscard]] int threadCount() const;

    // Performance timing counters (ms).
    static double _time_broadPhase;
    static double _time_colourGraph;
    static double _time_solve;
    static double _time_finish;

    // Performance timing getters (ms). Returns 0 if timing is not active.
    [[nodiscard]] static double get_broadPhase_time() { return 0.0; }
    [[nodiscard]] static double get_colourGraph_time() { return 0.0; }
    [[nodiscard]] static double get_solve_time() { return 0.0; }
    [[nodiscard]] static double get_finish_time() { return 0.0; }
    // Reset all timing counters to zero.
    static void reset_timing();
    // Get current CPU time in microseconds.
    static double get_time();

private:
    // Every body in list order. The inertial/warmstart and BDF1 phases run over this,
    // because they also have to keep the cached state of static bodies current: contact
    // and joint code reads `initialLin`/`initialAng` of *both* bodies of a constraint, and
    // the ground is usually one of them.
    std::vector<Rigid *> warmstartOrder;

    // Broad-phase: the body list flattened to index-addressable order (rebuilt every
    // step), the hierarchical spatial hash over those bodies, and the candidate
    // pair list (sorted unique (min index, other), the manifold creation order).
    std::vector<Rigid *> bodiesInOrder;
    SpatialHash broadHash;
    std::vector<std::pair<int, int>> sweepPairs;

    // Bodies that take part in the primal update (movable, awake ones), in list order.
    // The coloured Gauss-Seidel primal solves one block system per body of each
    // colour and applies it immediately; bodies the greedy colouring cannot fit
    // inside kMaxColours land in an overflow bucket that runs one block-Jacobi
    // round (solve every delta against the same snapshot, then combine), using
    // the delta buffers below.
    std::vector<Rigid *> updateOrder;
    std::vector<float3> bodyDeltaLin;
    std::vector<float3> bodyDeltaAng;
    // Previous linear iterate per overflow-bucket body (the Chebyshev combination's
    // x_{k-1}; only read when chebyshevRho > 0). Angular iterates are not mixed:
    // quats are not a vector space, so the angular delta only takes the gamma.
    std::vector<float3> bodyPrevLin;

    // Constraint-graph colouring, rebuilt by colourGraph() every step. Bodies of
    // one colour share no force, so a colour's whole set can solve and apply in
    // parallel (Gauss-Seidel between colours, one barrier each). bodyColour is
    // indexed by updateOrder position and holds either a colour < kMaxColours or
    // kMaxColours itself for the block-Jacobi overflow bucket.
    static constexpr int kMaxColours = 8;
    static constexpr uint8_t kOverflowColour = kMaxColours;
    std::vector<uint8_t> bodyColour;
    std::vector<int> colourBegins; // kMaxColours+1 offsets into colourOrder
    std::vector<int> colourOrder;  // update-set indices grouped by colour
    std::vector<int> overflowOrder; // update-set indices of the overflow bucket
    int colourCount = 0;           // colours in use (<= kMaxColours)

    // The force list rebuilt each step so the phases can address forces by index, plus the
    // per-force result of initialize() so that deletion can be done in one thread after the
    // parallel initialisation.
    std::vector<Force *> forceOrder;
    std::vector<uint8_t> forceActive;

    // Island scratch (rebuilt when sleep is active): body index base for the union-find,
    // per-force body slots, and the island id of every movable body (-1 = no island).
    std::vector<Rigid *> islandBodies;
    std::vector<std::pair<int, int>> islandForceSlots;
    std::vector<int> bodyIsland;

    std::unique_ptr<detail::JobPool> pool;

    // Current collision margin for this frame (dynamically adjusted)
    float currentCollisionMargin = COLLISION_MARGIN_BASE;

    // Adaptive iteration count parameters
    float maxPenetrationError = 1e-4f; // Target penetration error tolerance
    bool adaptiveIterations = true;    // Enable adaptive iteration count

    // Compute dynamic collision margin based on penetration depth and timestep
    inline float computeCollisionMargin(float penetrationDepth, float timestep)
    {
        // Dynamically reduce penetration margin as solver converges
        // Use penetration depth to avoid oscillation
        return std::min(COLLISION_MARGIN_BASE + penetrationDepth * 0.1f *timestep, COLLISION_MARGIN_MAX);
    }

    // Create the pool on first use so that a serial solver never spawns threads.
    void ensurePool();
    // Fill `forceOrder` from the solver's force list and return the count.
    int collectForces();

    // The phases of a step, in the order `step()` runs them.
    void broadPhase();
    int warmstartForces();
    void warmstartBodies();
    void solveIterations(int forceCount);
    // The whole iteration set as one persistent-worker loop: `targetIterations` rounds of
    // the block-Jacobi primal (one full block solve per body against the iteration-start
    // pose, then a combine pass) plus dual passes, fenced by LoopSync barriers. Declared
    // here so `iterate` stays testable.
    void iterate(int targetIterations, int forceCount, int newtonRounds);
    void finishVelocities();

    // Island-wide sleep bookkeeping (runs when sleeping is enabled): movable bodies are
    // grouped into islands through force connectivity, and an island whose every body
    // cleared idle detection this step puts all its members to sleep together. Prevents
    // the partially-asleep contact chain from thrashing between states.
    void islandSleep();
};

} // namespace avbd
