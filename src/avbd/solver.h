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

#include "avbd/maths.h"

namespace avbd {

// Defined in job_pool.hpp, which is only included by solver.cpp.
namespace detail {
class JobPool;
}

// Minimum penalty parameter
inline constexpr float PENALTY_MIN = 1.0f;
// Maximum penalty parameter
inline constexpr float PENALTY_MAX = 10000000000.0f;
// Margin for collision detection to avoid flickering contacts
inline constexpr float COLLISION_MARGIN = 0.01f;
// Position threshold for sticking contacts (ie static friction)
inline constexpr float STICK_THRESH = 0.00001f;

struct Rigid;
struct Force;
struct Manifold;
struct Solver;

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
    float3 size; // Full widths in each dimension
    float mass;
    float3 moment;
    float friction;
    float radius;

    Rigid(Solver *p_solver, float3 p_size, float p_density, float p_friction, float3 p_position, float3 p_velocity = float3{0, 0, 0});
    ~Rigid();

    // A body is a node in the solver's intrusive list; copying it would corrupt the links.
    Rigid(const Rigid &) = delete;
    Rigid &operator=(const Rigid &) = delete;
    Rigid(Rigid &&) = delete;
    Rigid &operator=(Rigid &&) = delete;

    // Graph-colouring label used by the parallel solver: bodies in one colour never share
    // a force, so they can be updated concurrently. Written only by Solver::colourGraph().
    int colour = -1;

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
          quat positionAng, quat inertialAng) noexcept :
            lhsLin(diagonal(mass, mass, mass) / dtSq),
            lhsAng(diagonal(moment.x, moment.y, moment.z) / dtSq),
            lhsCross{0, 0, 0, 0, 0, 0, 0, 0, 0},
            rhsLin(diagonal(mass, mass, mass) / dtSq * (positionLin - inertialLin)),
            rhsAng(diagonal(moment.x, moment.y, moment.z) / dtSq * (positionAng - inertialAng)) {}

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

    // Whether this constraint has broken and stopped acting on its bodies. Only
    // breakable joints ever report true; every other force is permanent.
    virtual bool isBroken() const { return false; }

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
    };

    Contact contacts[8];
    float3x3 basis; // Normal in the first row (pointing from B to A), and tangents in the second and third rows
    int numContacts;
    float friction;

    Manifold(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB);

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;
    int contactPointCount() const override { return numContacts; }

    static int collide(Rigid *bodyA, Rigid *bodyB, Contact *contacts, float3x3 &basis);
};

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

    // Worker threads for the per-body phases. 0 = one per hardware thread, 1 = run
    // everything inline on the calling thread. The result does not depend on this value:
    // the update order is fixed by the colouring, and threads only spread one colour's
    // independent bodies over more cores.
    int threads = 0;

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
    // Delete every body and force.
    void clear();
    void step();

    // Split the bodies into colours such that no two bodies in a colour share a force, and
    // store them grouped by colour in `updateOrder`. Rebuilt from the current constraint
    // graph at the start of every step, because contacts change as bodies move.
    void colourGraph();

    // Number of colours in the current ordering (1 means nothing can run in parallel).
    [[nodiscard]] int colourCount() const { return colours; }
    // Number of bodies in the largest colour: the most work one parallel phase can spread
    // over the worker threads. Reporting only.
    [[nodiscard]] int widestColour() const { return widest; }
    // Number of worker threads actually in use (1 when running inline).
    [[nodiscard]] int threadCount() const;

private:
    // Every body in list order. The inertial/warmstart and BDF1 phases run over this,
    // because they also have to keep the cached state of static bodies current: contact
    // and joint code reads `initialLin`/`initialAng` of *both* bodies of a constraint, and
    // the ground is usually one of them.
    std::vector<Rigid *> warmstartOrder;

    // Bodies that take part in the primal update (movable ones), grouped by colour:
    // colour c owns updateOrder[colourStart[c] .. colourStart[c + 1]).
    std::vector<Rigid *> updateOrder;
    std::vector<int> colourStart;
    int colours = 0;
    int widest = 0;

    // The force list rebuilt each step so the phases can address forces by index, plus the
    // per-force result of initialize() so that deletion can be done in one thread after the
    // parallel initialisation.
    std::vector<Force *> forceOrder;
    std::vector<uint8_t> forceActive;

    std::unique_ptr<detail::JobPool> pool;

    // Create the pool on first use so that a serial solver never spawns threads.
    void ensurePool();
    // Fill `forceOrder` from the solver's force list and return the count.
    int collectForces();

    // The phases of a step, in the order `step()` runs them.
    void broadPhase();
    int warmstartForces();
    void warmstartBodies();
    void solveIterations(int forceCount);
    void finishVelocities();
    void updatePrimal(Rigid *body);
};

} // namespace avbd
