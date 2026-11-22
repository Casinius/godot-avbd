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

#include "avbd/types.hpp"

namespace avbd {

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
    // contact order, for the convergence probe. No-op for non-contact forces.
    virtual void accumulatePenetration(float &r_total, int &r_count) const { (void)r_total; (void)r_count; }

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

    Joint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB, float p_stiffnessLin = std::numeric_limits<float>::infinity(), float p_stiffnessAng = 0.0f, float p_fracture = std::numeric_limits<float>::infinity());

    bool initialize() override;
    void updatePrimal(Rigid *body, float alpha, Block &block) override;
    void updateDual(float alpha) override;

    // Break the joint immediately, discarding its dual state. Used when a host
    // application wants the fracture decision to survive being re-created (the
    // solver has no way to remember a broken joint once it is deleted).
    void breakNow();

    bool isBroken() const override { return broken; }

private:
    // The same constraint value is read at step start (C0) and on every primal/dual pass.
    // Keeping the expression in one place stops the three callers from drifting apart.
    [[nodiscard]] float3 constraintLin() const noexcept;
    [[nodiscard]] float3 constraintAng() const noexcept;
    // Warm-start decay + material-stiffness cap, applied componentwise to a penalty vector.
    [[nodiscard]] static float3 warmPenalty(float3 p, float gamma, float stiffness) noexcept;
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


} // namespace avbd
