/*
 * General 6-DOF joint.
 *
 * Part of godot-avbd; the augmented-Lagrangian structure follows the reference
 * implementation (github.com/savant117/avbd-demo3d, MIT), as the rest of src/avbd does.
 *
 * Each of the six degrees of freedom is Free, Locked or Limited, and may carry a spring. A
 * degree of freedom is measured in the joint frame, which is body A's local frame anchored at
 * rA - or the world when bodyA is null, which is how the joint pins a body to a fixed pose.
 * Measurements are in physical units: metres for a linear axis, radians for an angular one,
 * and radians for the limits, so they mean the same thing whatever size the bodies are.
 *
 * Two things here are not a direct copy of anything in the reference implementation.
 *
 * The angular side has to survive a spinning axle. A locked axis cannot read its component of
 * the rotation vector between the two bodies: a wheel turning about its axle sweeps those
 * components through the full circle every half turn, so a hinge built that way tears itself
 * apart after a few revolutions. Instead the deviation from the joint's rest orientation is
 * split into the twist about each non-locked axis and the remaining swing, and only the swing
 * is handed to the locked axes.
 *
 * And a constraint anchored away from a body's centre acts through a moment arm, so it is
 * stamped through both halves of the system at once - see stampAxis.
 */

#include <cmath>

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {

namespace {

const quat kIdentity{0, 0, 0, 1};

const float3 &frameAxis(int p_axis)
{
    static const float3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    return axes[p_axis];
}

} // namespace

GenericJoint::GenericJoint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB)
    : Force(p_solver, p_bodyA, p_bodyB), rA(p_rA), rB(p_rB)
{
    // The arrangement the joint is created in is its zero. Held as B-relative-to-A in A's frame,
    // the same form angularValues() computes for the live pose, so the two compose directly.
    rest = normalize(conjugate(frameOrientation()) * p_bodyB->positionAng);
}

quat GenericJoint::frameOrientation() const
{
    return bodyA ? bodyA->positionAng : kIdentity;
}



float3 GenericJoint::anchorA() const
{
    return bodyA ? transform(bodyA->positionLin, bodyA->positionAng, rA) : rA;
}

float3 GenericJoint::anchorB() const
{
    return transform(bodyB->positionLin, bodyB->positionAng, rB);
}

// Remove the twist about axis `p_axis` from `p_dev`, writing its angle to `r_angle`.
//
// A quaternion's twist about n is (n * (v.n), w) normalised, so the angle is 2*atan2(v.n, w).
// That stays continuous through a full turn, which is what a spinning axle needs. Near a half
// turn the twist quaternion degenerates (both the projected vector part and w vanish), and
// there is then no twist to remove.
static void removeTwist(quat &p_dev, int p_axis, float &r_angle)
{
    const float3 &n = frameAxis(p_axis);
    const float proj = dot(float3{p_dev.x, p_dev.y, p_dev.z}, n);

    const quat twist{n.x * proj, n.y * proj, n.z * proj, p_dev.w};
    if (lengthSq(twist) <= 1.0e-12f)
    {
        r_angle = 0.0f;
        return;
    }

    r_angle = 2.0f * std::atan2(proj, p_dev.w);
    p_dev = normalize(p_dev * inverse(twist));
}

void GenericJoint::angularValues(float p_out[3]) const
{
    // Deviation from the rest configuration, in the joint frame (body A's).
    //
    // Both quantities below are B relative to A as seen *by A*: the current one, and the one
    // captured when the joint was created. The deviation is then the rest pose taken away from
    // the current one, which is identity exactly when the bodies are in their rest arrangement.
    //
    // The order matters and is easy to get wrong in a way nothing notices: putting the rest term
    // anywhere else happens to be harmless while the rest pose is the identity (every joint built
    // from upright bodies), and silently wrong as soon as a body starts out rotated - a wheel
    // turned onto its axle, say. Then the locked axes see an error proportional to the rest
    // rotation, which grows as the body spins and tears the joint apart.
    // B relative to A now, and the same quantity at rest. Taking the rest arrangement away from
    // the current one gives identity exactly when the bodies are back where they started, and its
    // axis-angle in A's frame is what each degree of freedom measures.
    const quat qA = frameOrientation();
    const quat relativeNow = normalize(conjugate(qA) * bodyB->positionAng);
    quat dev = normalize(rest * conjugate(relativeNow));

    // The twist about every non-locked axis belongs to that axis - a spring or a limit reads
    // its angle - and must not be visible to the locked ones. Removing one twist is exact; with
    // several non-locked axes it is applied one axis at a time, the usual approximation.
    for (int i = 0; i < 3; i++)
    {
        if (angular[i].mode != AxisMode::Locked)
            removeTwist(dev, i, p_out[i]);
    }

    // What remains is swing: a rotation about an axis perpendicular to every non-locked axis,
    // so its components along the locked axes are the errors those axes have to correct.
    for (int i = 0; i < 3; i++)
    {
        if (angular[i].mode == AxisMode::Locked)
            p_out[i] = dev[i] * 2.0f;
    }
}

float GenericJoint::linearValue(int p_axis) const
{
    const float3 offset = rotate(conjugate(frameOrientation()), anchorA() - anchorB());
    return offset[p_axis];
}

float GenericJoint::angularValue(int p_axis) const
{
    float values[3];
    angularValues(values);
    return values[p_axis];
}

bool GenericJoint::axisConstraint(const JointAxis &p_axis, float p_value, float &r_constraint)
{
    switch (p_axis.mode)
    {
        case AxisMode::Locked:
            r_constraint = p_value;
            return true;

        case AxisMode::Limited:
            // Whichever bound has been left. The sign of the constraint already points back
            // towards the allowed range, so no separate sign handling is needed.
            if (p_value > p_axis.upper)
            {
                r_constraint = p_value - p_axis.upper;
                return true;
            }
            if (p_value < p_axis.lower)
            {
                r_constraint = p_value - p_axis.lower;
                return true;
            }
            return false;

        case AxisMode::Free:
        default:
            return false;
    }
}

// The force one degree of freedom applies this pass, and the constraint value behind it.
//
// Both passes must agree on this value: the dual pass stores the force it produces as the next
// multiplier and grows the penalty from the same error. Computing it in one place is what keeps
// the augmented Lagrangian consistent - a lock stabilised in only one of the two would have the
// primal fighting the full error while the dual believed it was correcting a fraction of it,
// which converges to the wrong answer rather than merely converging slowly.
bool GenericJoint::axisForce(const JointAxis &p_axis, float p_value, float p_alpha, float &r_constraint,
        float &r_force) const
{
    if (!axisConstraint(p_axis, p_value, r_constraint))
        return false;

    // A lock is stabilised against where the step started, so its error is worked out over
    // several steps instead of being fought in one. A limit is not: biasing it would let it
    // settle short of its bound.
    if (p_axis.mode == AxisMode::Locked)
        r_constraint -= p_axis.initial * p_alpha;

    float F = p_axis.penalty * r_constraint + p_axis.lambda;

    // A limit may only push back inside its range, never pull further out.
    if (p_axis.mode == AxisMode::Limited)
        F = r_constraint > 0.0f ? max(F, 0.0f) : min(F, 0.0f);

    r_force = F;
    return true;
}

void GenericJoint::stampAxis(Block &block, const JointAxis &p_axis, float p_value, float p_alpha, float3 p_jLin,
        float3 p_jAng, float p_rate) const
{
    // Constraint term of the augmented Lagrangian, stamped through both Jacobian rows. The
    // cross block is what couples the body's translation and rotation, and it is not optional:
    // an anchor away from the centre of mass moves as much by turning the body as by moving it,
    // so dropping that term constrains the wrong point.
    float C;
    float F;
    if (axisForce(p_axis, p_value, p_alpha, C, F))
    {
        block.lhsLin += outer(p_jLin, p_jLin) * p_axis.penalty;
        block.lhsAng += outer(p_jAng, p_jAng) * p_axis.penalty;
        block.lhsCross += outer(p_jAng, p_jLin) * p_axis.penalty;
        block.rhsLin += p_jLin * F;
        block.rhsAng += p_jAng * F;
    }

    // Spring term: a plain penalty force with no dual variable, plus velocity damping. Written
    // in the same sign convention as the Spring force, so that it pulls towards equilibrium.
    if (p_axis.springStiffness > 0.0f)
    {
        const float Cs = p_value - p_axis.springEquilibrium;
        const float Fs = p_axis.springStiffness * Cs - p_axis.springDamping * p_rate;
        block.lhsLin += outer(p_jLin, p_jLin) * p_axis.springStiffness;
        block.lhsAng += outer(p_jAng, p_jAng) * p_axis.springStiffness;
        block.lhsCross += outer(p_jAng, p_jLin) * p_axis.springStiffness;
        block.rhsLin += p_jLin * Fs;
        block.rhsAng += p_jAng * Fs;
    }
}

// How stiff a degree of freedom may become before the solver can no longer resolve it.
//
// The augmented Lagrangian grows a constraint's penalty while it is violated, which is what makes
// it converge to zero error. Left unbounded that is a trap: the constraint's own frequency is
// sqrt(penalty / effective inertia), and once it approaches the step rate the iteration stops
// converging, so the error grows, so the penalty grows further. A degree of freedom that is
// *already* satisfied has nothing to gain from a large penalty and everything to lose, and a
// body spinning fast is exactly where it shows up first.
//
// The cap is the penalty at which the constraint's frequency is still a comfortable fraction of
// the step rate, using the smaller of the two bodies' effective inertia about the axis - the
// lighter side is the one that limits the timestep.
float GenericJoint::penaltyLimit(const Rigid *p_body, float3 p_axis, bool p_angular) const
{
    // Inertia (or mass) seen along the axis by this body.
    float effective;
    if (p_angular)
    {
        const float3 local = rotate(conjugate(p_body->positionAng), p_axis);
        const float3 moment = p_body->moment;
        // n^T I n for a diagonal inertia tensor.
        effective = local.x * local.x * moment.x + local.y * local.y * moment.y + local.z * local.z * moment.z;
    }
    else
    {
        effective = p_body->mass;
    }

    if (effective <= 0.0f)
        return PENALTY_MAX;

    // omega * dt <= kMaxFrequencyRatio, i.e. penalty <= effective * k^2 / dt^2.
    constexpr float kMaxFrequencyRatio = 0.25f;
    const float limit = effective * kMaxFrequencyRatio * kMaxFrequencyRatio / (solver->dt * solver->dt);
    return min(limit, PENALTY_MAX);
}

bool GenericJoint::initialize()
{
    // Record where each degree of freedom starts the step, and warm-start its dual state
    // (Eq. 19), exactly as Joint does.
    float angularValue_[3];
    angularValues(angularValue_);

    const float3 offset = rotate(conjugate(frameOrientation()), anchorA() - anchorB());

    for (int i = 0; i < 3; i++)
    {
        linear[i].initial = offset[i];
        angular[i].initial = angularValue_[i];

        JointAxis *both[2] = {&linear[i], &angular[i]};
        for (JointAxis *axis : both)
        {
            axis->lambda *= solver->alpha * solver->gamma;
            axis->penalty = clamp(axis->penalty * solver->gamma, PENALTY_MIN, PENALTY_MAX);
        }
    }

    return true;
}

void GenericJoint::updatePrimal(Rigid *body, float alpha, Block &block)
{
    const quat qA = frameOrientation();
    const float3 rAWorld = rotate(qA, rA);
    const float3 rBWorld = rotate(bodyB->positionAng, rB);

    // Measurements, in the joint frame.
    const float3 offset = rotate(conjugate(qA), anchorA() - anchorB());
    float angularValue_[3];
    angularValues(angularValue_);

    // Every Jacobian flips sign with the body it is stamped into.
    const bool isA = body == bodyA;
    const float sign = isA ? 1.0f : -1.0f;
    const float3 momentArm = isA ? rAWorld : rBWorld;

    // Relative motion of the two anchors, for spring damping: the rate of change of the linear
    // measurements, so that damping always opposes the spring.
    const float3 vA = bodyA ? bodyA->velocityLin + cross(bodyA->velocityAng, rAWorld) : float3{0, 0, 0};
    const float3 vB = bodyB->velocityLin + cross(bodyB->velocityAng, rBWorld);
    const float3 relativeVelocity = vA - vB;
    const float3 wA = bodyA ? bodyA->velocityAng : float3{0, 0, 0};
    const float3 wB = bodyB->velocityAng;

    for (int i = 0; i < 3; i++)
    {
        // Linear axes are directions in the world, taken from A's reference frame.
        const float3 axis = rotate(qA, frameAxis(i));

        // Linear: translating the body moves the anchor, and so does turning it - hence the
        // moment arm term.
        const float3 jLin = axis * sign;
        const float3 jAng = cross(momentArm, jLin);
        stampAxis(block, linear[i], offset[i], alpha, jLin, jAng, dot(relativeVelocity, axis));

        // Angular: a pure rotation about the joint axis.
        stampAxis(block, angular[i], angularValue_[i], alpha, float3{0, 0, 0}, axis * sign,
                dot(wA - wB, axis));
    }
}

void GenericJoint::updateDual(float alpha)
{
    float angularValue_[3];
    angularValues(angularValue_);

    const float3 offset = rotate(conjugate(frameOrientation()), anchorA() - anchorB());

    for (int i = 0; i < 3; i++)
    {
        JointAxis *both[2] = {&linear[i], &angular[i]};
        const float values[2] = {offset[i], angularValue_[i]};
        const float betas[2] = {solver->betaLin, solver->betaAng};

        for (int k = 0; k < 2; k++)
        {
            JointAxis &axis = *both[k];
            float C;
            float F;
            if (!axisForce(axis, values[k], alpha, C, F))
            {
                // Complementary slackness: a limit that is not being pushed against applies no
                // force, so its multiplier is zero. Keeping the previous one would re-apply it
                // the moment the axis came back to the bound, which reads as a violent kick -
                // and the stiffer the axis, the harder the kick, so a stiffer limit would
                // behave *worse* than a soft one.
                if (axis.mode == AxisMode::Limited)
                    axis.lambda = 0.0f;
                continue;
            }

            axis.lambda = F;

            // Grow the penalty only up to what the solver can actually resolve; see penaltyLimit.
            const float3 axisDir = rotate(frameOrientation(), frameAxis(i));
            float cap = penaltyLimit(bodyB, axisDir, k == 1);
            if (bodyA != 0) cap = min(cap, penaltyLimit(bodyA, axisDir, k == 1));
            axis.penalty = min(axis.penalty + abs(C) * betas[k], min(cap, PENALTY_MAX));
        }
    }
}

} // namespace avbd
