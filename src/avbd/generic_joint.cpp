// General 6-DOF joint; AL structure follows avbd-demo3d (MIT). DOFs Free/Locked/Limited (+spring),
// measured in A's frame (world if bodyA null). Swing-twist split survives spinning axles; off-centre anchors stamp both halves.

#include <algorithm>
#include <cmath>

#include "avbd/solver.h"
#include "avbd/maths.h"

namespace avbd {

namespace {

const quat kIdentity = quatIdentity;

const float3 &frameAxis(int p_axis)
{
    static const std::array<float3, 3> axes = {float3(1, 0, 0), float3(0, 1, 0), float3(0, 0, 1)};
    return axes[p_axis];
}

} // namespace

GenericJoint::GenericJoint(Solver *p_solver, Rigid *p_bodyA, Rigid *p_bodyB, float3 p_rA, float3 p_rB)
    : Force(p_solver, p_bodyA, p_bodyB), rA(p_rA), rB(p_rB)
{
    // The arrangement the joint is created in is its zero. Held as B-relative-to-A in A's frame,
    // the same form angularValues() computes for the live pose, so the two compose directly.
    rest = qmul(qconj(frameOrientation()), p_bodyB->positionAng);
}

quat GenericJoint::frameOrientation() const
{
    return bodyA ? bodyA->positionAng : kIdentity;
}



float3 GenericJoint::anchorA() const
{
    return bodyA ? bodyA->positionLin + qrot(bodyA->positionAng, rA) : rA;
}

float3 GenericJoint::anchorB() const
{
    return bodyB->positionLin + qrot(bodyB->positionAng, rB);
}

// Remove twist about `p_axis` from `p_dev` into `r_angle`: (n*(v.n), w) normalised, angle
// 2*atan2(v.n, w) - continuous through a full turn (spinning axles); degenerate near half turn means no twist.
static void removeTwist(quat &p_dev, int p_axis, float &r_angle)
{
    const float3 &n = frameAxis(p_axis);
    const float3 vec_part(p_dev.x, p_dev.y, p_dev.z);
    const float proj = dot(vec_part, n);

    // Pure twist quaternion: vector part is the projection onto n, w unchanged. linalg
    // stores quaternions as {x, y, z, w}.
    const quat twist{p_dev.x * proj, p_dev.y * proj, p_dev.z * proj, p_dev.w};
    if (length(twist) <= 1.0e-12f)
    {
        r_angle = 0.0f;
        return;
    }

    r_angle = 2.0f * std::atan2(proj, static_cast<float>(p_dev.w));
    p_dev = qmul(normalize(p_dev), qinv(twist));
}

void GenericJoint::angularValues(float p_out[3]) const
{
    // deviation = rest^-1 * (B-relative-to-A now), identity at rest; rest-term placement is CRITICAL:
    // Harmless only while rest is identity; once rotated, locked-axis error grows with spin and tears the joint apart.
    const quat qA = frameOrientation();
    const quat relativeNow = qmul(qconj(qA), normalize(bodyB->positionAng));
    quat dev = qmul(normalize(rest), qconj(relativeNow));

    // Each non-locked axis owns its twist (spring/limit reads it); locked axes must not see it.
    // Removing one twist is exact; several are removed one axis at a time (usual approximation).
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
            p_out[i] = static_cast<float>(dev[i] * 2.0f);
    }
}

float GenericJoint::linearValue(int p_axis) const
{
    const float3 offset = qrot(qconj(frameOrientation()), anchorA() - anchorB());
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

// One source of truth for the constraint value: the dual pass stores its force as the next multiplier
// and grows the penalty from the same error; disagreeing passes converge to the wrong answer, not just slowly.
bool GenericJoint::axisForce(const JointAxis &p_axis, float p_value, float p_alpha, float &r_constraint,
        float &r_force) const
{
    if (!axisConstraint(p_axis, p_value, r_constraint))
        return false;

    // A lock is stabilised against the step's start (error spread over steps); a limit is not,
    // since biasing would let it settle short of its bound.
    if (p_axis.mode == AxisMode::Locked)
        r_constraint -= p_axis.initial * p_alpha;

    float F = p_axis.penalty * r_constraint + p_axis.lambda;

    // A limit may only push back inside its range, never pull further out.
    if (p_axis.mode == AxisMode::Limited)
        F = r_constraint > 0.0f ? std::max(F, 0.0f) : std::min(F, 0.0f);

    r_force = F;
    return true;
}

void GenericJoint::stampAxis(Block &block, const JointAxis &p_axis, float p_value, float p_alpha, float3 p_jLin,
        float3 p_jAng, float p_rate) const
{
    // Constraint term stamped through both Jacobian rows; the cross block couples translation and
    // rotation and is not optional - an off-centre anchor turns as much as it translates.
    float C;
    float F;
    if (axisForce(p_axis, p_value, p_alpha, C, F))
    {
        block.lhsLin += outerprod(p_jLin, p_jLin) * p_axis.penalty;
        block.lhsAng += outerprod(p_jAng, p_jAng) * p_axis.penalty;
        block.lhsCross += outerprod(p_jAng, p_jLin) * p_axis.penalty;
        block.rhsLin += p_jLin * F;
        block.rhsAng += p_jAng * F;
    }

    // Spring term: a plain penalty force with no dual variable, plus velocity damping. Written
    // in the same sign convention as the Spring force, so that it pulls towards equilibrium.
    if (p_axis.springStiffness > 0.0f)
    {
        const float Cs = p_value - p_axis.springEquilibrium;
        const float Fs = p_axis.springStiffness * Cs - p_axis.springDamping * p_rate;
        block.lhsLin += outerprod(p_jLin, p_jLin) * p_axis.springStiffness;
        block.lhsAng += outerprod(p_jAng, p_jAng) * p_axis.springStiffness;
        block.lhsCross += outerprod(p_jAng, p_jLin) * p_axis.springStiffness;
        block.rhsLin += p_jLin * Fs;
        block.rhsAng += p_jAng * Fs;
    }
}

// Cap penalty before sqrt(k/I) nears the step rate: convergence stops, then error and penalty run away
// (a fast-spinning body first); use the smaller body's effective inertia - the lighter side limits dt.
float GenericJoint::penaltyLimit(const Rigid *p_body, float3 p_axis, bool p_angular) const
{
    // Inertia (or mass) seen along the axis by this body.
    float effective;
    if (p_angular)
    {
        const float3 local = qrot(qconj(p_body->positionAng), p_axis);
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
    return std::min(limit, PENALTY_MAX);
}

bool GenericJoint::initialize()
{
    // Record where each degree of freedom starts the step, and warm-start its dual state
    // (Eq. 19), exactly as Joint does.
    std::array<float, 3> angularValue_;
    angularValues(angularValue_.data());

    const float3 offset = qrot(qconj(frameOrientation()), anchorA() - anchorB());

    for (int i = 0; i < 3; i++)
    {
        linear[i].initial = offset[i];
        angular[i].initial = angularValue_[i];

        JointAxis *both[2] = {&linear[i], &angular[i]};
        for (JointAxis *axis : both)
        {
            axis->lambda *= solver->alpha * solver->gamma;
            axis->penalty = std::clamp(axis->penalty * solver->gamma, PENALTY_MIN, PENALTY_MAX);
        }
    }

    return true;
}

void GenericJoint::updatePrimal(Rigid *body, float alpha, Block &block)
{
    const quat qA = frameOrientation();
    const float3 rAWorld = qrot(qA, rA);
    const float3 rBWorld = qrot(bodyB->positionAng, rB);

    // Measurements, in the joint frame.
    const float3 offset = qrot(qconj(qA), anchorA() - anchorB());
    std::array<float, 3> angularValue_;
    angularValues(angularValue_.data());

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
        const float3 axis = qrot(qA, frameAxis(i));

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
    std::array<float, 3> angularValue_;
    angularValues(angularValue_.data());

    const float3 offset = qrot(qconj(frameOrientation()), anchorA() - anchorB());

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
                // Complementary slackness: an inactive limit applies no force, so lambda = 0;
                // keeping it re-applies as a violent kick on return, harsher for stiffer axes.
                if (axis.mode == AxisMode::Limited)
                    axis.lambda = 0.0f;
                continue;
            }

            axis.lambda = F;

            // Grow the penalty only up to what the solver can actually resolve; see penaltyLimit.
            const float3 axisDir = qrot(frameOrientation(), frameAxis(i));
            float cap = penaltyLimit(bodyB, axisDir, k == 1);
            if (bodyA != 0) cap = std::min(cap, penaltyLimit(bodyA, axisDir, k == 1));
            axis.penalty = std::clamp(axis.penalty + std::abs(C) * betas[k], 0.0f, cap);
        }
    }
}

} // namespace avbd
