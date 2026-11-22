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

#include <algorithm>
#include <cmath>
#include <array>

#include "avbd/solver.h"

namespace avbd {

uint32_t ManifoldStore::acquire(Rigid *p_bodyA, Rigid *p_bodyB)
{
    uint32_t slot;
    if (!freeList.empty())
    {
        slot = freeList.back();
        freeList.pop_back();
    }
    else
    {
        slot = static_cast<uint32_t>(bodyA.size());
        bodyA.push_back(nullptr);
        bodyB.push_back(nullptr);
        basis.push_back(float3x3{});
        numContacts.push_back(0);
        friction.push_back(0.0f);
        flags.push_back(0);
        dormantTicks.push_back(0);
        dormantSavedCount.push_back(0);
        rotACaptured.push_back(float3x3{});
        rotBCaptured.push_back(float3x3{});
        contacts.resize(contacts.size() + kMaxManifoldContacts);
        contactC0.resize(contactC0.size() + kMaxManifoldContacts);
        contactPenalty.resize(contactPenalty.size() + kMaxManifoldContacts);
        contactLambda.resize(contactLambda.size() + kMaxManifoldContacts);
        contactForceCache.resize(contactForceCache.size() + kMaxManifoldContacts);
        rAWorld.resize(rAWorld.size() + kMaxManifoldContacts);
        rBWorld.resize(rBWorld.size() + kMaxManifoldContacts);
        contactStick.resize(contactStick.size() + kMaxManifoldContacts);
        contactBetaScale.resize(contactBetaScale.size() + kMaxManifoldContacts);
        linHessian.resize(linHessian.size() + kMaxManifoldContacts);
    }

    // A recycled slot carries the previous manifold's state; every field the solver reads
    // before writing it again is reset here.
    bodyA[slot] = p_bodyA;
    bodyB[slot] = p_bodyB;
    basis[slot] = float3x3{};
    numContacts[slot] = 0;
    friction[slot] = 0.0f;
    flags[slot] = 0;
    dormantTicks[slot] = 0;
    dormantSavedCount[slot] = 0;
    rotACaptured[slot] = float3x3{};
    rotBCaptured[slot] = float3x3{};

    const int first = base(slot);
    for (int k = 0; k < kMaxManifoldContacts; k++)
    {
        contacts[first + k] = Manifold::Contact{};
        contactC0[first + k] = float3{};
        contactPenalty[first + k] = float3{};
        contactLambda[first + k] = float3{};
        contactForceCache[first + k] = float3{};
        rAWorld[first + k] = float3{};
        rBWorld[first + k] = float3{};
        contactStick[first + k] = 0;
        contactBetaScale[first + k] = 1.0f;
        linHessian[first + k] = Manifold::BasisHessian{};
    }

    live.push_back(slot);
    return slot;
}

void ManifoldStore::release(uint32_t slot) noexcept
{
    freeList.push_back(slot);
}

void ManifoldStore::clear()
{
    bodyA.clear();
    bodyB.clear();
    basis.clear();
    numContacts.clear();
    friction.clear();
    flags.clear();
    dormantTicks.clear();
    dormantSavedCount.clear();
    rotACaptured.clear();
    rotBCaptured.clear();
    contacts.clear();
    contactC0.clear();
    contactPenalty.clear();
    contactLambda.clear();
    contactForceCache.clear();
    rAWorld.clear();
    rBWorld.clear();
    contactStick.clear();
    contactBetaScale.clear();
    linHessian.clear();
    live.clear();
    freeList.clear();
}

void ForceAdjacency::clear()
{
    start.clear();
    refs.clear();
    cursor.clear();
}

void ForceAdjacency::build(std::span<Rigid *const> bodies, const Force *forceHead,
        const ManifoldStore &manifolds)
{
    const int count = static_cast<int>(bodies.size());
    start.assign(static_cast<size_t>(count) + 1, 0);

    const auto addBody = [this](Rigid *body) {
        if (body != nullptr && body->broadIndex >= 0)
            start[static_cast<size_t>(body->broadIndex) + 1]++;
    };
    for (const Force *f = forceHead; f != nullptr; f = f->next)
    {
        addBody(f->bodyA);
        addBody(f->bodyB);
    }
    for (uint32_t slot : manifolds.live)
    {
        addBody(manifolds.bodyA[slot]);
        addBody(manifolds.bodyB[slot]);
    }

    for (int i = 0; i < count; i++)
        start[static_cast<size_t>(i) + 1] += start[static_cast<size_t>(i)];

    refs.resize(start[static_cast<size_t>(count)]);
    cursor.assign(start.begin(), start.begin() + count);

    const auto push = [this](Rigid *body, ForceRef ref) {
        if (body != nullptr && body->broadIndex >= 0)
            refs[cursor[static_cast<size_t>(body->broadIndex)]++] = ref;
    };
    for (const Force *f = forceHead; f != nullptr; f = f->next)
    {
        const ForceRef ref = ForceRef::force(const_cast<Force *>(f));
        push(f->bodyA, ref);
        push(f->bodyB, ref);
    }
    for (uint32_t slot : manifolds.live)
    {
        const ForceRef ref = ForceRef::manifold(slot);
        push(manifolds.bodyA[slot], ref);
        push(manifolds.bodyB[slot], ref);
    }
}

bool ManifoldStore::frozen(uint32_t slot) const noexcept
{
    const auto awake = [](const Rigid *b) { return b != nullptr && b->mass > 0.0f && !b->sleeping; };
    return !awake(bodyA[slot]) && !awake(bodyB[slot]);
}

namespace {
// S = B^T K B (six symmetric entries) = weighted sum of basis-row outer products, K diagonal.
// Expression order matches the old per-stamp code exactly: caching does not change the result.
inline Manifold::BasisHessian basisHessian(const float3x3 &basis, const float3 &k) noexcept
{
    const float3 b0 = row(basis, 0);
    const float3 b1 = row(basis, 1);
    const float3 b2 = row(basis, 2);
    const float b0x2 = b0.x * b0.x, b0y2 = b0.y * b0.y, b0z2 = b0.z * b0.z;
    const float b0xy = b0.x * b0.y, b0xz = b0.x * b0.z, b0yz = b0.y * b0.z;
    const float b1x2 = b1.x * b1.x, b1y2 = b1.y * b1.y, b1z2 = b1.z * b1.z;
    const float b1xy = b1.x * b1.y, b1xz = b1.x * b1.z, b1yz = b1.y * b1.z;
    const float b2x2 = b2.x * b2.x, b2y2 = b2.y * b2.y, b2z2 = b2.z * b2.z;
    const float b2xy = b2.x * b2.y, b2xz = b2.x * b2.z, b2yz = b2.y * b2.z;
    return {k.x * b0x2 + k.y * b1x2 + k.z * b2x2,
            k.x * b0xy + k.y * b1xy + k.z * b2xy,
            k.x * b0xz + k.y * b1xz + k.z * b2xz,
            k.x * b0y2 + k.y * b1y2 + k.z * b2y2,
            k.x * b0yz + k.y * b1yz + k.z * b2yz,
            k.x * b0z2 + k.y * b1z2 + k.z * b2z2};
}

// C = (1-alpha)*step-start error + pose bracket; penetration damped to avoid deep-impact overshoot.
// Separated speculative contacts MUST use the full gap (1% damping deepens jams); margin off = bit-for-bit original.
inline float3 constraintError(const float3 &c0, const float3 &bracket, const float3x3 &basis,
        float oneMinusAlpha, bool speculative) noexcept
{
    const bool separated = speculative && c0.x > COLLISION_MARGIN + CONTACT_RELEASE_GAP;
    return c0 * (separated ? 1.0f : oneMinusAlpha) + mul(basis, bracket);
}

// Cap penetration on the NORMAL row only: scaling the whole vector by its norm cut 9 cm to 0.65 cm - the pair sank.
// Tangential rows are friction error (anchor offsets) needing no cap; the friction cone limits their force.
inline void clampPenetration(float3 &c0, float cap) noexcept
{
    if (c0.x < -cap)
        c0.x = -cap;
}

// What the friction-cone clamp decided, so callers that need the pre-clamp tangential
// magnitude (the dual's stick test) do not recompute it.
struct ForceClamp
{
    float bounds;
    float scaleBeforeClamp;
};

// Normal force one-sided, tangential inside the cone; keep expressions and their order verbatim - bit-exact.
inline ForceClamp clampContactForce(float3 &F, float mu) noexcept
{
    F[0] = std::min(F[0], 0.0f);

    const float bounds = std::fabs(F[0]) * mu;
    const float frictionScale = length(float2{F[1], F[2]});
    if (frictionScale > bounds && frictionScale > 0)
    {
        F[1] *= bounds / frictionScale;
        F[2] *= bounds / frictionScale;
    }
    return {bounds, frictionScale};
}
} // namespace

bool ManifoldStore::initialize(uint32_t slot, Solver &solver)
{
    Rigid *a = this->bodyA[slot];
    Rigid *b = this->bodyB[slot];
    const int first = base(slot);

    // Compute friction
    friction[slot] = std::sqrt(a->friction * b->friction);

    // New contact set: cached moment arms describe the old geometry, drop them.
    flags[slot] = static_cast<uint8_t>(flags[slot] & ~(kArmsCached | kPrimalCacheReady));

    // Grace revival: a dormant manifold re-enters with the contact set it preserved, so
    // the feature-key merge below can match the returning contacts against their history.
    if (numContacts[slot] == 0 && dormantSavedCount[slot] > 0)
        numContacts[slot] = static_cast<uint8_t>(dormantSavedCount[slot]);
    const int liveCount = numContacts[slot];

    // Compute new contacts
    std::array<Manifold::Contact, kMaxManifoldContacts> newContacts = {};
    int newNumContacts = Manifold::collide(a, b, std::span<Manifold::Contact>(newContacts),
            basis[slot]);

    if (solver.traceMatch(a, b))
    {
        const float3 d = (a->velocityLin - b->velocityLin) * solver.collisionDt;
        std::printf("TRACE narrow slot=%u new=%d prev=%d gap=%.4f relV=%.2f margin=%.3f\n", slot,
                newNumContacts, numContacts[slot],
                length(a->positionLin - b->positionLin), length(d), solver.speculativeMargin);
    }

    // Hot state carried from the matched previous contacts, indexed by the new contact.
    std::array<float3, kMaxManifoldContacts> newPenalty{};
    std::array<float3, kMaxManifoldContacts> newLambda{};
    std::array<uint8_t, kMaxManifoldContacts> newStick{};
    std::array<float, kMaxManifoldContacts> newBetaScale;
    newBetaScale.fill(1.0f);

    // Merge old contact data with new contacts. prevErr[i] carries the matched contact's
    // error from the previous step so the diagnostics below can count fresh contacts.
    std::array<float, kMaxManifoldContacts> prevErr{};
    prevErr.fill(-1.0f);
    for (int i = 0; i < newNumContacts; i++)
    {
        for (int j = 0; j < numContacts[slot]; j++)
        {
            if (newContacts[i].feature.key == contacts[first + j].feature.key)
            {
                const float3 newRA = newContacts[i].rA;
                const float3 newRB = newContacts[i].rB;
                newContacts[i] = contacts[first + j];
                // Residual penetration depth, not |C0|: C0 carries the collision margin,
                // so a resting contact would otherwise look like a constant error.
                prevErr[i] = std::max(0.0f, COLLISION_MARGIN - contactC0[first + j].x);
                newPenalty[i] = contactPenalty[first + j];
                newLambda[i] = contactLambda[first + j];
                newStick[i] = contactStick[first + j];
                newBetaScale[i] = contactBetaScale[first + j];

                // If no static friction in last frame, use the new contact point locations
                if (!contactStick[first + j])
                {
                    newContacts[i].rA = newRA;
                    newContacts[i].rB = newRB;
                }
                else
                {
                    // Stick heal: glued anchors hide sliding from STICK_THRESH; compare to fresh geometry.
                    // Past stickDriftHeal: phantom C0 -> real force; re-anchor, drop force, keep penalty ramp.
                    const float3 oldA = a->positionLin + qrot(a->positionAng,
                            contacts[first + j].rA);
                    const float3 oldB = b->positionLin + qrot(b->positionAng,
                            contacts[first + j].rB);
                    const float3 newA = a->positionLin + qrot(a->positionAng, newRA);
                    const float3 newB = b->positionLin + qrot(b->positionAng, newRB);
                    const float drift = std::max(length(oldA - newA), length(oldB - newB));
#ifdef AVBD_DIAG
                    if (drift > solver.stickDriftHeal)
                        solver.diag.stickHeals++;
#endif
                    if (drift > solver.stickDriftHeal)
                    {
                        newContacts[i].rA = newRA;
                        newContacts[i].rB = newRB;
                        newStick[i] = 0;
                        // Keep normal support: press-based, not anchor-based; dropping it made the pile re-sink.
                        // Only friction is tied to the anchors that moved.
                        newLambda[i] = float3{newLambda[i].x, 0.0f, 0.0f};
                    }
                }
                break;
            }
        }
    }

#ifdef AVBD_DIAG
    for (int i = 0; i < newNumContacts; i++)
        if (prevErr[i] < 0.0f)
            solver.diag.newContacts++;
#endif

    // Fresh contact (unmatched feature key) seeds penalty at a fraction of betaLin, so impacts don't ramp from k=1.
    // Gated on grace ticks; default 0 keeps the historical floor start.
    if (solver.contactGraceTicks > 0)
    {
        const float seed = std::max(PENALTY_MIN, 0.05f * solver.betaLin);
        for (int i = 0; i < newNumContacts; i++)
            if (prevErr[i] < 0.0f)
                newPenalty[i] = float3{seed, seed, seed};
    }

    // Copy new contacts and their hot state to the slot
    numContacts[slot] = static_cast<uint8_t>(newNumContacts);
    for (int i = 0; i < newNumContacts; i++)
    {
        contacts[first + i] = newContacts[i];
        contactPenalty[first + i] = newPenalty[i];
        contactLambda[first + i] = newLambda[i];
        contactStick[first + i] = newStick[i];
        contactBetaScale[first + i] = newBetaScale[i];
    }

    // Grace dormancy: parking (not reaping) keeps hot state so a return within grace ticks re-matches by feature key.
    // numContacts stays 0, so every solve loop over this manifold is a no-op.
    if (numContacts[slot] == 0)
    {
        if (solver.contactGraceTicks > 0 && dormantTicks[slot] < solver.contactGraceTicks)
        {
            // Park stiffness, not force: old lambda supports a dead pose; reviving with it shoves a box several m/s.
            // The live-path release rule is the same: a separated pair may carry no force.
            for (int i = 0; i < liveCount; i++)
            {
                contactLambda[first + i] = float3{};
                contactStick[first + i] = 0;
            }
            dormantTicks[slot]++;
            dormantSavedCount[slot] = liveCount;
            return false;
        }
        dormantSavedCount[slot] = 0;
    }
    dormantTicks[slot] = 0;

    // Compute error at q- and update penalty and lambdas
    for (int i = 0; i < numContacts[slot]; i++)
    {
        // Error at q-
        float3 xA = a->positionLin + qrot(a->positionAng, contacts[first + i].rA);
        float3 xB = b->positionLin + qrot(b->positionAng, contacts[first + i].rB);
        contactC0[first + i] = mul(basis[slot], xA - xB) + float3{COLLISION_MARGIN, 0, 0};

        // Cap depth: BDF1 v=dp/dt means an unclamped 10 cm overlap expels at ~24 m/s, re-ejecting neighbours.
        // Measured: double pyramid at 0.6 m overlaps, boxes left at 18 m/s; residual drains over later steps.
        clampPenetration(contactC0[first + i], solver.maxPenetration);

        // Adaptive beta: non-shrinking penetration doubles the ramp, converging halves toward 1; fresh starts neutral.
        // Residual = penetration depth, not margin; <= floor = solved so float noise can't ratchet to cap.
        if (solver.adaptiveBeta && prevErr[i] >= 0.0f)
        {
            const float err = std::max(0.0f, COLLISION_MARGIN - contactC0[first + i].x);
            const float floor = solver.adaptiveBetaErrorFloor;
            float scale = newBetaScale[i];
            if (err <= floor)
                scale = std::max(scale * 0.5f, 1.0f);
            else if (prevErr[i] <= floor)
                scale = std::min(scale * 2.0f, ADAPTIVE_BETA_MAX);
            else if (err > 0.5f * prevErr[i])
                scale = std::min(scale * 2.0f, ADAPTIVE_BETA_MAX);
            else if (err < 0.1f * prevErr[i])
                scale = std::max(scale * 0.5f, 1.0f);
            contactBetaScale[first + i] = scale;
        }
        else
        {
            contactBetaScale[first + i] = 1.0f;
        }

        // Contact release: apart means no force; stale impulse must not cross the gap - decay alone lasts ages.
        // Penalty restamps it, popping a box at 10+ m/s; stiffness survives, so support rebuilds in one round.
        if (contactC0[first + i][0] > COLLISION_MARGIN + CONTACT_RELEASE_GAP)
        {
            contactLambda[first + i] = float3{};
            contactStick[first + i] = 0;
        }

        // Initial penalty (Eq. 19): old log10/pow pen^2 always floored (pen<=0.02 -> 4e-4 < MIN), dead ~300 us/step.
        // Starts at PENALTY_MIN, ramped per updateDual; warmstart decay below kept bit-for-bit.
        contactLambda[first + i] = contactLambda[first + i] * solver.alpha * solver.gamma;
        // Warmstart penalty (Eq. 19): floor resets discarded the ramp, so resting bodies sank until the dual caught up.
        // Reference decays/clamps the carried value; fresh contacts start at floor (value-initialised).
        contactPenalty[first + i] = clamp(contactPenalty[first + i] * solver.gamma,
                float3{PENALTY_MIN, PENALTY_MIN, PENALTY_MIN},
                float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});

        linHessian[first + i] = basisHessian(basis[slot], contactPenalty[first + i]);
    }

    return numContacts[slot] > 0;
}

void ManifoldStore::decayPenalty(uint32_t slot, float decay) noexcept
{
    const int first = base(slot);
    const float3 floor{PENALTY_MIN, PENALTY_MIN, PENALTY_MIN};
    const float3 cap{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX};
    for (int i = 0; i < numContacts[slot]; i++)
        contactPenalty[first + i] = clamp(contactPenalty[first + i] * decay, floor, cap);
}

// Recompute C0 at the current pose, keeping the contact set (rA/rB/basis) from the collision pass.
// The force cache belongs to the old C0, so invalidate it; the next primal touch precomputes against the rebased error.
void ManifoldStore::rebase(uint32_t slot, float maxPenetration) noexcept
{
    Rigid *a = this->bodyA[slot];
    Rigid *b = this->bodyB[slot];
    const int first = base(slot);
    for (int i = 0; i < numContacts[slot]; i++)
    {
        const float3 xA = a->positionLin + qrot(a->positionAng, contacts[first + i].rA);
        const float3 xB = b->positionLin + qrot(b->positionAng, contacts[first + i].rB);
        contactC0[first + i] = mul(basis[slot], xA - xB) + float3{COLLISION_MARGIN, 0, 0};
        clampPenetration(contactC0[first + i], maxPenetration);
    }
    flags[slot] = static_cast<uint8_t>(flags[slot] & ~kPrimalCacheReady);
}

// Primal split: preCompute (per force, once/round) caches the shared half; stampCached (per body) stamps each side.
// Same inputs in the same order, so cached values are bit-identical to the old recompute.
void ManifoldStore::preCompute(uint32_t slot, const Solver &solver, float alpha) noexcept
{
    Rigid *a = this->bodyA[slot];
    Rigid *b = this->bodyB[slot];
    const int first = base(slot);

    const float3 dqALin = a->positionLin - a->initialLin;
    const float3 dqAAng = a->positionAng - a->initialAng;
    const float3 dqBLin = b->positionLin - b->initialLin;
    const float3 dqBAng = b->positionAng - b->initialAng;

    // Same rotation-matrix hoist as updateDual: body quats are loop-invariant here.
    const float3x3 rotA = qmat(a->positionAng);
    const float3x3 rotB = qmat(b->positionAng);
    const bool speculative = solver.speculativeMargin > 0.0f;

    // Constraint in basis frame: Jlin=±basis rows, Jang=rWorld x basis row; four old products collapse to one bracket:
    //   C = C0(1-alpha) + B*[(dqA_lin-dqB_lin) + (dqA_ang x rA) - (dqB_ang x rB)]; basis is the only frame.
    for (int i = 0; i < numContacts[slot]; i++)
    {
        const float3 rA = mul(rotA, contacts[first + i].rA);
        const float3 rB = mul(rotB, contacts[first + i].rB);
        rAWorld[first + i] = rA;
        rBWorld[first + i] = rB;

        const float3 bracket = (dqALin - dqBLin) + cross(dqAAng, rA) - cross(dqBAng, rB);
        const float3 C = constraintError(contactC0[first + i], bracket, basis[slot], 1 - alpha,
                speculative);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = cmul(contactPenalty[first + i], C) + contactLambda[first + i];

        clampContactForce(F, friction[slot]);
        contactForceCache[first + i] = F;
    }
    flags[slot] = static_cast<uint8_t>(flags[slot] | kPrimalCacheReady);
}

void ManifoldStore::stampCached(uint32_t slot, Rigid *body, Block &block) const noexcept
{
    const bool fromIsA = body == bodyA[slot];
    const float side = fromIsA ? 1.0f : -1.0f;
    const int first = base(slot);
    const float3x3 &frame = this->basis[slot];

    // S=B^T K B (K diagonal): sum of basis-row outer products; R=skew(rWorld), Jang=±B R^T -> closed form:
    //   lhsLin=S, lhsAng=R S R^T, lhsCross=R S, rhsLin=±B^T F, rhsAng=±R(B^T F); mirror symmetric for pivot fallback.
    const float3 b0 = row(frame, 0);
    const float3 b1 = row(frame, 1);
    const float3 b2 = row(frame, 2);

    // S is refreshed when the penalty changes, so the second body's stamp need not recompute it.
    // Accumulate in registers, then write back once.
    float3x3 accLin{};
    float3x3 accAng{};
    float3x3 accCross{};
    float3 accRhsLin{};
    float3 accRhsAng{};

    // A zero cached force is not inactive: normal clamped at zero while compressed, penalty still in Hessian.
    // Skipping it removes support stiffness, letting loaded stacks sink and then launch.
    for (int i = 0; i < numContacts[slot]; i++)
    {
        const float3 &F = contactForceCache[first + i];
        const float3 &r = fromIsA ? rAWorld[first + i] : rBWorld[first + i];

        const Manifold::BasisHessian &s = linHessian[first + i];
        const float s00 = s.s00, s01 = s.s01, s02 = s.s02;
        const float s11 = s.s11, s12 = s.s12, s22 = s.s22;

        // Column-major: column j of S is (S[0][j], S[1][j], S[2][j]).
        const float3x3 S{float3{s00, s01, s02}, float3{s01, s11, s12}, float3{s02, s12, s22}};

        // R S, one cross product per column of S.
        const float3x3 RS{cross(r, S.x), cross(r, S.y), cross(r, S.z)};

        // Precompute second-order products of r.
        const float rx2 = r.x * r.x, ry2 = r.y * r.y, rz2 = r.z * r.z;
        const float rxy = r.x * r.y, rxz = r.x * r.z, ryz = r.y * r.z;

        // R S R^T, symmetric. Closed form in r and S (verified against the matrix
        // product); h01 = H[0][1], h02 = H[0][2], h12 = H[1][2].
        const float h00 = s11 * rz2 - 2.0f * s12 * ryz + s22 * ry2;
        const float h11 = s22 * rx2 - 2.0f * s02 * rxz + s00 * rz2;
        const float h22 = s00 * ry2 - 2.0f * s01 * rxy + s11 * rx2;
        const float h01 = s12 * rxz + s02 * ryz - s22 * rxy - s01 * rz2;
        const float h02 = s12 * rxy + s01 * ryz - s11 * rxz - s02 * ry2;
        const float h12 = s02 * rxy + s01 * rxz - s00 * ryz - s12 * rx2;

        // B^T F and its moment about the contact arm: the only per-body part of the RHS.
        const float3 fw = b0 * F.x + b1 * F.y + b2 * F.z;
        const float3 fwAng = cross(r, fw);

        accLin += S;
        accAng += float3x3{float3{h00, h01, h02}, float3{h01, h11, h12}, float3{h02, h12, h22}};
        accCross += RS;
        accRhsLin += fw;
        accRhsAng += fwAng;
    }

    block.lhsLin += accLin;
    block.lhsAng += accAng;
    block.lhsCross += accCross;
    block.rhsLin += accRhsLin * side;
    block.rhsAng += accRhsAng * side;
}

void ManifoldStore::updatePrimal(uint32_t slot, const Solver &solver, Rigid *body, float alpha,
        Block &block)
{
    // First touch consumes the previous dual's cache: no pose moved yet, so recomputing gives the same floats.
    // The second body's colour runs after poses moved, so it recomputes - the round's one fresh evaluation.
    if ((flags[slot] & kPrimalCacheReady) == 0)
        preCompute(slot, solver, alpha);
    stampCached(slot, body, block);
    flags[slot] = static_cast<uint8_t>(flags[slot] & ~kPrimalCacheReady);
}

void ManifoldStore::updateDual(uint32_t slot, const Solver &solver, float alpha) noexcept
{
    Rigid *a = this->bodyA[slot];
    Rigid *b = this->bodyB[slot];
    const int first = base(slot);
    const float3x3 &frame = this->basis[slot];
    const float mu = this->friction[slot];

    const float3 dqALin = a->positionLin - a->initialLin;
    const float3 dqAAng = a->positionAng - a->initialAng;
    const float3 dqBLin = b->positionLin - b->initialLin;
    const float3 dqBAng = b->positionAng - b->initialAng;

    // Loop-invariant quats: one qmat (~20 ns) replaces two quat multiplies/contact (~16 ns each; 8 = ~256 ns/manifold).
    const float3x3 rotA = qmat(a->positionAng);
    const float3x3 rotB = qmat(b->positionAng);

    // Drift-gated arm reuse (jacobianRebuildDistance): under the rotation threshold, reuse cached world arms.
    // Only arms are cached now; angular jacobian matrices are gone.
    const bool cacheable = solver.jacobianRebuildDistance > 0.0f;
    const float rebuild2 = solver.jacobianRebuildDistance * solver.jacobianRebuildDistance;
    const float3 dA = rotA[0] - rotACaptured[slot][0];
    const float3 dB = rotB[0] - rotBCaptured[slot][0];
    const bool driftUsable = cacheable && (flags[slot] & kArmsCached) != 0 &&
            length2(dA) <= rebuild2 && length2(dB) <= rebuild2;
    const bool speculative = solver.speculativeMargin > 0.0f;

    for (int i = 0; i < numContacts[slot]; i++)
    {
        float3 rA, rB;
        if (driftUsable)
        {
            rA = rAWorld[first + i];
            rB = rBWorld[first + i];
        }
        else
        {
            rA = mul(rotA, contacts[first + i].rA);
            rB = mul(rotB, contacts[first + i].rB);
            rAWorld[first + i] = rA;
            rBWorld[first + i] = rB;
        }

        const float3 bracket = (dqALin - dqBLin) + cross(dqAAng, rA) - cross(dqBAng, rB);
        // Same C the primal solved against: the dual must not ramp the penalty on a
        // different error, or primal and dual disagree about the contact's state.
        const float3 C = constraintError(contactC0[first + i], bracket, frame, 1 - alpha,
                speculative);

        // Compute force - the penalty is a diagonal, so K*C is componentwise.
        float3 F = cmul(contactPenalty[first + i], C) + contactLambda[first + i];

        const ForceClamp cone = clampContactForce(F, mu);

        // Store updated force
        contactLambda[first + i] = F;

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        const float beta = solver.betaLin * contactBetaScale[first + i];
        if (F[0] < 0)
            contactPenalty[first + i][0] = std::min(contactPenalty[first + i][0] + beta * std::fabs(C[0]), PENALTY_MAX);
        if (cone.scaleBeforeClamp <= cone.bounds)
        {
            contactPenalty[first + i][1] = std::min(contactPenalty[first + i][1] + beta * std::fabs(C[1]), PENALTY_MAX);
            contactPenalty[first + i][2] = std::min(contactPenalty[first + i][2] + beta * std::fabs(C[2]), PENALTY_MAX);
            contactStick[first + i] = length(float2{C[1], C[2]}) < STICK_THRESH;
        }

        // The penalty just changed, so the shared linear Hessian is stale for the next
        // round's stamps. Refresh it here, where the new penalty is already in hand.
        linHessian[first + i] = basisHessian(frame, contactPenalty[first + i]);

        // Leave the primal cache ready: pose is frozen until stampCached, so preCompute repeats this bracket/force.
        // Doing it here removes the standalone preCompute after round one.
        float3 fNext = cmul(contactPenalty[first + i], C) + contactLambda[first + i];
        clampContactForce(fNext, mu);
        contactForceCache[first + i] = fNext;
    }

    if (!driftUsable && cacheable)
    {
        rotACaptured[slot] = rotA;
        rotBCaptured[slot] = rotB;
        flags[slot] = static_cast<uint8_t>(flags[slot] | kArmsCached);
    }

    // Ready for the next round's first primal touch: nothing moves between here and it.
    flags[slot] = static_cast<uint8_t>(flags[slot] | kPrimalCacheReady);
}

// For the convergence probe: with speculative contacts a positive normal is a gap, not penetration - don't count it.
void ManifoldStore::accumulatePenetration(uint32_t slot, const Solver &solver, float &r_total,
        int &r_count) const noexcept
{
    const bool speculative = solver.speculativeMargin > 0.0f;
    const int first = base(slot);
    for (int i = 0; i < numContacts[slot]; i++)
    {
        if (speculative && contactC0[first + i][0] > 0.0f)
            continue;
        // The normal row alone: the tangential rows are anchor offsets on the friction axes
        // and would dominate the probe without saying anything about penetration.
        const float penetration = -contactC0[first + i][0];
        if (penetration > 0.0f)
        {
            r_total += penetration;
            r_count++;
        }
    }
}

} // namespace avbd
