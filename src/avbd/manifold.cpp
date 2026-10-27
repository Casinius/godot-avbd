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
// S = B^T K B in the six independent entries of a symmetric 3x3. The basis rows are the
// linear jacobian rows (up to sign) and K is diagonal, so S is a weighted sum of the
// three basis-row outer products. The expression order matches the old per-stamp code
// exactly, so caching it does not change the result.
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

// The constraint value the primal and dual work against: the step-start error damped by
// (1-alpha), plus the pose bracket in the contact basis. Penetration is damped so a deep
// impact does not overshoot; a speculative contact that is still geometrically apart is
// NOT penetrating yet, and damping its gap to 1% makes any millimetre of jitter look like
// a violation - the pair starts pushing from up to `speculativeMargin` away, an invisible
// inflation layer that measurably deepens every jam in a dense pile. Use the full gap for
// those, so the pair only acts once the bracket predicts a real overlap. With the margin
// disabled the geometry never reports a separated pair, so the branch is unreachable and
// the expression is the original one bit for bit.
inline float3 constraintError(const float3 &c0, const float3 &bracket, const float3x3 &basis,
        float oneMinusAlpha, bool speculative) noexcept
{
    const bool separated = speculative && c0.x > COLLISION_MARGIN + CONTACT_RELEASE_GAP;
    return c0 * (separated ? 1.0f : oneMinusAlpha) + mul(basis, bracket);
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

    // Hot state carried from the matched previous contacts, indexed by the new contact.
    std::array<float3, kMaxManifoldContacts> newPenalty{};
    std::array<float3, kMaxManifoldContacts> newLambda{};
    std::array<uint8_t, kMaxManifoldContacts> newStick{};

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
                prevErr[i] = length(contactC0[first + j]);
                newPenalty[i] = contactPenalty[first + j];
                newLambda[i] = contactLambda[first + j];
                newStick[i] = contactStick[first + j];

                // If no static friction in last frame, use the new contact point locations
                if (!contactStick[first + j])
                {
                    newContacts[i].rA = newRA;
                    newContacts[i].rB = newRB;
                }
                else
                {
                    // Stick heal: a sticky pair's anchors are glued to each other, so their
                    // separation (what STICK_THRESH checks every dual update) stays ~zero no
                    // matter how far the pair has slid. Compare the frozen anchors against the
                    // fresh contact geometry at the current pose instead: past STICK_DRIFT_HEAL
                    // the "static" contact is really sliding, the old anchors no longer sit on
                    // the touching surfaces, and C0 measures anchor separation - phantom
                    // penetration the dual then answers with real force. Re-anchor from the
                    // current pass and drop the force history (the penalty ramp survives, so
                    // support rebuilds in-round).
                    const float3 oldA = a->positionLin + qrot(a->positionAng,
                            contacts[first + j].rA);
                    const float3 oldB = b->positionLin + qrot(b->positionAng,
                            contacts[first + j].rB);
                    const float3 newA = a->positionLin + qrot(a->positionAng, newRA);
                    const float3 newB = b->positionLin + qrot(b->positionAng, newRB);
                    const float drift = std::max(length(oldA - newA), length(oldB - newB));
#ifdef AVBD_DIAG
                    if (drift > STICK_DRIFT_HEAL)
                        solver.diag.stickHeals++;
#endif
                    if (drift > STICK_DRIFT_HEAL)
                    {
                        newContacts[i].rA = newRA;
                        newContacts[i].rB = newRB;
                        newStick[i] = 0;
                        // Keep the normal support impulse: it is a function of how hard the
                        // pair is pressed together, not of where the anchors sit. Dropping it
                        // made every sliding contact lose its support each tick - the pile
                        // never settled, it just kept re-sinking. Only friction is tied to the
                        // anchors that just moved.
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

    // Fresh-contact seed: a contact that matched no feature key starts its penalty at a
    // fraction of betaLin instead of the floor, so an impact-time contact does not spend its
    // first ticks ramping from k=1 while the load runs through it. Gated on the grace
    // setting (default 0 keeps the historical floor start).
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
    }

    // Grace dormancy: the pair separated and no contacts came back. Parking the manifold
    // (instead of reporting it dead for warmstartForces to reap) keeps the hot state of
    // `liveCount` contacts untouched in the arrays above, so a pair that returns within
    // Solver::contactGraceTicks re-matches by feature key instead of rebuilding from the
    // penalty floor. numContacts stays 0, which makes every solve loop over this manifold
    // a no-op.
    if (numContacts[slot] == 0)
    {
        if (solver.contactGraceTicks > 0 && dormantTicks[slot] < solver.contactGraceTicks)
        {
            // Park the stiffness, not the force. What a revived pair needs is the penalty
            // ramp it earned; the lambda from before the separation is support for a pose
            // that no longer exists. Kept alive, it is applied in full by the first primal
            // round after revival - a kick of lambda*dt that a settled box feels as a
            // several-m/s shove out of the pile. The release rule in the live path exists
            // for the same reason: a separated pair may carry no force.
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

        // The reference implementation reports the full penetration depth here:
        // a deep impact then produces a proportionally large error and force, and
        // the primal pushes the body back out within the step. Capping the depth
        // (as this port once did, at 2 cm) removes the restoring gradient past the
        // cap - bodies settle anywhere below it and rebound off the stored spring
        // energy, which showed up as endless bouncing and sink-through on the
        // round-shape scenes. No clamp.

        // Contact release: the anchors say the pair is geometrically apart, so no contact
        // force can act - whatever impulse the dual still holds came from an earlier pose
        // and must not be applied across the gap. The warm-start decay alone would keep a
        // big impact impulse alive for hundreds of ticks, and the penalty term keeps
        // stamping it into the solve every round; the pile then pops a box out at 10+ m/s
        // long after it was still. The penalty (stiffness) survives, so the support
        // rebuilds in one round when the pair closes again.
        if (contactC0[first + i][0] > COLLISION_MARGIN + CONTACT_RELEASE_GAP)
        {
            contactLambda[first + i] = float3{};
            contactStick[first + i] = 0;
        }

        // Initial penalty (Eq. 19): the original computed base * 10^(2*log10(pen)) —
        // algebraically pen^2 — and clamped to [PENALTY_MIN, PENALTY_MAX]. pen is capped
        // at 0.02 above, so the product is <= 4e-4 < PENALTY_MIN: the clamp always
        // saturates at the floor and the log10/pow chain was dead computation
        // (~300 us/step across a full manifold set). Penalty starts at PENALTY_MIN and is
        // ramped per iteration in updateDual (betaLin * |C|, capped by PENALTY_MAX).
        // The warmstart decay below is kept bit-for-bit.
        contactLambda[first + i] = contactLambda[first + i] * solver.alpha * solver.gamma;
        // Warmstart the penalty too (Eq. 19). Resetting it to the floor every frame
        // threw away the ramp the previous step had already built, so a resting body
        // had to re-grow its support from PENALTY_MIN each frame and sank until the
        // dual happened to catch up. The reference implementation decays and clamps
        // the carried value instead; a fresh contact starts from the floor because
        // its slots were value-initialised.
        contactPenalty[first + i] = clamp(contactPenalty[first + i] * solver.gamma,
                float3{PENALTY_MIN, PENALTY_MIN, PENALTY_MIN},
                float3{PENALTY_MAX, PENALTY_MAX, PENALTY_MAX});

        linHessian[first + i] = basisHessian(basis[slot], contactPenalty[first + i]);
    }

    return numContacts[slot] > 0;
}

// Recompute C0 at the current pose, keeping the contact set (rA/rB/basis) from the tick's
// collision pass. The prepared force cache belongs to the old C0, so invalidate it; the
// next primal touch will precompute against the rebased error.
void ManifoldStore::rebase(uint32_t slot) noexcept
{
    Rigid *a = this->bodyA[slot];
    Rigid *b = this->bodyB[slot];
    const int first = base(slot);
    for (int i = 0; i < numContacts[slot]; i++)
    {
        const float3 xA = a->positionLin + qrot(a->positionAng, contacts[first + i].rA);
        const float3 xB = b->positionLin + qrot(b->positionAng, contacts[first + i].rB);
        contactC0[first + i] = mul(basis[slot], xA - xB) + float3{COLLISION_MARGIN, 0, 0};
    }
    flags[slot] = static_cast<uint8_t>(flags[slot] & ~kPrimalCacheReady);
}

// The primal runs in two phases so the two bodies of a manifold share one assembly:
// preCompute (parallel over forces, once per round) evaluates everything both stamps
// need - world arms, angular jacobians, and the clamped contact forces - from the
// round's pose snapshot; stampCached (parallel over bodies) stamps each body's side
// from that cache. Before the split, each of the two per-round updatePrimal calls
// recomputed the whole shared half and used only its own side of the stamps. Same
// inputs in the same order, so the cached values are bit-identical to what the stamps
// used to recompute.
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

    // Constraint value in the basis frame. The linear jacobian rows are the basis rows
    // (±), and each angular row is rWorld x that basis row, so the four matrix-vector
    // products the old assembly built collapse into one bracket the basis then rotates:
    //   C = C0(1-alpha) + B * [ (dqA_lin - dqB_lin) + (dqA_ang x rA) - (dqB_ang x rB) ]
    // No angular jacobian matrices are needed: the basis is the only frame C lives in.
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

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * friction[slot];
        float frictionScale = length(float2{F[1], F[2]});
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

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

    // Both bodies share the linear jacobian up to sign, Jlin = ±basis, so the linear
    // Hessian Jlin^T K Jlin = B^T K B is the same for A and B. K is diagonal, so that
    // "basis sandwich" S is just a weighted sum of the three basis-row outer products.
    // With R = skew(rWorld) the angular jacobian is Jang = ±B R^T, which gives the whole
    // stamp in closed form:
    //   lhsLin   = S
    //   lhsAng   = Jang^T K Jang = R S R^T
    //   lhsCross = Jang^T K Jlin = R S
    //   rhsLin   = ±B^T F
    //   rhsAng   = ±R (B^T F)
    // This replaces the three generic 3x3 products (and their transposes) the old stamp
    // built per contact, per body. The symmetric results are written mirrored so the
    // pivot fallback in Block::solveDelta still sees a full matrix.
    const float3 b0 = row(frame, 0);
    const float3 b1 = row(frame, 1);
    const float3 b2 = row(frame, 2);

    // S is refreshed by initialize()/updateDual() whenever the penalty changes, so the
    // per-body stamp does not recompute the same sandwich for the second body.
    // Accumulate in registers, then write back once.
    float3x3 accLin{};
    float3x3 accAng{};
    float3x3 accCross{};
    float3 accRhsLin{};
    float3 accRhsAng{};

    // A zero cached force does not make the contact inactive: the normal force is clamped
    // at zero while the contact can still be compressed within the margin, and its penalty
    // still belongs in the Hessian. Skipping it removes the support stiffness, which lets
    // loaded stacks sink through each other and, once deep enough, launch.
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
    // The first body to touch this manifold in a round consumes the force/jacobian cache
    // the previous round's dual left behind: no pose in this computation has moved since
    // then (the other body has not been solved yet either), so recomputing would produce
    // the same floats. The second body's colour runs after poses have moved, so it
    // recomputes - that is the round's one fresh evaluation.
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

    // Both contact rotations share the two body quats, which never change inside this
    // loop: one rotation-matrix conversion (~20 ns) replaces two quat multiplies per
    // contact (~16 ns each at 8 contacts = ~256 ns/manifold).
    const float3x3 rotA = qmat(a->positionAng);
    const float3x3 rotB = qmat(b->positionAng);

    // Drift-gated arm reuse (the jacobianRebuildDistance semantics, now living in the
    // dual): when neither body's rotation moved further than the threshold since the last
    // capture, the cached world moment arms are reused for C and for the next round's
    // stamp. Only the arms are cached now - the angular jacobian matrices are gone.
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

        // Clamp normal force
        F[0] = std::min(F[0], 0.0f);

        // Clamp norm of friction forces to achieve a friction cone
        float bounds = std::fabs(F[0]) * mu;
        float frictionScale = length(float2{F[1], F[2]});
        if (frictionScale > bounds && frictionScale > 0)
        {
            F[1] *= bounds / frictionScale;
            F[2] *= bounds / frictionScale;
        }

        // Store updated force
        contactLambda[first + i] = F;

        // Update the penalty parameter and clamp to material stiffness if we are within the force bounds (Eq. 16)
        const float beta = solver.betaLin;
        if (F[0] < 0)
            contactPenalty[first + i][0] = std::min(contactPenalty[first + i][0] + beta * std::fabs(C[0]), PENALTY_MAX);
        if (frictionScale <= bounds)
        {
            contactPenalty[first + i][1] = std::min(contactPenalty[first + i][1] + beta * std::fabs(C[1]), PENALTY_MAX);
            contactPenalty[first + i][2] = std::min(contactPenalty[first + i][2] + beta * std::fabs(C[2]), PENALTY_MAX);
            contactStick[first + i] = length(float2{C[1], C[2]}) < STICK_THRESH;
        }

        // The penalty just changed, so the shared linear Hessian is stale for the next
        // round's stamps. Refresh it here, where the new penalty is already in hand.
        linHessian[first + i] = basisHessian(frame, contactPenalty[first + i]);

        // Leave the primal's cache ready for the next Newton round. The pose is frozen
        // between here and the next stampCached (the solve runs after the next
        // preCompute would), so the next preCompute would recompute exactly this bracket
        // - and, with the ramped penalty and the just-stored lambda, exactly this force.
        // Doing it here removes the standalone preCompute from every round after the first.
        float3 fNext = cmul(contactPenalty[first + i], C) + contactLambda[first + i];
        fNext[0] = std::min(fNext[0], 0.0f);
        const float boundsNext = std::fabs(fNext[0]) * mu;
        const float frictionScaleNext = length(float2{fNext[1], fNext[2]});
        if (frictionScaleNext > boundsNext && frictionScaleNext > 0)
        {
            fNext[1] *= boundsNext / frictionScaleNext;
            fNext[2] *= boundsNext / frictionScaleNext;
        }
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

// Accumulate this manifold's penetration error for the convergence probe. With
// speculative contacts on, a positive normal error is a gap the bodies have not closed
// yet, not penetration: counting it would skew the residual probe.
void ManifoldStore::accumulatePenetration(uint32_t slot, const Solver &solver, float &r_total,
        int &r_count) const noexcept
{
    const bool speculative = solver.speculativeMargin > 0.0f;
    const int first = base(slot);
    for (int i = 0; i < numContacts[slot]; i++)
    {
        if (speculative && contactC0[first + i][0] > 0.0f)
            continue;
        const float penetration = length(contactC0[first + i]);
        if (penetration > 0.0f)
        {
            r_total += penetration;
            r_count++;
        }
    }
}

} // namespace avbd
