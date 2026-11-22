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

#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>
#include <array>
#include <limits>
#include <span>

#include "avbd/broadphase_pair.hpp"
#include "avbd/bvh_broadphase.hpp"
#include "avbd/maths.h"


#include "avbd/types.hpp"
#include "avbd/forces.hpp"
#include "avbd/manifold.hpp"

namespace avbd {

// Defined in job_pool.hpp, which is only included by solver.cpp.
namespace detail {
class JobPool;
}

// Work per item differs by orders of magnitude between the phases below, so each one gets
// its own break-even batch size. Measured on a 12-core machine: one dispatch into the pool
// costs ~45 us, a body's primal update ~1 us, one manifold's contact solve ~85 ns, and a
// dual update ~200 ns. Break-even is therefore roughly 45 bodies, 500 contact solves and
// 250 dual updates - i.e. a scene has to be a few hundred bodies before splitting a phase
// across threads pays for the dispatch itself.
constexpr int kMinBodiesPerDispatch = 48;
constexpr int kMinForcesPerDispatch = 256;
// The broad phase's per-body preparation (AABB bound, level, cell key) costs tens of
// nanoseconds per body, so a dispatch only pays off at a much larger batch than a body's
// primal update does (~1 us). Below this the loop runs inline and saves the round trip.
constexpr int kMinBodiesPerCheapDispatch = 512;

// Core solver class which holds all the rigid bodies and forces, and has logic to step the simulation forward in time
struct Solver
{
    // Parameter defaults. Tuned for metre/kilogram/second scenes; `beta*` in particular is
    // unit dependent, which is why they are exposed rather than fixed.
    float dt = 1.0f / 60.0f; // Timestep
    // Set by step/stepInterleaved to the interval until the next contact refresh.
    float collisionDt = dt;
    float gravity = -10.0f;  // Gravity, along -Z
    int iterations = 8;      // Solver iterations per step

    float alpha = 0.99f;    // Stabilization parameter, in (0, 1]
    float betaLin = 10000.0f;  // Penalty ramping for linear constraints
    float betaAng = 100.0f;    // Penalty ramping for angular constraints
    float gamma = 0.999f;      // Warmstarting decay, < 1

    // Primal rounds the current step actually runs (set by step()/stepInterleaved before
    // any constraint is touched). A constraint penalty is only as stiff as the iteration
    // budget can resolve, so the joint penalty cap reads this: the augmented-Lagrangian
    // dual gains `penalty * C` every round, and a penalty that a long iteration set
    // converges with makes a short one diverge into a launcher. See
    // GenericJoint::penaltyLimit.
    int roundsThisStep = 8;

    // Round split: the first ceil(newtonRatio * rounds) rounds run the full primal (the
    // dense 6x6 block solve); the rest skip the primal and run dual-only relaxation on the
    // existing forces, with contact stiffness decayed by `stiffnessDecay` per round so the
    // cheap rounds finish what the Newton rounds started instead of fighting a stiff ramp.
    // 1.0 / 1.0 = every round is a Newton round (historical behaviour).
    float newtonRatio = 1.0f;
    float stiffnessDecay = 1.0f;

    // Adaptive beta: a per-contact multiplier on betaLin, carried across steps with the
    // warm start. A contact whose penetration did not shrink against last step doubles its
    // multiplier (capped at ADAPTIVE_BETA_MAX); a converging one halves it back toward 1.
    // This is the augmented-Lagrangian kappa update applied per constraint: a loaded stack
    // contact stiffens until it actually carries its share, instead of every contact
    // sharing one global ramp rate.
    bool adaptiveBeta = false;
    // Residual penetration (metres) at or below which a contact counts as solved: the
    // multiplier decays instead of doubling, so float noise on a resting contact cannot
    // ratchet it to the cap.
    float adaptiveBetaErrorFloor = 1e-3f;

#ifdef AVBD_DIAG
    // Root-cause instrumentation, compiled out of the shipped library. Reset at the top
    // of every stepInterleaved/step; the driver reads it per tick.
    struct Diag {
        int actualRounds = 0;      // rounds each substep actually ran
        int softContacts = 0;      // contacts with penalty < 0.1*betaLin after warm start (H1)
        int newContacts = 0;       // contact points that matched no previous feature key (H1 churn)
        int reapedManifolds = 0;   // manifolds deleted after a one-tick separation (H1 history loss)
        int stickHeals = 0;        // sticky contacts re-anchored after drifting off the surface (#13)
    };
    Diag diag;
#endif

    // Early-exit threshold on the mean penetration residual: when a Newton round
    // improves the average residual by less than this (or worsens it), the iteration set
    // ends early - remaining rounds would polish float noise at full price. 0 disables
    // the check. Realistic value ~1e-5.
    float convergenceThreshold = 0.0f;

    // Combine-pass relaxation for the overflow bucket of the coloured primal - the
    // bodies the greedy colouring could not fit inside kMaxColours, which run as one
    // block-Jacobi round against a pose snapshot. Such a snapshot solve can overshoot
    // where the coloured Gauss-Seidel updates do not; scaling the applied delta by
    // omega < 1 trades per-iteration progress for stability. The
    // bucket is empty in ordinary scenes, so these are tuning knobs for pathological
    // constraint cliques. Suite sweep: 0.5-0.7 hold every rigid bucket scene; 1.0
    // drifts a fully-cliqued stack.
    float jacobiRelaxation = 0.7f;

    // Route every body whose forces are all contact manifolds (Rigid::contactOnly) into the
    // block-Jacobi overflow bucket instead of the coloured Gauss-Seidel. True in a pile of
    // boxes - which makes the bucket the whole primal for a contact-only scene, not an
    // exception. False runs those bodies through the colouring.
    //
    // A/B on the 420-box double pyramid (`test_penetration.gd --pen-seconds=30`, identical
    // wall time, 626 ms / 600-frame pyramid run either way):
    //                        p50 pen   p99 pen   worst *settled*   max speed
    //   true  (Jacobi)        8.0 mm    11.3 mm        15.3 mm        17.9 m/s
    //   false (colour GS)     7.4 mm     8.8 mm        36.3 mm        18.4 m/s
    // Gauss-Seidel improves the median (it is the better solver) but leaves a systematic
    // tail: a pair that the colour order revisits last settles 3.6 cm deep, because the
    // update order is a property of the scene, not of the load. Jacobi averages those
    // order effects away and has the better worst case, so it stays the default until the
    // colouring gets an ordering pass (lowest body first, or a random sweep per round).
    // Exposed as `physics/avbd/contact_gs` on the Godot side.
    bool contactOnlyOverflow = true;

    // Anchor drift (metres) at which a sticky contact's frozen arms are re-adopted from the
    // current collision pass. STICK_DRIFT_HEAL is the default; lowering it re-anchors more
    // eagerly, which trades warm-start continuity for a C0 that always measures the surfaces
    // that are actually touching.
    float stickDriftHeal = STICK_DRIFT_HEAL;

    // Depth (metres) the constraint solve is allowed to see for one contact. A position-level
    // step turns the correction it applies into velocity (v = dp/dt), so an unbounded deep
    // overlap is expelled as a launch. This is the historical cap; see Solver docs.
    float maxPenetration = CONTACT_MAX_PENETRATION;

    // Incremental arm-rebuild distance (metres): later rounds reuse a contact's cached
    // world moment arms while the body rotations have drifted less than this from the
    // captured pose. 0 disables caching (fresh arms every round; default).
    // Realistic value 1e-4 (0.1 mm): in-iteration arm drift is millimetre-scale.
    float jacobianRebuildDistance = 0.0f;

    // Base gap (metres), extended by predicted relative motion until the next collision refresh; 0 disables.
    float speculativeMargin = 0.01f;

    // Position correction excludes resting contacts and does not feed displacement into velocity.
    bool positionProjection = true;
    // Jacobi passes measure all contacts before moving bodies, avoiding manifold-order pose bias.
    // `projectionIterations` is an upper bound: a pass whose worst remaining depth is below
    // 1% of projectionSlop stops the loop, and sleeping/parked pairs are not walked at all.
    int projectionIterations = 4;
    float projectionBeta = 1.0f;
    // Leave millimetre-scale resting penetration to the penalty solve.
    float projectionSlop = 0.008f;
    // Per-body, per-pass cap on summed translation (metres) and rotation (radians); 0 is uncapped.
    float projectionMaxCorrection = 0.0f;

    // Ticks a manifold whose pair slipped out of the swept AABBs stays dormant instead of
    // being deleted. A reaped manifold loses its whole hot state - penalty ramp, duals,
    // stick flags - so the impact-time feature churn (contacts flickering in and out as
    // faces toggle) restarts every returning contact from the penalty floor, and a loaded
    // stack pays that restart hundreds of times per tick. A dormant manifold keeps its
    // contacts and hot state (it contributes nothing while separated: every solve loop
    // walks contacts) and re-matches by feature key when the pair returns. Freshly formed
    // contacts inside a live manifold also seed their penalty from 5% of betaLin instead
    // of the floor while this is enabled. 0 = delete on first separation, as before.
    int contactGraceTicks = 0;

    // Worker threads for the per-body phases. 0 = one per hardware thread, 1 = run
    // everything inline on the calling thread. The result does not depend on this value:
    // every body solves its own block system against the same pose snapshot, and applies
    // its own delta.
    int threads = 0;

    // Idle detection: a movable body whose BDF1 velocities stay under both thresholds for
    // this many consecutive steps falls asleep (velocity zeroed, Godot semantics).
    // 0 disables automatic sleeping; the core tests rely on the default for bit-identical
    // digests. The horizon is a TIME (Godot's 0.5 s) converted to steps by the host, not a
    // fixed 60 frames. The server layer overrides all three from physics/avbd/sleep_*
    // settings (defaults 0.05 m/s, 0.05 rad/s, 1.0x horizon); the core defaults here stay
    // as they were so pinned digests are unaffected.
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
    // Manifold store accessors for hosts and diagnostics. `index` is a position in the
    // live list, not a store slot; the slot is an implementation detail.
    [[nodiscard]] int manifoldCount() const { return manifolds.count(); }
    // Primal partition sizes, for diagnostics: colour groups in use and the number of
    // bodies solved by the block-Jacobi overflow bucket instead of coloured Gauss-Seidel.
    [[nodiscard]] int colourGroups() const { return colourCount; }
    [[nodiscard]] int overflowBodies() const { return static_cast<int>(overflowOrder.size()); }
    // Diagnostics only: did the last broad phase emit this pair as a narrow-phase candidate?
    // A pair that genuinely overlaps and answers no here is a broad-phase miss, not a solve
    // failure. Linear scan - the callers run at diagnostic rates.
        [[nodiscard]] bool sweepHasPair(const Rigid *a, const Rigid *b) const
    {
        const uint64_t key = packPair(a->broadIndex, b->broadIndex);
        for (const uint64_t p : sweepPairs)
            if (p == key)
                return true;
        return false;
    }
    // Diagnostics for a probe run: when a pair matches these two bodies (either order), the
    // broad phase, the manifold pass and the reaper print every lifecycle event for it. A
    // pair with a real overlap and no constraint needs exactly this to be attributed.
    Rigid *tracePair[2] = { nullptr, nullptr };
    [[nodiscard]] bool traceMatch(const Rigid *a, const Rigid *b) const
    {
        return tracePair[0] != nullptr && ((a == tracePair[0] && b == tracePair[1]) ||
                (a == tracePair[1] && b == tracePair[0]));
    }
    [[nodiscard]] Rigid *manifoldBodyA(int index) const
    {
        return manifolds.bodyA[manifolds.live[static_cast<size_t>(index)]];
    }
    [[nodiscard]] Rigid *manifoldBodyB(int index) const
    {
        return manifolds.bodyB[manifolds.live[static_cast<size_t>(index)]];
    }
    [[nodiscard]] int manifoldContactCount(int index) const
    {
        return manifolds.contactCount(manifolds.live[static_cast<size_t>(index)]);
    }
    [[nodiscard]] const Manifold::Contact &manifoldContact(int index, int k) const
    {
        const uint32_t slot = manifolds.live[static_cast<size_t>(index)];
        return manifolds.contacts[static_cast<size_t>(ManifoldStore::base(slot) + k)];
    }
    // Constraint value at the step-start pose: x is the normal error (positive = overlap).
    [[nodiscard]] float3 manifoldContactError(int index, int k) const
    {
        const uint32_t slot = manifolds.live[static_cast<size_t>(index)];
        return manifolds.contactC0[static_cast<size_t>(ManifoldStore::base(slot) + k)];
    }
    // First basis row: the contact normal, pointing from B to A.
    [[nodiscard]] float3 manifoldNormal(int index) const
    {
        return row(manifolds.basis[manifolds.live[static_cast<size_t>(index)]], 0);
    }
    // Candidate pairs the last broad phase emitted (diagnostics: the two broad phases are
    // compared on equal inputs, so their candidate counts must be close).
    [[nodiscard]] int broadPhasePairCount() const { return static_cast<int>(sweepPairs.size()); }
    // Total contact points across live manifolds (statistics only).
    [[nodiscard]] int manifoldContactPoints() const
    {
        int total = 0;
        for (uint32_t slot : manifolds.live)
            total += manifolds.numContacts[slot];
        return total;
    }
    // Delete every body and force.
    void clear();
    void step();

    // Decide which bodies the primal phases may move this step (wake-up, island sleep and
    // the sleep-mode gate), size the per-body delta buffers, and colour the constraint
    // graph for the Gauss-Seidel primal.
    void selectUpdateBodies();
    // Interleaved substeps: one broad phase / force warmstart / update-set selection per
    // tick, then `substeps` rounds of (body warmstart, primal+dual, velocity finish) with
    // the contact error rebased between them. `roundsPerSubstep` 0 splits the tick's
    // `iterations` evenly. `dt` must already be the substep dt.
    void stepInterleaved(int substeps, int roundsPerSubstep);
    // Greedy kMaxColours-bounded colouring of the update set's constraint graph;
    // bodies that do not fit join the block-Jacobi overflow bucket.
    void colourGraph();

    // Number of worker threads actually in use (1 when running inline).
    [[nodiscard]] int threadCount() const;

    // Per-phase wall time of the last step, milliseconds. Diagnostics only - nothing in the
    // solver reads it - and it is what the host's get_*_time() accessors report, so a profile
    // can attribute a step without a second copy of the schedule. A host that reads it from
    // another thread while a step runs gets the previous step's numbers.
    struct PhaseTimes
    {
        double broadPhase = 0.0;  // broad phase: bounds, LBVH build, pair filter, manifold slots
        // The algorithm under the broad phase alone: the LBVH build/refit plus candidate
        // emission, excluding the bounds/AABB pass, the pair sort, the exact filter and
        // manifold slot creation. This separates the broad-phase index cost from the
        // shared pipeline around it.
        double broadIndex = 0.0;
        double graph = 0.0;       // force warmstart + adjacency rebuild + update-set selection/colouring
        double solve = 0.0;       // the iteration region (primal rounds + dual passes)
        double finish = 0.0;      // velocity finish, position projection, contact rebase
    };
    PhaseTimes lastStepTimes;
    // Set by broadPhase() so step()/stepInterleaved can fold it into lastStepTimes.
    double lastBroadIndexMs = 0.0;

    // Per-phase timings of the last step (ms).
    [[nodiscard]] double get_broadPhase_time() const { return lastStepTimes.broadPhase; }
    [[nodiscard]] double get_colourGraph_time() const { return lastStepTimes.graph; }
    [[nodiscard]] double get_solve_time() const { return lastStepTimes.solve; }
    [[nodiscard]] double get_finish_time() const { return lastStepTimes.finish; }

private:
    // Every body in list order. The inertial/warmstart and BDF1 phases run over this,
    // because they also have to keep the cached state of static bodies current: contact
    // and joint code reads `initialLin`/`initialAng` of *both* bodies of a constraint, and
    // the ground is usually one of them.
    std::vector<Rigid *> warmstartOrder;

    // Broad-phase: the body list flattened to index-addressable order (rebuilt every
    // step), the Morton LBVH over those bodies, and the candidate pair list (sorted
    // unique (min index, other), the manifold creation order).
    std::vector<Rigid *> bodiesInOrder;
    // Swept AABBs per bodiesInOrder: the tight bounds grown by this step's linear
    // displacement, the rotation sweep and the speculative margin. center is the
    // current position plus half the displacement; half is the matching extent.
    std::vector<float3> bodyAabbHalf;
    std::vector<float3> bodyAabbCenter;
    LbvhBroadPhase broadLbvh;
    std::vector<uint64_t> sweepPairs;
    // Per-step scratch (reused capacity): exact-filter results for the candidate pairs and
    // the radix sort's ping-pong buffer.
    std::vector<uint8_t> broadPairValid;
    std::vector<uint64_t> pairSortScratch;

    // Bodies that take part in the primal update (movable, awake ones), in list order.
    // The coloured Gauss-Seidel primal solves one block system per body of each
    // colour and applies it immediately; bodies the greedy colouring cannot fit
    // inside kMaxColours land in an overflow bucket that runs one block-Jacobi
    // round (solve every delta against the same snapshot, then combine), using
    // the delta buffers below.
    std::vector<Rigid *> updateOrder;
    std::vector<float3> bodyDeltaLin;
    std::vector<float3> bodyDeltaAng;
    std::vector<float3> projDeltaLin;
    std::vector<float3> projDeltaAng;
    std::vector<int> projManifoldCount;
    // Slots whose pair can still move this step, rebuilt by projectContacts. The brake
    // pass and every projection iteration walk this instead of the full live list.
    std::vector<uint32_t> projectionActive;

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

    // Every collision manifold of this solver, as stable store slots (no per-step
    // allocation, no global pool).
    ManifoldStore manifolds;
    // Per-body constraint adjacency, keyed by Rigid::broadIndex. Rebuilt at the start of
    // the broad phase (pair filter) and again after warm-starting (colour/primal/sleep).
    ForceAdjacency forceAdjacency;
    // The non-manifold force list rebuilt each step so the phases can address joints,
    // springs and ignore pairs by index, plus the per-force result of initialize() so the
    // deletion can be done in one thread after the parallel initialisation.
    std::vector<Force *> forceOrder;
    std::vector<uint8_t> forceActive;
    // Per-slot initialize() result for the step's manifolds (parallel fill, serial reap).
    std::vector<uint8_t> manifoldActive;

    // Island scratch (rebuilt when sleep is active): the movable bodies that may sleep and
    // the island id of every body, indexed by Rigid::broadIndex (-1 = not an island body).
    std::vector<Rigid *> islandBodies;
    std::vector<int> bodyIsland;

    std::unique_ptr<detail::JobPool> pool;
    // Requested thread count the pool was built for; a change rebuilds it (a solver whose
    // `threads` is set after the first step otherwise keeps the old width forever).
    int poolThreads = -1;

    // Create the pool on first use so that a serial solver never spawns threads.
    void ensurePool();
    // Fill `forceOrder` from the solver's force list and return the count.
    int collectForces();

    // The phases of a step (see `step()`, which sequences them).
    void broadPhase();
    int warmstartForces();
    // One body's inertial/warm-start update (independent across the list).
    void warmstartBody(int i);
    // The whole iteration set as one persistent-worker region: `targetIterations` rounds of
    // the block-Jacobi primal (one full block solve per body against the iteration-start
    // pose, then a combine pass) plus dual passes, fenced by LoopSync barriers. Both sync
    // types run the identical body; NullSync is the inline/reference path.
    template <typename Sync>
    void iterateBody(int workerIndex, Sync &sync, int targetIterations, int forceCount,
            float &r_prevResidual, std::atomic<bool> &r_converged);
    // One body's BDF1 velocity update (independent across the list). Idle detection only
    // runs on the final substep, so substeps do not sleep bodies faster than ticks.
    void finishBody(int i, bool idleDetection);
    // One substep's body warmstart, iteration region, and velocity finish.
    void runSubstep(int forceCount, int targetIterations, bool idleDetection);

    // Rebase every manifold's contact error at the current pose (parallel over slots).
    void rebaseContacts();
    // Brake approach velocities against the gap that is still open (`brakeWindow`, default:
    // the whole collision interval), then optionally correct overlapping positions.
    void projectContacts(bool projectPositions, float brakeWindow = 0.0f);
    // Rebuild the per-body constraint adjacency from the current force list and manifold
    // store. Serial: the lists are mutated only on the calling thread at this point.
    void buildForceAdjacency();

    // Island-wide sleep bookkeeping (runs when sleeping is enabled): movable bodies are
    // grouped into islands through force connectivity, and an island whose every body
    // cleared idle detection this step puts all its members to sleep together. Prevents
    // the partially-asleep contact chain from thrashing between states.
    void islandSleep();
};






} // namespace avbd
