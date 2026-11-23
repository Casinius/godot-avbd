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

// Each phase has its own break-even batch size. Measured on 12 cores: dispatch ~45 us, primal ~1 us,
// contact solve ~85 ns, dual ~200 ns; break-even ~45 bodies / 500 solves / 250 duals.
constexpr int kMinBodiesPerDispatch = 48;
constexpr int kMinForcesPerDispatch = 256;
// Broad-phase per-body preparation (AABB, level, cell key) costs tens of ns vs ~1 us for a primal update,
// so dispatch pays only at a much larger batch; below this the loop runs inline and saves the round trip.
constexpr int kMinBodiesPerCheapDispatch = 1024;
// Pair filtering is cheaper still (~10-30 ns/pair), so the task-submission round trip only
// pays off when the pair list is genuinely large. Below this the filter runs inline.
constexpr int kMinBroadPairsPerDispatch = 2048;

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

    // Rounds actually run this step (set before any constraint is touched). The joint penalty cap reads it:
    // a penalty stable over many rounds launches a short run. See GenericJoint::penaltyLimit.
    int roundsThisStep = 8;

    // First ceil(newtonRatio * rounds) rounds run the full 6x6 primal; the rest are dual-only with contact
    // stiffness decayed by stiffnessDecay. 1.0/1.0 = all-Newton (historical behaviour).
    float newtonRatio = 1.0f;
    float stiffnessDecay = 1.0f;

    // Per-contact betaLin multiplier carried in the warm start: doubles if penetration did not shrink
    // (capped at ADAPTIVE_BETA_MAX), halves if converging. Per-constraint augmented-Lagrangian kappa update.
    bool adaptiveBeta = false;
    // Residual penetration (m) at or below which a contact counts as solved: multiplier decays rather than
    // doubles, so float noise cannot ratchet a resting contact to the cap.
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

    // Early-exit when a Newton round improves mean penetration by less than this (or worsens it): further
    // rounds would polish float noise. 0 disables; realistic ~1e-5.
    float convergenceThreshold = 0.0f;

    // Relaxation for the block-Jacobi overflow bucket's combine pass: snapshot solves can overshoot, so
    // scale the delta. Sweep: 0.5-0.7 hold all rigid bucket scenes; 1.0 drifts a cliqued stack.
    float jacobiRelaxation = 0.7f;

    // All-contact bodies use block-Jacobi, not colour GS: 420-box A/B worst 15.3 vs 36.3 mm
    // (GS p50/p99 7.4/8.8 vs 8.0/11.3 mm); GS's 3.6 cm colour-order tail keeps Jacobi the default.
    bool contactOnlyOverflow = true;

    // Anchor drift (m) at which sticky frozen arms are re-adopted from the collision pass. Defaults to
    // STICK_DRIFT_HEAL; lowering trades warm-start continuity for a C0 on actually-touching surfaces.
    float stickDriftHeal = STICK_DRIFT_HEAL;

    // Depth (m) the solve sees per contact. Position correction becomes velocity (v = dp/dt), so an
    // unbounded overlap launches bodies. Historical cap; see Solver docs.
    float maxPenetration = CONTACT_MAX_PENETRATION;

    // Arm-rebuild threshold (m): later rounds reuse cached world arms until rotations drift this far.
    // 0 disables caching (default); realistic 1e-4 (0.1 mm), as in-iteration arm drift is millimetre-scale.
    float jacobianRebuildDistance = 0.0f;

    // Base gap (metres), extended by predicted relative motion until the next collision refresh; 0 disables.
    float speculativeMargin = 0.01f;

    // Position correction excludes resting contacts and does not feed displacement into velocity.
    bool positionProjection = true;
    // Jacobi passes measure before moving, avoiding manifold-order bias. projectionIterations is an upper
    // bound: stop when worst depth < 1% of projectionSlop; sleeping/parked pairs are skipped.
    int projectionIterations = 4;
    float projectionBeta = 1.0f;
    // Leave millimetre-scale resting penetration to the penalty solve.
    float projectionSlop = 0.008f;
    // Per-body, per-pass cap on summed translation (metres) and rotation (radians); 0 is uncapped.
    float projectionMaxCorrection = 0.0f;

    // Separated pairs stay dormant (hot state kept, re-matched by feature key): reaping would restart
    // churning contacts from the penalty floor hundreds of times/tick. New contacts seed at 5% betaLin; 0 = old delete.
    int contactGraceTicks = 0;

    // Worker threads for per-body phases: 0 = per hardware thread, 1 = inline. Results are thread-count
    // independent: each body solves its block against the same pose snapshot.
    int threads = 0;

    // Steps under both BDF1 velocity thresholds before sleep (Godot semantics, velocity zeroed); 0 disables.
    // Horizon is a TIME (Godot 0.5 s -> steps, not 60 frames). Keep defaults for bit-identical digests; host overrides.
    float sleepLinearThreshold = 0.1f;       // m/s
    float sleepAngularThreshold = 0.139626f; // rad/s (Godot default: 8 deg/s)
    int sleepFrames = 0;                     // steps below thresholds before sleeping

    // Heads of the two intrusive lists. Owned: `clear()` deletes through them.
    Rigid *bodies = nullptr;
    Force *forces = nullptr;

    Solver();
    ~Solver();

    // Raw-pointer intrusive lists are shared and a step deletes forces, so copying or moving would share list
    // heads and double-delete; hence the deleted special members.
    Solver(const Solver &) = delete;
    Solver &operator=(const Solver &) = delete;
    Solver(Solver &&) = delete;
    Solver &operator=(Solver &&) = delete;

    // Closest body along a ray, and the hit point in that body's local space (null if none).
    [[nodiscard]] Rigid *pick(float3 origin, float3 dir, float3 &local);
    // Masked pick: candidate when (body.layer & p_mask) != 0 (2-arg form forwards an all-accept mask). Both
    // walk the last broad-phase tree; before one exists, they fall back to a body-list walk.
    [[nodiscard]] Rigid *pick(float3 origin, float3 dir, float3 &local, uint32_t p_mask);
    // Bodies whose last broad-phase bounds overlap the AABB and pass (layer & p_mask) != 0, for space
    // shape/point queries; same no-tree fallback.
    void collectAABB(float3 min, float3 max, uint32_t p_mask, std::vector<Rigid *> &out) const;
    // Manifold store accessors for hosts and diagnostics. `index` is a position in the
    // live list, not a store slot; the slot is an implementation detail.
    [[nodiscard]] int manifoldCount() const { return manifolds.count(); }
    // Primal partition sizes, for diagnostics: colour groups in use and the number of
    // bodies solved by the block-Jacobi overflow bucket instead of coloured Gauss-Seidel.
    [[nodiscard]] int colourGroups() const { return colourCount; }
    [[nodiscard]] int overflowBodies() const { return static_cast<int>(overflowOrder.size()); }
    // Diagnostics only: did the last broad phase emit this pair? Overlap + "no" means a broad-phase
    // miss, not a solve failure. Linear scan (callers run at diagnostic rates).
        [[nodiscard]] bool sweepHasPair(const Rigid *a, const Rigid *b) const
    {
        const uint64_t key = packPair(a->broadIndex, b->broadIndex);
        for (const uint64_t p : sweepPairs)
            if (p == key)
                return true;
        return false;
    }
    // Probe diagnostics: a pair matching these bodies (either order) prints every lifecycle event in
    // broad phase, manifold pass and reaper - needed to attribute a real overlap with no constraint.
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

    // Choose movable bodies this step (wake-up, island sleep, sleep-mode gate), size delta buffers and
    // colour the constraint graph for the Gauss-Seidel primal.
    void selectUpdateBodies();
    // One broad phase / force warmstart / selection per tick, then substeps of (warmstart, primal+dual,
    // finish) with C0 rebased between. roundsPerSubstep 0 splits `iterations`; `dt` must be substep dt.
    void stepInterleaved(int substeps, int roundsPerSubstep);
    // Greedy kMaxColours-bounded colouring of the update set's constraint graph;
    // bodies that do not fit join the block-Jacobi overflow bucket.
    void colourGraph();

    // Number of worker threads actually in use (1 when running inline).
    [[nodiscard]] int threadCount() const;

    // Last step's per-phase wall time (ms), diagnostics only - host get_*_time() reports it. A cross-thread
    // read during a step sees the previous step's numbers.
    struct PhaseTimes
    {
        double broadPhase = 0.0;  // broad phase: bounds, LBVH build, pair filter, manifold slots
        // LBVH build/refit + candidate emission only; excludes bounds, pair sort, exact filter and slot
        // creation, separating index cost from the surrounding pipeline. Equals broadTree + broadEmit.
        double broadIndex = 0.0;
        double broadTree = 0.0;   // LBVH topology build + refit (serial, amortised rebuild every 16 steps)
        double broadEmit = 0.0;   // tree-vs-tree candidate pair emission (the parallelisation target)
        double graph = 0.0;       // force warmstart + adjacency rebuild + update-set selection/colouring
        double solve = 0.0;       // the iteration region (primal rounds + dual passes)
        double finish = 0.0;      // velocity finish, position projection, contact rebase
    };
    PhaseTimes lastStepTimes;
    // Set by broadPhase() so step()/stepInterleaved can fold it into lastStepTimes.
    double lastBroadIndexMs = 0.0;
    double lastBroadTreeMs = 0.0;
    double lastBroadEmitMs = 0.0;

    // Per-phase timings of the last step (ms).
    [[nodiscard]] double get_broadPhase_time() const { return lastStepTimes.broadPhase; }
    [[nodiscard]] double get_broadTree_time() const { return lastStepTimes.broadTree; }
    [[nodiscard]] double get_broadEmit_time() const { return lastStepTimes.broadEmit; }
    [[nodiscard]] double get_colourGraph_time() const { return lastStepTimes.graph; }
    [[nodiscard]] double get_solve_time() const { return lastStepTimes.solve; }
    [[nodiscard]] double get_finish_time() const { return lastStepTimes.finish; }

private:
    // All bodies in list order: warmstart/BDF1 must refresh cached state even for static bodies, since
    // constraints read initialLin/initialAng of both ends (often the ground).
    std::vector<Rigid *> warmstartOrder;

    // Flattened index-addressable body list (rebuilt each step), Morton LBVH over it, and the sorted
    // unique candidate pair list (also the manifold creation order).
    std::vector<Rigid *> bodiesInOrder;
    // Swept AABBs: tight bounds grown by linear displacement, rotation sweep and speculative margin.
    // center = current position + half displacement; half is the matching extent.
    std::vector<float3> bodyAabbHalf;
    std::vector<float3> bodyAabbCenter;
    LbvhBroadPhase broadLbvh;
    std::vector<uint64_t> sweepPairs;
    // Per-step scratch (reused capacity): exact-filter results for the candidate pairs and
    // the radix sort's ping-pong buffer.
    std::vector<uint8_t> broadPairValid;
    std::vector<uint64_t> pairSortScratch;

    // Movable, awake bodies in list order. Colour GS solves each body's block and applies it at once;
    // those not fitting kMaxColours form an overflow block-Jacobi round using the delta buffers below.
    std::vector<Rigid *> updateOrder;
    std::vector<float3> bodyDeltaLin;
    std::vector<float3> bodyDeltaAng;
    std::vector<float3> projDeltaLin;
    std::vector<float3> projDeltaAng;
    std::vector<int> projManifoldCount;
    // Slots whose pair can still move this step, rebuilt by projectContacts. The brake
    // pass and every projection iteration walk this instead of the full live list.
    std::vector<uint32_t> projectionActive;

    // Per-step colouring: same-colour bodies share no force, so solve/apply in parallel (GS between
    // colours). bodyColour[updateOrder pos] is a colour or kMaxColours for the Jacobi overflow bucket.
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
    // Non-manifold forces rebuilt each step for index access to joints/springs/ignores, plus each
    // initialize() result so deletion runs in one thread after parallel init.
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
    // Persistent-worker region: targetIterations block-Jacobi primal rounds (solve vs iteration-start
    // pose, then combine) + dual passes, fenced by LoopSync. Both sync types run the identical body.
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

    // With sleeping enabled, group movable bodies into force-connected islands; if all pass idle
    // detection, sleep the whole island, so partially-asleep contact chains cannot thrash.
    void islandSleep();
};






} // namespace avbd
