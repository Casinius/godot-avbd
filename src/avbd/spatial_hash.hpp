/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 *
 * Hierarchical hashed-grid broad phase ("big-block screen + layered hash").
 *
 * The default broad phase (`Solver::broadPhaseMode = 0`); the Morton LBVH is the
 * alternative selected by `broadPhaseMode = 1`.
 *
 * Bodies are bucketed by the largest full extent of their world AABB (the solver passes
 * swept bounds) into dyadic levels: level 0 has cells of `baseCell` (the big blocks),
 * level l has cells of baseCell * 2^l, and a body lives at the smallest level whose cell
 * fits its AABB. Each level is hashed as packed 20-bit integer coordinates; the level
 * index rides in the sort key's top four bits, so one radix sort produces both the
 * per-level CSR ranges and the bucketed cells.
 *
 * Pair generation: a body walks the 27-cell neighbourhood of its AABB centre at every
 * non-empty level from its own up to the coarsest. Two bodies whose AABBs overlap are at
 * most one cell apart (Chebyshev) at any level that fits both AABBs, so every overlapping
 * pair is found from the smaller body's walk. A same-level pair is emitted only by the
 * smaller body index, and the swept AABBs are tested before the pair is written, so the
 * output is the exact AABB-overlap set rather than the 27-cell superset. The solver still
 * applies its layer / constraint / both-static filter on top.
 *
 * Determinism: keys, bodies and pairs are sorted integers - the output order does not
 * depend on hash iteration or thread count. The per-body work (extent, level, cell key)
 * is written by index and runs on the job pool; the sequential parts are the two O(n)
 * scans that pick the big-block size and the top level.
 */

#pragma once
#include "avbd/broadphase_pair.hpp"
#include "avbd/job_pool.hpp"
#include <array>
#include <cstdint>
#include <span>
#include <vector>
#include "maths.h"
namespace avbd {

struct Rigid;

class SpatialHash
{
public:
    // Largest level index the hierarchy will use. The coarsest cell must fit
    // the largest body's AABB extent; build() grows maxLevel until it does.
    static constexpr int kMaxLevels = 16;

    // Rebuild the hierarchy over `bodies` and emit every pair whose swept world AABBs
    // overlap into `out`, packed by `packPair` (see above). The solver's exact filter
    // re-checks layers and constraints but the geometric test is already done here.
    // `baseCell` <= 0 selects it automatically: the median largest AABB extent
    // of the scene, so half the bodies land in the level-0 big blocks.
    // `aabbHalfExt` / `aabbCenter` are the precomputed world AABB bounds per body
    // (the solver passes swept bounds, so a fast body is not missed). They are also
    // the overlap test's inputs, so the bound is computed once per body per step.
    void buildAndCollectPairs(std::span<Rigid *const> bodies, float baseCell,
                              std::vector<uint64_t> &out,
                              std::span<const float3> aabbHalfExt,
                              std::span<const float3> aabbCenter,
                              detail::JobPool &pool);

    // Statistics of the last build (diagnostics only).
    int levelCount = 0;
    int cellCount = 0;
    int top = 0; // coarsest level index in use

private:
    // Per-level CSR bucket table over the sorted key list.
    std::vector<uint64_t> keys;    // sorted, 20-bit-coords key per occupied cell
    std::vector<int> cellStart;    // keys.size()+1 offsets
    std::vector<int> cellBodies;   // body indices grouped by cell
    // Level l occupies keys[cellStart offsets levelStart[l] .. levelStart[l+1]).
    std::array<int, kMaxLevels + 1> levelStart{};

    std::vector<float> extents;  // per body: largest AABB full extent, scratch
    std::vector<int> levels;       // per body level assignment
    // Scratch grown once and reused every step: the (level, cell, body) sort inputs, the
    // radix ping-pong buffers, and the median filter for the automatic cell size.
    std::vector<uint64_t> keyedKeys;
    std::vector<int> keyedBodies;
    std::vector<uint64_t> radixKeysTmp;
    std::vector<int> radixBodiesTmp;
    std::vector<float> nonEmpty;
    // Per-chunk pair buffers for the parallel walk (capacity survives across steps, so a
    // steady scene allocates nothing here).
    std::vector<std::vector<uint64_t>> pairChunks;
};

} // namespace avbd
