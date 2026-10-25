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
 * Replaces the Morton LBVH. Bodies are bucketed by the largest full extent of
 * their tight world AABB into dyadic levels: level 0 has cells of `baseCell`
 * (the big blocks), level l has cells of baseCell * 2^l, and a body lives at the
 * smallest level whose cell fits its AABB. Each level is hashed as packed 20-bit
 * integer coordinates; the level index rides in the sort key's top four bits, so
 * one radix sort produces both the per-level CSR ranges and the bucketed cells.
 *
 * Pair generation: a body walks the 27-cell neighbourhood of its own centre at
 * every level from its own up to the coarsest. Two bodies whose AABBs overlap
 * are at most one cell apart (Chebyshev) at any level that fits both AABBs, so
 * every overlapping pair is found from the smaller body's walk; same-level pairs
 * are found from both sides and de-duplicated by the caller's sort+unique.
 *
 * Determinism: keys, bodies and pairs are sorted integers - the output order
 * does not depend on hash iteration or thread count. The per-body work (extent,
 * level, cell key) is written by index and runs on the job pool; the sequential
 * parts are the two O(n) scans that pick the big-block size and the top level.
 */

#pragma once
#include "avbd/job_pool.hpp"
#include <array>
#include <cstdint>
#include <span>
#include <vector>
#include "maths.h"
namespace avbd {

struct Rigid;

// A broad-phase candidate pair packed into one integer for the sort: the smaller body
// index in the high 32 bits, the larger in the low 32. Lexicographic order over
// `std::pair<int,int>` is therefore plain integer order, which is what lets the unlink
// step use a real radix sort instead of a comparator sort.
[[nodiscard]] inline constexpr uint64_t packPair(int p_a, int p_b) noexcept
{
    const uint64_t lo = static_cast<uint64_t>(static_cast<uint32_t>(p_a < p_b ? p_a : p_b));
    const uint64_t hi = static_cast<uint64_t>(static_cast<uint32_t>(p_a < p_b ? p_b : p_a));
    return (lo << 32) | hi;
}

[[nodiscard]] inline constexpr int pairFirst(uint64_t p_key) noexcept
{
    return static_cast<int>(p_key >> 32);
}

[[nodiscard]] inline constexpr int pairSecond(uint64_t p_key) noexcept
{
    return static_cast<int>(p_key & 0xFFFFFFFFu);
}

class SpatialHash
{
public:
    // Largest level index the hierarchy will use. The coarsest cell must fit
    // the largest body's AABB extent; build() grows maxLevel until it does.
    static constexpr int kMaxLevels = 16;

    // Rebuild the hierarchy over `bodies` and emit every candidate pair whose
    // world AABBs can overlap into `out`, packed by `packPair` (see above).
    // `baseCell` <= 0 selects it automatically: the median largest AABB extent
    // of the scene, so half the bodies land in the level-0 big blocks.
    // `aabbHalfExt` is filled with each body's tight world AABB half extents (the
    // exact filter in the solver reuses them), so the per-body bound is computed once.
    void buildAndCollectPairs(std::span<Rigid *const> bodies, float baseCell,
                              std::vector<uint64_t> &out, std::span<float3> aabbHalfExt,
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
    // Per-chunk pair buffers for the parallel walk (capacity survives across steps, so a
    // steady scene allocates nothing here).
    std::vector<std::vector<uint64_t>> pairChunks;
};

} // namespace avbd
