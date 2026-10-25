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
 * Replaces the Morton LBVH. Bodies are bucketed by bounding-sphere size into
 * dyadic levels: level 0 has cells of `baseCell` (the big blocks), level l has
 * cells of baseCell * 2^l, and a body lives at the smallest level whose cell
 * fits its diameter. Each level is hashed as packed 21-bit integer coordinates
 * into a sorted key array, so a level's occupied cells form a CSR-indexed
 * bucket table without any per-step hash-map allocation.
 *
 * Pair generation: a body walks the 27-cell neighbourhood of its own centre at
 * every level from its own up to the coarsest. Two bodies whose spheres overlap
 * are at most one cell apart (Chebyshev) at any level that fits both diameters,
 * so every overlapping pair is found from the smaller body's walk; same-level
 * pairs are found from both sides and de-duplicated by the caller's sort+unique.
 *
 * Determinism: keys, bodies and pairs are sorted integers - the output order
 * does not depend on hash iteration or thread count.
 */

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

namespace avbd {

struct Rigid;

class SpatialHash
{
public:
    // Largest level index the hierarchy will use. The coarsest cell must fit
    // the largest body's diameter; build() grows maxLevel until it does.
    static constexpr int kMaxLevels = 16;

    // Rebuild the hierarchy over `bodies` and emit every candidate pair whose
    // bounding spheres can overlap into `out` as (min index, other index).
    // `baseCell` <= 0 selects it automatically: the median bounding-sphere
    // diameter of the scene, so half the bodies land in the level-0 big blocks.
    void buildAndCollectPairs(std::span<Rigid *const> bodies, float baseCell,
                              std::vector<std::pair<int, int>> &out);

    // Statistics of the last build (diagnostics only).
    int levelCount = 0;
    int cellCount = 0;
    int top = 0; // coarsest level index in use

private:
    // Per-level CSR bucket table over the sorted key list.
    std::vector<uint64_t> keys;    // sorted, 21-bit-coords key per occupied cell
    std::vector<int> cellStart;    // keys.size()+1 offsets
    std::vector<int> cellBodies;   // body indices grouped by cell
    // Level l occupies keys[cellStart offsets levelStart[l] .. levelStart[l+1]).
    std::array<int, kMaxLevels + 1> levelStart{};

    std::vector<float> diameters;  // per body, scratch
    std::vector<int> levels;       // per body level assignment
};

} // namespace avbd
