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

#include "spatial_hash.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "solver.h"

namespace avbd {

// Pack one cell coordinate triple into a 63-bit key. Coordinates are masked
// into [0, 2^21): a scene wider than 2^21 cells at a level wraps keys, which
// can only merge distant cells and add candidate pairs the exact filter in the
// solver then rejects - true overlaps can never be missed by it.
[[nodiscard]] static inline uint64_t cellKey(int x, int y, int z) noexcept
{
    const uint64_t mx = static_cast<uint64_t>(static_cast<uint32_t>(x) & 0x1FFFFFu);
    const uint64_t my = static_cast<uint64_t>(static_cast<uint32_t>(y) & 0x1FFFFFu);
    const uint64_t mz = static_cast<uint64_t>(static_cast<uint32_t>(z) & 0x1FFFFFu);
    return (mx << 42) | (my << 21) | mz;
}

void SpatialHash::buildAndCollectPairs(std::span<Rigid *const> bodies, float baseCell,
                                        std::vector<std::pair<int, int>> &out)
{
    const int count = static_cast<int>(bodies.size());
    out.clear();
    if (count < 2)
        return;

    // Automatic big-block size: the median largest AABB extent, so about half
    // the bodies live in the level-0 big blocks and one level up covers each
    // further doubling of size.
    extents.resize(count);
    for (int i = 0; i < count; i++) {
        const float3 half = aabbHalfExtent(*bodies[i]);
        extents[i] = 2.0f * std::max(half.x, std::max(half.y, half.z));
    }
    if (!(baseCell > 0.0f)) {
        std::vector<float> nonEmpty;
        nonEmpty.reserve(count);
        for (int i = 0; i < count; i++)
            if (extents[i] > 0.0f)
                nonEmpty.push_back(extents[i]);
        if (nonEmpty.empty())
            return;
        std::nth_element(nonEmpty.begin(), nonEmpty.begin() + nonEmpty.size() / 2, nonEmpty.end());
        baseCell = nonEmpty[nonEmpty.size() / 2];
        if (!(baseCell > 0.0f))
            baseCell = 1.0f;
    }

    // Level assignment: the smallest level whose cell fits the AABB's largest
    // extent. The coarsest level is grown to fit the largest body, so every
    // body's AABB fits its own level's cell and the 27-cell walk below is exact.
    float largest = 0.0f;
    for (int i = 0; i < count; i++)
        largest = std::max(largest, extents[i]);
    top = 0;
    while (baseCell * (float)(1 << top) < largest && top < kMaxLevels - 1)
        ++top;
    levelCount = top + 1;

    levels.resize(count);
    for (int i = 0; i < count; i++) {
        int l = 0;
        while (baseCell * (float)(1 << l) < extents[i] && l < top)
            ++l;
        levels[i] = l;
    }

    // Hash: (level, key) pairs sorted, then run-lengthed into a CSR bucket
    // table. The level is a separate sort key: packing it into the key's top
    // bits would overlap the 21-bit X coordinate and merge/split level ranges
    // for scenes whose cell coordinates are large.
    struct Keyed
    {
        uint64_t key;
        int level;
        int body;
    };
    std::vector<Keyed> keyed;
    keyed.reserve(count);
    for (int i = 0; i < count; i++) {
        const int l = levels[i];
        const float invCell = 1.0f / (baseCell * (float)(1 << l));
        const float3 &p = bodies[i]->positionLin;
        const uint64_t key = cellKey((int)std::floor(p.x * invCell),
                (int)std::floor(p.y * invCell), (int)std::floor(p.z * invCell));
        keyed.push_back({key, l, i});
    }
    std::sort(keyed.begin(), keyed.end(), [](const Keyed &a, const Keyed &b) {
        if (a.level != b.level)
            return a.level < b.level;
        return a.key != b.key ? a.key < b.key : a.body < b.body;
    });

    keys.clear();
    cellStart.clear();
    cellBodies.clear();
    cellStart.push_back(0);
    std::array<int, kMaxLevels + 1> levelCells{};
    for (size_t k = 0; k < keyed.size();) {
        size_t e = k;
        while (e < keyed.size() && keyed[e].level == keyed[k].level && keyed[e].key == keyed[k].key)
            ++e;
        keys.push_back(keyed[k].key);
        ++levelCells[keyed[k].level];
        for (size_t j = k; j < e; j++)
            cellBodies.push_back(keyed[j].body);
        cellStart.push_back(static_cast<int>(cellBodies.size()));
        k = e;
    }
    cellCount = static_cast<int>(keys.size());

    // Per-level ranges over the cell table, including levels with no cells.
    levelStart[0] = 0;
    for (int l = 0; l < levelCount; l++)
        levelStart[l + 1] = levelStart[l] + levelCells[l];

    // Pair walk. For body i, every level k from its own to the top: the
    // 27-cell neighbourhood at level k, each cell's bodies emitted as
    // (min(i,j), max(i,j)). A pair of overlapping AABBs is found when the
    // smaller-at-coarser body walks their common level (same-level pairs are
    // found from both sides); duplicates collapse in the caller's sort+unique.
    // The caller applies the exact AABB overlap / layer / constraint filter.
    for (int i = 0; i < count; i++) {
        const float3 &pi = bodies[i]->positionLin;
        for (int l = levels[i]; l < levelCount; l++) {
            const float invCell = 1.0f / (baseCell * (float)(1 << l));
            const int cx = (int)std::floor(pi.x * invCell);
            const int cy = (int)std::floor(pi.y * invCell);
            const int cz = (int)std::floor(pi.z * invCell);
            const int first = levelStart[l];
            const int last = levelStart[l + 1];
            for (int dz = -1; dz <= 1; dz++) {
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        const uint64_t key = cellKey(cx + dx, cy + dy, cz + dz);
                        // Binary search the key within this level's range.
                        int lo = first, hi = last;
                        while (lo < hi) {
                            const int mid = (lo + hi) / 2;
                            if (keys[mid] < key)
                                lo = mid + 1;
                            else
                                hi = mid;
                        }
                        if (lo >= last || keys[lo] != key)
                            continue;
                        for (int s = cellStart[lo]; s < cellStart[lo + 1]; s++) {
                            const int j = cellBodies[s];
                            if (j != i)
                                out.emplace_back(i < j ? i : j, i < j ? j : i);
                        }
                    }
                }
            }
        }
    }
}

} // namespace avbd
