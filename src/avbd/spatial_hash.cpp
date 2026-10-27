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
#include <array>
#include <cmath>
#include <utility>
#include <memory>

#include "solver.h"
#include "avbd/job_pool.hpp"
#include "radixsort.hpp"
namespace avbd {

// Cell coordinates are packed into 20 bits per axis (60 bits total). The remaining top
// four bits carry the level index, so `(level << 60) | cellKey` is a single sortable
// integer and the per-level CSR ranges fall out of one radix sort. A scene wider than
// 2^20 cells at a level wraps keys, which can only merge distant cells and add
// candidate pairs the exact filter in the solver then rejects - true overlaps can
// never be missed by it.
inline constexpr int kCellCoordBits = 20;
inline constexpr uint64_t kCellCoordMask = (uint64_t{1} << kCellCoordBits) - 1;

[[nodiscard]] static inline uint64_t cellKey(int x, int y, int z) noexcept
{
    const uint64_t mx = static_cast<uint64_t>(static_cast<uint32_t>(x)) & kCellCoordMask;
    const uint64_t my = static_cast<uint64_t>(static_cast<uint32_t>(y)) & kCellCoordMask;
    const uint64_t mz = static_cast<uint64_t>(static_cast<uint32_t>(z)) & kCellCoordMask;
    return (mx << (2 * kCellCoordBits)) | (my << kCellCoordBits) | mz;
}

[[nodiscard]] static inline uint64_t cellSortKey(int level, uint64_t cell) noexcept
{
    return (static_cast<uint64_t>(level) << (3 * kCellCoordBits)) | cell;
}

namespace {

// Stable LSD radix sort of packed cell keys with a parallel body-index payload. Sorting
// `(level, cell)` pairs this way replaces the former comparator over a 24-byte struct:
// the keys are unsigned integers, the payload rides along, and equal keys may come out
// in any order (a cell's bodies are re-sorted by the pair sort downstream anyway).
void radixSortKeyed(std::vector<uint64_t> &keys, std::vector<int> &payload,
        std::vector<uint64_t> &tmpKeys, std::vector<int> &tmpPayload)
{
    const std::size_t n = keys.size();
    if (n < 2)
        return;

    if (tmpKeys.size() < n)
        tmpKeys.resize(n);
    if (tmpPayload.size() < n)
        tmpPayload.resize(n);
    uint64_t *srcKeys = keys.data();
    uint64_t *dstKeys = tmpKeys.data();
    int *srcPayload = payload.data();
    int *dstPayload = tmpPayload.data();

    std::array<std::size_t, 256> count{};
    for (unsigned shift = 0; shift < 64; shift += 8)
    {
        count.fill(0);
        for (std::size_t i = 0; i < n; ++i)
            ++count[(srcKeys[i] >> shift) & 0xFFu];

        std::size_t total = 0;
        for (auto &c : count)
        {
            const std::size_t old = c;
            c = total;
            total += old;
        }

        for (std::size_t i = 0; i < n; ++i)
        {
            const unsigned bucket = (srcKeys[i] >> shift) & 0xFFu;
            const std::size_t dst = count[bucket]++;
            dstKeys[dst] = srcKeys[i];
            dstPayload[dst] = srcPayload[i];
        }
        std::swap(srcKeys, dstKeys);
        std::swap(srcPayload, dstPayload);
    }

    if (srcKeys != keys.data())
    {
        std::copy(srcKeys, srcKeys + n, keys.data());
        std::copy(srcPayload, srcPayload + n, payload.data());
    }
}

} // namespace

void SpatialHash::buildAndCollectPairs(std::span<Rigid *const> bodies, float baseCell,
                                        std::vector<uint64_t> &out,
                                        std::span<const float3> aabbHalfExt,
                                        std::span<const float3> aabbCenter,
                                        detail::JobPool &pool)
{
    const int count = static_cast<int>(bodies.size());
    out.clear();
    if (count < 2)
        return;

    // Automatic big-block size: the median largest AABB extent, so about half
    // the bodies live in the level-0 big blocks and one level up covers each
    // further doubling of size. The extents come from the caller's bounds (swept,
    // so a fast body's level follows the ground it covers this step).
    extents.resize(count);
    pool.forCount(count, [&](int i) {
        const float3 h = aabbHalfExt[i];
        extents[i] = 2.0f * std::max(h.x, std::max(h.y, h.z));
    }, kMinBodiesPerCheapDispatch);

    if (!(baseCell > 0.0f)) {
        nonEmpty.clear();
        nonEmpty.reserve(static_cast<size_t>(count));
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
    //
    // `top` and the median above need the whole extent array, so they stay serial
    // (two O(n) scans); the per-body level decision is independent and runs on the
    // pool. Nothing here reduces across threads: each body writes its own slot.
    float largest = 0.0f;
    for (int i = 0; i < count; i++)
        largest = std::max(largest, extents[i]);
    top = 0;
    while (baseCell * (float)(1 << top) < largest && top < kMaxLevels - 1)
        ++top;
    levelCount = top + 1;

    levels.resize(count);
    pool.forCount(count, [&](int i) {
        int l = 0;
        while (baseCell * (float)(1 << l) < extents[i] && l < top)
            ++l;
        levels[i] = l;
    }, kMinBodiesPerCheapDispatch);

    // Packed (level, cell) key plus the body index, filled by index so the dispatch is
    // race-free (this used to be a serial push_back). The radix sort below is stable but
    // need not preserve body order for equal keys: cell contents have no defined order,
    // and the caller's pair sort removes every trace of it.
    keyedKeys.resize(static_cast<size_t>(count));
    keyedBodies.resize(static_cast<size_t>(count));
    pool.forCount(count, [&](int i) {
        const int l = levels[i];
        const float invCell = 1.0f / (baseCell * (float)(1 << l));
        const float3 &p = aabbCenter[i];
        keyedKeys[i] = cellSortKey(l, cellKey((int)std::floor(p.x * invCell),
                (int)std::floor(p.y * invCell), (int)std::floor(p.z * invCell)));
        keyedBodies[i] = i;
    }, kMinBodiesPerCheapDispatch);
    radixSortKeyed(keyedKeys, keyedBodies, radixKeysTmp, radixBodiesTmp);

    keys.clear();
    cellStart.clear();
    cellBodies.clear();
    cellStart.push_back(0);
    std::array<int, kMaxLevels + 1> levelCells{};
    for (std::size_t k = 0; k < keyedKeys.size();) {
        const uint64_t key = keyedKeys[k];
        std::size_t e = k;
        while (e < keyedKeys.size() && keyedKeys[e] == key)
            ++e;
        keys.push_back(key & ((uint64_t{1} << (3 * kCellCoordBits)) - 1));
        ++levelCells[(key >> (3 * kCellCoordBits)) & 0xFu];
        for (std::size_t j = k; j < e; j++)
            cellBodies.push_back(keyedBodies[j]);
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
    //
    // Emission is the phase's serial bottleneck (the per-body loops above and the
    // exact filter are already on the pool), so bodies are split into contiguous
    // chunks and each chunk appends to its own buffer. Pair order after the walk is
    // therefore chunk-dependent, which is harmless: the caller radix-sorts and
    // de-duplicates the packed pairs, so the final order is a pure function of the
    // candidate set.
    //
    // The cell key is (x << 40) | (y << 20) | z, so inside one level the table is
    // sorted by x, then y, then z. A column of constant (x, y) is therefore a
    // contiguous run sorted by z. One binary search lands on the first cell with
    // z >= the window's low end; the at-most-three occupied z values in the
    // window are then adjacent entries, so the walk needs 9 searches, not 27.
    const auto walk = [&](int chunk, int firstBody, int lastBody) {
        std::vector<uint64_t> &buf = pairChunks[chunk];
        buf.clear();
        for (int i = firstBody; i < lastBody; i++) {
            const float3 &pi = aabbCenter[i];
            for (int l = levels[i]; l < levelCount; l++) {
                const float invCell = 1.0f / (baseCell * (float)(1 << l));
                const int cx = (int)std::floor(pi.x * invCell);
                const int cy = (int)std::floor(pi.y * invCell);
                const int cz = (int)std::floor(pi.z * invCell);
                const int first = levelStart[l];
                const int last = levelStart[l + 1];

                // One z interval per column, or two when the window wraps through 0.
                int zLo[2];
                int zHi[2];
                int zRanges = 0;
                const int czMasked = static_cast<int>(static_cast<uint32_t>(cz) & kCellCoordMask);
                if (czMasked == 0) {
                    zLo[zRanges] = 0;
                    zHi[zRanges++] = 1;
                    zLo[zRanges] = static_cast<int>(kCellCoordMask);
                    zHi[zRanges++] = static_cast<int>(kCellCoordMask);
                } else if (czMasked == static_cast<int>(kCellCoordMask)) {
                    zLo[zRanges] = czMasked - 1;
                    zHi[zRanges++] = czMasked;
                    zLo[zRanges] = 0;
                    zHi[zRanges++] = 0;
                } else {
                    zLo[zRanges] = czMasked - 1;
                    zHi[zRanges++] = czMasked + 1;
                }

                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        const uint64_t columnStart = cellKey(cx + dx, cy + dy, 0);
                        const uint64_t columnPrefix = columnStart >> kCellCoordBits;
                        for (int r = 0; r < zRanges; r++) {
                            // First cell in the column at or after the window's low end.
                            const uint64_t startKey = cellKey(cx + dx, cy + dy, zLo[r]);
                            int lo = first, hi = last;
                            while (lo < hi) {
                                const int mid = (lo + hi) / 2;
                                if (keys[mid] < startKey)
                                    lo = mid + 1;
                                else
                                    hi = mid;
                            }
                            for (int p = lo; p < last; p++) {
                                const uint64_t key = keys[p];
                                if ((key >> kCellCoordBits) != columnPrefix)
                                    break;
                                const int z = static_cast<int>(key & kCellCoordMask);
                                if (z > zHi[r])
                                    break;
                                for (int s = cellStart[p]; s < cellStart[p + 1]; s++) {
                                    const int j = cellBodies[s];
                                    if (j != i)
                                        buf.push_back(packPair(i, j));
                                }
                            }
                        }
                    }
                }
            }
        }
    };

    if (count >= kMinBodiesPerCheapDispatch) {
        const int blocks = pool.blockCount(count);
        pairChunks.resize(static_cast<size_t>(blocks));
        pool.forBlocks(count, walk);
    } else {
        // Small scene: one chunk on the calling thread, no dispatch round trip.
        pairChunks.resize(1);
        walk(0, 0, count);
    }

    std::size_t total = 0;
    for (const std::vector<uint64_t> &chunk : pairChunks)
        total += chunk.size();
    out.reserve(total);
    for (const std::vector<uint64_t> &chunk : pairChunks)
        out.insert(out.end(), chunk.begin(), chunk.end());
}

} // namespace avbd
