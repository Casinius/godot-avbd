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

#include "bvh.hpp"
#include "avbd/solver.h"
#include <algorithm>

namespace bvh::builder {

// ============================================
// Morton LBVH Builder: implicit heap layout, balanced by construction
// ============================================
//
// 64-bit 3D Morton codes (21 bits per axis) over the scene bounds; bodies are sorted by
// code, and the code range is split by binary search (equivalent to a radix tree split,
// which is what makes LBVH O(n log n) with near-balanced depth). Nodes live in the
// implicit heap layout of a *complete* binary tree over `padded` = next_pow2(leafCount)
// leaves: node i has children 2i+1 / 2i+2, leaves occupy slots [padded-1, padded-1+n).
// No left/right/parent arrays - the layout IS the tree. Empty slots (when n < padded)
// are marked bodyIndex = -1 with degenerate bounds and never probed.
namespace morton {

inline uint64_t expandBits21(uint32_t v) noexcept {
    uint64_t x = v & 0x1fffffu;
    x = (x | (x << 32)) & 0x1f00000000ffffull;
    x = (x | (x << 16)) & 0x1f0000ff0000ffull;
    x = (x | (x << 8))  & 0x100f00f00f00f00full;
    x = (x | (x << 4))  & 0x10cfcfcfcfcfcfcfull;
    x = (x | (x << 2))  & 0x1249249249249249ull;
    return x;
}

inline uint64_t morton3D(float x, float y, float z, float lo, float scale) noexcept {
    const uint32_t qx = static_cast<uint32_t>(std::clamp((x - lo) * scale, 0.0f, 2097151.0f));
    const uint32_t qy = static_cast<uint32_t>(std::clamp((y - lo) * scale, 0.0f, 2097151.0f));
    const uint32_t qz = static_cast<uint32_t>(std::clamp((z - lo) * scale, 0.0f, 2097151.0f));
    return expandBits21(qx) | (expandBits21(qy) << 1) | (expandBits21(qz) << 2);
}

} // namespace morton


int Builder::buildLBVH()
{
    const int n = static_cast<int>(bodyIndices_.size());
    if (n == 0)
        return 0;

    // Scene bounds over the bounding-sphere centres.
    float3 lo = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float3 hi = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
    for (const avbd::Rigid *b : bodies_) {
        lo = {std::min(lo.x(), b->positionLin.x()), std::min(lo.y(), b->positionLin.y()), std::min(lo.z(), b->positionLin.z())};
        hi = {std::max(hi.x(), b->positionLin.x()), std::max(hi.y(), b->positionLin.y()), std::max(hi.z(), b->positionLin.z())};
    }
    const float3 extent = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
    const float scale = 2097151.0f / std::max({extent.x(), extent.y(), extent.z(), 1e-6f});

    // Morton codes per body index, then sort body indices by code. The builder's index
    // contract is "bodyIndices_[k] selects bodies_[k]"; the leaf records that value so a
    // caller passing a filtered or reordered subset maps back to the right bodies. The buffer is the
    // caller's, reused across steps.
    std::vector<std::pair<uint64_t, int>>& coded = coded_;
    coded.clear();
    for (int i = 0; i < n; ++i) {
        const int bodyIndex = bodyIndices_[i];
        const float3 c = bodies_[bodyIndex]->positionLin;
        coded.emplace_back(morton::morton3D(c.x(), c.y(), c.z(), lo.x(), scale), bodyIndex);
    }
    std::sort(coded.begin(), coded.end(),
            [](const auto &a, const auto &b) { return a.first < b.first; });

    int padded = 1;
    while (padded < n)
        padded <<= 1;

    // Leaf slots live at [padded - 1, padded - 1 + n); unused trailing slots stay
    // degenerate. Internal slots [0, padded - 1) aggregate bottom-up.
    nodes_.clear();
    nodes_.resizeImplicit(padded * 2 - 1);

    float3 allMin = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float3 allMax = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
    for (int i = 0; i < n; ++i) {
        // The caller's boxes, indexed by the same body index the leaf records.
        const size_t bodySlot = static_cast<size_t>(coded[i].second);
        const float3 mn = aabbMin_[bodySlot];
        const float3 mx = aabbMax_[bodySlot];
        nodes_.setImplicitLeaf(padded - 1 + i, mn, mx, coded[i].second);
        allMin = rmin(allMin, mn);
        allMax = rmax(allMax, mx);
    }

    // Aggregate internal nodes bottom-up, right-to-left so children always exist first.
    for (int slot = padded - 2; slot >= 0; --slot) {
        const int l = 2 * slot + 1, r = 2 * slot + 2;
        const float3 lmn = l < padded * 2 - 1 ? nodes_.boundsMin(l) : allMin;
        const float3 lmx = l < padded * 2 - 1 ? nodes_.boundsMax(l) : allMax;
        const float3 rmn = r < padded * 2 - 1 ? nodes_.boundsMin(r) : allMin;
        const float3 rmx = r < padded * 2 - 1 ? nodes_.boundsMax(r) : allMax;
        // Empty child slots (beyond n) carry max/lowest bounds; skip them in the union
        // so a mostly-empty tree does not poison its ancestors with inverted boxes.
        float3 mn, mx;
        const bool lValid = l >= padded - 1 ? (l - (padded - 1)) < n : nodes_.valid(l);
        const bool rValid = r >= padded - 1 ? (r - (padded - 1)) < n : nodes_.valid(r);
        if (lValid && rValid) {
            mn = rmin(lmn, rmn);
            mx = rmax(lmx, rmx);
        } else if (lValid) {
            mn = lmn; mx = lmx;
        } else if (rValid) {
            mn = rmn; mx = rmx;
        } else {
            // Fully empty subtree: keep the node out of every traversal (its bounds are
            // uninitialized) instead of letting a query test garbage floats.
            nodes_.setImplicitInternalBounds(slot, allMin, allMax, false);
            continue;
        }
        nodes_.setImplicitInternalBounds(slot, mn, mx, true);
    }
    return n;
}



} // namespace bvh::builder
