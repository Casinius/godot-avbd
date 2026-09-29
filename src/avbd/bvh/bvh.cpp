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
// cppcheck-suppress missingIncludeSystem
#include <algorithm>

namespace bvh::builder {

// ============================================
// Helper: Calculate bounding-sphere center and AABB for a body
// ============================================
static inline float3 computeAABBCenter(const avbd::Rigid* body) noexcept {
    return body->positionLin;
}
static inline float computeAABBRadius(const avbd::Rigid* body) noexcept {
    return body->radius;
}
static inline void computeAABB(const avbd::Rigid* body, float3& min, float3& max) noexcept {
    const float r = body->radius;
    min = {body->positionLin.x() - r, body->positionLin.y() - r, body->positionLin.z() - r};
    max = {body->positionLin.x() + r, body->positionLin.y() + r, body->positionLin.z() + r};
}

// ============================================
// Build Recursive Implementation
// ============================================
// Returns the node index of the subtree root for this range of body indices.
int Builder::buildRange(std::vector<int> local, int depth) noexcept {
    const int count = static_cast<int>(local.size());

    if (count == 1) {
        const int bodyIndex = local[0];
        const avbd::Rigid* body = bodies_[bodyIndex];
        float3 min, max;
        computeAABB(body, min, max);
        return nodes_.createLeaf(min, max, bodyIndex);
    }

    // Split on the axis with largest span of bounding-sphere centres.
    float3 min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float3 max = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

    for (int i = 0; i < count; ++i) {
        const float3 c = computeAABBCenter(bodies_[local[i]]);
        min = {std::min(min.x(), c.x()), std::min(min.y(), c.y()), std::min(min.z(), c.z())};
        max = {std::max(max.x(), c.x()), std::max(max.y(), c.y()), std::max(max.z(), c.z())};
    }

    const float3 span = {max.x() - min.x(), max.y() - min.y(), max.z() - min.z()};
    const int splitAxis = span.x() >= span.y() ? (span.x() >= span.z() ? 0 : 2) : (span.y() >= span.z() ? 1 : 2);
    const float splitPos = min[splitAxis] + span[splitAxis] * 0.5f;

    int splitPosIdx = 0;
    for (int i = 0; i < count; ++i) {
        if (computeAABBCenter(bodies_[local[i]])[splitAxis] < splitPos) {
            std::swap(local[i], local[splitPosIdx]);
            ++splitPosIdx;
        }
    }
    // Degenerate partition (all centres on one side): median split keeps the tree
    // balanced and the recursion terminating.
    if (splitPosIdx == 0 || splitPosIdx == count)
        splitPosIdx = count / 2;

    // Sub-ranges recurse on their own copies; the (const) span is never written.
    std::vector<int> leftLocal(local.begin(), local.begin() + splitPosIdx);
    std::vector<int> rightLocal(local.begin() + splitPosIdx, local.end());

    const int leftNode = buildRange(std::move(leftLocal), depth + 1);
    const int rightNode = buildRange(std::move(rightLocal), depth + 1);

    const float3 mn = rmin(nodes_.boundsMin()[leftNode], nodes_.boundsMin()[rightNode]);
    const float3 mx = rmax(nodes_.boundsMax()[leftNode], nodes_.boundsMax()[rightNode]);
    return nodes_.createInternal(mn, mx, leftNode, rightNode);
}

// Morton LBVH: compact binary nodes built over contiguous ranges of Morton-sorted leaves.
// Every stored node is valid and carries explicit child indices; no power-of-two padding.
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
        return -1;

    LBVHScratch localScratch;
    LBVHScratch &scratch = scratch_ ? *scratch_ : localScratch;
    auto &coded = scratch.coded;
    auto &sortedBodies = scratch.sortedBodies;
    coded.clear();
    sortedBodies.clear();

    // Scene bounds over the bounding-sphere centres.
    float3 lo = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float3 hi = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};
    for (const int bodyIndex : bodyIndices_) {
        const avbd::Rigid *b = bodies_[bodyIndex];
        lo = {std::min(lo.x(), b->positionLin.x()), std::min(lo.y(), b->positionLin.y()), std::min(lo.z(), b->positionLin.z())};
        hi = {std::max(hi.x(), b->positionLin.x()), std::max(hi.y(), b->positionLin.y()), std::max(hi.z(), b->positionLin.z())};
    }
    const float3 extent = {hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z()};
    const float scale = 2097151.0f / std::max({extent.x(), extent.y(), extent.z(), 1e-6f});

    // Morton codes per body index, then sort body indices by code.
    coded.reserve(n);
    for (const int bodyIndex : bodyIndices_) {
        const float3 c = bodies_[bodyIndex]->positionLin;
        coded.emplace_back(morton::morton3D(c.x(), c.y(), c.z(), lo.x(), scale), bodyIndex);
    }
    std::sort(coded.begin(), coded.end(),
            [](const auto &a, const auto &b) {
                if (a.first != b.first)
                    return a.first < b.first;
                return a.second < b.second;
            });

    nodes_.clear();
    nodes_.reserve(2 * n - 1);
    sortedBodies.reserve(n);
    for (const auto &[code, bodyIndex] : coded)
        sortedBodies.push_back(bodyIndex);

    const auto build = [&](auto &&self, int begin, int end) -> int {
        if (end - begin == 1) {
            const avbd::Rigid *body = bodies_[sortedBodies[begin]];
            float3 mn, mx;
            computeAABB(body, mn, mx);
            return nodes_.createLeaf(mn, mx, sortedBodies[begin]);
        }
        const int mid = begin + (end - begin) / 2;
        const int left = self(self, begin, mid);
        const int right = self(self, mid, end);
        const float3 mn = rmin(nodes_.boundsMin(left), nodes_.boundsMin(right));
        const float3 mx = rmax(nodes_.boundsMax(left), nodes_.boundsMax(right));
        return nodes_.createInternal(mn, mx, left, right);
    };
    return build(build, 0, n);
}



} // namespace bvh::builder
