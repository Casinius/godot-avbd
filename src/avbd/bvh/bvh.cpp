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

} // namespace bvh::builder
