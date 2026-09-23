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
// Helper: Calculate AABB for a range of bodies
// ============================================
static inline float3 computeAABB(const avbd::Rigid* body) noexcept {
    // Simplified AABB - in a real implementation, this would use the actual shape bounds
    float3 pos = body->positionLin;
    return {pos.x(), pos.y(), pos.z()};
}

// ============================================
// Build Recursive Implementation
// ============================================
void Builder::buildRecursive(int start, int end, int depth) noexcept {
    int count = end - start;

    // If we have a single body, create a leaf node
    if (count == 1) {
        int bodyIndex = bodyIndices_[start];
        const avbd::Rigid* body = bodies_[bodyIndex];
        float3 pos = computeAABB(body);
        float3 min = pos;
        float3 max = pos;
        nodes_.createLeaf(min, max, bodyIndex);
        return;
    }

    // Otherwise, find the best split
    // TODO: Implement SAH-based split
    // For now, use a simple axis-aligned split
    float3 min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
    float3 max = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

    // Create a copy of the indices for sorting
    std::vector<int> sortedIndices(bodyIndices_.begin(), bodyIndices_.end());

    for (int i = 0; i < count; ++i) {
        int bodyIndex = sortedIndices[i];
        const avbd::Rigid* body = bodies_[bodyIndex];
        float3 pos = computeAABB(body);
        min = {std::min(min.x(), pos.x()), std::min(min.y(), pos.y()), std::min(min.z(), pos.z())};
        max = {std::max(max.x(), pos.x()), std::max(max.y(), pos.y()), std::max(max.z(), pos.z())};
    }

    // Simple split: divide on the axis with largest span
    float3 span = {max.x() - min.x(), max.y() - min.y(), max.z() - min.z()};
    int splitAxis = span.x() >= span.y() ? (span.x() >= span.z() ? 0 : 2) : (span.y() >= span.z() ? 1 : 2);
    float splitPos = min[splitAxis] + span[splitAxis] * 0.5f;

    // Find split position
    int splitPosIdx = 0;
    for (int i = 0; i < count; ++i) {
        int bodyIndex = sortedIndices[i];
        const avbd::Rigid* body = bodies_[bodyIndex];
        float3 pos = computeAABB(body);
        if (pos[splitAxis] < splitPos) {
            std::swap(sortedIndices[i], sortedIndices[splitPosIdx]);
            ++splitPosIdx;
        }
    }

    // If split didn't partition (all bodies on one side), make a leaf
    if (splitPosIdx == 0 || splitPosIdx == count) {
        for (int i = 0; i < count; ++i) {
            int bodyIndex = sortedIndices[i];
            const avbd::Rigid* body = bodies_[bodyIndex];
            float3 pos = computeAABB(body);
            float3 nodeMin = pos;
            float3 nodeMax = pos;
            nodes_.createLeaf(nodeMin, nodeMax, bodyIndex);
        }
        return;
    }

    // Recursively build children
    buildRecursive(start, start + splitPosIdx, depth + 1);
    buildRecursive(start + splitPosIdx, end, depth + 1);
}

} // namespace bvh::builder
