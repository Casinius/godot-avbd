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

#include <vector>
#include <span>
#include <limits>
#include <algorithm>

#include "avbd/maths.h"
#include "node_storage.hpp"

// Forward declarations
namespace avbd { struct Rigid; }

namespace bvh {

// Helper functions
// using avbd::min;
// using avbd::abs;
using avbd::quat;
using avbd::float3x3;
// using avbd::diagonal;
using std::max;
static inline float3 rmin(float3 a, float3 b) { return {std::min(a.x(), b.x()), std::min(a.y(), b.y()), std::min(a.z(), b.z())}; }
static inline float3 rmax(float3 a, float3 b) { return {std::max(a.x(), b.x()), std::max(a.y(), b.y()), std::max(a.z(), b.z())}; }

// ============================================
// SAH (Surface Area Heuristic) Utilities
// ============================================
namespace sah {

/// Compute AABB surface area (constexpr for optimization)
[[nodiscard]] constexpr inline float surfaceArea(const float3& min, const float3& max) noexcept {
    float dx = max.x() - min.x();
    float dy = max.y() - min.y();
    float dz = max.z() - min.z();
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}

/// Split axis and position using SAH
struct SplitResult {
    int axis;        // 0=x, 1=y, 2=z
    float pos;       // Split position
    float cost;      // SAH cost
};

[[nodiscard]] inline SplitResult computeBestSplit(const float3& min, const float3& max,
                                                   std::span<const int> bodyIndices,
                                                   std::span<const avbd::Rigid* const> bodies) noexcept;

} // namespace sah

// ============================================
// BVH Builder (Top-down, SAH-based)
// ============================================
namespace builder {

// Forward declaration
class Builder;

} // namespace builder


// ============================================
// BVH Builder Implementation
// ============================================
namespace builder {

class Builder {
public:
    Builder(bvh::nodes::NodeStorage& nodes, std::span<const int> bodyIndices,
            std::span<const avbd::Rigid* const> bodies, int maxDepth = 64) noexcept
        : nodes_(nodes), bodyIndices_(bodyIndices), bodies_(bodies), maxDepth_(maxDepth) {}

    // Builds the tree and returns the root node index.
    int build() noexcept {
        if (bodyIndices_.empty()) return -1;
        return buildRange(std::vector<int>(bodyIndices_.begin(), bodyIndices_.end()), 0);
    }

private:
    bvh::nodes::NodeStorage& nodes_;
    std::span<const int> bodyIndices_;
    std::span<const avbd::Rigid* const> bodies_;
    int maxDepth_;

    // Returns the node index of the subtree root for this index range.
    int buildRange(std::vector<int> local, int depth) noexcept;
};

} // namespace builder

} // namespace bvh
