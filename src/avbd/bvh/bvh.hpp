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
// BVH Builder (Morton LBVH, implicit heap layout)
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
    // The bounding boxes are the caller's: the broad phase has to compute them anyway to decide which
    // bodies moved, and a second traversal of the same shapes inside the builder is pure duplication.
    // `scratch` is the caller's (code, body) buffer, reused across steps: the leaves are sorted by
    // Morton code every tick, so a per-tick allocation for ~1600 entries is a measurable slice of a
    // build that has to stay well under a millisecond.
    Builder(bvh::nodes::NodeStorage& nodes, std::span<const int> bodyIndices,
            std::span<const avbd::Rigid* const> bodies, std::span<const float3> aabbMin,
            std::span<const float3> aabbMax, std::vector<std::pair<uint64_t, int>>& scratch,
            int maxDepth = 64) noexcept
        : nodes_(nodes), bodyIndices_(bodyIndices), bodies_(bodies), aabbMin_(aabbMin), aabbMax_(aabbMax),
          coded_(scratch), maxDepth_(maxDepth) {}

    // Morton LBVH: implicit heap layout in `nodes` (children of i are 2i+1 / 2i+2),
    // balanced by construction. Returns the number of leaves stored (also the offset of
    // the first leaf slot: leaves live at [padded-1 .. padded-2+n] of the node arrays).
    int buildLBVH();

private:
    bvh::nodes::NodeStorage& nodes_;
    std::span<const int> bodyIndices_;
    std::span<const avbd::Rigid* const> bodies_;
    std::span<const float3> aabbMin_;
    std::span<const float3> aabbMax_;
    std::vector<std::pair<uint64_t, int>>& coded_;
    int maxDepth_;
};

} // namespace builder

} // namespace bvh
