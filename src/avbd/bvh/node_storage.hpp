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

// cppcheck-suppress missingIncludeSystem
#include <vector>
// cppcheck-suppress missingIncludeSystem
#include <span>

#include "avbd/maths.h"

namespace bvh {

using avbd::float3;

// Forward declarations
namespace nodes { class NodeStorage; }

namespace nodes {

class NodeStorage {
private:
    // Contiguous arrays for cache locality (SafeC++: no raw pointers)
    std::vector<float3> boundsMin_;
    std::vector<float3> boundsMax_;
    std::vector<int> left_;
    std::vector<int> right_;
    std::vector<int> parent_;
    std::vector<int> bodyIndex_;

public:
    // Clear and reserve
    void clear() noexcept {
        boundsMin_.clear();
        boundsMax_.clear();
        left_.clear();
        right_.clear();
        parent_.clear();
        bodyIndex_.clear();
    }

    void reserve(int capacity) noexcept {
        boundsMin_.reserve(capacity);
        boundsMax_.reserve(capacity);
        left_.reserve(capacity);
        right_.reserve(capacity);
        parent_.reserve(capacity);
        bodyIndex_.reserve(capacity);
    }

    // Getters using std::span (zero-copy, SafeC++)
    [[nodiscard]] std::span<const float3> boundsMin() const noexcept { return boundsMin_; }
    [[nodiscard]] std::span<const float3> boundsMax() const noexcept { return boundsMax_; }
    [[nodiscard]] std::span<const int> left() const noexcept { return left_; }
    [[nodiscard]] std::span<const int> right() const noexcept { return right_; }
    [[nodiscard]] std::span<const int> parent() const noexcept { return parent_; }
    [[nodiscard]] std::span<const int> bodyIndex() const noexcept { return bodyIndex_; }

    // Node operations
    [[nodiscard]] int size() const noexcept { return static_cast<int>(boundsMin_.size()); }
    [[nodiscard]] const float3& boundsMin(int index) const noexcept { return boundsMin_[index]; }
    [[nodiscard]] const float3& boundsMax(int index) const noexcept { return boundsMax_[index]; }
    [[nodiscard]] const int& bodyIndex(int index) const noexcept { return bodyIndex_[index]; }
    [[nodiscard]] int countLeaves() const noexcept;
    [[nodiscard]] int root() const noexcept { return 0; }

    // Implicit (heap) layout support: resize to an exact node count and write nodes by
    // index. Children of node i are 2i+1 / 2i+2 by construction - no left/right arrays.
    void resizeImplicit(int nodeCount) noexcept {
        boundsMin_.resize(nodeCount);
        boundsMax_.resize(nodeCount);
        bodyIndex_.resize(nodeCount, -1);
    }
    void setImplicitLeaf(int slot, const float3& min, const float3& max, int bodyIdx) noexcept {
        boundsMin_[slot] = min;
        boundsMax_[slot] = max;
        bodyIndex_[slot] = bodyIdx;
    }
    void setImplicitInternalBounds(int slot, const float3& min, const float3& max) noexcept {
        boundsMin_[slot] = min;
        boundsMax_[slot] = max;
    }

    // Node creation (SafeC++: returns int, no nullptr)
    int createLeaf(const float3& min, const float3& max, int bodyIdx) noexcept;
    int createInternal(const float3& min, const float3& max, int leftIdx, int rightIdx, int parentIdx = -1) noexcept;
    void updateBounds(int nodeId, const float3& min, const float3& max) noexcept;
};

} // namespace nodes

} // namespace bvh
