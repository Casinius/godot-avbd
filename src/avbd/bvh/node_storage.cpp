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

#include "node_storage.hpp"

namespace bvh::nodes {

int NodeStorage::countLeaves() const noexcept {
    // Count how many internal nodes vs leaf nodes
    // For a simple implementation, we can just count based on bodyIndex size
    return static_cast<int>(bodyIndex_.size());
}

int NodeStorage::createLeaf(const float3& min, const float3& max, int bodyIdx) noexcept {
    // Allocate a new node
    int nodeId = static_cast<int>(bodyIndex_.size());

    boundsMin_.push_back(min);
    boundsMax_.push_back(max);
    left_.push_back(-1);   // No children for leaf
    right_.push_back(-1);
    parent_.push_back(-1);
    bodyIndex_.push_back(bodyIdx);

    return nodeId;
}

int NodeStorage::createInternal(const float3& min, const float3& max, int leftIdx, int rightIdx, int parentIdx) noexcept {
    // Allocate a new node
    int nodeId = static_cast<int>(bodyIndex_.size());

    boundsMin_.push_back(min);
    boundsMax_.push_back(max);
    left_.push_back(leftIdx);
    right_.push_back(rightIdx);
    parent_.push_back(parentIdx);
    bodyIndex_.push_back(-1);  // No body index for internal node

    return nodeId;
}

void NodeStorage::updateBounds(int nodeId, const float3& min, const float3& max) noexcept {
    if (nodeId >= 0 && nodeId < static_cast<int>(bodyIndex_.size())) {
        boundsMin_[nodeId] = min;
        boundsMax_[nodeId] = max;
    }
}

} // namespace bvh::nodes
