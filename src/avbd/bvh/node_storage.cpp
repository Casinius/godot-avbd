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

int NodeStorage::createLeaf(const float3& min, const float3& max, int bodyIndex) noexcept {
    // Allocate a new node
    int nodeId = static_cast<int>(bodyIndex_.size());

    boundsMin_.push_back(min);
    boundsMax_.push_back(max);
    left_.push_back(-1);   // No children for leaf
    right_.push_back(-1);
    parent_.push_back(-1);
    bodyIndex_.push_back(bodyIndex);

    return nodeId;
}

int NodeStorage::createInternal(const float3& min, const float3& max, int left, int right, int parent) noexcept {
    // Allocate a new node
    int nodeId = static_cast<int>(bodyIndex_.size());

    boundsMin_.push_back(min);
    boundsMax_.push_back(max);
    left_.push_back(left);
    right_.push_back(right);
    parent_.push_back(parent);
    bodyIndex_.push_back(-1);  // No body index for internal node

    return nodeId;
}

void NodeStorage::updateBounds(int nodeId, const float3& min, const float3& max) noexcept {
    if (nodeId >= 0 && nodeId < static_cast<int>(bodyIndex_.size())) {
        boundsMin_[nodeId] = min;
        boundsMax_[nodeId] = max;
    }
}

int NodeStorage::addNode(const float3& min, const float3& max, int left, int right, int parent, int bodyIndex) noexcept {
    // Either create a new leaf or internal node based on bodyIndex
    if (bodyIndex >= 0) {
        return createLeaf(min, max, bodyIndex);
    } else {
        return createInternal(min, max, left, right, parent);
    }
}

} // namespace bvh::nodes
