/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, and distribute this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#pragma once

#include <limits>
#include <span>
#include <vector>

#include "avbd/maths.h"

namespace bvh {

using avbd::float3;

namespace nodes {

// One node of the implicit-heap LBVH.
//
// The bounds and the liveness of the subtree live together on purpose. The pair traversal issues tens
// of thousands of node-pair tests per step and it is the single hottest loop outside the solver, and
// the old layout spread each node over three parallel arrays (boundsMin, boundsMax, a validity byte)
// - three cache streams per test, where this is one line per node. An empty subtree (the right tail
// when the leaf count is not a power of two) is encoded in the bounds themselves, min > max on every
// axis, so there is no fourth field and no separate array to keep in step with it.
struct Node
{
    float3 min{emptyMin()};
    float3 max{emptyMax()};

    // One compare, not six: an empty subtree is written with min = +inf / max = -inf on *every* axis,
    // so the x pair decides it, and a real box always satisfies min.x <= max.x (degenerate boxes are
    // equal, never inverted). The traversal runs this twice per test, tens of thousands of times a
    // step, which is exactly why it is one compare.
    [[nodiscard]] bool valid() const noexcept { return min.x() <= max.x(); }

    // +inf / -inf, so `valid()` is false and any union with a live node is the live node's bounds.
    [[nodiscard]] static float3 emptyMin() noexcept
    {
        constexpr float inf = std::numeric_limits<float>::infinity();
        return float3{inf, inf, inf};
    }
    [[nodiscard]] static float3 emptyMax() noexcept
    {
        constexpr float inf = std::numeric_limits<float>::infinity();
        return float3{-inf, -inf, -inf};
    }
};

class NodeStorage {
public:
    void clear() noexcept {
        nodes_.clear();
        bodyIndex_.clear();
    }

    void reserve(int capacity) noexcept {
        nodes_.reserve(static_cast<size_t>(capacity));
        bodyIndex_.reserve(static_cast<size_t>(capacity));
    }

    [[nodiscard]] int size() const noexcept { return static_cast<int>(nodes_.size()); }

    // The whole node array, for a traversal that hoists the pointer out of its inner loop.
    [[nodiscard]] std::span<const Node> nodes() const noexcept { return nodes_; }
    [[nodiscard]] const Node &node(int index) const noexcept { return nodes_[static_cast<size_t>(index)]; }
    [[nodiscard]] const float3 &boundsMin(int index) const noexcept { return nodes_[static_cast<size_t>(index)].min; }
    [[nodiscard]] const float3 &boundsMax(int index) const noexcept { return nodes_[static_cast<size_t>(index)].max; }
    [[nodiscard]] bool valid(int index) const noexcept { return nodes_[static_cast<size_t>(index)].valid(); }
    [[nodiscard]] const int &bodyIndex(int index) const noexcept { return bodyIndex_[static_cast<size_t>(index)]; }
    [[nodiscard]] int root() const noexcept { return 0; }

    // Implicit (heap) layout support: resize to an exact node count and write nodes by
    // index. Children of node i are 2i+1 / 2i+2 by construction - no left/right arrays.
    void resizeImplicit(int nodeCount) noexcept {
        nodes_.resize(static_cast<size_t>(nodeCount)); // value-initialised: every node starts invalid
        bodyIndex_.assign(static_cast<size_t>(nodeCount), -1);
    }
    void setImplicitLeaf(int slot, const float3 &min, const float3 &max, int bodyIndex) noexcept {
        const size_t index = static_cast<size_t>(slot);
        nodes_[index].min = min;
        nodes_[index].max = max;
        bodyIndex_[index] = bodyIndex;
    }
    void setImplicitInternalBounds(int slot, const float3 &min, const float3 &max, bool subtreeValid) noexcept {
        const size_t index = static_cast<size_t>(slot);
        nodes_[index].min = subtreeValid ? min : Node::emptyMin();
        nodes_[index].max = subtreeValid ? max : Node::emptyMax();
    }

private:
    std::vector<Node> nodes_;
    // Leaf payload, read once per emitted pair rather than per test, so it stays in its own array.
    std::vector<int> bodyIndex_;
};

} // namespace nodes

} // namespace bvh
