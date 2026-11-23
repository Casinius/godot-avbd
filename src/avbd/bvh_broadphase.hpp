// Morton LBVH, sole broad phase: ported to linalg, swept-AABB leaves, folded NodeStorage.
// Serial/deterministic ((code,index) sort, fixed child order; `pool` reserved). Output: packPair(i,j) once each, unsorted; caller sorts/filters.

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "avbd/maths.h"

namespace avbd {

struct Rigid;
namespace detail {
class JobPool;
}

class LbvhBroadPhase
{
public:
    // Build or refit, then emit every candidate pair whose swept AABBs can touch.
    void buildAndCollectPairs(std::span<Rigid *const> bodies, std::vector<uint64_t> &out,
            std::span<const float3> aabbHalfExt, std::span<const float3> aabbCenter,
            detail::JobPool &pool);

    // Statistics of the last build (diagnostics only).
    int nodeCount = 0;
    int builds = 0;
    int refits = 0;

    // Last call split (milliseconds, diagnostics only): topology build + refit vs the
    // tree-vs-tree candidate emission that follows. Their sum is the caller's broadIndex.
    double lastTreeMs = 0.0;
    double lastEmitMs = 0.0;

    // Refits allowed before rebuilding from fresh Morton codes. A refit is cheap (one post-order
    // walk) but the topology ages: far-moved bodies unbalance the tree and worsen overlap tests.
    int refitBudget = 16;

private:
    // One 32-byte node: min(12) + max(12) + two ints(8). A leaf sets left=-1 and keeps
    // the body index in right; an internal node keeps both child indices. Siblings are
    // allocated consecutively, so one 64B cache line holds a child pair during descent.
    struct Node
    {
        float3 min;
        float3 max;
        int left;  // < 0: leaf (body index in right); >= 0: left child
        int right; // leaf: body index; internal: right child
    };
    static_assert(sizeof(Node) == 32, "LBVH node must stay half a cache line");
    std::vector<Node> nodes;

    // The exact body list the topology was built from, to know when a refit is still legal.
    std::vector<Rigid *> treeBodies;
    int root = -1;
    int age = 0;

    // Traversal stack, reused every step.
    std::vector<std::array<int, 2>> pending;
    // Scratch for the build, reused every step.
    std::vector<uint64_t> codes;
    std::vector<int> order;

    [[nodiscard]] int leaf(int bodyIndex, const float3 &mn, const float3 &mx);
    [[nodiscard]] int internal(int left, int right);
    // Recompute every leaf's bounds from `aabbHalfExt`/`aabbCenter` and then every internal
    // node's from its children (post-order). Returns the root.
    int refit(std::span<const float3> aabbHalfExt, std::span<const float3> aabbCenter);
    int buildTopology(std::span<const float3> aabbCenter);
};

} // namespace avbd
