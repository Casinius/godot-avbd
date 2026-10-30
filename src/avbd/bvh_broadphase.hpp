/*
 * Morton LBVH broad phase - the alternative to the hierarchical hashed grid
 * (`SpatialHash`), ported from `main:src/avbd/bvh/`. Selected by `broadPhaseMode = 1`;
 * both implementations are fed the same swept AABBs and feed the same downstream.
 *
 * The port keeps the algorithm and changes three things, all of which are needed for the
 * comparison to mean anything:
 *
 *   1. Eigen -> linalg (the `main` tree is the pre-port branch).
 *   2. Leaves are sized from the solver's **swept world AABBs** (`half`/`centre`, the same
 *      spans `SpatialHash::buildAndCollectPairs` takes) instead of `body->radius` bounding
 *      spheres. A sphere-vs-AABB difference is a bigger effect than the tree-vs-grid
 *      difference, so sharing the input is the whole point.
 *   3. `NodeStorage` is folded into one class. The original spread the arrays over a
 *      separate translation unit with accessors per node; that costs nothing measurable but
 *      makes the file count four instead of two.
 *
 * Output contract matches the hash: packed `packPair(i, j)` candidates, unsorted and with
 * each overlapping pair emitted once. The caller sorts, dedupes and applies the exact
 * AABB/layer/constraint filter, so both broad phases feed the identical downstream work.
 *
 * Determinism: the build sorts by (morton code, body index), the traversal pushes children
 * in a fixed order, and the pair list is integer-sorted by the caller. Unlike the hash's
 * walk, this one is serial - that is a real property of the two designs, not an oversight.
 */

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

    // Refits allowed before the topology is rebuilt from fresh Morton codes. A refit is
    // cheap (one post-order walk) but the topology ages: bodies that move far from where
    // their code said they were leave the tree unbalanced and the overlap tests worse.
    int refitBudget = 16;

private:
    // Node arrays, index-addressed. `body` >= 0 marks a leaf; the internal nodes keep
    // explicit child indices, so there is no padding and every node is reachable.
    std::vector<float3> nodeMin;
    std::vector<float3> nodeMax;
    std::vector<int> nodeLeft;
    std::vector<int> nodeRight;
    std::vector<int> nodeBody;

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
