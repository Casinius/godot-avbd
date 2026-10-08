#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "avbd/bvh/bvh.hpp"
#include "avbd/bvh/node_storage.hpp"

namespace avbd {
struct Rigid;

// Incremental BVH broad phase. Body indices in returned pairs are stable for the
// duration of a call and follow the existing Morton/tree traversal contract.
class BroadPhase {
public:
    void clear() noexcept;
    void rebuild(const std::vector<Rigid *> &bodies);
    void refit(const std::vector<Rigid *> &bodies);
    [[nodiscard]] const std::vector<std::pair<int, int>> &collectPairs(
            const std::vector<Rigid *> &bodies);

private:
    void build(const std::vector<Rigid *> &bodies);
    void refitTree();

    std::size_t rebuildInterval = 30;
    float rebuildDisplacement = 1.0f;
    int root = -1;
    bool valid = false;
    bool needsRebuild = true;
    std::size_t framesSinceRebuild = 0;
    std::vector<Rigid *> builtBodies;
    std::vector<float3> builtPositions;
    std::vector<const Rigid *> bodyPointers;
    std::vector<int> indices;
    bvh::builder::LBVHScratch scratch;
    bvh::nodes::NodeStorage nodes;
    std::vector<std::pair<int, int>> pairs;
    std::vector<std::pair<int, int>> pending;
};
} // namespace avbd
