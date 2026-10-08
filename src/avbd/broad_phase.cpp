#include "avbd/broad_phase.hpp"
#include "avbd/solver.h"

#include <algorithm>
#include <numeric>

namespace avbd {

void BroadPhase::clear() noexcept {
    root = -1;
    valid = false;
    needsRebuild = true;
    framesSinceRebuild = 0;
    builtBodies.clear();
    builtPositions.clear();
    bodyPointers.clear();
    indices.clear();
    nodes.clear();
    pairs.clear();
    pending.clear();
}

static void computeBounds(const Rigid *body, float3 &min, float3 &max) noexcept {
    const float r = body->radius;
    min = {body->positionLin.x() - r, body->positionLin.y() - r, body->positionLin.z() - r};
    max = {body->positionLin.x() + r, body->positionLin.y() + r, body->positionLin.z() + r};
}

void BroadPhase::build(const std::vector<Rigid *> &bodies) {
    indices.resize(bodies.size());
    std::iota(indices.begin(), indices.end(), 0);
    bodyPointers.assign(bodies.begin(), bodies.end());
    bvh::builder::Builder builder(nodes, indices, bodyPointers, &scratch);
    root = builder.buildLBVH().value_or(-1);
    builtBodies = bodies;
    builtPositions.resize(bodies.size());
    for (std::size_t i = 0; i < bodies.size(); ++i)
        builtPositions[i] = bodies[i]->positionLin;
    needsRebuild = false;
    valid = true;
    framesSinceRebuild = 0;
}

void BroadPhase::rebuild(const std::vector<Rigid *> &bodies) {
    nodes.clear();
    build(bodies);
}

void BroadPhase::refitTree() {
    const int count = nodes.size();
    const auto refitNode = [&](auto &&self, int node) -> void {
        if (node < 0 || node >= count) return;
        if (nodes.isLeaf(node)) {
            const int index = nodes.bodyIndex(node);
            if (index >= 0 && index < static_cast<int>(builtBodies.size())) {
                float3 min, max;
                computeBounds(builtBodies[index], min, max);
                nodes.updateBounds(node, min, max);
            }
            return;
        }
        const int left = nodes.left()[node], right = nodes.right()[node];
        self(self, left); self(self, right);
        const float3 min = comp3<float3, float3>(nodes.boundsMin(left), nodes.boundsMin(right), [](auto x, auto y) { return std::min(x, y); });
        const float3 max = comp3<float3, float3>(nodes.boundsMax(left), nodes.boundsMax(right), [](auto x, auto y) { return std::max(x, y); });
        nodes.updateBounds(node, min, max);
    };
    refitNode(refitNode, root);
}

void BroadPhase::refit(const std::vector<Rigid *> &bodies) {
    if (builtBodies.size() != bodies.size()) { needsRebuild = true; return; }
    const float limit2 = rebuildDisplacement * rebuildDisplacement;
    for (std::size_t i = 0; i < bodies.size(); ++i)
        if (builtBodies[i] != bodies[i] || (bodies[i]->positionLin - builtPositions[i]).squaredNorm() > limit2) { needsRebuild = true; return; }
    refitTree();
    ++framesSinceRebuild;
    if (framesSinceRebuild >= rebuildInterval) needsRebuild = true;
}

const std::vector<std::pair<int, int>> &BroadPhase::collectPairs(const std::vector<Rigid *> &bodies) {
    if (bodies.size() < 2) { pairs.clear(); return pairs; }
    if (!valid || needsRebuild) rebuild(bodies); else refit(bodies);
    pairs.clear(); pending.clear();
    if (root < 0 || nodes.size() < 2) return pairs;
    const auto overlap = [this](int a, int b) {
        const auto amin = nodes.boundsMin(a), amax = nodes.boundsMax(a), bmin = nodes.boundsMin(b), bmax = nodes.boundsMax(b);
        return !(bmin.x() > amax.x() || bmin.y() > amax.y() || bmin.z() > amax.z() || bmax.x() < amin.x() || bmax.y() < amin.y() || bmax.z() < amin.z());
    };
    pending.emplace_back(root, root);
    while (!pending.empty()) {
        auto [a, b] = pending.back(); pending.pop_back();
        if (a == b) { if (!nodes.isLeaf(a)) { int l=nodes.left()[a], r=nodes.right()[a]; pending.emplace_back(l,l); pending.emplace_back(l,r); pending.emplace_back(r,r); } continue; }
        if (!overlap(a,b)) continue;
        const bool al=nodes.isLeaf(a), bl=nodes.isLeaf(b);
        if (al && bl) { int i=nodes.bodyIndex(a), j=nodes.bodyIndex(b); if (i > j) std::swap(i,j); pairs.emplace_back(i,j); }
        else if (al) { pending.emplace_back(a,nodes.left()[b]); pending.emplace_back(a,nodes.right()[b]); }
        else if (bl) { pending.emplace_back(nodes.left()[a],b); pending.emplace_back(nodes.right()[a],b); }
        else { int la=nodes.left()[a],ra=nodes.right()[a],lb=nodes.left()[b],rb=nodes.right()[b]; pending.emplace_back(la,lb); pending.emplace_back(la,rb); pending.emplace_back(ra,lb); pending.emplace_back(ra,rb); }
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}
} // namespace avbd
