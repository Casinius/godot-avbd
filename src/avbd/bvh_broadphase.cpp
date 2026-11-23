/* Morton LBVH broad phase; bvh_broadphase.hpp documents what the port changed and why. */

#include "bvh_broadphase.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>

#include "avbd/solver.h"
#include "avbd/broadphase_pair.hpp"

namespace avbd {

namespace {

// 21 bits per axis, interleaved; the same encoding the reference builder uses.
inline uint64_t expandBits21(uint32_t v) noexcept
{
    uint64_t x = v & 0x1fffffu;
    x = (x | (x << 32)) & 0x1f00000000ffffull;
    x = (x | (x << 16)) & 0x1f0000ff0000ffull;
    x = (x | (x << 8)) & 0x100f00f00f00f00full;
    x = (x | (x << 4)) & 0x10cfcfcfcfcfcfcfull;
    x = (x | (x << 2)) & 0x1249249249249249ull;
    return x;
}

inline uint64_t morton3D(float x, float y, float z, float lo, float scale) noexcept
{
    const auto q = [lo, scale](float v) {
        return static_cast<uint32_t>(std::clamp((v - lo) * scale, 0.0f, 2097151.0f));
    };
    return expandBits21(q(x)) | (expandBits21(q(y)) << 1) | (expandBits21(q(z)) << 2);
}

} // namespace

int LbvhBroadPhase::leaf(int bodyIndex, const float3 &mn, const float3 &mx)
{
    const int id = static_cast<int>(nodes.size());
    nodes.push_back(Node{mn, mx, -1, bodyIndex});
    return id;
}

int LbvhBroadPhase::internal(int left, int right)
{
    const int id = static_cast<int>(nodes.size());
    const Node &a = nodes[left];
    const Node &b = nodes[right];
    nodes.push_back(Node{
            float3{std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y), std::min(a.min.z, b.min.z)},
            float3{std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y), std::max(a.max.z, b.max.z)},
            left, right});
    return id;
}

// Topology only: leaves are created from the morton-sorted order, then merged bottom-up.
// Bounds are filled by the refit that follows, so this stays a pure shape decision.
int LbvhBroadPhase::buildTopology(std::span<const float3> aabbCenter)
{
    // The body count comes from the caller's spans, not from `order`: `order` is scratch
    // that the previous build left sized, and on the first call it is empty.
    const int n = static_cast<int>(aabbCenter.size());
    if (n == 0)
        return -1;

    // Scene bounds over the centres, then one morton code per body.
    float3 lo{std::numeric_limits<float>::max(), std::numeric_limits<float>::max(),
            std::numeric_limits<float>::max()};
    float3 hi{std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(),
            std::numeric_limits<float>::lowest()};
    for (int i = 0; i < n; ++i)
    {
        const float3 &c = aabbCenter[i];
        lo = float3{std::min(lo.x, c.x), std::min(lo.y, c.y), std::min(lo.z, c.z)};
        hi = float3{std::max(hi.x, c.x), std::max(hi.y, c.y), std::max(hi.z, c.z)};
    }
    const float3 extent{hi.x - lo.x, hi.y - lo.y, hi.z - lo.z};
    const float scale = 2097151.0f / std::max({extent.x, extent.y, extent.z, 1.0e-6f});

    codes.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        codes[i] = morton3D(aabbCenter[i].x, aabbCenter[i].y, aabbCenter[i].z, lo.x, scale);

    // Sort indices by (code, index): equal codes must keep a fixed order or the topology
    // would depend on the sort's internals.
    order.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [this](int a, int b) {
        if (codes[a] != codes[b])
            return codes[a] < codes[b];
        return a < b;
    });

    nodes.clear();
    nodes.reserve(static_cast<size_t>(2 * n - 1));

    // Leaves in sorted order, merged bottom-up as a range walk. The open-range stack is
    // iterative: recursion would go log2(n) deep per level with a heap allocation each step.
    std::vector<int> level;
    level.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
        level.push_back(leaf(order[i], float3{}, float3{}));

    while (level.size() > 1)
    {
        std::vector<int> next;
        next.reserve(level.size() / 2 + 1);
        for (size_t i = 0; i + 1 < level.size(); i += 2)
            next.push_back(internal(level[i], level[i + 1]));
        if (level.size() % 2 == 1)
            next.push_back(level.back());
        level.swap(next);
    }
    return level[0];
}

int LbvhBroadPhase::refit(std::span<const float3> aabbHalfExt,
        std::span<const float3> aabbCenter)
{
    // Build order puts children before parents only within a pair, so refit is a post-order
    // walk over stored child indices, iterative (same reason as the build).
    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(root);
    // `visited` marks a node whose children are already updated.
    std::vector<uint8_t> done(nodes.size(), 0);

    while (!stack.empty())
    {
        const int id = stack.back();
        const int left = nodes[id].left;
        const int right = nodes[id].right;
        if (left < 0)
        {
            // Leaf: bounds come straight from the body's swept AABB.
            const int body = nodes[id].right;
            const float3 c = aabbCenter[body];
            const float3 h = aabbHalfExt[body];
            nodes[id].min = float3{c.x - h.x, c.y - h.y, c.z - h.z};
            nodes[id].max = float3{c.x + h.x, c.y + h.y, c.z + h.z};
            // Mark it done: a leaf has no children to wait for, and without the flag the
            // parent would re-push it forever.
            done[static_cast<size_t>(id)] = 1;
            stack.pop_back();
            continue;
        }
        if (done[static_cast<size_t>(left)] && done[static_cast<size_t>(right)])
        {
            const Node &a = nodes[left];
            const Node &b = nodes[right];
            nodes[id].min = float3{std::min(a.min.x, b.min.x), std::min(a.min.y, b.min.y),
                    std::min(a.min.z, b.min.z)};
            nodes[id].max = float3{std::max(a.max.x, b.max.x), std::max(a.max.y, b.max.y),
                    std::max(a.max.z, b.max.z)};
            done[static_cast<size_t>(id)] = 1;
            stack.pop_back();
            continue;
        }
        if (!done[static_cast<size_t>(left)])
            stack.push_back(left);
        if (!done[static_cast<size_t>(right)])
            stack.push_back(right);
    }
    return root;
}

void LbvhBroadPhase::buildAndCollectPairs(std::span<Rigid *const> bodies,
        std::vector<uint64_t> &out, std::span<const float3> aabbHalfExt,
        std::span<const float3> aabbCenter, detail::JobPool &pool)
{
    using Clock = std::chrono::steady_clock;
    const auto tStart = Clock::now();
    (void)pool; // the tree build and the traversal are serial by construction
    const int count = static_cast<int>(bodies.size());
    out.clear();
    if (count < 2)
    {
        root = -1;
        treeBodies.clear();
        nodeCount = 0;
        lastTreeMs = 0.0;
        lastEmitMs = 0.0;
        return;
    }

    // A refit is legal only while every leaf still addresses the same body.
    const bool sameBodySet = root >= 0 && treeBodies.size() == bodies.size() &&
            std::equal(bodies.begin(), bodies.end(), treeBodies.begin());
    if (!sameBodySet || age >= refitBudget || nodes.size() < 2)
    {
        treeBodies.assign(bodies.begin(), bodies.end());
        root = buildTopology(aabbCenter);
        age = 0;
        ++builds;
        refit(aabbHalfExt, aabbCenter);
    }
    else
    {
        refit(aabbHalfExt, aabbCenter);
        ++age;
        ++refits;
    }
    nodeCount = static_cast<int>(nodes.size());
#ifdef AVBD_LBVH_DEBUG
    {
        int leaves = 0, badLeaves = 0;
        for (size_t i = 0; i < nodes.size(); ++i)
            if (nodes[i].left < 0) {
                ++leaves;
                if (nodes[i].min.x == 0.0f && nodes[i].max.x == 0.0f) ++badLeaves;
            }
        std::fprintf(stderr, "[lbvh] nodes=%zu leaves=%d badLeaves=%d root=%d rootMin=(%.2f %.2f %.2f) rootMax=(%.2f %.2f %.2f) builds=%d refits=%d\n",
                nodes.size(), leaves, badLeaves, root, nodes[root].min.x, nodes[root].min.y,
                nodes[root].min.z, nodes[root].max.x, nodes[root].max.y, nodes[root].max.z, builds, refits);
    }
#endif

    const auto tEmit = Clock::now();
    lastTreeMs = std::chrono::duration<double, std::milli>(tEmit - tStart).count();

    // Tree-vs-tree self query. A node pair enters the stack only if the two bounds
    // overlap, so every entry is a live pair; a leaf pair is the answer.
    pending.clear();
    pending.reserve(128);
    pending.push_back({root, root});

    const auto overlaps = [this](int a, int b) {
        const Node &x = nodes[a];
        const Node &y = nodes[b];
        return !(y.min.x > x.max.x || y.max.x < x.min.x ||
                y.min.y > x.max.y || y.max.y < x.min.y ||
                y.min.z > x.max.z || y.max.z < x.min.z);
    };
    const auto push = [this, &overlaps](int a, int b) {
        if (a == b)
        {
            if (nodes[a].left >= 0) // a leaf self-pair carries no information
                pending.push_back({a, b});
            return;
        }
        if (overlaps(a, b))
            pending.push_back({a, b});
    };

    while (!pending.empty())
    {
        const std::array<int, 2> ab = pending.back();
        pending.pop_back();
        const int a = ab[0];
        const int b = ab[1];

        if (a == b)
        {
            const int left = nodes[a].left;
            const int right = nodes[a].right;
            push(left, left);
            push(left, right);
            push(right, right);
            continue;
        }

        const bool aLeaf = nodes[a].left < 0;
        const bool bLeaf = nodes[b].left < 0;
        if (aLeaf && bLeaf)
        {
            out.push_back(packPair(nodes[a].right, nodes[b].right));
            continue;
        }
        if (aLeaf)
        {
            push(a, nodes[b].left);
            push(a, nodes[b].right);
        }
        else if (bLeaf)
        {
            push(nodes[a].left, b);
            push(nodes[a].right, b);
        }
        else
        {
            push(nodes[a].left, nodes[b].left);
            push(nodes[a].left, nodes[b].right);
            push(nodes[a].right, nodes[b].left);
            push(nodes[a].right, nodes[b].right);
        }
    }
    lastEmitMs = std::chrono::duration<double, std::milli>(Clock::now() - tEmit).count();
}

} // namespace avbd
