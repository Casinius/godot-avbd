/*
 * Morton LBVH broad phase. See bvh_broadphase.hpp for what the port changed and why.
 */

#include "bvh_broadphase.hpp"

#include <algorithm>
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
    const int id = static_cast<int>(nodeMin.size());
    nodeMin.push_back(mn);
    nodeMax.push_back(mx);
    nodeLeft.push_back(-1);
    nodeRight.push_back(-1);
    nodeBody.push_back(bodyIndex);
    return id;
}

int LbvhBroadPhase::internal(int left, int right)
{
    const int id = static_cast<int>(nodeMin.size());
    nodeMin.push_back(float3{std::min(nodeMin[left].x, nodeMin[right].x),
            std::min(nodeMin[left].y, nodeMin[right].y),
            std::min(nodeMin[left].z, nodeMin[right].z)});
    nodeMax.push_back(float3{std::max(nodeMax[left].x, nodeMax[right].x),
            std::max(nodeMax[left].y, nodeMax[right].y),
            std::max(nodeMax[left].z, nodeMax[right].z)});
    nodeLeft.push_back(left);
    nodeRight.push_back(right);
    nodeBody.push_back(-1);
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

    nodeMin.clear();
    nodeMax.clear();
    nodeLeft.clear();
    nodeRight.clear();
    nodeBody.clear();
    nodeMin.reserve(static_cast<size_t>(2 * n - 1));
    nodeMax.reserve(static_cast<size_t>(2 * n - 1));
    nodeLeft.reserve(static_cast<size_t>(2 * n - 1));
    nodeRight.reserve(static_cast<size_t>(2 * n - 1));
    nodeBody.reserve(static_cast<size_t>(2 * n - 1));

    // Leaves first (in sorted order), so the merge below is a simple range walk. A range
    // of leaves [begin, end) is replaced by the internal node that spans it; the stack of
    // open ranges is rebuilt without recursion (a 4k-body scene would otherwise recurse
    // log2(n) deep per level with a heap allocation at every step).
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
    // Every node was created in a build order that puts children before parents inside a
    // pair but not globally, so the walk is a post-order over the stored child indices.
    // Iterative with an explicit stack: same reason as the build.
    std::vector<int> stack;
    stack.reserve(64);
    stack.push_back(root);
    // `visited` marks a node whose children are already updated.
    std::vector<uint8_t> done(nodeMin.size(), 0);

    while (!stack.empty())
    {
        const int id = stack.back();
        const int left = nodeLeft[id];
        const int right = nodeRight[id];
        if (left < 0)
        {
            // Leaf: bounds come straight from the body's swept AABB.
            const int body = nodeBody[id];
            const float3 c = aabbCenter[body];
            const float3 h = aabbHalfExt[body];
            nodeMin[id] = float3{c.x - h.x, c.y - h.y, c.z - h.z};
            nodeMax[id] = float3{c.x + h.x, c.y + h.y, c.z + h.z};
            // Mark it done: a leaf has no children to wait for, and without the flag the
            // parent would re-push it forever.
            done[static_cast<size_t>(id)] = 1;
            stack.pop_back();
            continue;
        }
        if (done[static_cast<size_t>(left)] && done[static_cast<size_t>(right)])
        {
            nodeMin[id] = float3{std::min(nodeMin[left].x, nodeMin[right].x),
                    std::min(nodeMin[left].y, nodeMin[right].y),
                    std::min(nodeMin[left].z, nodeMin[right].z)};
            nodeMax[id] = float3{std::max(nodeMax[left].x, nodeMax[right].x),
                    std::max(nodeMax[left].y, nodeMax[right].y),
                    std::max(nodeMax[left].z, nodeMax[right].z)};
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
    (void)pool; // the tree build and the traversal are serial by construction
    const int count = static_cast<int>(bodies.size());
    out.clear();
    if (count < 2)
    {
        root = -1;
        treeBodies.clear();
        nodeCount = 0;
        return;
    }

    // A refit is legal only while every leaf still addresses the same body.
    const bool sameBodySet = root >= 0 && treeBodies.size() == bodies.size() &&
            std::equal(bodies.begin(), bodies.end(), treeBodies.begin());
    if (!sameBodySet || age >= refitBudget || nodeMin.size() < 2)
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
    nodeCount = static_cast<int>(nodeMin.size());
#ifdef AVBD_LBVH_DEBUG
    {
        int leaves = 0, badLeaves = 0;
        for (size_t i = 0; i < nodeMin.size(); ++i)
            if (nodeBody[i] >= 0) {
                ++leaves;
                if (nodeMin[i].x == 0.0f && nodeMax[i].x == 0.0f) ++badLeaves;
            }
        std::fprintf(stderr, "[lbvh] nodes=%zu leaves=%d badLeaves=%d root=%d rootMin=(%.2f %.2f %.2f) rootMax=(%.2f %.2f %.2f) builds=%d refits=%d\n",
                nodeMin.size(), leaves, badLeaves, root, nodeMin[root].x, nodeMin[root].y,
                nodeMin[root].z, nodeMax[root].x, nodeMax[root].y, nodeMax[root].z, builds, refits);
    }
#endif

    // Tree-vs-tree self query. A node pair enters the stack only if the two bounds
    // overlap, so every entry is a live pair; a leaf pair is the answer.
    pending.clear();
    pending.reserve(128);
    pending.push_back({root, root});

    const auto overlaps = [this](int a, int b) {
        return !(nodeMin[b].x > nodeMax[a].x || nodeMax[b].x < nodeMin[a].x ||
                nodeMin[b].y > nodeMax[a].y || nodeMax[b].y < nodeMin[a].y ||
                nodeMin[b].z > nodeMax[a].z || nodeMax[b].z < nodeMin[a].z);
    };
    const auto push = [this, &overlaps](int a, int b) {
        if (a == b)
        {
            if (nodeLeft[a] >= 0) // a leaf self-pair carries no information
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
            const int left = nodeLeft[a];
            const int right = nodeRight[a];
            push(left, left);
            push(left, right);
            push(right, right);
            continue;
        }

        const bool aLeaf = nodeLeft[a] < 0;
        const bool bLeaf = nodeLeft[b] < 0;
        if (aLeaf && bLeaf)
        {
            out.push_back(packPair(nodeBody[a], nodeBody[b]));
            continue;
        }
        if (aLeaf)
        {
            push(a, nodeLeft[b]);
            push(a, nodeRight[b]);
        }
        else if (bLeaf)
        {
            push(nodeLeft[a], b);
            push(nodeRight[a], b);
        }
        else
        {
            push(nodeLeft[a], nodeLeft[b]);
            push(nodeLeft[a], nodeRight[b]);
            push(nodeRight[a], nodeLeft[b]);
            push(nodeRight[a], nodeRight[b]);
        }
    }
}

} // namespace avbd
