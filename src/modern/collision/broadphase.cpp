#include "collision/broadphase.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace collision {

namespace {
[[nodiscard]] inline bool aabb_overlap(const bvh::Node &a, const bvh::Node &b) noexcept {
  return !(b.lo.x > a.hi.x || b.lo.y > a.hi.y || b.lo.z > a.hi.z || a.lo.x > b.hi.x ||
           a.lo.y > b.hi.y || a.lo.z > b.hi.z);
}
} // namespace

std::uint64_t Broadphase::expand_bits21(std::uint32_t v) noexcept {
  std::uint64_t x = v & 0x1fffffu;
  x               = (x | (x << 32)) & 0x1f00000000ffffull;
  x               = (x | (x << 16)) & 0x1f0000ff0000ffull;
  x               = (x | (x << 8)) & 0x100f00f00f00f00full;
  x               = (x | (x << 4)) & 0x10cfcfcfcfcfcfcfull;
  x               = (x | (x << 2)) & 0x1249249249249249ull;
  return x;
}

std::uint64_t Broadphase::morton3d(Scalar_t x, Scalar_t y, Scalar_t z) noexcept {
  const auto qx = static_cast<std::uint32_t>(std::clamp(x, Scalar_t(0), Scalar_t(2097151)));
  const auto qy = static_cast<std::uint32_t>(std::clamp(y, Scalar_t(0), Scalar_t(2097151)));
  const auto qz = static_cast<std::uint32_t>(std::clamp(z, Scalar_t(0), Scalar_t(2097151)));
  return expand_bits21(qx) | (expand_bits21(qy) << 1) | (expand_bits21(qz) << 2);
}

void Broadphase::update(std::span<const Aabb> aabbs) {
  const std::size_t count = aabbs.size();
  nodes_.clear();
  pairs_.clear();
  leaf_count_ = 0;
  first_leaf_ = 0;
  if (count < 2)
    return;

  // 场景中心范围（不是 AABB 范围）用于 21-bit 量化。
  Vec3 lo{std::numeric_limits<Scalar_t>::max()};
  Vec3 hi{std::numeric_limits<Scalar_t>::lowest()};
  for (const Aabb &box : aabbs) {
    if (!box.valid())
      continue;
    const Vec3 c = box.center();
    lo           = Vec3{std::min(lo.x, c.x), std::min(lo.y, c.y), std::min(lo.z, c.z)};
    hi           = Vec3{std::max(hi.x, c.x), std::max(hi.y, c.y), std::max(hi.z, c.z)};
  }
  if (hi.x < lo.x)
    return; // 没有任何有效盒

  const Vec3     ext     = hi - lo;
  const Scalar_t max_ext = std::max({ext.x, ext.y, ext.z, Scalar_t(1e-6)});
  const Scalar_t scale   = Scalar_t(2097151) / max_ext;

  coded_.clear();
  coded_.reserve(count);
  for (std::size_t b = 0; b < count; ++b) {
    const Aabb &box = aabbs[b];
    if (!box.valid())
      continue;
    const Vec3 c = box.center();
    coded_.push_back(
        MortonEntry{morton3d((c.x - lo.x) * scale, (c.y - lo.y) * scale, (c.z - lo.z) * scale), b});
  }
  // 同码用 body index 破平局：确定性不依赖排序实现。
  std::sort(coded_.begin(), coded_.end(), [](const MortonEntry &x, const MortonEntry &y) {
    return x.code < y.code || (x.code == y.code && x.body < y.body);
  });

  const int n      = static_cast<int>(coded_.size());
  int       padded = 1;
  while (padded < n)
    padded <<= 1;
  leaf_count_ = n;
  first_leaf_ = padded - 1;

  nodes_.resize_implicit(padded * 2 - 1);

  // 空树/单叶时全部节点无效；但单叶也需要写叶子供查询。
  Aabb all;
  for (int i = 0; i < n; ++i) {
    const Aabb &box = aabbs[coded_[static_cast<std::size_t>(i)].body];
    nodes_.set_leaf(first_leaf_ + i, box.lo, box.hi,
                    static_cast<std::int32_t>(coded_[static_cast<std::size_t>(i)].body));
    all.merge(box);
  }

  const auto subtree_valid = [&](int slot) noexcept -> bool {
    if (slot >= first_leaf_)
      return (slot - first_leaf_) < leaf_count_;
    return nodes_.valid(slot);
  };
  const auto subtree_bounds = [&](int slot, Vec3 &r_lo, Vec3 &r_hi) noexcept {
    if (slot >= first_leaf_) {
      const bvh::Node &node = nodes_.node(slot);
      r_lo                  = node.lo;
      r_hi                  = node.hi;
    } else {
      const bvh::Node &node = nodes_.node(slot);
      r_lo                  = node.lo;
      r_hi                  = node.hi;
    }
  };

  // 自底向上聚合。右到左保证孩子先算好。
  for (int slot = padded - 2; slot >= 0; --slot) {
    const int  l  = 2 * slot + 1;
    const int  r  = 2 * slot + 2;
    const bool lv = (l < padded * 2 - 1) && subtree_valid(l);
    const bool rv = (r < padded * 2 - 1) && subtree_valid(r);
    if (!lv && !rv) {
      nodes_.set_internal(slot, all.lo, all.hi, false);
      continue;
    }
    Vec3 llo{}, lhi{}, rlo{}, rhi{};
    if (lv)
      subtree_bounds(l, llo, lhi);
    if (rv)
      subtree_bounds(r, rlo, rhi);
    Vec3 mn, mx;
    if (lv && rv) {
      mn = Vec3{std::min(llo.x, rlo.x), std::min(llo.y, rlo.y), std::min(llo.z, rlo.z)};
      mx = Vec3{std::max(lhi.x, rhi.x), std::max(lhi.y, rhi.y), std::max(lhi.z, rhi.z)};
    } else if (lv) {
      mn = llo;
      mx = lhi;
    } else {
      mn = rlo;
      mx = rhi;
    }
    nodes_.set_internal(slot, mn, mx, true);
  }

  // 自碰撞：每个内部节点只测一次左右孩子，叶子对收进 (min,max) 升序。
  const std::span<const bvh::Node> nodes            = nodes_.nodes();
  const auto                       collide_subtrees = [&](auto &&self, int sa, int sb) -> void {
    const bvh::Node &na = nodes[static_cast<std::size_t>(sa)];
    const bvh::Node &nb = nodes[static_cast<std::size_t>(sb)];
    if (!na.valid() || !nb.valid())
      return;
    if (!aabb_overlap(na, nb))
      return;

    const bool a_leaf = sa >= first_leaf_;
    const bool b_leaf = sb >= first_leaf_;
    if (a_leaf && b_leaf) {
      const std::int32_t ia = nodes_.body_index(sa);
      const std::int32_t ib = nodes_.body_index(sb);
      if (ia >= 0 && ib >= 0 && ia != ib) {
        const std::size_t a = static_cast<std::size_t>(std::min(ia, ib));
        const std::size_t b = static_cast<std::size_t>(std::max(ia, ib));
        pairs_.push_back(BodyPair{a, b});
      }
      return;
    }
    if (!a_leaf) {
      self(self, 2 * sa + 1, sb);
      self(self, 2 * sa + 2, sb);
    } else {
      self(self, sa, 2 * sb + 1);
      self(self, sa, 2 * sb + 2);
    }
  };
  const auto self_traverse = [&](auto &&self, int slot) -> void {
    if (slot >= first_leaf_)
      return;
    const int l = 2 * slot + 1;
    const int r = 2 * slot + 2;
    collide_subtrees(collide_subtrees, l, r);
    self(self, l);
    self(self, r);
  };
  self_traverse(self_traverse, 0);

  std::sort(pairs_.begin(), pairs_.end());
  pairs_.erase(std::unique(pairs_.begin(), pairs_.end()), pairs_.end());
}

} // namespace collision
