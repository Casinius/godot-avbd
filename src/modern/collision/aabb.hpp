#pragma once

// ============================================================================
//  collision/aabb.hpp —— 轴对齐包围盒（纯值类型，无拥有关系）
// ----------------------------------------------------------------------------
//  空盒用 lo=+inf / hi=-inf 表示，valid() 只比较 x 一对分量就够了：真实盒永远
//  满足 lo.x <= hi.x，空盒则相反。这个编码让 BVH 的空子树不需要额外的标志位。
// ============================================================================

#include <algorithm>
#include <cmath>
#include <limits>

#include "phys_core/state.hpp"

namespace collision {

struct Aabb {
  Vec3 lo{empty_lo()};
  Vec3 hi{empty_hi()};

  [[nodiscard]] static constexpr Vec3 empty_lo() noexcept {
    constexpr Scalar_t inf = std::numeric_limits<Scalar_t>::infinity();
    return Vec3{inf, inf, inf};
  }
  [[nodiscard]] static constexpr Vec3 empty_hi() noexcept {
    constexpr Scalar_t inf = std::numeric_limits<Scalar_t>::infinity();
    return Vec3{-inf, -inf, -inf};
  }

  [[nodiscard]] bool valid() const noexcept { return lo.x <= hi.x; }

  void reset() noexcept {
    lo = empty_lo();
    hi = empty_hi();
  }

  void expand(Vec3 p) noexcept {
    lo = Vec3{std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = Vec3{std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }

  void merge(const Aabb &o) noexcept {
    if (!o.valid())
      return;
    lo = Vec3{std::min(lo.x, o.lo.x), std::min(lo.y, o.lo.y), std::min(lo.z, o.lo.z)};
    hi = Vec3{std::max(hi.x, o.hi.x), std::max(hi.y, o.hi.y), std::max(hi.z, o.hi.z)};
  }

  [[nodiscard]] Aabb expanded(Scalar_t margin) const noexcept {
    return Aabb{lo - Vec3{margin}, hi + Vec3{margin}};
  }

  [[nodiscard]] Vec3 center() const noexcept { return (lo + hi) * Scalar_t(0.5); }
  [[nodiscard]] Vec3 extent() const noexcept { return hi - lo; }
};

// 六分离轴测试。接触是允许的（<=），与参考实现的 broadphase 语义一致。
[[nodiscard]] inline bool overlaps(const Aabb &a, const Aabb &b) noexcept {
  return !(b.lo.x > a.hi.x || b.lo.y > a.hi.y || b.lo.z > a.hi.z || a.lo.x > b.hi.x ||
           a.lo.y > b.hi.y || a.lo.z > b.hi.z);
}

} // namespace collision
