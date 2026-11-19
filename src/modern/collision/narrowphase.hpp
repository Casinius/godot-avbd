#pragma once

// ============================================================================
//  collision/narrowphase.hpp —— 形状对 → 接触点 + 约束基
// ----------------------------------------------------------------------------
//  【法线/基的约定（硬契约，别改）】
//   ContactPatch::basis.normal 是**约束基**的法线，方向从 B 指向 A：
//     ContactConstraint::evaluate = n·(xA - xB) - target
//   穿透时 n·(xA - xB) < 0，于是 C<0，梯度 {+n,-n} 把 A 沿 +n 推开。
//   `target` 由 contact_constraints.hpp 设为 -margin。
//
//  【输出契约】
//   1) 每个 ContactPoint 的 xA/xB 是世界空间见证点，xA 属于 A，xB 属于 B；
//   2) rA/rB 是相对各自质心的局部锚点；
//   3) feature_key 跨帧稳定，用于 manifold warm start；
//   4) 最多 kMaxContactsPerManifold 个点，顺序即 narrowphase 输出顺序；
//   5) 本函数不分配堆内存，不写任何全局状态。
// ============================================================================

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

#include "collision/shape.hpp"

namespace collision {

inline constexpr std::size_t kMaxContactsPerManifold = 8;

struct BasisFrame {
  Vec3 normal{0, 1, 0}; // B -> A（约束推开 A 的方向）
  Vec3 tangent0{1, 0, 0};
  Vec3 tangent1{0, 0, 1};
};

struct ContactPoint {
  std::uint64_t feature_key = 0;
  Vec3          rA{0}; // A 的局部锚点
  Vec3          rB{0}; // B 的局部锚点
  Vec3          xA{0}; // A 的世界见证点
  Vec3          xB{0}; // B 的世界见证点
  Scalar_t      depth = 0;
};

struct ContactPatch {
  std::size_t                                       body_a  = 0;
  std::size_t                                       body_b  = 0;
  std::size_t                                       proxy_a = 0;
  std::size_t                                       proxy_b = 0;
  BasisFrame                                        basis{};
  std::array<ContactPoint, kMaxContactsPerManifold> points{};
  std::size_t                                       count    = 0;
  Scalar_t                                          friction = Scalar_t(0.5);
  // 与 ContactEmitOptions::margin 相同（法线 target = -margin）。ManifoldStore 只
  // 用它换算“新鲜接触的几何间隙”，做释放判定，不参与窄相。
  Scalar_t margin = Scalar_t(0.005);
};

// 从“B->A 的约束法线”构造右手正交基：normal 为第一行，两条切线为另两行。
[[nodiscard]] BasisFrame make_basis(Vec3 normal_b_to_a) noexcept;

// 离散形状对检测。out 至少有 8 个槽位；返回实际接触点数。
// 只有 count>0 时 basis 才有意义。
std::size_t collide_shapes(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                           BasisFrame &basis);

// NarrowphaseModel 的具体实现：把自由函数包成一个可替换的策略对象。
// 以后加 CCD / GPU narrowphase 时，换这个对象即可，CollisionWorld 的流程不变。
struct DiscreteNarrowphase {
  std::size_t collide(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                      BasisFrame &basis) const {
    return collide_shapes(a, b, out, basis);
  }
};

} // namespace collision
