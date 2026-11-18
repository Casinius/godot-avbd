#pragma once

// ============================================================================
//  collision/collision_digest.hpp —— 碰撞层的“可观测摘要”
// ----------------------------------------------------------------------------
//  【用途】给 collision 层单独建立一份可记录的指纹，和 State/ConstraintList 的
//  主基线解耦。它回答的是：“碰撞层的输出（manifold / 接触约束 / 统计）变了没有？”
//  而不是“求解器最终位姿变了没有？”——后者由 tests/baseline.txt 负责。
//
//  【为什么单独一份】
//  接下来会把 Box/Cylinder/Convex 统一到 PolytopeView + SAT 管线，接触点生成方式
//  会**有意**变化。主基线不经过 CollisionWorld，不应该被牵连；而碰撞层需要一份
//  能对比“变化了多少、变化是否合理”的观察值。所以这里只做确定性哈希，不在
//  tests/run.sh 里做强制的逐位门禁。
//
//  【确定性契约】
//   1) 只按固定字段顺序 hash，不对结构体整体 memcpy（避免 padding）；
//   2) 浮点按 float 字节喂入，整数按 64 位喂入；
//   3) 遍历顺序即 ABI：Manifold 按 store 顺序、contact 按数组顺序、约束按列表顺序；
//   4) 不 hash 指针/容量/未初始化字段。
// ============================================================================

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <type_traits>
#include <variant>

#include "collision/contact_constraints.hpp"
#include "collision/manifold.hpp"
#include "collision/narrowphase.hpp"
#include "collision/shape.hpp"
#include "phys_core/constraint.hpp"

namespace collision::collision_digest {

inline constexpr std::uint64_t kFnvOffset = 1469598103934665603ull;
inline constexpr std::uint64_t kFnvPrime  = 1099511628211ull;

inline std::uint64_t mix_bytes(std::uint64_t h, const unsigned char *b, std::size_t n) noexcept {
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<std::uint64_t>(b[i]);
    h *= kFnvPrime;
  }
  return h;
}

inline std::uint64_t mix_u64(std::uint64_t h, std::uint64_t v) noexcept {
  return mix_bytes(h, reinterpret_cast<const unsigned char *>(&v), sizeof(v));
}

inline std::uint64_t mix_scalar(std::uint64_t h, Scalar_t v) noexcept {
  return mix_bytes(h, reinterpret_cast<const unsigned char *>(&v), sizeof(v));
}

inline std::uint64_t mix_vec3(std::uint64_t h, Vec3 v) noexcept {
  h = mix_scalar(h, v.x);
  h = mix_scalar(h, v.y);
  return mix_scalar(h, v.z);
}

inline std::uint64_t mix_quat(std::uint64_t h, Quat q) noexcept {
  h = mix_scalar(h, q.x);
  h = mix_scalar(h, q.y);
  h = mix_scalar(h, q.z);
  return mix_scalar(h, q.w);
}

inline std::uint64_t mix_mat3(std::uint64_t h, const Mat3 &m) noexcept {
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r)
      h = mix_scalar(h, m[c][r]);
  return h;
}

inline std::uint64_t mix_basis(std::uint64_t h, const BasisFrame &b) noexcept {
  h = mix_vec3(h, b.normal);
  h = mix_vec3(h, b.tangent0);
  return mix_vec3(h, b.tangent1);
}

// ShapeView 本身也进摘要：这样“改了一个 shape 参数但接触点恰好没变”的情况
// 仍然可见。
inline std::uint64_t shape_view(std::uint64_t h, const ShapeView &v) noexcept {
  h = mix_u64(h, static_cast<std::uint64_t>(v.type));
  h = mix_vec3(h, v.center);
  h = mix_quat(h, v.rotation);
  h = mix_mat3(h, v.rot);
  h = mix_vec3(h, v.body_center);
  h = mix_quat(h, v.body_rotation);
  h = mix_vec3(h, v.local_half);
  h = mix_vec3(h, v.half);
  h = mix_scalar(h, v.radius);
  h = mix_scalar(h, v.half_height);
  h = mix_vec3(h, v.axis);
  return mix_scalar(h, v.friction);
}

inline std::uint64_t contact_point(std::uint64_t h, const ContactPoint &p) noexcept {
  h = mix_u64(h, static_cast<std::uint64_t>(p.feature_key));
  h = mix_vec3(h, p.rA);
  h = mix_vec3(h, p.rB);
  h = mix_vec3(h, p.xA);
  h = mix_vec3(h, p.xB);
  return mix_scalar(h, p.depth);
}

inline std::uint64_t contact_patch(std::uint64_t h, const ContactPatch &p) noexcept {
  h = mix_u64(h, static_cast<std::uint64_t>(p.body_a));
  h = mix_u64(h, static_cast<std::uint64_t>(p.body_b));
  h = mix_u64(h, static_cast<std::uint64_t>(p.proxy_a));
  h = mix_u64(h, static_cast<std::uint64_t>(p.proxy_b));
  h = mix_basis(h, p.basis);
  h = mix_scalar(h, p.friction);
  h = mix_u64(h, static_cast<std::uint64_t>(p.count));
  for (std::size_t i = 0; i < p.count; ++i)
    h = contact_point(h, p.points[i]);
  return h;
}

inline std::uint64_t manifold(uint64_t h, const Manifold &m) noexcept {
  h = mix_u64(h, static_cast<std::uint64_t>(m.body_a));
  h = mix_u64(h, static_cast<std::uint64_t>(m.body_b));
  h = mix_u64(h, static_cast<std::uint64_t>(m.proxy_a));
  h = mix_u64(h, static_cast<std::uint64_t>(m.proxy_b));
  h = mix_basis(h, m.basis);
  h = mix_scalar(h, m.friction);
  h = mix_u64(h, static_cast<std::uint64_t>(m.count));
  for (std::size_t i = 0; i < m.count; ++i) {
    const ManifoldContact &c = m.contacts[i];
    h                        = mix_u64(h, static_cast<std::uint64_t>(c.feature_key));
    h                        = mix_vec3(h, c.rA);
    h                        = mix_vec3(h, c.rB);
    h                        = mix_scalar(h, c.depth);
    h                        = mix_vec3(h, c.lambda);
    h                        = mix_vec3(h, c.penalty);
    h                        = mix_scalar(h, c.beta_scale);
    h                        = mix_u64(h, static_cast<std::uint64_t>(c.stick));
  }
  return h;
}

inline std::uint64_t manifolds(std::span<const Manifold> ms) noexcept {
  std::uint64_t h = kFnvOffset;
  h               = mix_u64(h, static_cast<std::uint64_t>(ms.size()));
  for (const Manifold &m : ms)
    h = manifold(h, m);
  return h;
}

inline std::uint64_t contact_constraint(std::uint64_t h, const ContactConstraint &c) noexcept {
  h = mix_u64(h, static_cast<std::uint64_t>(c.bodies_[0]));
  h = mix_u64(h, static_cast<std::uint64_t>(c.bodies_[1]));
  h = mix_u64(h, static_cast<std::uint64_t>(c.manifold_id));
  h = mix_u64(h, static_cast<std::uint64_t>(c.contact_index));
  h = mix_u64(h, static_cast<std::uint64_t>(c.kind));
  for (int i = 0; i < 3; ++i)
    h = mix_vec3(h, c.basis[i]);
  h = mix_vec3(h, c.rA);
  h = mix_vec3(h, c.rB);
  h = mix_scalar(h, c.target);
  h = mix_scalar(h, c.lambda);
  h = mix_scalar(h, c.compliance);
  h = mix_scalar(h, c.penalty);
  h = mix_scalar(h, c.beta_scale);
  return mix_u64(h, static_cast<std::uint64_t>(c.stick));
}

inline std::uint64_t constraints(const ConstraintList &list) noexcept {
  std::uint64_t h = kFnvOffset;
  h               = mix_u64(h, static_cast<std::uint64_t>(list.size()));
  for (const ConstraintVariant &v : list) {
    h = mix_u64(h, static_cast<std::uint64_t>(v.index()));
    std::visit(
        [&h](const auto &c) {
          using C = std::decay_t<decltype(c)>;
          if constexpr (std::is_same_v<C, ContactConstraint>) {
            h = contact_constraint(h, c);
          } else if constexpr (std::is_same_v<C, DistanceConstraint>) {
            h = mix_u64(h, static_cast<std::uint64_t>(c.bodies_[0]));
            h = mix_u64(h, static_cast<std::uint64_t>(c.bodies_[1]));
            h = mix_scalar(h, c.rest_length);
            h = mix_scalar(h, c.lambda);
            h = mix_scalar(h, c.compliance);
          }
        },
        v);
  }
  return h;
}

inline std::string to_hex(std::uint64_t v) {
  static const char *kDigits = "0123456789abcdef";
  std::string        out(16, '0');
  for (int i = 15; i >= 0; --i) {
    out[static_cast<std::size_t>(i)] = kDigits[v & 0xFull];
    v >>= 4;
  }
  return out;
}

} // namespace collision::collision_digest
