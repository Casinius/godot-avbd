#pragma once

// ============================================================================
//  collision/convex_hull.hpp —— 可选的点云 → ConvexMesh builder
// ----------------------------------------------------------------------------
//  【定位】
//  客户端通常已经能提供 vertices + face loops（Unity/PhysX/自研资产管线）。
//  这个 builder 只是资产导入便利路径，**不在 CollisionWorld::step() 热路径**，
//  也不被 collision world 自动调用。实现放在独立 .cpp，避免污染窄相。
//
//  【确定性】
//   输入先排序 + epsilon 去重；初始四面体、外侧点处理、horizon 边都按排序后的
//   索引顺序推进。不使用随机数、不使用哈希容器。
// ============================================================================

#include <cstddef>
#include <cstdint>
#include <span>

#include "collision/convex_mesh.hpp"

namespace collision {

enum class HullStatus : std::uint8_t {
  Ok = 0,
  TooFewPoints,
  Degenerate,
  NonFiniteInput,
  TooManyVertices,
  TooManyFaces,
};

struct HullOptions {
  Scalar_t    epsilon      = Scalar_t(1e-5);
  std::size_t max_vertices = kMaxConvexVertices;
  std::size_t max_faces    = kMaxConvexFaces;
};

// 成功后 out 已经过 finalize_convex_mesh() 校验；失败时 out 内容未定义。
[[nodiscard]] HullStatus build_convex_hull(std::span<const Vec3> points, const HullOptions &opt,
                                           ConvexMesh &out);

[[nodiscard]] inline const char *to_string(HullStatus s) noexcept {
  switch (s) {
  case HullStatus::Ok:
    return "Ok";
  case HullStatus::TooFewPoints:
    return "TooFewPoints";
  case HullStatus::Degenerate:
    return "Degenerate";
  case HullStatus::NonFiniteInput:
    return "NonFiniteInput";
  case HullStatus::TooManyVertices:
    return "TooManyVertices";
  case HullStatus::TooManyFaces:
    return "TooManyFaces";
  }
  return "Unknown";
}

} // namespace collision
