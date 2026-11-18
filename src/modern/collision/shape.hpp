#pragma once

// ============================================================================
//  collision/shape.hpp —— 碰撞形状的局部定义 + 世界视图
// ----------------------------------------------------------------------------
//  【不变量】
//   1) Shape 只描述**局部**几何与材质，不含位置/旋转。世界位姿由 ShapeView 携带。
//   2) 圆柱的轴固定为局部 +Y，half_height 沿 +Y，半径在局部 XZ 平面内。
//      参考实现里 makeShape 用 +Z、closestPointOnCylinder 用 +Y、pointInCylinder
//      又用 +Z，三处不一致；这里只允许 +Y，避免把那个 bug 带进来。
//   3) ShapeView 是每步由 State 的 pose 派生的只读快照，不拥有任何东西。
//   4) 本文件不 include 求解器/线程/容器拥有关系，只依赖 state.hpp 的数学类型。
// ============================================================================

#include <cmath>
#include <cstdint>

#include "collision/aabb.hpp"
#include "collision/convex_mesh.hpp"
#include "phys_core/state.hpp"

namespace collision {

// convex 数据由 ConvexMeshStore 拥有，ShapeView 只借看。
enum class ShapeType : std::uint8_t { Box = 0, Sphere = 1, Cylinder = 2, Convex = 3 };

struct Material {
  Scalar_t friction    = Scalar_t(0.5);
  Scalar_t restitution = Scalar_t(0);
};

struct Shape {
  ShapeType type = ShapeType::Box;
  Vec3      half{0.5f, 0.5f, 0.5f}; // Box: 半尺寸；其它类型只读 x 作为兼容值
  Scalar_t  radius      = 0.5f;     // Sphere / Cylinder
  Scalar_t  half_height = 0.5f;     // Cylinder，沿局部 +Y
  ConvexMeshId mesh     = kInvalidConvexMesh; // Convex：注册表 id
  Material  material{};

  [[nodiscard]] static Shape box(Vec3 half_extents, Material m = {}) noexcept {
    Shape s;
    s.type     = ShapeType::Box;
    s.half     = half_extents;
    s.material = m;
    return s;
  }
  [[nodiscard]] static Shape sphere(Scalar_t r, Material m = {}) noexcept {
    Shape s;
    s.type     = ShapeType::Sphere;
    s.radius   = r;
    s.half     = Vec3{r, r, r};
    s.material = m;
    return s;
  }
  [[nodiscard]] static Shape convex(ConvexMeshId id, Material m = {}) noexcept {
    Shape s;
    s.type     = ShapeType::Convex;
    s.mesh     = id;
    s.material = m;
    return s;
  }
  [[nodiscard]] static Shape cylinder(Scalar_t r, Scalar_t height, Material m = {}) noexcept {
    Shape s;
    s.type        = ShapeType::Cylinder;
    s.radius      = r;
    s.half_height = height * Scalar_t(0.5);
    s.half        = Vec3{r, s.half_height, r};
    s.material    = m;
    return s;
  }

  // 包围球半径：用于快速拒绝与 broadphase 的保守 AABB。
  [[nodiscard]] Scalar_t bounding_radius() const noexcept {
    switch (type) {
    case ShapeType::Sphere:
      return radius;
    case ShapeType::Cylinder:
      return std::sqrt(radius * radius + half_height * half_height);
    case ShapeType::Convex:
      return glm::length(half); // 保守 fallback；权威值来自 ConvexMesh.bounding_radius
    case ShapeType::Box:
    default:
      return std::sqrt(glm::dot(half, half));
    }
  }

  // 局部 AABB 半尺寸。圆柱按局部 +Y 轴。
  [[nodiscard]] Vec3 local_half_extents() const noexcept {
    switch (type) {
    case ShapeType::Sphere:
      return Vec3{radius, radius, radius};
    case ShapeType::Cylinder:
      return Vec3{radius, half_height, radius};
    case ShapeType::Convex:
      return half; // 保守 fallback；CollisionWorld 会从 mesh 填 ShapeView.local_half
    case ShapeType::Box:
    default:
      return half;
    }
  }
};

// 世界空间只读视图。center/rotation 是 **shape** 的世界位姿；body_center/body_rotation
// 是所属刚体的位姿。compound 下两者不同：几何计算用 shape 位姿，接触锚点 rA/rB
// 必须用 body 位姿，否则 solver 会把 shape-local 当成 body-local。
struct ShapeView {
  ShapeType type = ShapeType::Box;
  Vec3      center{0};
  Quat      rotation = identity_quat();
  Mat3      rot      = Mat3(1);
  Vec3      half{0.5f};
  Scalar_t  radius      = 0.5f;
  Scalar_t  half_height = 0.5f;
  Vec3      axis{0, 1, 0}; // Cylinder 局部 +Y 的世界方向
  Scalar_t  friction = 0.5f;

  // —— body 位姿：只用于 world → body-local 锚点，不参与 shape 几何 ——
  Vec3 body_center{0};
  Quat body_rotation = identity_quat();

  // 解析后的局部 AABB 半尺寸。Box/Sphere/Cylinder 来自 Shape；Convex 来自 mesh。
  Vec3 local_half{0.5f};

  // 非拥有指针。生命周期由 CollisionWorld 的 ConvexMeshStore 保证。
  const ConvexMesh *convex = nullptr;

  [[nodiscard]] Vec3 local_to_world(Vec3 local) const noexcept { return center + rot * local; }
  // shape-local：几何算法（point_in_box / point_in_cylinder / convex support）用。
  [[nodiscard]] Vec3 world_to_local(Vec3 world) const noexcept {
    return glm::conjugate(rotation) * (world - center);
  }
  // body-local：只有接触锚点 rA/rB 用。compound 下这是正确性关键。
  [[nodiscard]] Vec3 world_to_body_local(Vec3 world) const noexcept {
    return glm::conjugate(body_rotation) * (world - body_center);
  }
};

[[nodiscard]] inline bool quat_is_identity(const Quat &q) noexcept {
  return q.x == 0.0f && q.y == 0.0f && q.z == 0.0f && q.w == 1.0f;
}

// 完整构造：body pose + instance local transform → shape pose。
// 对 identity local transform 做短路，保证旧单 shape 路径逐位不变。
[[nodiscard]] inline ShapeView make_shape_view(const Shape &shape, Vec3 body_center,
                                               Quat body_rotation, Vec3 local_position,
                                               Quat local_rotation,
                                               const ConvexMesh *convex = nullptr) noexcept {
  ShapeView v;
  v.type          = shape.type;
  v.body_center   = body_center;
  v.body_rotation = body_rotation;

  const bool identity_local = local_position.x == 0.0f && local_position.y == 0.0f &&
                              local_position.z == 0.0f && quat_is_identity(local_rotation);
  if (identity_local) {
    v.center   = body_center;
    v.rotation = body_rotation;
  } else {
    v.center   = body_center + rotate(body_rotation, local_position);
    v.rotation = glm::normalize(body_rotation * local_rotation);
  }

  v.rot         = rot_matrix(v.rotation);
  v.half        = shape.half;
  v.radius      = shape.radius;
  v.half_height = shape.half_height;
  v.axis        = v.rot * Vec3{0, 1, 0};
  v.friction    = shape.material.friction;
  v.local_half  = (shape.type == ShapeType::Convex && convex != nullptr)
                      ? convex->local_half
                      : shape.local_half_extents();
  v.convex      = convex;
  return v;
}

// 兼容旧调用：body 位姿 == shape 位姿，local transform 为零。
[[nodiscard]] inline ShapeView make_shape_view(const Shape &shape, Vec3 center,
                                               Quat rotation) noexcept {
  return make_shape_view(shape, center, rotation, Vec3{0}, identity_quat(), nullptr);
}

// 紧致世界 AABB：|R| * local_half 是 OBB 在三个世界轴上的投影半径。
[[nodiscard]] inline Aabb shape_aabb(const ShapeView &v) noexcept {
  const Vec3  h = v.local_half;
  const Mat3 &r = v.rot;
  const Vec3  wx{std::abs(r[0][0]) * h.x + std::abs(r[1][0]) * h.y + std::abs(r[2][0]) * h.z,
                std::abs(r[0][1]) * h.x + std::abs(r[1][1]) * h.y + std::abs(r[2][1]) * h.z,
                std::abs(r[0][2]) * h.x + std::abs(r[1][2]) * h.y + std::abs(r[2][2]) * h.z};
  return Aabb{v.center - wx, v.center + wx};
}

} // namespace collision
