#pragma once

// ============================================================================
//  collision/mass_properties.hpp —— compound 的质量/质心/惯量 helper
// ----------------------------------------------------------------------------
//  solver 约定 body pose 的位置是 **COM**。客户端从引擎导入 compound 时，body
//  原点经常不是 COM；如果直接把各 part 的局部位姿塞进 ShapeInstance，旋转动力学
//  会错。这个 helper 做两件事：
//
//    1) 由 ShapeInstance 列表 + density 汇总 mass/com/inertia（相对 COM）；
//    2) recenter_instances() 把 local_position 平移到新 COM，方便客户端直接写
//       State::pose.set_position(body, com)。
//
//  只在 setup-time 调用，允许分配；不进入 step() 热路径。
// ============================================================================

#include <cmath>
#include <cstddef>
#include <span>

#include "collision/convex_mesh.hpp"
#include "collision/shape_instance.hpp"

namespace collision {

struct MassProperties {
  Scalar_t mass = 0;
  Vec3     com{0};
  Mat3     inertia{0}; // 相对 com
};

enum class MassStatus : std::uint8_t {
  Ok = 0,
  Empty,
  InvalidShape,
  Degenerate,
};

[[nodiscard]] inline Mat3 diag3(Vec3 d) noexcept {
  Mat3 m(0);
  m[0][0] = d.x;
  m[1][1] = d.y;
  m[2][2] = d.z;
  return m;
}

[[nodiscard]] inline Mat3 inertia_about_origin_from_com(Mat3 inertia_com, Scalar_t mass,
                                                        Vec3 com) noexcept {
  return inertia_com + (Mat3(1) * glm::dot(com, com) - outer(com, com)) * mass;
}

[[nodiscard]] inline Mat3 inertia_about_com_from_origin(Mat3 inertia_origin, Scalar_t mass,
                                                        Vec3 com) noexcept {
  return inertia_origin - (Mat3(1) * glm::dot(com, com) - outer(com, com)) * mass;
}

// 单个 shape 的质量属性。convex 必须传 finalized mesh；density 是质量密度。
[[nodiscard]] inline MassStatus shape_mass_properties(const Shape &shape, const ConvexMesh *mesh,
                                                      Scalar_t        density,
                                                      MassProperties &out) noexcept {
  if (!(density > 0))
    return MassStatus::Degenerate;

  switch (shape.type) {
  case ShapeType::Box: {
    const Scalar_t hx = shape.half.x, hy = shape.half.y, hz = shape.half.z;
    if (!(hx > 0 && hy > 0 && hz > 0))
      return MassStatus::Degenerate;
    const Scalar_t m = density * Scalar_t(8) * hx * hy * hz;
    out.mass         = m;
    out.com          = Vec3{0};
    out.inertia =
        diag3(Vec3{m * (hy * hy + hz * hz) / Scalar_t(3), m * (hx * hx + hz * hz) / Scalar_t(3),
                   m * (hx * hx + hy * hy) / Scalar_t(3)});
    return MassStatus::Ok;
  }
  case ShapeType::Sphere: {
    const Scalar_t r = shape.radius;
    if (!(r > 0))
      return MassStatus::Degenerate;
    const Scalar_t m = density * Scalar_t(4.188790204786391) * r * r * r; // 4/3 π r³
    const Scalar_t i = Scalar_t(0.4) * m * r * r;
    out.mass         = m;
    out.com          = Vec3{0};
    out.inertia      = diag3(Vec3{i, i, i});
    return MassStatus::Ok;
  }
  case ShapeType::Cylinder: {
    const Scalar_t r = shape.radius;
    const Scalar_t h = shape.half_height * Scalar_t(2);
    if (!(r > 0 && h > 0))
      return MassStatus::Degenerate;
    const Scalar_t m  = density * Scalar_t(3.141592653589793) * r * r * h;
    const Scalar_t ix = m * (Scalar_t(3) * r * r + h * h) / Scalar_t(12);
    const Scalar_t iy = m * r * r * Scalar_t(0.5);
    out.mass          = m;
    out.com           = Vec3{0};
    out.inertia       = diag3(Vec3{ix, iy, ix});
    return MassStatus::Ok;
  }
  case ShapeType::Convex: {
    if (mesh == nullptr || mesh->vertices.empty() || mesh->faces.empty())
      return MassStatus::InvalidShape;

    Scalar_t volume = 0;
    Vec3     first_moment{0};
    Mat3     covariance(0);
    for (const ConvexFace &face : mesh->faces) {
      const std::uint32_t i0 = mesh->face_indices[face.first];
      const Vec3          a  = mesh->vertices[i0];
      for (std::uint32_t j = 1; j + 1 < face.count; ++j) {
        const Vec3     b   = mesh->vertices[mesh->face_indices[face.first + j]];
        const Vec3     c   = mesh->vertices[mesh->face_indices[face.first + j + 1]];
        const Scalar_t det = glm::dot(a, glm::cross(b, c));
        volume += det / Scalar_t(6); // 四面体体积
        first_moment += (a + b + c) * (det / Scalar_t(24));
        covariance += (outer(a, a) + outer(b, b) + outer(c, c) + outer(a + b + c, a + b + c)) *
                      (det / Scalar_t(120));
      }
    }
    if (!(volume > Scalar_t(1e-9)))
      return MassStatus::Degenerate;

    const Vec3 com = first_moment / volume;
    const Mat3 inertia_origin =
        diag3(Vec3{covariance[0][0] + covariance[1][1] + covariance[2][2],
                   covariance[0][0] + covariance[1][1] + covariance[2][2],
                   covariance[0][0] + covariance[1][1] + covariance[2][2]}) -
        covariance;
    const Scalar_t mass = density * volume;

    out.mass    = mass;
    out.com     = com;
    out.inertia = inertia_about_com_from_origin(inertia_origin, mass, com);
    return MassStatus::Ok;
  }
  }
  return MassStatus::InvalidShape;
}

// compound 汇总：所有 local transform 相对 body 原点；结果 com/inertia 相对 COM。
[[nodiscard]] inline MassStatus body_mass_properties(std::span<const ShapeInstance> instances,
                                                     const ConvexMeshStore         &meshes,
                                                     Scalar_t                       density,
                                                     MassProperties                &out) noexcept {
  if (instances.empty())
    return MassStatus::Empty;

  Scalar_t total_mass = 0;
  Vec3     first_moment{0};
  Mat3     inertia_origin(0);

  for (const ShapeInstance &inst : instances) {
    const ConvexMesh *mesh =
        inst.shape.type == ShapeType::Convex ? meshes.get(inst.shape.mesh) : nullptr;
    MassProperties   mp;
    const MassStatus st = shape_mass_properties(inst.shape, mesh, density, mp);
    if (st != MassStatus::Ok)
      return st;

    const Mat3 R        = rot_matrix(inst.local_rotation);
    const Vec3 com_body = inst.local_position + rotate(inst.local_rotation, mp.com);
    const Mat3 I_com    = R * mp.inertia * transpose(R);
    inertia_origin += inertia_about_origin_from_com(I_com, mp.mass, com_body);
    first_moment += com_body * mp.mass;
    total_mass += mp.mass;
  }

  if (!(total_mass > 0))
    return MassStatus::Degenerate;

  const Vec3 com = first_moment / total_mass;
  out.mass       = total_mass;
  out.com        = com;
  out.inertia    = inertia_about_com_from_origin(inertia_origin, total_mass, com);
  return MassStatus::Ok;
}

// 把实例表平移到新 COM：调用后 body pose 应设为 com。
inline void recenter_instances(std::span<ShapeInstance> instances, Vec3 com) noexcept {
  for (ShapeInstance &inst : instances)
    inst.local_position -= com;
}

} // namespace collision
