#pragma once

// ============================================================================
//  collision_selftest.hpp —— 碰撞层的结构自检
// ----------------------------------------------------------------------------
//  这些检查不进入基线，只在 --selftest 里跑。它们回答“碰撞层自己是否自洽”：
//  宽相与暴力枚举一致、窄相法线符号正确、manifold warm start 不丢、约束布局连续、
//  以及 CollisionWorld 的拓扑签名/版本号行为。端到端解算只断言有限性与静态冻结。
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <format>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "collision/collision_digest.hpp"
#include "collision/collision_world.hpp"
#include "collision/convex_hull.hpp"
#include "collision/mass_properties.hpp"
#include "collision/contact_constraints.hpp"
#include "collision/manifold.hpp"
#include "collision/narrowphase.hpp"
#include "collision/shape.hpp"
#include "phys_core/policy.hpp"
#include "phys_core/sched.hpp"
#include "phys_core/solver_phase.hpp"
#include "task.hpp"

namespace collision_selftest {

using namespace collision;

// 单 shape 的测试糖。正式 API 已经删掉 set_shape：一个 body 必须显式给
// ShapeInstance 列表。这里只是让旧用例保持一行，不代表 collision 层还有单 shape 特例。
inline void set_single_shape(CollisionWorld &world, std::size_t body, Shape shape) {
  const std::array<ShapeInstance, 1> one{ShapeInstance{shape}};
  world.set_shapes(body, one);
}

// 本头文件不在 main 的匿名命名空间里（否则 namespace collision 会被嵌进匿名
// 命名空间，和 .cpp 里的定义对不上），所以自带一个同形结果类型，由 main 转换。
struct CollisionCheck {
  std::string name;
  bool        ok = true;
  std::string detail;
};

inline ConvexMesh make_cube_convex_mesh() {
  ConvexMesh m;
  m.vertices = {
      Vec3{-0.5f, -0.5f, -0.5f}, Vec3{0.5f, -0.5f, -0.5f}, Vec3{0.5f, 0.5f, -0.5f},
      Vec3{-0.5f, 0.5f, -0.5f},  Vec3{-0.5f, -0.5f, 0.5f}, Vec3{0.5f, -0.5f, 0.5f},
      Vec3{0.5f, 0.5f, 0.5f},    Vec3{-0.5f, 0.5f, 0.5f},
  };
  m.face_indices = {
      0, 3, 2, 1, // -z
      4, 5, 6, 7, // +z
      0, 1, 5, 4, // -y
      3, 7, 6, 2, // +y
      0, 4, 7, 3, // -x
      1, 2, 6, 5, // +x
  };
  m.faces = {
      ConvexFace{0, 4}, ConvexFace{4, 4}, ConvexFace{8, 4},
      ConvexFace{12, 4}, ConvexFace{16, 4}, ConvexFace{20, 4},
  };
  return m;
}

inline CollisionCheck check_convex_mesh_and_hull() {
  CollisionCheck c{"collision convex mesh validation + deterministic hull", true, ""};

  ConvexMesh cube = make_cube_convex_mesh();
  const ConvexMeshStatus st = finalize_convex_mesh(cube);
  if (st != ConvexMeshStatus::Ok) {
    c.ok     = false;
    c.detail = std::format("finalize_convex_mesh: {}", to_string(st));
    return c;
  }
  if (cube.vertices.size() != 8 || cube.faces.size() != 6 || cube.edges.size() != 12) {
    c.ok     = false;
    c.detail = std::format("cube: vertices={} faces={} edges={}", cube.vertices.size(),
                           cube.faces.size(), cube.edges.size());
    return c;
  }
  if (std::abs(cube.local_half.x - 0.5f) > 1e-5f || cube.bounding_radius < 0.86f) {
    c.ok     = false;
    c.detail = std::format("cube extents: half=({},{},{}) radius={}", cube.local_half.x,
                           cube.local_half.y, cube.local_half.z, cube.bounding_radius);
    return c;
  }

  ConvexMeshStore store;
  const ConvexMeshId id = store.add(cube);
  if (id == kInvalidConvexMesh || store.get(id) == nullptr ||
      store.get(id)->faces.size() != 6) {
    c.ok     = false;
    c.detail = "ConvexMeshStore add/get failed";
    return c;
  }

  // 点云 helper：8 个立方体角点 → 至少 6 个面、8 个顶点，且必须通过 finalize。
  const std::array<Vec3, 8> pts{
      Vec3{-0.5f, -0.5f, -0.5f}, Vec3{0.5f, -0.5f, -0.5f}, Vec3{0.5f, 0.5f, -0.5f},
      Vec3{-0.5f, 0.5f, -0.5f},  Vec3{-0.5f, -0.5f, 0.5f}, Vec3{0.5f, -0.5f, 0.5f},
      Vec3{0.5f, 0.5f, 0.5f},    Vec3{-0.5f, 0.5f, 0.5f},
  };
  ConvexMesh hull;
  const HullStatus hs = build_convex_hull(pts, HullOptions{}, hull);
  if (hs != HullStatus::Ok) {
    c.ok     = false;
    c.detail = std::format("build_convex_hull: {}", to_string(hs));
    return c;
  }
  if (hull.vertices.size() != 8 || hull.faces.size() < 6 || hull.edges.size() < 12) {
    c.ok     = false;
    c.detail = std::format("hull: vertices={} faces={} edges={}", hull.vertices.size(),
                           hull.faces.size(), hull.edges.size());
    return c;
  }

  const std::array<Vec3, 3> degenerate{
      Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}};
  ConvexMesh bad;
  if (build_convex_hull(degenerate, HullOptions{}, bad) != HullStatus::TooFewPoints) {
    c.ok     = false;
    c.detail = "degenerate point set did not return TooFewPoints";
    return c;
  }

  const Shape convex_shape = Shape::convex(id);
  if (convex_shape.type != ShapeType::Convex || convex_shape.mesh != id) {
    c.ok     = false;
    c.detail = "Shape::convex did not store the mesh id";
  }
  return c;
}

inline CollisionCheck check_shape_table_layout() {
  CollisionCheck c{"collision shape table body ranges + revision", true, ""};

  ShapeTable table;
  table.resize_bodies(3);
  if (table.body_count() != 3 || table.instance_count() != 0) {
    c.ok     = false;
    c.detail = std::format("resize_bodies(3): bodies={} instances={}", table.body_count(),
                           table.instance_count());
    return c;
  }

  const std::uint64_t r0 = table.revision();
  const std::array<ShapeInstance, 2> two{
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f})},
      ShapeInstance{Shape::sphere(0.5f), Vec3{1, 0, 0}, Quat{2, 0, 0, 0}, 7},
  };
  table.set_shapes(1, two);
  const std::span<const ShapeInstance> b1 = table.shapes(1);
  if (b1.size() != 2 || table.instance_count() != 2 || table.revision() == r0) {
    c.ok     = false;
    c.detail = std::format("set_shapes: size={} instances={} revision={}", b1.size(),
                           table.instance_count(), table.revision());
    return c;
  }
  if (b1[1].user_id != 7 || std::abs(glm::dot(b1[1].local_rotation, b1[1].local_rotation) - 1) >
                                 Scalar_t(1e-5)) {
    c.ok     = false;
    c.detail = "ShapeInstance local rotation was not normalized / user_id lost";
    return c;
  }

  const std::array<ShapeInstance, 1> one{ShapeInstance{Shape::sphere(0.25f)}};
  table.add_shapes(0, one);
  if (table.shapes(0).size() != 1 || table.shapes(1).size() != 2 ||
      table.shapes(2).size() != 0 || table.instance_count() != 3) {
    c.ok     = false;
    c.detail = std::format("add_shapes order: b0={} b1={} b2={} total={}", table.shapes(0).size(),
                           table.shapes(1).size(), table.shapes(2).size(), table.instance_count());
    return c;
  }

  table.clear_shapes(1);
  if (table.shapes(1).size() != 0 || table.instance_count() != 1) {
    c.ok     = false;
    c.detail = "clear_shapes did not remove the body range";
    return c;
  }

  table.resize_bodies(5);
  if (table.body_count() != 5 || table.shapes(0).size() != 1 || table.shapes(4).size() != 0) {
    c.ok     = false;
    c.detail = "resize_bodies did not preserve existing ranges";
    return c;
  }
  table.resize_bodies(2);
  if (table.body_count() != 2 || table.shapes(0).size() != 1) {
    c.ok     = false;
    c.detail = "resize_bodies(shrink) lost existing data";
  }
  return c;
}

inline CollisionCheck check_convex_polytope_collision() {
  CollisionCheck c{"collision convex/polytope narrowphase", true, ""};

  ConvexMesh cube = make_cube_convex_mesh();
  if (finalize_convex_mesh(cube) != ConvexMeshStatus::Ok) {
    c.ok     = false;
    c.detail = "cube finalize failed";
    return c;
  }
  ConvexMeshStore store;
  const ConvexMeshId id = store.add(cube);
  const ConvexMesh  *mesh = store.get(id);
  if (mesh == nullptr) {
    c.ok     = false;
    c.detail = "cube store add failed";
    return c;
  }

  const ShapeView convex_a =
      make_shape_view(Shape::convex(id), Vec3{0, 0, 0}, identity_quat(), Vec3{0},
                      identity_quat(), mesh);
  const ShapeView convex_b =
      make_shape_view(Shape::convex(id), Vec3{0, 0.9f, 0}, identity_quat(), Vec3{0},
                      identity_quat(), mesh);
  const ShapeView box_b = make_shape_view(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0, 0.9f, 0},
                                          identity_quat());
  const ShapeView sphere_b = make_shape_view(Shape::sphere(0.5f), Vec3{0, 0.8f, 0},
                                             identity_quat());

  auto collide = [](const ShapeView &a, const ShapeView &b, BasisFrame &basis) {
    std::array<ContactPoint, kMaxContactsPerManifold> points{};
    const std::size_t n = collide_shapes(a, b, points, basis);
    return std::pair{n, points};
  };

  BasisFrame basis{};
  auto [n_box, pts_box] = collide(convex_a, box_b, basis);
  if (n_box == 0 || !(basis.normal.y < -0.9f)) {
    c.ok     = false;
    c.detail = std::format("convex/box: n={} normal.y={}", n_box, basis.normal.y);
    return c;
  }

  auto [n_convex, pts_convex] = collide(convex_a, convex_b, basis);
  if (n_convex == 0 || !(basis.normal.y < -0.9f)) {
    c.ok     = false;
    c.detail = std::format("convex/convex: n={} normal.y={}", n_convex, basis.normal.y);
    return c;
  }

  auto [n_sphere, pts_sphere] = collide(convex_a, sphere_b, basis);
  if (n_sphere == 0 || !(basis.normal.y < -0.9f)) {
    c.ok     = false;
    c.detail = std::format("convex/sphere: n={} normal.y={}", n_sphere, basis.normal.y);
    return c;
  }
  if (sphere_b.center.y - convex_a.center.y > 1.0f) {
    c.ok     = false;
    c.detail = "test setup invalid";
  }
  return c;
}

inline CollisionCheck check_mass_properties() {
  CollisionCheck c{"collision compound mass properties", true, ""};

  ConvexMesh cube = make_cube_convex_mesh();
  if (finalize_convex_mesh(cube) != ConvexMeshStatus::Ok) {
    c.ok     = false;
    c.detail = "cube finalize failed";
    return c;
  }
  ConvexMeshStore store;
  const ConvexMeshId id = store.add(cube);
  if (id == kInvalidConvexMesh) {
    c.ok     = false;
    c.detail = "cube store add failed";
    return c;
  }

  MassProperties box_mp;
  if (shape_mass_properties(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), nullptr, Scalar_t(1),
                            box_mp) != MassStatus::Ok) {
    c.ok     = false;
    c.detail = "box mass properties failed";
    return c;
  }
  MassProperties convex_mp;
  if (shape_mass_properties(Shape::convex(id), store.get(id), Scalar_t(1), convex_mp) !=
      MassStatus::Ok) {
    c.ok     = false;
    c.detail = "convex mass properties failed";
    return c;
  }
  const Scalar_t want_i = Scalar_t(1) / Scalar_t(6);
  if (std::abs(box_mp.mass - 1.0f) > 1e-5f || std::abs(convex_mp.mass - 1.0f) > 1e-5f ||
      std::abs(box_mp.inertia[0][0] - want_i) > 1e-5f ||
      std::abs(convex_mp.inertia[0][0] - want_i) > 1e-4f ||
      glm::length(convex_mp.com) > 1e-5f) {
    c.ok     = false;
    c.detail = std::format("cube mismatch: box m={} Ixx={} convex m={} Ixx={} com=({},{},{})",
                           box_mp.mass, box_mp.inertia[0][0], convex_mp.mass,
                           convex_mp.inertia[0][0], convex_mp.com.x, convex_mp.com.y,
                           convex_mp.com.z);
    return c;
  }

  std::array<ShapeInstance, 2> parts{
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{-0.6f, 0, 0}},
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0.6f, 0, 0}},
  };
  MassProperties compound;
  if (body_mass_properties(parts, store, Scalar_t(1), compound) != MassStatus::Ok) {
    c.ok     = false;
    c.detail = "compound mass properties failed";
    return c;
  }
  const Scalar_t want_yy = Scalar_t(2) * (want_i + Scalar_t(0.36));
  if (std::abs(compound.mass - 2.0f) > 1e-5f || glm::length(compound.com) > 1e-5f ||
      std::abs(compound.inertia[0][0] - Scalar_t(2) * want_i) > 1e-4f ||
      std::abs(compound.inertia[1][1] - want_yy) > 1e-4f ||
      std::abs(compound.inertia[2][2] - want_yy) > 1e-4f) {
    c.ok     = false;
    c.detail = std::format("compound mismatch: m={} I=({},{},{}) com=({},{},{})", compound.mass,
                           compound.inertia[0][0], compound.inertia[1][1], compound.inertia[2][2],
                           compound.com.x, compound.com.y, compound.com.z);
    return c;
  }

  recenter_instances(parts, Vec3{0.2f, 0, 0});
  if (std::abs(parts[0].local_position.x - (-0.8f)) > 1e-6f ||
      std::abs(parts[1].local_position.x - 0.4f) > 1e-6f) {
    c.ok     = false;
    c.detail = "recenter_instances did not shift local positions";
  }
  return c;
}

inline CollisionCheck check_broadphase_matches_bruteforce() {
  CollisionCheck c{"collision broadphase == brute force pair set", true, ""};

  // 9 个盒：一个 3x3 网格，相邻重叠、对角线分开；故意让几个盒完全同 AABB。
  std::vector<Aabb> boxes;
  for (int y = 0; y < 3; ++y)
    for (int x = 0; x < 3; ++x) {
      const Scalar_t px = Scalar_t(x) * Scalar_t(0.9);
      const Scalar_t py = Scalar_t(y) * Scalar_t(0.9);
      boxes.push_back(Aabb{Vec3{px - 0.5f, py - 0.5f, -0.5f}, Vec3{px + 0.5f, py + 0.5f, 0.5f}});
    }
  boxes.push_back(boxes[4]); // 完全重合的盒，pair 必须去重

  std::vector<BodyPair> want;
  for (std::size_t i = 0; i < boxes.size(); ++i)
    for (std::size_t j = i + 1; j < boxes.size(); ++j)
      if (overlaps(boxes[i], boxes[j]))
        want.push_back(BodyPair{i, j});
  std::sort(want.begin(), want.end());

  Broadphase bp;
  bp.update(boxes);
  const std::span<const BodyPair> got = bp.pairs();
  if (got.size() != want.size()) {
    c.ok     = false;
    c.detail = std::format("pair count {} != {}", got.size(), want.size());
    return c;
  }
  for (std::size_t i = 0; i < got.size(); ++i) {
    if (!(got[i] == want[i])) {
      c.ok     = false;
      c.detail = std::format("pair {} = ({},{}) want ({},{})", i, got[i].a, got[i].b, want[i].a,
                             want[i].b);
      return c;
    }
    if (got[i].a >= got[i].b) {
      c.ok     = false;
      c.detail = std::format("pair {} not canonical", i);
      return c;
    }
  }
  return c;
}

inline CollisionCheck check_narrowphase_sphere_sphere() {
  CollisionCheck c{"collision narrowphase sphere/sphere normal and depth", true, ""};

  const ShapeView a = make_shape_view(Shape::sphere(0.5f), Vec3{0, 0, 0}, identity_quat());
  const ShapeView b = make_shape_view(Shape::sphere(0.5f), Vec3{0.8f, 0, 0}, identity_quat());
  std::array<ContactPoint, kMaxContactsPerManifold> points{};
  BasisFrame                                        basis;
  const std::size_t                                 n = collide_shapes(a, b, points, basis);
  if (n != 1) {
    c.ok     = false;
    c.detail = std::format("contact count {}", n);
    return c;
  }
  const Vec3 want_n{-1, 0, 0}; // B->A
  if (glm::dot(basis.normal, want_n) < 0.99f) {
    c.ok     = false;
    c.detail = std::format("normal ({},{},{})", basis.normal.x, basis.normal.y, basis.normal.z);
    return c;
  }
  const Scalar_t C = glm::dot(basis.normal, points[0].xA - points[0].xB);
  if (!(C < 0) || std::abs(points[0].depth - 0.2f) > 1e-4f) {
    c.ok     = false;
    c.detail = std::format("C={} depth={}", C, points[0].depth);
  }
  return c;
}

inline CollisionCheck check_narrowphase_box_stack_normal() {
  CollisionCheck c{"collision narrowphase box/box A below B -> normal B->A", true, ""};

  const ShapeView a =
      make_shape_view(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0, 0, 0}, identity_quat());
  const ShapeView b =
      make_shape_view(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0, 0.9f, 0}, identity_quat());
  std::array<ContactPoint, kMaxContactsPerManifold> points{};
  BasisFrame                                        basis;
  const std::size_t                                 n = collide_shapes(a, b, points, basis);
  if (n == 0) {
    c.ok     = false;
    c.detail = "no contact for a clearly overlapping stack";
    return c;
  }
  if (!(basis.normal.y < -0.9f)) {
    c.ok     = false;
    c.detail = std::format("normal.y={} (want B->A = -y)", basis.normal.y);
    return c;
  }
  for (std::size_t i = 0; i < n; ++i) {
    const Scalar_t C = glm::dot(basis.normal, points[i].xA - points[i].xB);
    if (!(C < 0)) {
      c.ok     = false;
      c.detail = std::format("contact {} C={} (want <0)", i, C);
      return c;
    }
  }
  return c;
}

inline CollisionCheck check_narrowphase_sphere_box() {
  CollisionCheck c{"collision narrowphase box/sphere normal and depth", true, ""};

  const ShapeView a =
      make_shape_view(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0, 0, 0}, identity_quat());
  const ShapeView b = make_shape_view(Shape::sphere(0.5f), Vec3{0, 0.8f, 0}, identity_quat());
  std::array<ContactPoint, kMaxContactsPerManifold> points{};
  BasisFrame                                        basis;
  const std::size_t                                 n = collide_shapes(a, b, points, basis);
  if (n == 0) {
    c.ok     = false;
    c.detail = "no contact for overlapping sphere/box";
    return c;
  }
  if (!(basis.normal.y < -0.9f)) {
    c.ok     = false;
    c.detail = std::format("normal.y={} (want B->A = -y)", basis.normal.y);
    return c;
  }
  if (!(glm::dot(basis.normal, points[0].xA - points[0].xB) < 0)) {
    c.ok     = false;
    c.detail = "sphere/box C is not negative";
  }
  return c;
}

inline CollisionCheck check_narrowphase_cylinder_box() {
  CollisionCheck c{"collision narrowphase cylinder/box normal and depth", true, ""};

  const ShapeView a = make_shape_view(Shape::cylinder(0.5f, 1.0f), Vec3{0, 0, 0}, identity_quat());
  const ShapeView b =
      make_shape_view(Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0, 0.9f, 0}, identity_quat());
  std::array<ContactPoint, kMaxContactsPerManifold> points{};
  BasisFrame                                        basis;
  const std::size_t                                 n = collide_shapes(a, b, points, basis);
  if (n == 0) {
    c.ok     = false;
    c.detail = "no contact for overlapping cylinder/box";
    return c;
  }
  // 采样型 manifold 允许有 C==0 的边缘点（它们不推、也不该被双端更新）。要求至少
  // 一个点沿全局 basis 是穿透的即可。
  bool any_penetrating = false;
  for (std::size_t i = 0; i < n; ++i)
    any_penetrating = any_penetrating || (glm::dot(basis.normal, points[i].xA - points[i].xB) < 0);
  if (!any_penetrating) {
    c.ok     = false;
    c.detail = "cylinder/box manifold has no penetrating point";
  }
  return c;
}

inline CollisionCheck check_manifold_warm_start_by_feature() {
  CollisionCheck c{"collision manifold warm start survives by feature key", true, ""};

  ContactPatch p{};
  p.body_a                = 0;
  p.body_b                = 1;
  p.proxy_a               = 0;
  p.proxy_b               = 1;
  p.count                 = 1;
  p.points[0].feature_key = 42;
  p.points[0].rA          = Vec3{0.1f, 0, 0};
  p.points[0].rB          = Vec3{-0.1f, 0, 0};
  p.points[0].xA          = Vec3{0.1f, 0, 0};
  p.points[0].xB          = Vec3{0.1f, 0, 0};

  ManifoldStore store;
  store.begin({});
  store.match(p);
  store.commit();
  if (store.manifolds().size() != 1) {
    c.ok     = false;
    c.detail = "first manifold missing";
    return c;
  }
  // 测试专用：直接给上一帧的 contact 写入 warm start 状态。
  Manifold &seed              = const_cast<Manifold &>(store.manifolds()[0]);
  seed.contacts[0].lambda     = Vec3{-3, 0.25f, -0.5f};
  seed.contacts[0].penalty    = Vec3{7, 8, 9};
  seed.contacts[0].beta_scale = 4;
  seed.contacts[0].stick      = true;

  store.begin({});
  store.match(p);
  store.commit();
  if (store.manifolds().size() != 1) {
    c.ok     = false;
    c.detail = "second manifold missing";
    return c;
  }
  const Manifold &m = store.manifolds()[0];
  if (m.contacts[0].lambda.x != -3 || m.contacts[0].penalty.y != 8 ||
      m.contacts[0].beta_scale != 4 || !m.contacts[0].stick) {
    c.ok     = false;
    c.detail = "feature-matched warm start was not carried over";
  }
  return c;
}

inline CollisionCheck check_manifold_store_is_sorted_and_drops_old() {
  CollisionCheck c{"collision manifold store sorted + stale pair removed", true, ""};

  ContactPatch p1{};
  p1.body_a                = 0;
  p1.body_b                = 1;
  p1.proxy_a               = 0;
  p1.proxy_b               = 1;
  p1.count                 = 1;
  p1.points[0].feature_key = 1;
  ContactPatch p2{};
  p2.body_a                = 0;
  p2.body_b                = 2;
  p2.proxy_a               = 0;
  p2.proxy_b               = 2;
  p2.count                 = 1;
  p2.points[0].feature_key = 2;

  ManifoldStore store;
  store.begin({});
  store.match(p1);
  store.match(p2);
  store.commit();
  if (store.size() != 2 || store.manifolds()[0].body_b != 1 || store.manifolds()[1].body_b != 2) {
    c.ok     = false;
    c.detail = "manifolds are not in ascending pair order";
    return c;
  }

  // 下一帧只有 (0,2)：旧 (0,1) 必须消失。
  store.begin({});
  store.match(p2);
  store.commit();
  if (store.size() != 1 || store.manifolds()[0].body_b != 2) {
    c.ok     = false;
    c.detail = std::format("stale manifold not dropped (size={})", store.size());
  }
  return c;
}

inline CollisionCheck check_contact_constraint_layout() {
  CollisionCheck c{"collision contact constraints contiguous + targets", true, ""};

  ContactPatch p1{};
  p1.body_a                = 0;
  p1.body_b                = 1;
  p1.proxy_a               = 0;
  p1.proxy_b               = 1;
  p1.count                 = 2;
  p1.points[0].feature_key = 10;
  p1.points[1].feature_key = 11;
  p1.basis.normal          = Vec3{0, 1, 0};
  p1.basis.tangent0        = Vec3{1, 0, 0};
  p1.basis.tangent1        = Vec3{0, 0, 1};
  ContactPatch p2{};
  p2.body_a                = 0;
  p2.body_b                = 2;
  p2.proxy_a               = 0;
  p2.proxy_b               = 2;
  p2.count                 = 1;
  p2.points[0].feature_key = 20;

  ManifoldStore store;
  store.begin({});
  store.match(p1);
  store.match(p2);
  store.commit();

  ConstraintList list;
  list.emplace_back(DistanceConstraint{{0, 2}, Scalar_t(1)});
  ContactEmitOptions opt;
  opt.margin                = Scalar_t(0.005);
  const std::size_t emitted = emit_contact_constraints(store.manifolds(), opt, list);
  if (emitted != 9 || list.size() != 10) {
    c.ok     = false;
    c.detail = std::format("emitted={} list={}", emitted, list.size());
    return c;
  }
  if (!std::holds_alternative<DistanceConstraint>(list[0])) {
    c.ok     = false;
    c.detail = "non-contact constraint was rewritten";
    return c;
  }

  const std::uint32_t want_id[9] = {0, 0, 0, 0, 0, 0, 1, 1, 1};
  const std::uint8_t  want_ci[9] = {0, 0, 0, 1, 1, 1, 0, 0, 0};
  for (std::size_t i = 0; i < 9; ++i) {
    if (!std::holds_alternative<ContactConstraint>(list[i + 1])) {
      c.ok     = false;
      c.detail = std::format("slot {} is not a ContactConstraint", i + 1);
      return c;
    }
    const auto &cc = std::get<ContactConstraint>(list[i + 1]);
    if (cc.manifold_id != want_id[i] || cc.contact_index != want_ci[i] ||
        static_cast<int>(cc.kind) != static_cast<int>(i % 3)) {
      c.ok     = false;
      c.detail = std::format("slot {} id={} ci={} kind={}", i + 1, cc.manifold_id, cc.contact_index,
                             static_cast<int>(cc.kind));
      return c;
    }
    const Scalar_t want_target = (cc.kind == ContactKind::Normal) ? -opt.margin : Scalar_t(0);
    if (cc.target != want_target) {
      c.ok     = false;
      c.detail = std::format("slot {} target={} want {}", i + 1, cc.target, want_target);
      return c;
    }
  }
  return c;
}

inline CollisionCheck check_collision_world_step_and_version() {
  CollisionCheck c{"collision world step emits contacts and bumps topology once", true, ""};

  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList list;
  std::uint64_t  version  = 0;
  const bool     changed1 = world.sync(s, list, version);
  if (!changed1 || version != 1 || world.manifolds().empty()) {
    c.ok     = false;
    c.detail = std::format("changed={} version={} manifolds={}", changed1, version,
                           world.manifolds().size());
    return c;
  }
  std::size_t contact_constraints = 0;
  for (const auto &v : list) {
    if (std::holds_alternative<ContactConstraint>(v))
      ++contact_constraints;
  }
  if (contact_constraints == 0) {
    c.ok     = false;
    c.detail = "no ContactConstraint emitted";
    return c;
  }

  // dual 状态必须能从 ContactConstraint 回写到 manifold，否则 warm start 是死的。
  for (auto &v : list) {
    if (!std::holds_alternative<ContactConstraint>(v))
      continue;
    auto &cc  = std::get<ContactConstraint>(v);
    cc.lambda = Scalar_t(-7);
    world.absorb_warm_start(list);
    const Manifold &m = world.manifolds()[cc.manifold_id];
    if (m.contacts[cc.contact_index].lambda[static_cast<int>(cc.kind)] != Scalar_t(-7)) {
      c.ok     = false;
      c.detail = "absorb_warm_start did not copy lambda back";
      return c;
    }
    break;
  }

  // 同一 pose 再 sync：拓扑没变，version 不许再涨。
  const std::uint64_t v_before = version;
  const bool          changed2 = world.sync(s, list, version);
  if (changed2 || version != v_before) {
    c.ok     = false;
    c.detail = std::format("static pose bumped version {} -> {}", v_before, version);
    return c;
  }

  // 分离：签名变化，version 必须 bump，接触约束必须清空。
  s.pose.set_position(1, Vec3{0, 5.0f, 0});
  const bool changed3 = world.sync(s, list, version);
  if (!changed3 || version <= v_before) {
    c.ok     = false;
    c.detail = "separating pose did not bump topology version";
    return c;
  }
  contact_constraints = 0;
  for (const auto &v : list) {
    if (std::holds_alternative<ContactConstraint>(v))
      ++contact_constraints;
  }
  if (contact_constraints != 0) {
    c.ok     = false;
    c.detail = "contact constraints survived separation";
  }
  return c;
}

inline CollisionCheck check_collision_end_to_end_resting_box() {
  CollisionCheck c{"collision end-to-end: dynamic box pushed out, static box frozen", true, ""};

  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList       list;
  std::uint64_t        version = 0;
  AVBDStep             policy;
  task::SerialExecutor serial;
  const SolveSchedule  sched{};
  const Vec3           static_before = s.pose.position(0);
  const Scalar_t       start_y       = s.pose.position(1).y;

  for (int step = 0; step < 30; ++step) {
    s.reset_external();
    s.begin_step();
    world.sync(s, list, version);
    solve(s, list, policy, 20, serial, sched, static_cast<NullPhaseObserver *>(nullptr), version);
    s.finish_step();
    s.cache.commit_prev(s.velocity);
  }

  const Vec3 static_after = s.pose.position(0);
  if (std::memcmp(&static_before, &static_after, sizeof(Vec3)) != 0) {
    c.ok     = false;
    c.detail = "static box moved";
    return c;
  }
  const Vec3 p1 = s.pose.position(1);
  if (!all_finite(p1) || !all_finite(s.velocity.velocity(1)) || std::abs(p1.y) > 10.0f) {
    c.ok     = false;
    c.detail = std::format("dynamic box diverged: y={}", p1.y);
    return c;
  }
  // 初始 0.1 m 穿透，接触应把它推回到接近 1.0（允许 margin）附近。
  if (!(p1.y > start_y + 0.01f)) {
    c.ok     = false;
    c.detail = std::format("dynamic box was not pushed out: y {} -> {}", start_y, p1.y);
    return c;
  }
  return c;
}


// ==================== 可观测摘要（step 1）====================
// 这些值不进入 tests/baseline.txt，也不要求跨算法版本逐位一致。它们是 collision
// 层的“体检表”：统一 PolytopeView/SAT 管线之后，用同一组场景对比 manifold / 接触
// 约束到底变了哪些部分。主基线继续只盯求解器状态。
//
// 保存/对比入口：
//   avbd-mo --save-collision tests/collision_observable.txt
//   avbd-mo --check-collision tests/collision_observable.txt
//   avbd-mo --collision-digest            # 直接打到 stdout

struct NarrowDigest {
  std::uint64_t all   = 0;
  std::size_t   count = 0;
};

inline NarrowDigest narrowphase_digest(const ShapeView &a, const ShapeView &b) {
  ContactPatch patch{};
  patch.body_a   = 0;
  patch.body_b   = 1;
  patch.friction = std::sqrt(std::max(Scalar_t(0), a.friction) *
                             std::max(Scalar_t(0), b.friction));
  patch.count    = collide_shapes(a, b, std::span<ContactPoint>{patch.points}, patch.basis);

  std::uint64_t h = collision_digest::kFnvOffset;
  h               = collision_digest::shape_view(h, a);
  h               = collision_digest::shape_view(h, b);
  h               = collision_digest::contact_patch(h, patch);
  return NarrowDigest{h, patch.count};
}

inline std::string narrow_line(std::string_view name, const ShapeView &a, const ShapeView &b) {
  const NarrowDigest d = narrowphase_digest(a, b);
  return std::format("{} all={} contacts={}", name, collision_digest::to_hex(d.all), d.count);
}

inline std::string world_line(std::string_view name, const CollisionWorld &world,
                              const ConstraintList &list, std::uint64_t version) {
  std::uint64_t h = collision_digest::kFnvOffset;
  h = collision_digest::mix_u64(h, collision_digest::manifolds(world.manifolds()));
  h = collision_digest::mix_u64(h, collision_digest::constraints(list));

  const CollisionWorldStats &st = world.stats();
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.bodies));
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.proxies));
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.broad_pairs));
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.narrow_pairs));
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.contact_points));
  h = collision_digest::mix_u64(h, static_cast<std::uint64_t>(st.manifolds));
  h = collision_digest::mix_u64(h, world.layout_signature());
  h = collision_digest::mix_u64(h, version);

  std::size_t contacts = 0;
  for (const Manifold &m : world.manifolds())
    contacts += m.count;

  return std::format("{} all={} manifolds={} contacts={} constraints={}", name,
                     collision_digest::to_hex(h), world.manifolds().size(), contacts,
                     list.size());
}

inline std::string observable_world_box_pair() {
  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);
  return world_line("world/box-pair", world, list, version);
}

inline std::string observable_world_box_stack() {
  State s;
  s.resize(3);
  s.dt = Scalar_t(1) / 60;
  for (std::size_t b = 0; b < 3; ++b) {
    s.pose.set_position(b, Vec3{0, Scalar_t(0.9) * Scalar_t(b), 0});
    s.inertial.set_inv_mass(b, b == 0 ? 0.0f : 1.0f);
    if (b != 0)
      s.inertial.set_local_inv_inertia(b, Mat3(1));
  }

  CollisionWorld world;
  world.resize(3);
  for (std::size_t b = 0; b < 3; ++b)
    set_single_shape(world, b, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);
  return world_line("world/box-stack", world, list, version);
}

inline std::string observable_world_mixed() {
  State s;
  s.resize(3);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.pose.set_position(2, Vec3{0.6f, 1.2f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_inv_mass(2, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));
  s.inertial.set_local_inv_inertia(2, Mat3(1));

  CollisionWorld world;
  world.resize(3);
  set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 2, Shape::sphere(0.5f));

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);
  return world_line("world/mixed", world, list, version);
}

inline std::string observable_world_compound() {
  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  const std::array<ShapeInstance, 2> parts{
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{-0.6f, 0, 0}},
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0.6f, 0, 0}},
  };
  world.set_shapes(0, parts);
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);
  return world_line("world/compound", world, list, version);
}

inline std::string observable_world_convex() {
  ConvexMesh cube = make_cube_convex_mesh();
  if (finalize_convex_mesh(cube) != ConvexMeshStatus::Ok)
    return "world/convex all=finalize-failed";

  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  const ConvexMeshId id = world.add_convex_mesh(std::move(cube));
  if (id == kInvalidConvexMesh)
    return "world/convex all=add-failed";
  const std::array<ShapeInstance, 1> a{ShapeInstance{Shape::convex(id)}};
  const std::array<ShapeInstance, 1> b{ShapeInstance{Shape::convex(id)}};
  world.set_shapes(0, a);
  world.set_shapes(1, b);

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);
  return world_line("world/convex", world, list, version);
}

inline std::vector<std::string> observable_lines() {
  const Shape box      = Shape::box(Vec3{0.5f, 0.5f, 0.5f});
  const Shape sphere   = Shape::sphere(0.5f);
  const Shape cylinder = Shape::cylinder(0.5f, 1.0f);

  const ShapeView box_a = make_shape_view(box, Vec3{0, 0, 0}, identity_quat());
  const ShapeView box_b = make_shape_view(box, Vec3{0, 0.9f, 0}, identity_quat());
  const ShapeView box_rot =
      make_shape_view(box, Vec3{0, 0.9f, 0}, quat_from_rotvec(Vec3{0, 0, 0.5235988f}));
  const ShapeView sphere_a = make_shape_view(sphere, Vec3{0, 0, 0}, identity_quat());
  const ShapeView sphere_b = make_shape_view(sphere, Vec3{0.8f, 0, 0}, identity_quat());
  const ShapeView sphere_c = make_shape_view(sphere, Vec3{0.6f, 0, 0}, identity_quat());
  const ShapeView cyl_a    = make_shape_view(cylinder, Vec3{0, 0, 0}, identity_quat());
  const ShapeView cyl_b    = make_shape_view(cylinder, Vec3{0.4f, 0, 0}, identity_quat());

  std::vector<std::string> out;
  out.push_back(narrow_line("narrow/box-box", box_a, box_b));
  out.push_back(narrow_line("narrow/box-box-rot", box_a, box_rot));
  out.push_back(narrow_line("narrow/sphere-sphere", sphere_a, sphere_b));
  out.push_back(narrow_line("narrow/box-sphere", box_a, sphere_b));
  out.push_back(narrow_line("narrow/cylinder-box", cyl_a, box_b));
  out.push_back(narrow_line("narrow/cylinder-cylinder", cyl_a, cyl_b));
  out.push_back(narrow_line("narrow/sphere-cylinder", sphere_c, cyl_a));
  out.push_back(observable_world_box_pair());
  out.push_back(observable_world_box_stack());
  out.push_back(observable_world_mixed());
  out.push_back(observable_world_compound());
  out.push_back(observable_world_convex());
  return out;
}

inline CollisionCheck check_compound_proxy_pipeline() {
  CollisionCheck c{"collision compound produces one manifold per shape pair", true, ""};

  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  const std::array<ShapeInstance, 2> parts{
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{-0.6f, 0, 0}},
      ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0.6f, 0, 0}},
  };
  world.set_shapes(0, parts);
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  ConstraintList list;
  std::uint64_t  version = 0;
  world.sync(s, list, version);

  const std::span<const Manifold> ms = world.manifolds();
  if (ms.size() != 2) {
    c.ok     = false;
    c.detail = std::format("expected 2 manifolds (one per shape pair), got {}", ms.size());
    return c;
  }
  if (ms[0].body_a != 0 || ms[0].body_b != 1 || ms[1].body_a != 0 || ms[1].body_b != 1) {
    c.ok     = false;
    c.detail = "manifold body ids are wrong";
    return c;
  }
  if (!(ms[0].proxy_a < ms[0].proxy_b && ms[1].proxy_a < ms[1].proxy_b) ||
      ms[0].proxy_a == ms[1].proxy_a) {
    c.ok     = false;
    c.detail = std::format("proxy keys not distinct/sorted: ({},{}) ({},{})", ms[0].proxy_a,
                           ms[0].proxy_b, ms[1].proxy_a, ms[1].proxy_b);
    return c;
  }

  // rA 必须是 body-local：两个 part 的接触锚点应该分居 body 原点的两侧。
  // 如果还是 shape-local，两个 manifold 的 rA.x 会几乎相同。
  Scalar_t rA_x[2] = {0, 0};
  for (std::size_t i = 0; i < 2; ++i) {
    rA_x[i] = ms[i].contacts[0].rA.x;
    if (ms[i].count == 0) {
      c.ok     = false;
      c.detail = std::format("manifold {} has no contacts", i);
      return c;
    }
  }
  if (!(rA_x[0] < 0.0f && rA_x[1] > 0.0f)) {
    c.ok     = false;
    c.detail = std::format("rA is not body-local: manifold0.x={} manifold1.x={}", rA_x[0],
                           rA_x[1]);
    return c;
  }

  // 同一 body pair 的两个 manifold 必须映射到不同的 manifold_id，约束连续。
  std::uint32_t seen_id[2] = {kStandaloneManifold, kStandaloneManifold};
  std::size_t   seen       = 0;
  for (const ConstraintVariant &v : list) {
    if (!std::holds_alternative<ContactConstraint>(v))
      continue;
    const ContactConstraint &cc = std::get<ContactConstraint>(v);
    if (cc.bodies_[0] != 0 || cc.bodies_[1] != 1)
      continue;
    if (seen == 0) {
      seen_id[0] = cc.manifold_id;
      seen       = 1;
    } else if (seen == 1 && cc.manifold_id != seen_id[0]) {
      seen_id[1] = cc.manifold_id;
      seen       = 2;
      break;
    }
  }
  if (seen < 2 || seen_id[0] == seen_id[1]) {
    c.ok     = false;
    c.detail = "contact constraints did not cover both manifolds";
    return c;
  }

  // warm start：同一 proxy pair + 同一 feature 必须把 lambda 带过下一帧。
  for (ConstraintVariant &v : list) {
    if (std::holds_alternative<ContactConstraint>(v))
      std::get<ContactConstraint>(v).lambda = Scalar_t(-3);
  }
  world.absorb_warm_start(list);
  world.sync(s, list, version);
  bool warm_started = false;
  for (const Manifold &m : world.manifolds())
    for (std::size_t i = 0; i < m.count; ++i)
      warm_started = warm_started || m.contacts[i].lambda[0] == Scalar_t(-3);
  if (!warm_started) {
    c.ok     = false;
    c.detail = "compound manifold warm start did not preserve lambda";
  }
  return c;
}

inline CollisionCheck check_compound_filters() {
  CollisionCheck c{"collision compound self/static filtering + proxy culling", true, ""};

  // 1) 同一 body 的两个重叠 part 不得自碰撞。
  {
    State s;
    s.resize(1);
    s.dt = Scalar_t(1) / 60;
    s.inertial.set_inv_mass(0, 1);
    s.inertial.set_local_inv_inertia(0, Mat3(1));
    CollisionWorld world;
    world.resize(1);
    const std::array<ShapeInstance, 2> parts{
        ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f})},
        ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{0.2f, 0, 0}},
    };
    world.set_shapes(0, parts);
    ConstraintList list;
    std::uint64_t  version = 0;
    world.sync(s, list, version);
    if (!world.manifolds().empty()) {
      c.ok     = false;
      c.detail = std::format("self-collision produced {} manifolds", world.manifolds().size());
      return c;
    }
  }

  // 2) 两个静态 body 不产生接触。
  {
    State s;
    s.resize(2);
    s.dt = Scalar_t(1) / 60;
    s.pose.set_position(1, Vec3{0, 0.9f, 0});
    s.inertial.set_inv_mass(0, 0);
    s.inertial.set_inv_mass(1, 0);
    CollisionWorld world;
    world.resize(2);
    set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
    set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
    ConstraintList list;
    std::uint64_t  version = 0;
    world.sync(s, list, version);
    if (!world.manifolds().empty()) {
      c.ok     = false;
      c.detail = "static-static pair produced contacts";
      return c;
    }
  }

  // 3) proxy 级宽相：body A 的两个 part 分别位于 x=0 和 x=100，body B 在 x=50。
  // body 级 AABB 会把 A 合并成 [−0.5,100.5] 并与 B 相交；proxy 级则完全不交。
  {
    State s;
    s.resize(2);
    s.dt = Scalar_t(1) / 60;
    s.pose.set_position(0, Vec3{0, 0, 0});
    s.pose.set_position(1, Vec3{50, 0, 0});
    s.inertial.set_inv_mass(0, 0);
    s.inertial.set_inv_mass(1, 1);
    s.inertial.set_local_inv_inertia(1, Mat3(1));
    CollisionWorld world;
    world.resize(2);
    const std::array<ShapeInstance, 2> parts{
        ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f})},
        ShapeInstance{Shape::box(Vec3{0.5f, 0.5f, 0.5f}), Vec3{100, 0, 0}},
    };
    world.set_shapes(0, parts);
    set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
    ConstraintList list;
    std::uint64_t  version = 0;
    world.sync(s, list, version);
    if (world.stats().broad_pairs != 0 || !world.manifolds().empty()) {
      c.ok     = false;
      c.detail = std::format("proxy broadphase did not cull: broad_pairs={} manifolds={}",
                             world.stats().broad_pairs, world.manifolds().size());
    }
  }
  return c;
}

inline CollisionCheck check_collision_world_pair_filter() {
  CollisionCheck c{"collision world pair filter gates layer/mask pairs", true, ""};

  // 两个动态盒在几何上重叠；只有过滤规则决定它们是否产生接触。
  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 1);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(0, Mat3(1));
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  set_single_shape(world, 0, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));
  set_single_shape(world, 1, Shape::box(Vec3{0.5f, 0.5f, 0.5f}));

  // 规则表：Godot 语义 —— 一方 layer 命中另一方 mask 即可（双向 OR）。
  struct LayerMask {
    std::uint32_t layer[2];
    std::uint32_t mask[2];
  } rule{{1, 2}, {1, 2}};
  const auto filter = [](void *user, std::size_t a, std::size_t b) noexcept {
    const auto *r = static_cast<const LayerMask *>(user);
    if (a >= 2 || b >= 2)
      return true;
    return (r->layer[a] & r->mask[b]) != 0 || (r->layer[b] & r->mask[a]) != 0;
  };

  ConstraintList list;
  std::uint64_t  version = 0;

  // 1) layer(1)/(2) 与对方的 mask 都不相交 ⇒ 不产生接触。
  world.set_pair_filter(filter, &rule);
  world.sync(s, list, version);
  if (!world.manifolds().empty()) {
    c.ok     = false;
    c.detail = std::format("filtered pair produced {} manifolds", world.manifolds().size());
    return c;
  }

  // 2) 让 body 0 的 mask 含 body 1 的 layer ⇒ 接触出现，版本号必须 bump。
  rule.mask[0] = 3;
  world.set_pair_filter(filter, &rule);
  const std::uint64_t before = version;
  world.sync(s, list, version);
  if (world.manifolds().empty() || version <= before) {
    c.ok     = false;
    c.detail = std::format("allowed pair manifolds={} version {} -> {}", world.manifolds().size(),
                           before, version);
    return c;
  }

  // 3) 清除过滤器后行为与从未设置过一致（除内建规则外全部允许）。
  world.set_pair_filter(nullptr);
  world.sync(s, list, version);
  if (world.manifolds().empty()) {
    c.ok     = false;
    c.detail = "clearing the filter removed contacts";
  }
  return c;
}

inline CollisionCheck check_convex_world_end_to_end() {
  CollisionCheck c{"collision convex world end-to-end", true, ""};

  ConvexMesh cube = make_cube_convex_mesh();
  if (finalize_convex_mesh(cube) != ConvexMeshStatus::Ok) {
    c.ok     = false;
    c.detail = "cube finalize failed";
    return c;
  }

  State s;
  s.resize(2);
  s.dt = Scalar_t(1) / 60;
  s.pose.set_position(0, Vec3{0, 0, 0});
  s.pose.set_position(1, Vec3{0, 0.9f, 0});
  s.inertial.set_inv_mass(0, 0);
  s.inertial.set_inv_mass(1, 1);
  s.inertial.set_local_inv_inertia(1, Mat3(1));

  CollisionWorld world;
  world.resize(2);
  const ConvexMeshId id = world.add_convex_mesh(std::move(cube));
  const std::array<ShapeInstance, 1> shape{ShapeInstance{Shape::convex(id)}};
  world.set_shapes(0, shape);
  world.set_shapes(1, shape);

  ConstraintList       list;
  std::uint64_t        version = 0;
  AVBDStep             policy;
  task::SerialExecutor serial;
  const SolveSchedule  sched{};
  const Vec3           static_before = s.pose.position(0);
  const Scalar_t       start_y       = s.pose.position(1).y;

  for (int step = 0; step < 30; ++step) {
    s.reset_external();
    s.begin_step();
    world.sync(s, list, version);
    solve(s, list, policy, 20, serial, sched, static_cast<NullPhaseObserver *>(nullptr), version);
    s.finish_step();
    s.cache.commit_prev(s.velocity);
  }

  const Vec3 static_after = s.pose.position(0);
  if (std::memcmp(&static_before, &static_after, sizeof(Vec3)) != 0) {
    c.ok     = false;
    c.detail = "static convex moved";
    return c;
  }
  const Vec3 p1 = s.pose.position(1);
  if (!all_finite(p1) || !all_finite(s.velocity.velocity(1)) || std::abs(p1.y) > 10.0f) {
    c.ok     = false;
    c.detail = std::format("dynamic convex diverged: y={}", p1.y);
    return c;
  }
  if (!(p1.y > start_y + 0.001f)) {
    c.ok     = false;
    c.detail = std::format("dynamic convex was not pushed out: y {} -> {}", start_y, p1.y);
  }
  return c;
}

inline CollisionCheck check_collision_digest_is_repeatable() {
  CollisionCheck c{"collision observable digest is repeatable", true, ""};
  const std::vector<std::string> a = observable_lines();
  const std::vector<std::string> b = observable_lines();
  if (a.size() != b.size()) {
    c.ok     = false;
    c.detail = "case count changed between runs";
    return c;
  }
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) {
      c.ok     = false;
      c.detail = std::format("case {} differs:\n     {} \n     {}", i + 1, a[i], b[i]);
      return c;
    }
    if (a[i].find("all=") == std::string::npos) {
      c.ok     = false;
      c.detail = std::format("case {} has no digest field: {}", i + 1, a[i]);
      return c;
    }
  }
  return c;
}

inline CollisionCheck check_collision_digest_tracks_input() {
  CollisionCheck c{"collision observable digest is sensitive to geometry", true, ""};
  const Shape     box = Shape::box(Vec3{0.5f, 0.5f, 0.5f});
  const ShapeView a   = make_shape_view(box, Vec3{0, 0, 0}, identity_quat());
  const ShapeView b1  = make_shape_view(box, Vec3{0, 0.9f, 0}, identity_quat());
  const ShapeView b2  = make_shape_view(box, Vec3{0, 0.8f, 0}, identity_quat());

  const NarrowDigest d1 = narrowphase_digest(a, b1);
  const NarrowDigest d2 = narrowphase_digest(a, b2);
  if (d1.count == 0 || d2.count == 0) {
    c.ok     = false;
    c.detail = "expected contacts for both box poses";
    return c;
  }
  if (d1.all == d2.all) {
    c.ok     = false;
    c.detail = "moving the second box did not change the digest";
  }
  return c;
}

inline std::vector<CollisionCheck> run() {
  return {
      check_shape_table_layout(),
      check_convex_mesh_and_hull(),
      check_convex_polytope_collision(),
      check_mass_properties(),
      check_broadphase_matches_bruteforce(),
      check_narrowphase_sphere_sphere(),
      check_narrowphase_box_stack_normal(),
      check_narrowphase_sphere_box(),
      check_narrowphase_cylinder_box(),
      check_manifold_warm_start_by_feature(),
      check_manifold_store_is_sorted_and_drops_old(),
      check_contact_constraint_layout(),
      check_collision_world_step_and_version(),
      check_collision_end_to_end_resting_box(),
      check_convex_world_end_to_end(),
      check_collision_digest_is_repeatable(),
      check_collision_digest_tracks_input(),
      check_compound_proxy_pipeline(),
      check_compound_filters(),
      check_collision_world_pair_filter(),
  };
}

} // namespace collision_selftest
