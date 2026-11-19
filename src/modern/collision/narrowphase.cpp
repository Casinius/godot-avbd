#include "collision/narrowphase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace collision {

namespace {

constexpr Scalar_t kSatAxisEpsilon    = Scalar_t(1e-6);
constexpr Scalar_t kPlaneEpsilon      = Scalar_t(1e-5);
constexpr Scalar_t kContactMergeDist2 = Scalar_t(1e-6);
constexpr int      kMaxPolyVerts      = 16;
// 预接触必须严格记录一个正间隙：depth = -gap < 0 是"这是 speculative 点"的
// 唯一标记。gap == 0（正好接触）交给 exact 路径；否则 -0.0 不是负数，
// 会被当成普通接触发射出去，而 C = margin > 0 会产生吸引。
constexpr Scalar_t kSpecGapEpsilon = Scalar_t(1e-6);

[[nodiscard]] inline int clip_polygon_against_plane(const Vec3 *in, int in_count, Vec3 plane_normal,
                                                    Scalar_t plane_offset, Vec3 *out) noexcept {
  if (in_count <= 0)
    return 0;
  int      out_count = 0;
  Vec3     a         = in[in_count - 1];
  Scalar_t da        = glm::dot(plane_normal, a) - plane_offset;
  for (int i = 0; i < in_count; ++i) {
    const Vec3     b        = in[i];
    const Scalar_t db       = glm::dot(plane_normal, b) - plane_offset;
    const bool     a_inside = da <= kPlaneEpsilon;
    const bool     b_inside = db <= kPlaneEpsilon;
    if (a_inside != b_inside) {
      Scalar_t       t     = 0;
      const Scalar_t denom = da - db;
      if (std::abs(denom) > kSatAxisEpsilon)
        t = std::clamp(da / denom, Scalar_t(0), Scalar_t(1));
      if (out_count < kMaxPolyVerts)
        out[out_count++] = a + (b - a) * t;
    }
    if (b_inside && out_count < kMaxPolyVerts)
      out[out_count++] = b;
    a  = b;
    da = db;
  }
  return out_count;
}

[[nodiscard]] inline Scalar_t contact_depth(Vec3 normal_b_to_a, Vec3 xA, Vec3 xB) noexcept {
  return std::max(Scalar_t(0), -glm::dot(normal_b_to_a, xA - xB));
}

inline bool add_contact(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                        std::size_t &count, Vec3 *midpoints, Vec3 xA, Vec3 xB,
                        std::uint64_t feature_key, Vec3 normal_b_to_a) noexcept {
  const Vec3 mid = (xA + xB) * Scalar_t(0.5);
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3 d = mid - midpoints[i];
    if (glm::dot(d, d) < kContactMergeDist2)
      return false;
  }
  if (count >= out.size() || count >= kMaxContactsPerManifold)
    return false;

  const Vec3    localA = a.world_to_body_local(xA);
  const Vec3    localB = b.world_to_body_local(xB);
  ContactPoint &c      = out[count];
  c.feature_key        = feature_key;
  c.rA                 = localA;
  c.rB                 = localB;
  c.xA                 = xA;
  c.xB                 = xB;
  c.depth              = contact_depth(normal_b_to_a, xA, xB);
  midpoints[count]     = mid;
  ++count;
  return true;
}

// speculative 版本：接触点表示的是**间隙**而不是穿透，depth 记成 -gap（负值）。
// 合并规则与 add_contact 相同，保证同一 manifold 内的点序仍然只由几何决定。
inline bool add_contact_gap(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                            std::size_t &count, Vec3 *midpoints, Vec3 xA, Vec3 xB,
                            std::uint64_t feature_key, Scalar_t gap) noexcept {
  if (!(gap > Scalar_t(0)))
    return false;
  const Vec3 mid = (xA + xB) * Scalar_t(0.5);
  for (std::size_t i = 0; i < count; ++i) {
    const Vec3 d = mid - midpoints[i];
    if (glm::dot(d, d) < kContactMergeDist2)
      return false;
  }
  if (count >= out.size() || count >= kMaxContactsPerManifold)
    return false;

  ContactPoint &c = out[count];
  c.feature_key   = feature_key;
  c.rA            = a.world_to_body_local(xA);
  c.rB            = b.world_to_body_local(xB);
  c.xA            = xA;
  c.xB            = xB;
  c.depth         = -gap; // gap > 0 ⇒ 严格负值
  midpoints[count] = mid;
  ++count;
  return true;
}

inline void closest_points_on_segments(Vec3 p0, Vec3 p1, Vec3 q0, Vec3 q1, Vec3 &c0,
                                       Vec3 &c1) noexcept {
  const Vec3     d1 = p1 - p0;
  const Vec3     d2 = q1 - q0;
  const Vec3     r  = p0 - q0;
  const Scalar_t aa = glm::dot(d1, d1);
  const Scalar_t ee = glm::dot(d2, d2);
  const Scalar_t ff = glm::dot(d2, r);
  Scalar_t       s = 0, t = 0;
  if (aa <= kSatAxisEpsilon && ee <= kSatAxisEpsilon) {
    c0 = p0;
    c1 = q0;
    return;
  }
  if (aa <= kSatAxisEpsilon) {
    t = std::clamp(ff / ee, Scalar_t(0), Scalar_t(1));
  } else {
    const Scalar_t cc = glm::dot(d1, r);
    if (ee <= kSatAxisEpsilon) {
      s = std::clamp(-cc / aa, Scalar_t(0), Scalar_t(1));
    } else {
      const Scalar_t bb    = glm::dot(d1, d2);
      const Scalar_t denom = aa * ee - bb * bb;
      if (std::abs(denom) > kSatAxisEpsilon)
        s = std::clamp((bb * ff - cc * ee) / denom, Scalar_t(0), Scalar_t(1));
      t = (bb * s + ff) / ee;
      if (t < 0) {
        t = 0;
        s = std::clamp(-cc / aa, Scalar_t(0), Scalar_t(1));
      } else if (t > 1) {
        t = 1;
        s = std::clamp((bb - cc) / aa, Scalar_t(0), Scalar_t(1));
      }
    }
  }
  c0 = p0 + d1 * s;
  c1 = q0 + d2 * t;
}

std::size_t collide_sphere_sphere(const ShapeView &a, const ShapeView &b,
                                  std::span<ContactPoint> out, BasisFrame &basis) noexcept {
  const Vec3     d   = b.center - a.center;
  const Scalar_t d2  = glm::dot(d, d);
  const Scalar_t sum = a.radius + b.radius;
  if (d2 >= sum * sum)
    return 0;
  const Scalar_t dist                             = std::sqrt(d2);
  const Vec3     normal_ab                        = dist > 1e-6f ? d / dist : Vec3{0, 0, 1};
  basis                                           = make_basis(-normal_ab);
  std::size_t                               count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact(a, b, out, count, midpoints.data(), a.center + normal_ab * a.radius,
              b.center - normal_ab * b.radius, 1u, basis.normal);
  return count;
}

// 间隙 ≤ max_gap 的球/球预接触。法线/见证点约定与重叠路径完全一致；
// 只有 depth 变成负的 gap。exact 路径已经返回 0 才会走到这里，所以 gap > 0。
std::size_t collide_sphere_sphere_speculative(const ShapeView &a, const ShapeView &b,
                                              std::span<ContactPoint> out, BasisFrame &basis,
                                              Scalar_t max_gap) noexcept {
  const Vec3     d   = b.center - a.center;
  const Scalar_t d2  = glm::dot(d, d);
  const Scalar_t sum = a.radius + b.radius;
  const Scalar_t reach = sum + max_gap;
  if (d2 >= reach * reach)
    return 0;
  const Scalar_t dist = std::sqrt(d2);
  if (!(dist > kSatAxisEpsilon) || !(dist >= sum))
    return 0;
  const Vec3 normal_ab = d / dist;
  const Scalar_t gap   = dist - sum;
  if (!(gap > kSpecGapEpsilon) || gap > max_gap)
    return 0;
  basis = make_basis(-normal_ab);
  std::size_t                               count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact_gap(a, b, out, count, midpoints.data(), a.center + normal_ab * a.radius,
                  b.center - normal_ab * b.radius, 1u, gap);
  return count;
}

// ============================================================================
//  统一 polytope 窄相：Box / Cylinder(棱柱化) / Convex
// ----------------------------------------------------------------------------
//  这里不再为每种形状写一个 pair 函数。所有非球形状都通过 PolytopeView 暴露
//  vertices/faces/edges，然后共用 SAT + 参考面裁剪；Sphere 走最接近点路径。
// ============================================================================

inline constexpr int         kCylinderSides   = 16;
inline constexpr Scalar_t    kPolytopeEpsilon = Scalar_t(1e-5);
inline constexpr std::size_t kMaxSatEdgeDirs  = 64;

enum class PolytopeKind : std::uint8_t { Box = 0, Cylinder = 1, Convex = 2 };

[[nodiscard]] inline std::uint64_t make_feature_key(std::uint64_t tag, std::uint64_t a,
                                                    std::uint64_t b, std::uint64_t c,
                                                    std::uint64_t d = 0) noexcept {
  std::uint64_t       h    = 1469598103934665603ull;
  const std::uint64_t v[5] = {tag, a, b, c, d};
  for (std::uint64_t x : v) {
    h ^= x;
    h *= 1099511628211ull;
    h ^= h >> 32;
  }
  return h;
}

struct PolytopeView {
  const ShapeView *shape = nullptr;
  PolytopeKind     kind  = PolytopeKind::Box;

  explicit PolytopeView(const ShapeView &v) noexcept : shape(&v) {
    switch (v.type) {
    case ShapeType::Box:
      kind = PolytopeKind::Box;
      break;
    case ShapeType::Cylinder:
      kind = PolytopeKind::Cylinder;
      break;
    case ShapeType::Convex:
      kind = PolytopeKind::Convex;
      break;
    default:
      kind = PolytopeKind::Box;
      break;
    }
  }

  [[nodiscard]] Vec3 center() const noexcept { return shape->center; }

  [[nodiscard]] std::uint32_t vertex_count() const noexcept {
    switch (kind) {
    case PolytopeKind::Box:
      return 8;
    case PolytopeKind::Cylinder:
      return static_cast<std::uint32_t>(2 * kCylinderSides);
    case PolytopeKind::Convex:
      return shape->convex != nullptr ? static_cast<std::uint32_t>(shape->convex->vertices.size())
                                      : 0;
    }
    return 0;
  }

  [[nodiscard]] Vec3 vertex(std::uint32_t i) const noexcept {
    if (kind == PolytopeKind::Box) {
      const Scalar_t sx = (i & 1u) ? 1.0f : -1.0f;
      const Scalar_t sy = (i & 2u) ? 1.0f : -1.0f;
      const Scalar_t sz = (i & 4u) ? 1.0f : -1.0f;
      return shape->center +
             shape->rot * Vec3{shape->half.x * sx, shape->half.y * sy, shape->half.z * sz};
    }
    if (kind == PolytopeKind::Cylinder) {
      const std::uint32_t ring = i % static_cast<std::uint32_t>(kCylinderSides);
      const bool          top  = i >= static_cast<std::uint32_t>(kCylinderSides);
      const Scalar_t a = Scalar_t(6.283185307179586) * Scalar_t(ring) / Scalar_t(kCylinderSides);
      return shape->center + shape->rot * Vec3{shape->radius * std::cos(a),
                                               top ? shape->half_height : -shape->half_height,
                                               shape->radius * std::sin(a)};
    }
    if (kind == PolytopeKind::Convex && shape->convex != nullptr)
      return shape->center + shape->rot * shape->convex->vertices[i];
    return shape->center;
  }

  [[nodiscard]] std::uint32_t face_count() const noexcept {
    switch (kind) {
    case PolytopeKind::Box:
      return 6;
    case PolytopeKind::Cylinder:
      return static_cast<std::uint32_t>(kCylinderSides + 2);
    case PolytopeKind::Convex:
      return shape->convex != nullptr ? static_cast<std::uint32_t>(shape->convex->faces.size()) : 0;
    }
    return 0;
  }

  [[nodiscard]] std::uint32_t face_vertex_count(std::uint32_t f) const noexcept {
    if (kind == PolytopeKind::Box)
      return 4;
    if (kind == PolytopeKind::Cylinder)
      return f < 2 ? static_cast<std::uint32_t>(kCylinderSides) : 4;
    if (kind == PolytopeKind::Convex && shape->convex != nullptr)
      return shape->convex->faces[f].count;
    return 0;
  }

  [[nodiscard]] Vec3 face_vertex(std::uint32_t f, std::uint32_t j) const noexcept {
    if (kind == PolytopeKind::Box) {
      static constexpr std::uint8_t kBoxFaceVerts[6][4] = {
          {1, 3, 7, 5}, // +X
          {0, 4, 6, 2}, // -X
          {2, 6, 7, 3}, // +Y
          {0, 1, 5, 4}, // -Y
          {4, 5, 7, 6}, // +Z
          {0, 2, 3, 1}, // -Z
      };
      return vertex(kBoxFaceVerts[f][j]);
    }
    if (kind == PolytopeKind::Cylinder) {
      const std::uint32_t n = static_cast<std::uint32_t>(kCylinderSides);
      if (f == 0)
        return vertex(j); // 底盖：角度递增 = -Y 外向
      if (f == 1)
        return vertex(n + (n - 1 - j)); // 顶盖反向
      const std::uint32_t i   = f - 2;
      const std::uint32_t nxt = (i + 1) % n;
      switch (j) {
      case 0:
        return vertex(i); // bottom i
      case 1:
        return vertex(n + i); // top i
      case 2:
        return vertex(n + nxt); // top next
      default:
        return vertex(nxt); // bottom next
      }
    }
    if (kind == PolytopeKind::Convex && shape->convex != nullptr) {
      const ConvexFace   &face = shape->convex->faces[f];
      const std::uint32_t idx  = shape->convex->face_indices[face.first + j];
      return shape->center + shape->rot * shape->convex->vertices[idx];
    }
    return shape->center;
  }

  [[nodiscard]] Vec3 face_normal(std::uint32_t f) const noexcept {
    if (kind == PolytopeKind::Box) {
      const int      axis = static_cast<int>(f) / 2;
      const Scalar_t sign = (f % 2 == 0) ? 1.0f : -1.0f;
      Vec3           local{0};
      local[axis] = sign;
      return shape->rot * local;
    }
    if (kind == PolytopeKind::Cylinder) {
      if (f == 0)
        return shape->rot * Vec3{0, -1, 0};
      if (f == 1)
        return shape->rot * Vec3{0, 1, 0};
      const std::uint32_t i = f - 2;
      const Scalar_t      a =
          Scalar_t(6.283185307179586) * (Scalar_t(i) + Scalar_t(0.5)) / Scalar_t(kCylinderSides);
      return shape->rot * Vec3{std::cos(a), 0, std::sin(a)};
    }
    if (kind == PolytopeKind::Convex && shape->convex != nullptr)
      return shape->rot * shape->convex->faces[f].normal;
    return Vec3{0, 1, 0};
  }

  [[nodiscard]] Scalar_t face_offset(std::uint32_t f) const noexcept {
    return glm::dot(face_normal(f), face_vertex(f, 0));
  }

  [[nodiscard]] std::uint32_t edge_count() const noexcept {
    switch (kind) {
    case PolytopeKind::Box:
      return 12;
    case PolytopeKind::Cylinder:
      return static_cast<std::uint32_t>(3 * kCylinderSides);
    case PolytopeKind::Convex:
      return shape->convex != nullptr ? static_cast<std::uint32_t>(shape->convex->edges.size()) : 0;
    }
    return 0;
  }

  [[nodiscard]] Vec3 edge_a(std::uint32_t e) const noexcept {
    if (kind == PolytopeKind::Box) {
      static constexpr std::uint8_t kBoxEdges[12][2] = {
          {0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3},
          {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7},
      };
      return vertex(kBoxEdges[e][0]);
    }
    if (kind == PolytopeKind::Cylinder) {
      const std::uint32_t n = static_cast<std::uint32_t>(kCylinderSides);
      if (e < n)
        return vertex(e);
      if (e < 2 * n)
        return vertex(n + (e - n));
      return vertex(e - 2 * n);
    }
    if (kind == PolytopeKind::Convex && shape->convex != nullptr)
      return vertex(shape->convex->edges[e].a);
    return shape->center;
  }

  [[nodiscard]] Vec3 edge_b(std::uint32_t e) const noexcept {
    if (kind == PolytopeKind::Box) {
      static constexpr std::uint8_t kBoxEdges[12][2] = {
          {0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3},
          {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7},
      };
      return vertex(kBoxEdges[e][1]);
    }
    if (kind == PolytopeKind::Cylinder) {
      const std::uint32_t n = static_cast<std::uint32_t>(kCylinderSides);
      if (e < n)
        return vertex((e + 1) % n);
      if (e < 2 * n)
        return vertex(n + ((e - n + 1) % n));
      return vertex(n + (e - 2 * n));
    }
    if (kind == PolytopeKind::Convex && shape->convex != nullptr)
      return vertex(shape->convex->edges[e].b);
    return shape->center;
  }

  [[nodiscard]] Vec3 edge_dir(std::uint32_t e) const noexcept {
    const Vec3     d   = edge_b(e) - edge_a(e);
    const Scalar_t len = glm::length(d);
    return len > kPolytopeEpsilon ? d / len : Vec3{1, 0, 0};
  }

  [[nodiscard]] Vec3 support(Vec3 dir) const noexcept {
    std::uint32_t       best   = 0;
    Scalar_t            best_d = -std::numeric_limits<Scalar_t>::infinity();
    const std::uint32_t n      = vertex_count();
    for (std::uint32_t i = 0; i < n; ++i) {
      const Scalar_t d = glm::dot(vertex(i), dir);
      if (d > best_d) {
        best_d = d;
        best   = i;
      }
    }
    return vertex(best);
  }

  // 在给定边方向上选“朝向 normal_ab（A->B）那一侧”的支撑边。
  void support_edge(Vec3 axis_dir, Vec3 normal_ab, Vec3 &a, Vec3 &b) const noexcept {
    const Vec3    axis       = glm::normalize(axis_dir);
    Scalar_t      best_score = -std::numeric_limits<Scalar_t>::infinity();
    std::uint32_t best       = kInvalidConvexMesh;
    for (std::uint32_t e = 0; e < edge_count(); ++e) {
      if (std::abs(glm::dot(edge_dir(e), axis)) < Scalar_t(1) - Scalar_t(1e-4))
        continue;
      const Vec3     pa    = edge_a(e);
      const Vec3     pb    = edge_b(e);
      const Scalar_t score = glm::dot((pa + pb) * Scalar_t(0.5), normal_ab);
      if (score > best_score) {
        best_score = score;
        best       = e;
      }
    }
    if (best == kInvalidConvexMesh) {
      a = support(normal_ab);
      b = a;
      return;
    }
    a = edge_a(best);
    b = edge_b(best);
  }
};

inline void polytope_interval(const PolytopeView &p, Vec3 dir, Scalar_t &lo,
                              Scalar_t &hi) noexcept {
  lo = std::numeric_limits<Scalar_t>::infinity();
  hi = -std::numeric_limits<Scalar_t>::infinity();
  for (std::uint32_t i = 0; i < p.vertex_count(); ++i) {
    const Scalar_t d = glm::dot(p.vertex(i), dir);
    lo               = std::min(lo, d);
    hi               = std::max(hi, d);
  }
}

enum class PolytopeAxisKind : std::uint8_t { FaceA = 0, FaceB = 1, Edge = 2 };

struct PolytopeSat {
  bool             separated = false;
  PolytopeAxisKind kind      = PolytopeAxisKind::FaceA;
  std::uint32_t    face_a    = 0;
  std::uint32_t    face_b    = 0;
  std::uint32_t    edge_a    = 0;
  std::uint32_t    edge_b    = 0;
  Vec3             edge_dir_a{1, 0, 0};
  Vec3             edge_dir_b{1, 0, 0};
  Vec3             normal_ab{0, 1, 0};
  Scalar_t         penetration = std::numeric_limits<Scalar_t>::infinity();
};

[[nodiscard]] bool sat_polytope(const PolytopeView &a, const PolytopeView &b,
                                PolytopeSat &result) noexcept {
  const Vec3 delta = b.center() - a.center();

  PolytopeSat best_face{};
  PolytopeSat best_edge{};
  bool        have_face = false;
  bool        have_edge = false;

  auto test_face = [&](bool is_a, std::uint32_t f) {
    const Vec3 n = is_a ? a.face_normal(f) : b.face_normal(f);
    if (is_a) {
      if (glm::dot(delta, n) < -kPolytopeEpsilon)
        return true; // B 在 A 的这个面的背面，它不是朝 B 的参考面
      Scalar_t min_a, max_a, min_b, max_b;
      polytope_interval(a, n, min_a, max_a);
      polytope_interval(b, n, min_b, max_b);
      const Scalar_t pen = max_a - min_b;
      if (pen < 0)
        return false;
      if (!have_face || pen < best_face.penetration) {
        best_face             = PolytopeSat{};
        best_face.kind        = PolytopeAxisKind::FaceA;
        best_face.face_a      = f;
        best_face.normal_ab   = n;
        best_face.penetration = pen;
        have_face             = true;
      }
    } else {
      if (glm::dot(-delta, n) < -kPolytopeEpsilon)
        return true; // A 在 B 的这个面的背面
      Scalar_t min_a, max_a, min_b, max_b;
      polytope_interval(a, n, min_a, max_a);
      polytope_interval(b, n, min_b, max_b);
      const Scalar_t pen = max_b - min_a;
      if (pen < 0)
        return false;
      if (!have_face || pen < best_face.penetration) {
        best_face             = PolytopeSat{};
        best_face.kind        = PolytopeAxisKind::FaceB;
        best_face.face_b      = f;
        best_face.normal_ab   = -n;
        best_face.penetration = pen;
        have_face             = true;
      }
    }
    return true;
  };

  for (std::uint32_t f = 0; f < a.face_count(); ++f)
    if (!test_face(true, f)) {
      result.separated = true;
      return false;
    }
  for (std::uint32_t f = 0; f < b.face_count(); ++f)
    if (!test_face(false, f)) {
      result.separated = true;
      return false;
    }

  // 唯一边方向（按边顺序 first-win，最多 kMaxSatEdgeDirs 个）。
  auto collect_dirs = [](const PolytopeView &p, Vec3 *dirs, std::uint32_t &count) {
    count = 0;
    for (std::uint32_t e = 0; e < p.edge_count() && count < kMaxSatEdgeDirs; ++e) {
      const Vec3 d   = p.edge_dir(e);
      bool       dup = false;
      for (std::uint32_t i = 0; i < count; ++i) {
        if (std::abs(glm::dot(dirs[i], d)) > Scalar_t(1) - Scalar_t(1e-5)) {
          dup = true;
          break;
        }
      }
      if (!dup)
        dirs[count++] = d;
    }
  };

  std::array<Vec3, kMaxSatEdgeDirs> da{}, db{};
  std::uint32_t                     na = 0, nb = 0;
  collect_dirs(a, da.data(), na);
  collect_dirs(b, db.data(), nb);

  for (std::uint32_t i = 0; i < na; ++i) {
    for (std::uint32_t j = 0; j < nb; ++j) {
      Vec3           axis = glm::cross(da[i], db[j]);
      const Scalar_t len  = glm::length(axis);
      if (!(len > kSatAxisEpsilon))
        continue;
      axis /= len;
      if (glm::dot(delta, axis) < 0)
        axis = -axis;

      Scalar_t min_a, max_a, min_b, max_b;
      polytope_interval(a, axis, min_a, max_a);
      polytope_interval(b, axis, min_b, max_b);
      const Scalar_t overlap = std::min(max_a, max_b) - std::max(min_a, min_b);
      if (overlap < 0) {
        result.separated = true;
        return false;
      }
      if (!have_edge || overlap < best_edge.penetration) {
        best_edge             = PolytopeSat{};
        best_edge.kind        = PolytopeAxisKind::Edge;
        best_edge.edge_a      = i;
        best_edge.edge_b      = j;
        best_edge.edge_dir_a  = da[i];
        best_edge.edge_dir_b  = db[j];
        best_edge.normal_ab   = axis;
        best_edge.penetration = overlap;
        have_edge             = true;
      }
    }
  }

  if (!have_face && !have_edge) {
    result.separated = true;
    return false;
  }
  if (!have_face) {
    result = best_edge;
    return true;
  }
  if (!have_edge) {
    result = best_face;
    return true;
  }

  // 和 box-box 一致：face 优先，edge 只有显著更浅时才切换。
  if (best_edge.penetration * Scalar_t(0.95) < best_face.penetration - Scalar_t(0.01))
    result = best_edge;
  else
    result = best_face;
  return true;
}

// ============================================================================
//  speculative SAT：允许轴上有正间隙，取"最小正间隙"的分离轴
// ----------------------------------------------------------------------------
//  与 sat_polytope 的区别只有两点：
//    1) 分离轴不再直接拒绝，而是记录 gap = min_b - max_a > 0 的最小者；
//    2) 任一轴 gap > max_gap 立即拒绝（投影间隙是真实距离的下界）。
//  面轴的朝向规则与 sat_polytope::test_face 完全一致：A 只测朝向 B 的面，
//  B 只测朝向 A 的面，normal_ab 统一指 A->B。这样生成器可以复用同一套
//  reference/incident face 约定，不会把"远面"当成参考面。
// ============================================================================
struct PolytopeGapSat {
  bool             found   = false;
  PolytopeAxisKind kind    = PolytopeAxisKind::FaceA;
  std::uint32_t    face_a  = 0;
  std::uint32_t    face_b  = 0;
  std::uint32_t    edge_a  = 0;
  std::uint32_t    edge_b  = 0;
  Vec3             edge_dir_a{1, 0, 0};
  Vec3             edge_dir_b{1, 0, 0};
  Vec3             normal_ab{0, 1, 0};
  Scalar_t         gap = std::numeric_limits<Scalar_t>::infinity();
};

[[nodiscard]] bool sat_polytope_speculative(const PolytopeView &a, const PolytopeView &b,
                                            Scalar_t max_gap, PolytopeGapSat &result) noexcept {
  const Vec3 delta = b.center() - a.center();

  auto test_face = [&](bool is_a, std::uint32_t f) -> bool {
    const Vec3 n = is_a ? a.face_normal(f) : b.face_normal(f);
    if (is_a) {
      if (glm::dot(delta, n) < -kPolytopeEpsilon)
        return true; // B 在 A 的这个面的背面：不是朝向 B 的参考面
    } else {
      if (glm::dot(-delta, n) < -kPolytopeEpsilon)
        return true; // A 在 B 的这个面的背面
    }
    const Vec3 axis = is_a ? n : -n; // A -> B
    Scalar_t   min_a, max_a, min_b, max_b;
    polytope_interval(a, axis, min_a, max_a);
    polytope_interval(b, axis, min_b, max_b);
    const Scalar_t gap = min_b - max_a;
    if (gap > max_gap)
      return false; // 真实距离 ≥ 该投影间隙
    if (gap > Scalar_t(0) && (!result.found || gap < result.gap)) {
      result          = PolytopeGapSat{};
      result.kind     = is_a ? PolytopeAxisKind::FaceA : PolytopeAxisKind::FaceB;
      result.face_a   = is_a ? f : 0;
      result.face_b   = is_a ? 0 : f;
      result.normal_ab = axis;
      result.gap      = gap;
      result.found    = true;
    }
    return true;
  };

  for (std::uint32_t f = 0; f < a.face_count(); ++f)
    if (!test_face(true, f))
      return false;
  for (std::uint32_t f = 0; f < b.face_count(); ++f)
    if (!test_face(false, f))
      return false;

  auto collect_dirs = [](const PolytopeView &p, Vec3 *dirs, std::uint32_t &count) {
    count = 0;
    for (std::uint32_t e = 0; e < p.edge_count() && count < kMaxSatEdgeDirs; ++e) {
      const Vec3 d   = p.edge_dir(e);
      bool       dup = false;
      for (std::uint32_t i = 0; i < count; ++i) {
        if (std::abs(glm::dot(dirs[i], d)) > Scalar_t(1) - Scalar_t(1e-5)) {
          dup = true;
          break;
        }
      }
      if (!dup)
        dirs[count++] = d;
    }
  };

  std::array<Vec3, kMaxSatEdgeDirs> da{}, db{};
  std::uint32_t                     na = 0, nb = 0;
  collect_dirs(a, da.data(), na);
  collect_dirs(b, db.data(), nb);

  for (std::uint32_t i = 0; i < na; ++i) {
    for (std::uint32_t j = 0; j < nb; ++j) {
      Vec3           axis = glm::cross(da[i], db[j]);
      const Scalar_t len  = glm::length(axis);
      if (!(len > kSatAxisEpsilon))
        continue;
      axis /= len;
      if (glm::dot(delta, axis) < 0)
        axis = -axis;

      Scalar_t min_a, max_a, min_b, max_b;
      polytope_interval(a, axis, min_a, max_a);
      polytope_interval(b, axis, min_b, max_b);
      const Scalar_t gap = min_b - max_a;
      if (gap > max_gap)
        return false;
      if (gap > Scalar_t(0) && (!result.found || gap < result.gap)) {
        result           = PolytopeGapSat{};
        result.kind      = PolytopeAxisKind::Edge;
        result.edge_a    = i;
        result.edge_b    = j;
        result.edge_dir_a = da[i];
        result.edge_dir_b = db[j];
        result.normal_ab  = axis;
        result.gap        = gap;
        result.found      = true;
      }
    }
  }
  return result.found;
}

[[nodiscard]] inline bool point_in_face(const PolytopeView &p, std::uint32_t f, Vec3 q,
                                        Scalar_t eps) noexcept {
  const Vec3          n     = p.face_normal(f);
  const std::uint32_t count = p.face_vertex_count(f);
  for (std::uint32_t i = 0; i < count; ++i) {
    const Vec3     a      = p.face_vertex(f, i);
    const Vec3     b      = p.face_vertex(f, (i + 1) % count);
    const Vec3     edge   = b - a;
    Vec3           inward = glm::cross(n, edge);
    const Scalar_t len    = glm::length(inward);
    if (!(len > kPolytopeEpsilon))
      continue;
    inward /= len;
    if (glm::dot(inward, q - a) < -eps)
      return false;
  }
  return true;
}

[[nodiscard]] inline bool point_in_polytope(const PolytopeView &p, Vec3 q, Vec3 &push,
                                            Scalar_t &depth) noexcept {
  Scalar_t best = std::numeric_limits<Scalar_t>::infinity();
  Vec3     best_n{0};
  for (std::uint32_t f = 0; f < p.face_count(); ++f) {
    const Vec3     n = p.face_normal(f);
    const Scalar_t d = p.face_offset(f) - glm::dot(n, q);
    if (d < -kPolytopeEpsilon)
      return false;
    if (d < best) {
      best   = d;
      best_n = n;
    }
  }
  if (!std::isfinite(best))
    return false;
  push  = best_n;
  depth = best;
  return true;
}

[[nodiscard]] inline Vec3 closest_point_on_segment(Vec3 p, Vec3 a, Vec3 b) noexcept {
  const Vec3     ab   = b - a;
  const Scalar_t len2 = glm::dot(ab, ab);
  if (len2 <= kPolytopeEpsilon * kPolytopeEpsilon)
    return a;
  const Scalar_t t = std::clamp(glm::dot(p - a, ab) / len2, Scalar_t(0), Scalar_t(1));
  return a + ab * t;
}

struct PolytopeClosest {
  Vec3     point{0};
  Scalar_t dist2 = std::numeric_limits<Scalar_t>::infinity();
};

[[nodiscard]] PolytopeClosest closest_point_on_polytope(const PolytopeView &p, Vec3 q) noexcept {
  PolytopeClosest best;

  auto consider = [&](Vec3 candidate) {
    const Vec3     d  = q - candidate;
    const Scalar_t d2 = glm::dot(d, d);
    if (d2 < best.dist2) {
      best.dist2 = d2;
      best.point = candidate;
    }
  };

  for (std::uint32_t f = 0; f < p.face_count(); ++f) {
    const Vec3 n    = p.face_normal(f);
    const Vec3 proj = q - n * (glm::dot(n, q) - p.face_offset(f));
    if (point_in_face(p, f, proj, kPolytopeEpsilon))
      consider(proj);
  }
  for (std::uint32_t e = 0; e < p.edge_count(); ++e)
    consider(closest_point_on_segment(q, p.edge_a(e), p.edge_b(e)));
  for (std::uint32_t v = 0; v < p.vertex_count(); ++v)
    consider(p.vertex(v));
  return best;
}

[[nodiscard]] std::size_t build_polytope_face_manifold(const ShapeView &sa, const ShapeView &sb,
                                                       const PolytopeView &pa,
                                                       const PolytopeView &pb, bool reference_is_a,
                                                       std::uint32_t reference_face, Vec3 normal_ab,
                                                       std::span<ContactPoint> out,
                                                       const BasisFrame       &basis) noexcept {
  const PolytopeView &reference        = reference_is_a ? pa : pb;
  const PolytopeView &incident         = reference_is_a ? pb : pa;
  const Vec3          reference_normal = reference_is_a ? normal_ab : -normal_ab;

  if (glm::dot(reference.face_normal(reference_face), reference_normal) < 0)
    return 0;

  // incident face：取与参考法线最反向的面。
  std::uint32_t incident_face = kInvalidConvexMesh;
  Scalar_t      best_dot      = std::numeric_limits<Scalar_t>::infinity();
  for (std::uint32_t f = 0; f < incident.face_count(); ++f) {
    const Scalar_t d = glm::dot(incident.face_normal(f), reference_normal);
    if (d < best_dot) {
      best_dot      = d;
      incident_face = f;
    }
  }
  if (incident_face == kInvalidConvexMesh || best_dot > kPolytopeEpsilon)
    return 0;

  const std::uint32_t inc_count = incident.face_vertex_count(incident_face);
  if (inc_count < 3 || inc_count > kMaxPolyVerts)
    return 0;

  std::array<Vec3, kMaxPolyVerts> poly_a{}, poly_b{};
  for (std::uint32_t i = 0; i < inc_count; ++i)
    poly_a[i] = incident.face_vertex(incident_face, i);

  Vec3 *in      = poly_a.data();
  Vec3 *out_buf = poly_b.data();
  int   count   = static_cast<int>(inc_count);

  const std::uint32_t ref_count = reference.face_vertex_count(reference_face);
  for (std::uint32_t i = 0; i < ref_count && count > 0; ++i) {
    const Vec3 a = reference.face_vertex(reference_face, i);
    const Vec3 b = reference.face_vertex(reference_face, (i + 1) % ref_count);
    // clip 的约定是“保留 dot(outward, p) <= offset”，所以这里必须用面外法线。
    Vec3           outward = glm::cross(b - a, reference_normal);
    const Scalar_t len     = glm::length(outward);
    if (!(len > kPolytopeEpsilon))
      continue;
    outward /= len;
    const Scalar_t offset = glm::dot(outward, a);
    count                 = clip_polygon_against_plane(in, count, outward, offset, out_buf);
    std::swap(in, out_buf);
  }
  std::size_t                               contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  const Vec3                                ref_point = reference.face_vertex(reference_face, 0);
  for (int i = 0; i < count && contact_count < kMaxContactsPerManifold; ++i) {
    const Vec3     p_incident = in[i];
    const Scalar_t distance   = glm::dot(p_incident - ref_point, reference_normal);
    if (distance > kPlaneEpsilon)
      continue;
    const Vec3          p_reference = p_incident - reference_normal * distance;
    const Vec3          xA          = reference_is_a ? p_reference : p_incident;
    const Vec3          xB          = reference_is_a ? p_incident : p_reference;
    const std::uint64_t key         = make_feature_key(1, reference_is_a ? 0 : 1, reference_face,
                                                       incident_face, static_cast<std::uint64_t>(i));
    add_contact(sa, sb, out, contact_count, midpoints.data(), xA, xB, key, basis.normal);
  }

  if (contact_count == 0) {
    const Vec3 xA = pa.support(normal_ab);
    const Vec3 xB = pb.support(-normal_ab);
    add_contact(sa, sb, out, contact_count, midpoints.data(), xA, xB,
                make_feature_key(2, reference_is_a ? 0 : 1, reference_face, incident_face),
                basis.normal);
  }
  return contact_count;
}

[[nodiscard]] std::size_t
build_polytope_edge_manifold(const ShapeView &sa, const ShapeView &sb, const PolytopeView &pa,
                             const PolytopeView &pb, const PolytopeSat &sat,
                             std::span<ContactPoint> out, const BasisFrame &basis) noexcept {
  Vec3 a0, a1, b0, b1;
  pa.support_edge(sat.edge_dir_a, sat.normal_ab, a0, a1);
  pb.support_edge(sat.edge_dir_b, -sat.normal_ab, b0, b1);
  Vec3 xA, xB;
  closest_points_on_segments(a0, a1, b0, b1, xA, xB);
  std::size_t                               contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact(sa, sb, out, contact_count, midpoints.data(), xA, xB,
              make_feature_key(3, sat.edge_a, sat.edge_b, 0), basis.normal);
  if (contact_count == 0) {
    const Vec3 sxA = pa.support(sat.normal_ab);
    const Vec3 sxB = pb.support(-sat.normal_ab);
    add_contact(sa, sb, out, contact_count, midpoints.data(), sxA, sxB,
                make_feature_key(2, sat.edge_a, sat.edge_b, 0), basis.normal);
  }
  return contact_count;
}

[[nodiscard]] std::size_t collide_polytope_polytope(const ShapeView &a, const ShapeView &b,
                                                    std::span<ContactPoint> out,
                                                    BasisFrame             &basis) noexcept {
  const PolytopeView pa(a);
  const PolytopeView pb(b);
  if (pa.vertex_count() == 0 || pb.vertex_count() == 0)
    return 0;

  PolytopeSat sat;
  if (!sat_polytope(pa, pb, sat))
    return 0;

  basis = make_basis(-sat.normal_ab);
  if (sat.kind == PolytopeAxisKind::Edge)
    return build_polytope_edge_manifold(a, b, pa, pb, sat, out, basis);
  return build_polytope_face_manifold(a, b, pa, pb, sat.kind == PolytopeAxisKind::FaceA,
                                      sat.kind == PolytopeAxisKind::FaceA ? sat.face_a : sat.face_b,
                                      sat.normal_ab, out, basis);
}

// speculative 面 manifold：与 build_polytope_face_manifold 的裁剪/投影过程逐行
// 相同，区别只有：入射点允许落在参考面外侧最多 max_gap 处，depth 记 -距离。
[[nodiscard]] std::size_t build_polytope_face_manifold_gap(
        const ShapeView &sa, const ShapeView &sb, const PolytopeView &pa,
        const PolytopeView &pb, bool reference_is_a, std::uint32_t reference_face,
        Vec3 normal_ab, Scalar_t max_gap, std::span<ContactPoint> out,
        const BasisFrame &basis) noexcept {
  const PolytopeView &reference        = reference_is_a ? pa : pb;
  const PolytopeView &incident         = reference_is_a ? pb : pa;
  const Vec3          reference_normal = reference_is_a ? normal_ab : -normal_ab;

  if (glm::dot(reference.face_normal(reference_face), reference_normal) < 0)
    return 0;

  std::uint32_t incident_face = kInvalidConvexMesh;
  Scalar_t      best_dot      = std::numeric_limits<Scalar_t>::infinity();
  for (std::uint32_t f = 0; f < incident.face_count(); ++f) {
    const Scalar_t d = glm::dot(incident.face_normal(f), reference_normal);
    if (d < best_dot) {
      best_dot      = d;
      incident_face = f;
    }
  }
  if (incident_face == kInvalidConvexMesh || best_dot > kPolytopeEpsilon)
    return 0;

  const std::uint32_t inc_count = incident.face_vertex_count(incident_face);
  if (inc_count < 3 || inc_count > kMaxPolyVerts)
    return 0;

  std::array<Vec3, kMaxPolyVerts> poly_a{}, poly_b{};
  for (std::uint32_t i = 0; i < inc_count; ++i)
    poly_a[i] = incident.face_vertex(incident_face, i);

  Vec3 *in      = poly_a.data();
  Vec3 *out_buf = poly_b.data();
  int   count   = static_cast<int>(inc_count);

  const std::uint32_t ref_count = reference.face_vertex_count(reference_face);
  for (std::uint32_t i = 0; i < ref_count && count > 0; ++i) {
    const Vec3 a = reference.face_vertex(reference_face, i);
    const Vec3 b = reference.face_vertex(reference_face, (i + 1) % ref_count);
    Vec3       outward = glm::cross(b - a, reference_normal);
    const Scalar_t len = glm::length(outward);
    if (!(len > kPolytopeEpsilon))
      continue;
    outward /= len;
    const Scalar_t offset = glm::dot(outward, a);
    count = clip_polygon_against_plane(in, count, outward, offset, out_buf);
    std::swap(in, out_buf);
  }

  std::size_t                               contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  const Vec3                                ref_point = reference.face_vertex(reference_face, 0);
  for (int i = 0; i < count && contact_count < kMaxContactsPerManifold; ++i) {
    const Vec3     p_incident = in[i];
    const Scalar_t distance   = glm::dot(p_incident - ref_point, reference_normal);
    if (distance > max_gap)
      continue;
    if (!(distance > kSpecGapEpsilon))
      continue; // exact 路径已返回 0：这里只生成严格正的间隙
    const Vec3          p_reference = p_incident - reference_normal * distance;
    const Vec3          xA          = reference_is_a ? p_reference : p_incident;
    const Vec3          xB          = reference_is_a ? p_incident : p_reference;
    const std::uint64_t key         = make_feature_key(1, reference_is_a ? 0 : 1, reference_face,
                                                       incident_face, static_cast<std::uint64_t>(i));
    add_contact_gap(sa, sb, out, contact_count, midpoints.data(), xA, xB, key, distance);
  }
  return contact_count;
}

[[nodiscard]] std::size_t collide_polytope_polytope_speculative(const ShapeView &a,
                                                                const ShapeView &b,
                                                                std::span<ContactPoint> out,
                                                                BasisFrame &basis,
                                                                Scalar_t max_gap) noexcept {
  const PolytopeView pa(a);
  const PolytopeView pb(b);
  if (pa.vertex_count() == 0 || pb.vertex_count() == 0)
    return 0;

  PolytopeGapSat sat;
  if (!sat_polytope_speculative(pa, pb, max_gap, sat))
    return 0;

  basis = make_basis(-sat.normal_ab);
  if (sat.kind == PolytopeAxisKind::Edge) {
    Vec3 a0, a1, b0, b1;
    pa.support_edge(sat.edge_dir_a, sat.normal_ab, a0, a1);
    pb.support_edge(sat.edge_dir_b, -sat.normal_ab, b0, b1);
    Vec3 xA, xB;
    closest_points_on_segments(a0, a1, b0, b1, xA, xB);
    const Scalar_t gap = glm::dot(basis.normal, xA - xB);
    if (!(gap > kSpecGapEpsilon) || gap > max_gap)
      return 0;
    std::size_t                               contact_count = 0;
    std::array<Vec3, kMaxContactsPerManifold> midpoints{};
    add_contact_gap(a, b, out, contact_count, midpoints.data(), xA, xB,
                    make_feature_key(3, sat.edge_a, sat.edge_b, 0), gap);
    return contact_count;
  }

  const bool          reference_is_a = sat.kind == PolytopeAxisKind::FaceA;
  const std::uint32_t reference_face = reference_is_a ? sat.face_a : sat.face_b;
  const std::size_t   n =
          build_polytope_face_manifold_gap(a, b, pa, pb, reference_is_a, reference_face,
                                           sat.normal_ab, max_gap, out, basis);
  if (n > 0)
    return n;

  // 裁剪退化（例如旋转后投影几乎平行）：退回沿分离轴的支撑点。
  // 支撑点在轴向上的距离正好是 sat.gap，仍是合法的最接近特征。
  const Vec3 xA = pa.support(sat.normal_ab);
  const Vec3 xB = pb.support(-sat.normal_ab);
  const Scalar_t gap = glm::dot(basis.normal, xA - xB);
  if (!(gap > kSpecGapEpsilon) || gap > max_gap)
    return 0;
  std::size_t                               contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact_gap(a, b, out, contact_count, midpoints.data(), xA, xB,
                  make_feature_key(2, reference_is_a ? 0 : 1, reference_face, 0), gap);
  return contact_count;
}

[[nodiscard]] std::size_t collide_sphere_polytope(const ShapeView &a, const ShapeView &b,
                                                  bool sphere_is_a, std::span<ContactPoint> out,
                                                  BasisFrame &basis) noexcept {
  const ShapeView   &sphere = sphere_is_a ? a : b;
  const ShapeView   &other  = sphere_is_a ? b : a;
  const PolytopeView poly(other);
  if (poly.vertex_count() == 0)
    return 0;

  Vec3     normal_other_to_sphere{0, 1, 0};
  Vec3     closest{sphere.center};
  Vec3     push;
  Scalar_t face_depth = 0;

  if (point_in_polytope(poly, sphere.center, push, face_depth)) {
    // 球心在 polytope 内：按最近面推出；接触点沿用原 sphere-box 的约定。
    normal_other_to_sphere = push;
    closest                = sphere.center;
  } else {
    const PolytopeClosest c  = closest_point_on_polytope(poly, sphere.center);
    const Scalar_t        r2 = sphere.radius * sphere.radius;
    if (!(c.dist2 < r2))
      return 0;
    const Scalar_t dist = std::sqrt(std::max(c.dist2, Scalar_t(0)));
    if (!(dist > kPolytopeEpsilon))
      return 0;
    normal_other_to_sphere = (sphere.center - c.point) / dist;
    closest                = c.point;
  }

  const Vec3 normal_ab = sphere_is_a ? -normal_other_to_sphere : normal_other_to_sphere;
  basis                = make_basis(-normal_ab);

  const Vec3  on_sphere     = sphere.center - normal_other_to_sphere * sphere.radius;
  const Vec3  xA            = sphere_is_a ? on_sphere : closest;
  const Vec3  xB            = sphere_is_a ? closest : on_sphere;
  std::size_t contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact(a, b, out, contact_count, midpoints.data(), xA, xB,
              make_feature_key(4, sphere_is_a ? 0 : 1, 0, 0), basis.normal);
  return contact_count;
}

[[nodiscard]] std::size_t collide_sphere_polytope_speculative(const ShapeView &a,
                                                              const ShapeView &b,
                                                              bool sphere_is_a,
                                                              std::span<ContactPoint> out,
                                                              BasisFrame &basis,
                                                              Scalar_t max_gap) noexcept {
  const ShapeView   &sphere = sphere_is_a ? a : b;
  const ShapeView   &other  = sphere_is_a ? b : a;
  const PolytopeView poly(other);
  if (poly.vertex_count() == 0)
    return 0;

  Vec3     push;
  Scalar_t face_depth = 0;
  if (point_in_polytope(poly, sphere.center, push, face_depth))
    return 0; // 球心在内：exact 路径已经处理

  const PolytopeClosest c    = closest_point_on_polytope(poly, sphere.center);
  const Scalar_t        dist = std::sqrt(std::max(c.dist2, Scalar_t(0)));
  if (!(dist > kPolytopeEpsilon))
    return 0;
  const Scalar_t gap = dist - sphere.radius;
  if (!(gap > kSpecGapEpsilon) || gap > max_gap)
    return 0;

  const Vec3 normal_other_to_sphere = (sphere.center - c.point) / dist;
  const Vec3 normal_ab = sphere_is_a ? -normal_other_to_sphere : normal_other_to_sphere;
  basis                = make_basis(-normal_ab);

  const Vec3  on_sphere     = sphere.center - normal_other_to_sphere * sphere.radius;
  const Vec3  xA            = sphere_is_a ? on_sphere : c.point;
  const Vec3  xB            = sphere_is_a ? c.point : on_sphere;
  std::size_t contact_count = 0;
  std::array<Vec3, kMaxContactsPerManifold> midpoints{};
  add_contact_gap(a, b, out, contact_count, midpoints.data(), xA, xB,
                  make_feature_key(4, sphere_is_a ? 0 : 1, 0, 0), gap);
  return contact_count;
}

} // namespace

BasisFrame make_basis(Vec3 normal_b_to_a) noexcept {
  normal_b_to_a = glm::normalize(normal_b_to_a);
  BasisFrame f;
  f.normal   = normal_b_to_a;
  Vec3 t0    = std::abs(normal_b_to_a.x) > std::abs(normal_b_to_a.z)
                   ? Vec3{-normal_b_to_a.y, normal_b_to_a.x, 0}
                   : Vec3{0, -normal_b_to_a.z, normal_b_to_a.y};
  t0         = glm::normalize(t0);
  f.tangent0 = t0;
  f.tangent1 = glm::cross(normal_b_to_a, t0);
  return f;
}

namespace {
std::size_t collide_shapes_exact(const ShapeView &a, const ShapeView &b,
                                 std::span<ContactPoint> out, BasisFrame &basis) {
  if (a.type == ShapeType::Sphere && b.type == ShapeType::Sphere)
    return collide_sphere_sphere(a, b, out, basis);
  if (a.type == ShapeType::Sphere || b.type == ShapeType::Sphere)
    return collide_sphere_polytope(a, b, a.type == ShapeType::Sphere, out, basis);
  return collide_polytope_polytope(a, b, out, basis);
}

std::size_t collide_shapes_speculative(const ShapeView &a, const ShapeView &b,
                                       std::span<ContactPoint> out, BasisFrame &basis,
                                       Scalar_t speculative_distance) {
  if (a.type == ShapeType::Sphere && b.type == ShapeType::Sphere)
    return collide_sphere_sphere_speculative(a, b, out, basis, speculative_distance);
  if (a.type == ShapeType::Sphere || b.type == ShapeType::Sphere)
    return collide_sphere_polytope_speculative(a, b, a.type == ShapeType::Sphere, out, basis,
                                               speculative_distance);
  return collide_polytope_polytope_speculative(a, b, out, basis, speculative_distance);
}
} // namespace

bool sweep_shapes_toi(const ShapeView &a, const ShapeView &b, Vec3 relative_motion,
                      Scalar_t &t) {
  const Scalar_t motion_len = glm::length(relative_motion);
  if (!(motion_len > kSatAxisEpsilon))
    return false;
  const Scalar_t max_gap = motion_len + Scalar_t(1e-6);

  if (a.type == ShapeType::Sphere && b.type == ShapeType::Sphere) {
    const Vec3     d   = b.center - a.center;
    const Scalar_t d2  = glm::dot(d, d);
    const Scalar_t sum = a.radius + b.radius;
    const Scalar_t dist = std::sqrt(d2);
    if (!(dist > kSatAxisEpsilon))
      return false;
    const Scalar_t gap = dist - sum;
    if (!(gap > kSpecGapEpsilon) || gap > max_gap)
      return false;
    const Vec3     normal_ab = d / dist;
    const Scalar_t closing   = glm::dot(relative_motion, normal_ab);
    if (!(closing > Scalar_t(0)))
      return false;
    t = gap / closing;
    return t > Scalar_t(0) && t <= Scalar_t(1);
  }

  if (a.type == ShapeType::Sphere || b.type == ShapeType::Sphere) {
    const bool        sphere_is_a = a.type == ShapeType::Sphere;
    const ShapeView  &sphere      = sphere_is_a ? a : b;
    const ShapeView  &other       = sphere_is_a ? b : a;
    const PolytopeView poly(other);
    if (poly.vertex_count() == 0)
      return false;
    Vec3     push;
    Scalar_t face_depth = 0;
    if (point_in_polytope(poly, sphere.center, push, face_depth))
      return false;
    const PolytopeClosest c    = closest_point_on_polytope(poly, sphere.center);
    const Scalar_t        dist = std::sqrt(std::max(c.dist2, Scalar_t(0)));
    if (!(dist > kPolytopeEpsilon))
      return false;
    const Scalar_t gap = dist - sphere.radius;
    if (!(gap > kSpecGapEpsilon) || gap > max_gap)
      return false;
    const Vec3 normal_other_to_sphere = (sphere.center - c.point) / dist;
    const Vec3 normal_ab = sphere_is_a ? -normal_other_to_sphere : normal_other_to_sphere;
    const Scalar_t closing = glm::dot(relative_motion, normal_ab);
    if (!(closing > Scalar_t(0)))
      return false;
    t = gap / closing;
    return t > Scalar_t(0) && t <= Scalar_t(1);
  }

  const PolytopeView pa(a);
  const PolytopeView pb(b);
  if (pa.vertex_count() == 0 || pb.vertex_count() == 0)
    return false;
  PolytopeGapSat sat;
  if (!sat_polytope_speculative(pa, pb, max_gap, sat))
    return false;
  const Scalar_t closing = glm::dot(relative_motion, sat.normal_ab);
  if (!(closing > Scalar_t(0)))
    return false;
  t = sat.gap / closing;
  return t > Scalar_t(0) && t <= Scalar_t(1);
}

std::size_t collide_shapes(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                           BasisFrame &basis) {
  return collide_shapes_exact(a, b, out, basis);
}

std::size_t collide_shapes(const ShapeView &a, const ShapeView &b, std::span<ContactPoint> out,
                           BasisFrame &basis, Scalar_t speculative_distance) {
  // 重叠路径逐位不变：只有 exact 明确返回 0、且允许预接触时，才走 speculative。
  const std::size_t n = collide_shapes_exact(a, b, out, basis);
  if (n > 0 || !(speculative_distance > Scalar_t(0)))
    return n;
  return collide_shapes_speculative(a, b, out, basis, speculative_distance);
}

} // namespace collision
