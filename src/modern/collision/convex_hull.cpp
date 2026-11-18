#include "collision/convex_hull.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace collision {
namespace {

struct BuildFace {
  std::uint32_t a = 0, b = 0, c = 0;
  Vec3          normal{0, 1, 0};
  Scalar_t      offset = 0;
};

[[nodiscard]] bool finite_vec(Vec3 v) noexcept {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

[[nodiscard]] bool make_face(std::uint32_t a, std::uint32_t b, std::uint32_t c, Vec3 interior,
                             const std::vector<Vec3> &pts, Scalar_t eps, BuildFace &out) {
  Vec3           n   = glm::cross(pts[b] - pts[a], pts[c] - pts[a]);
  const Scalar_t len = glm::length(n);
  if (!(len > eps))
    return false;
  n /= len;
  if (glm::dot(n, interior - pts[a]) > 0) {
    std::swap(b, c);
    n = -n;
  }
  out.a      = a;
  out.b      = b;
  out.c      = c;
  out.normal = n;
  out.offset = glm::dot(n, pts[a]);
  return true;
}

} // namespace

HullStatus build_convex_hull(std::span<const Vec3> points, const HullOptions &opt,
                             ConvexMesh &out) {
  if (points.size() < 4)
    return HullStatus::TooFewPoints;

  const Scalar_t eps  = std::max(opt.epsilon, Scalar_t(1e-7));
  const Scalar_t eps2 = eps * eps;

  std::vector<Vec3> pts(points.begin(), points.end());
  for (Vec3 p : pts)
    if (!finite_vec(p))
      return HullStatus::NonFiniteInput;

  std::sort(pts.begin(), pts.end(), [](Vec3 a, Vec3 b) {
    if (a.x != b.x)
      return a.x < b.x;
    if (a.y != b.y)
      return a.y < b.y;
    return a.z < b.z;
  });

  std::vector<Vec3> unique;
  unique.reserve(pts.size());
  for (Vec3 p : pts) {
    bool duplicate = false;
    for (Vec3 q : unique) {
      const Vec3 d = p - q;
      if (glm::dot(d, d) <= eps2) {
        duplicate = true;
        break;
      }
    }
    if (!duplicate)
      unique.push_back(p);
  }
  pts.swap(unique);
  if (pts.size() < 4)
    return HullStatus::TooFewPoints;
  if (pts.size() > opt.max_vertices)
    return HullStatus::TooManyVertices;

  // —— 确定性初始四面体：先找最长轴上的两个极值点，再找最远线外点/面外点 ——
  std::size_t p0 = 0, p1 = 0;
  Scalar_t    best_d2 = -1;
  for (int axis = 0; axis < 3; ++axis) {
    std::size_t lo = 0, hi = 0;
    for (std::size_t i = 1; i < pts.size(); ++i) {
      if (pts[i][axis] < pts[lo][axis])
        lo = i;
      if (pts[i][axis] > pts[hi][axis])
        hi = i;
    }
    const Vec3     d  = pts[hi] - pts[lo];
    const Scalar_t d2 = glm::dot(d, d);
    if (d2 > best_d2) {
      best_d2 = d2;
      p0      = lo;
      p1      = hi;
    }
  }
  if (!(best_d2 > eps2))
    return HullStatus::Degenerate;

  Vec3     line_dir = pts[p1] - pts[p0];
  Scalar_t line_len = glm::length(line_dir);
  line_dir /= line_len;

  std::size_t p2        = kInvalidConvexMesh;
  Scalar_t    best_area = -1;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Vec3     rel = pts[i] - pts[p0];
    const Vec3     c   = glm::cross(line_dir, rel);
    const Scalar_t d2  = glm::dot(c, c);
    if (d2 > best_area) {
      best_area = d2;
      p2        = i;
    }
  }
  if (p2 == kInvalidConvexMesh || !(best_area > eps2))
    return HullStatus::Degenerate;

  Vec3     plane_n   = glm::cross(pts[p1] - pts[p0], pts[p2] - pts[p0]);
  Scalar_t plane_len = glm::length(plane_n);
  plane_n /= plane_len;

  std::size_t p3       = kInvalidConvexMesh;
  Scalar_t    best_off = -1;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Scalar_t d = std::abs(glm::dot(plane_n, pts[i] - pts[p0]));
    if (d > best_off) {
      best_off = d;
      p3       = i;
    }
  }
  if (p3 == kInvalidConvexMesh || !(best_off > eps))
    return HullStatus::Degenerate;

  const Vec3 interior = (pts[p0] + pts[p1] + pts[p2] + pts[p3]) * Scalar_t(0.25);

  std::vector<BuildFace> faces;
  faces.reserve(64);
  auto add_face = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
    BuildFace f;
    if (make_face(a, b, c, interior, pts, eps, f))
      faces.push_back(f);
  };

  add_face(static_cast<std::uint32_t>(p0), static_cast<std::uint32_t>(p1),
           static_cast<std::uint32_t>(p2));
  add_face(static_cast<std::uint32_t>(p0), static_cast<std::uint32_t>(p3),
           static_cast<std::uint32_t>(p1));
  add_face(static_cast<std::uint32_t>(p1), static_cast<std::uint32_t>(p3),
           static_cast<std::uint32_t>(p2));
  add_face(static_cast<std::uint32_t>(p2), static_cast<std::uint32_t>(p3),
           static_cast<std::uint32_t>(p0));
  if (faces.size() < 4)
    return HullStatus::Degenerate;

  std::vector<bool> in_tetra(pts.size(), false);
  in_tetra[p0] = in_tetra[p1] = in_tetra[p2] = in_tetra[p3] = true;

  // —— 逐点增量：可见面 → horizon 边 → 新面 ——
  for (std::uint32_t q = 0; q < pts.size(); ++q) {
    if (in_tetra[q])
      continue;

    std::vector<std::uint32_t> visible;
    visible.reserve(8);
    for (std::uint32_t f = 0; f < faces.size(); ++f) {
      if (glm::dot(faces[f].normal, pts[q]) - faces[f].offset > eps)
        visible.push_back(f);
    }
    if (visible.empty())
      continue; // q 在凸包内部或面上

    std::vector<std::pair<std::uint32_t, std::uint32_t>> horizon;
    for (std::uint32_t fi : visible) {
      const BuildFace    &f           = faces[fi];
      const std::uint32_t edges[3][2] = {{f.a, f.b}, {f.b, f.c}, {f.c, f.a}};
      for (const auto &edge : edges) {
        const std::pair<std::uint32_t, std::uint32_t> rev{edge[1], edge[0]};
        const auto it = std::find(horizon.begin(), horizon.end(), rev);
        if (it != horizon.end())
          horizon.erase(it);
        else
          horizon.emplace_back(edge[0], edge[1]);
      }
    }

    std::vector<BuildFace> kept;
    kept.reserve(faces.size());
    for (std::uint32_t f = 0; f < faces.size(); ++f) {
      if (std::find(visible.begin(), visible.end(), f) == visible.end())
        kept.push_back(faces[f]);
    }

    for (const auto &[a, b] : horizon) {
      BuildFace f;
      if (make_face(a, b, q, interior, pts, eps, f))
        kept.push_back(f);
    }
    if (kept.size() < 4)
      return HullStatus::Degenerate;
    if (kept.size() > opt.max_faces)
      return HullStatus::TooManyFaces;
    faces.swap(kept);
  }

  // —— 压缩到被 face 引用的顶点，输出三角面；finalize 会再做闭合/凸性校验 ——
  std::vector<std::uint32_t> remap(pts.size(), kInvalidConvexMesh);
  std::vector<Vec3>          out_vertices;
  std::vector<std::uint32_t> out_indices;
  std::vector<ConvexFace>    out_faces;
  out_vertices.reserve(pts.size());
  out_indices.reserve(faces.size() * 3);
  out_faces.reserve(faces.size());

  for (const BuildFace &f : faces) {
    ConvexFace face;
    face.first                = static_cast<std::uint32_t>(out_indices.size());
    face.count                = 3;
    const std::uint32_t vs[3] = {f.a, f.b, f.c};
    for (std::uint32_t v : vs) {
      if (remap[v] == kInvalidConvexMesh) {
        remap[v] = static_cast<std::uint32_t>(out_vertices.size());
        out_vertices.push_back(pts[v]);
      }
      out_indices.push_back(remap[v]);
    }
    out_faces.push_back(face);
  }

  ConvexMesh mesh;
  mesh.vertices     = std::move(out_vertices);
  mesh.face_indices = std::move(out_indices);
  mesh.faces        = std::move(out_faces);

  const ConvexMeshStatus status = finalize_convex_mesh(mesh, eps);
  if (status != ConvexMeshStatus::Ok)
    return status == ConvexMeshStatus::TooLarge ? HullStatus::TooManyFaces : HullStatus::Degenerate;

  out = std::move(mesh);
  return HullStatus::Ok;
}

} // namespace collision
