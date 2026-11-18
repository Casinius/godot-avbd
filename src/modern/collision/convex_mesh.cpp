#include "collision/convex_mesh.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace collision {
namespace {

[[nodiscard]] bool finite_vec(Vec3 v) noexcept {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

[[nodiscard]] Vec3 newell_normal(const std::vector<Vec3>          &vertices,
                                 const std::vector<std::uint32_t> &loop) noexcept {
  Vec3 n{0};
  for (std::size_t i = 0; i < loop.size(); ++i) {
    const Vec3 a = vertices[loop[i]];
    const Vec3 b = vertices[loop[(i + 1) % loop.size()]];
    n.x += (a.y - b.y) * (a.z + b.z);
    n.y += (a.z - b.z) * (a.x + b.x);
    n.z += (a.x - b.x) * (a.y + b.y);
  }
  return n;
}

} // namespace

ConvexMeshStatus finalize_convex_mesh(ConvexMesh &mesh, Scalar_t epsilon) {
  if (mesh.vertices.size() < 4)
    return ConvexMeshStatus::TooFewVertices;
  if (mesh.vertices.size() > kMaxConvexVertices)
    return ConvexMeshStatus::TooLarge;
  if (mesh.faces.size() < 4)
    return ConvexMeshStatus::TooFewFaces;
  if (mesh.faces.size() > kMaxConvexFaces)
    return ConvexMeshStatus::TooLarge;

  const Scalar_t eps  = std::max(epsilon, Scalar_t(1e-7));
  const Scalar_t eps2 = eps * eps;

  // —— 1) 去重顶点（保持首次出现顺序，确定性）——
  std::vector<Vec3>          unique;
  std::vector<std::uint32_t> remap(mesh.vertices.size(), 0);
  unique.reserve(mesh.vertices.size());
  for (std::size_t i = 0; i < mesh.vertices.size(); ++i) {
    const Vec3 v = mesh.vertices[i];
    if (!finite_vec(v))
      return ConvexMeshStatus::NonFinite;

    std::uint32_t found = kInvalidConvexMesh;
    for (std::uint32_t j = 0; j < unique.size(); ++j) {
      const Vec3 d = unique[j] - v;
      if (glm::dot(d, d) <= eps2) {
        found = j;
        break;
      }
    }
    if (found == kInvalidConvexMesh) {
      found = static_cast<std::uint32_t>(unique.size());
      unique.push_back(v);
    }
    remap[i] = found;
  }
  if (unique.size() < 4)
    return ConvexMeshStatus::TooFewVertices;
  if (unique.size() > kMaxConvexVertices)
    return ConvexMeshStatus::TooLarge;

  Vec3 centroid{0};
  for (Vec3 v : unique)
    centroid += v;
  centroid /= Scalar_t(unique.size());

  // —— 2) 规范化每个 face 环：remap → 去连续重复 → 外向 → rotate 到最小索引 ——
  std::vector<std::uint32_t> face_indices;
  std::vector<ConvexFace>    faces;
  face_indices.reserve(mesh.face_indices.size());
  faces.reserve(mesh.faces.size());

  for (const ConvexFace &src_face : mesh.faces) {
    if (src_face.first > mesh.face_indices.size() ||
        src_face.count > mesh.face_indices.size() - src_face.first)
      return ConvexMeshStatus::BadFaceIndex;

    std::vector<std::uint32_t> loop;
    loop.reserve(src_face.count);
    for (std::uint32_t i = 0; i < src_face.count; ++i) {
      const std::uint32_t raw = mesh.face_indices[src_face.first + i];
      if (raw >= remap.size())
        return ConvexMeshStatus::BadFaceIndex;
      const std::uint32_t v = remap[raw];
      if (loop.empty() || loop.back() != v)
        loop.push_back(v);
    }
    if (loop.size() > 1 && loop.front() == loop.back())
      loop.pop_back();
    if (loop.size() < 3)
      return ConvexMeshStatus::DegenerateFace;

    Vec3           normal = newell_normal(unique, loop);
    const Scalar_t len    = glm::length(normal);
    if (!(len > eps))
      return ConvexMeshStatus::DegenerateFace;
    normal /= len;

    Vec3 face_center{0};
    for (std::uint32_t v : loop)
      face_center += unique[v];
    face_center /= Scalar_t(loop.size());

    if (glm::dot(normal, face_center - centroid) < 0) {
      std::reverse(loop.begin(), loop.end());
      normal = -normal;
    }

    const auto min_it = std::min_element(loop.begin(), loop.end());
    std::rotate(loop.begin(), min_it, loop.end());

    ConvexFace face;
    face.first  = static_cast<std::uint32_t>(face_indices.size());
    face.count  = static_cast<std::uint32_t>(loop.size());
    face.normal = normal;
    face.offset = glm::dot(normal, unique[loop.front()]);
    faces.push_back(face);
    face_indices.insert(face_indices.end(), loop.begin(), loop.end());
  }

  if (faces.size() < 4)
    return ConvexMeshStatus::TooFewFaces;
  if (faces.size() > kMaxConvexFaces)
    return ConvexMeshStatus::TooLarge;

  // —— 3) 凸性校验：所有顶点必须在该 face 的背面（外向法线）——
  for (const ConvexFace &face : faces) {
    for (Vec3 v : unique) {
      const Scalar_t d = glm::dot(face.normal, v) - face.offset;
      if (d > eps * (Scalar_t(1) + std::abs(face.offset)))
        return ConvexMeshStatus::NonConvex;
    }
  }

  // —— 4) 唯一边 + 闭合校验（每条无向边恰好出现两次）——
  struct EdgeKey {
    std::uint32_t a     = 0;
    std::uint32_t b     = 0;
    std::uint32_t count = 0;
    std::uint32_t index = 0;
  };
  std::vector<EdgeKey> edge_keys;
  edge_keys.reserve(faces.size() * 4);
  for (const ConvexFace &face : faces) {
    for (std::uint32_t i = 0; i < face.count; ++i) {
      const std::uint32_t a = face_indices[face.first + i];
      const std::uint32_t b = face_indices[face.first + (i + 1) % face.count];
      if (a == b)
        return ConvexMeshStatus::DegenerateFace;
      const std::uint32_t lo    = std::min(a, b);
      const std::uint32_t hi    = std::max(a, b);
      bool                found = false;
      for (EdgeKey &e : edge_keys) {
        if (e.a == lo && e.b == hi) {
          ++e.count;
          found = true;
          break;
        }
      }
      if (!found)
        edge_keys.push_back(EdgeKey{lo, hi, 1, static_cast<std::uint32_t>(edge_keys.size())});
    }
  }
  if (edge_keys.size() > kMaxConvexEdges)
    return ConvexMeshStatus::TooLarge;
  for (const EdgeKey &e : edge_keys)
    if (e.count != 2)
      return ConvexMeshStatus::NonClosed;

  std::vector<ConvexEdge> edges;
  edges.reserve(edge_keys.size());
  for (const EdgeKey &e : edge_keys) {
    const Vec3     d   = unique[e.b] - unique[e.a];
    const Scalar_t len = glm::length(d);
    if (!(len > eps))
      return ConvexMeshStatus::DegenerateFace;
    edges.push_back(ConvexEdge{e.a, e.b, d / len});
  }

  Vec3     local_half{0};
  Scalar_t bounding_radius = 0;
  for (Vec3 v : unique) {
    local_half.x    = std::max(local_half.x, std::abs(v.x));
    local_half.y    = std::max(local_half.y, std::abs(v.y));
    local_half.z    = std::max(local_half.z, std::abs(v.z));
    bounding_radius = std::max(bounding_radius, glm::length(v));
  }

  mesh.vertices        = std::move(unique);
  mesh.face_indices    = std::move(face_indices);
  mesh.faces           = std::move(faces);
  mesh.edges           = std::move(edges);
  mesh.local_half      = local_half;
  mesh.centroid        = centroid;
  mesh.bounding_radius = bounding_radius;
  return ConvexMeshStatus::Ok;
}

ConvexMeshId ConvexMeshStore::add(ConvexMesh mesh) {
  if (finalize_convex_mesh(mesh) != ConvexMeshStatus::Ok)
    return kInvalidConvexMesh;
  if (meshes_.size() >= kInvalidConvexMesh - 1)
    return kInvalidConvexMesh;
  meshes_.push_back(std::move(mesh));
  return static_cast<ConvexMeshId>(meshes_.size() - 1);
}

bool ConvexMeshStore::replace(ConvexMeshId id, ConvexMesh mesh) {
  if (id >= meshes_.size())
    return false;
  if (finalize_convex_mesh(mesh) != ConvexMeshStatus::Ok)
    return false;
  meshes_[id] = std::move(mesh);
  return true;
}

const ConvexMesh *ConvexMeshStore::get(ConvexMeshId id) const noexcept {
  return id < meshes_.size() ? &meshes_[id] : nullptr;
}

} // namespace collision
