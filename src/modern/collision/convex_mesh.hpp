#pragma once

// ============================================================================
//  collision/convex_mesh.hpp —— 客户端预构建凸包 + 注册表
// ----------------------------------------------------------------------------
//  【输入契约】
//  客户端提供的是**已构建好的凸包**：vertices + 每个 face 的顶点索引环。
//  环必须按外向 CCW 排列；这里会做一次规范化（去重顶点、去退化面、统一朝向、
//  面环 rotate 到最小索引），但不会尝试把非凸网格“修成”凸包。
//
//  【不变量】
//   1) ConvexMesh 由 CollisionWorld 的 ConvexMeshStore 拥有；Shape 只存 id；
//   2) vertices/faces/edges 在注册后不可变，避免 ShapeView 指针失效；
//   3) face.normal 预归一化，face.offset = dot(normal, 面上任一点)；
//   4) edges 是唯一边，dir 预归一化；SAT 用 edges 的方向，manifold 用端点；
//   5) 所有校验在 setup-time 完成，失败返回状态码；不使用异常。
// ============================================================================

#include <cstddef>
#include <cstdint>
#include <vector>

#include "phys_core/state.hpp"

namespace collision {

inline constexpr std::uint32_t kInvalidConvexMesh = 0xFFFFFFFFu;
inline constexpr std::size_t   kMaxConvexVertices = 256;
inline constexpr std::size_t   kMaxConvexFaces    = 512;
inline constexpr std::size_t   kMaxConvexEdges    = 1536;

using ConvexMeshId = std::uint32_t;

enum class ConvexMeshStatus : std::uint8_t {
  Ok = 0,
  TooFewVertices,
  TooFewFaces,
  NonFinite,
  BadFaceIndex,
  DegenerateFace,
  NonConvex,
  NonClosed,
  TooLarge,
};

struct ConvexFace {
  std::uint32_t first = 0; // face_indices 起点
  std::uint32_t count = 0; // 顶点数
  Vec3          normal{0, 1, 0};
  Scalar_t      offset = 0; // dot(normal, face 上一点)
};

struct ConvexEdge {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
  Vec3          dir{1, 0, 0}; // normalize(vertices[b] - vertices[a])
};

struct ConvexMesh {
  std::vector<Vec3>          vertices;
  std::vector<std::uint32_t> face_indices; // 每个 face 的 CCW 顶点环
  std::vector<ConvexFace>    faces;
  std::vector<ConvexEdge>    edges;
  Vec3                       local_half{0};
  Vec3                       centroid{0};
  Scalar_t                   bounding_radius = 0;
};

// 校验并规范化 mesh。成功时 mesh 的 faces/edges/local_half/centroid 被填好。
[[nodiscard]] ConvexMeshStatus finalize_convex_mesh(ConvexMesh &mesh,
                                                    Scalar_t    epsilon = Scalar_t(1e-5));

[[nodiscard]] inline const char *to_string(ConvexMeshStatus s) noexcept {
  switch (s) {
  case ConvexMeshStatus::Ok:
    return "Ok";
  case ConvexMeshStatus::TooFewVertices:
    return "TooFewVertices";
  case ConvexMeshStatus::TooFewFaces:
    return "TooFewFaces";
  case ConvexMeshStatus::NonFinite:
    return "NonFinite";
  case ConvexMeshStatus::BadFaceIndex:
    return "BadFaceIndex";
  case ConvexMeshStatus::DegenerateFace:
    return "DegenerateFace";
  case ConvexMeshStatus::NonConvex:
    return "NonConvex";
  case ConvexMeshStatus::NonClosed:
    return "NonClosed";
  case ConvexMeshStatus::TooLarge:
    return "TooLarge";
  }
  return "Unknown";
}

// 拥有型注册表。id 顺序分配，确定性；不支持删除（热重载用 replace）。
class ConvexMeshStore {
public:
  [[nodiscard]] ConvexMeshId      add(ConvexMesh mesh);
  [[nodiscard]] bool              replace(ConvexMeshId id, ConvexMesh mesh);
  [[nodiscard]] const ConvexMesh *get(ConvexMeshId id) const noexcept;
  [[nodiscard]] std::size_t       size() const noexcept { return meshes_.size(); }
  void                            clear() noexcept { meshes_.clear(); }

private:
  std::vector<ConvexMesh> meshes_{};
};

} // namespace collision
