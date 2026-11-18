#pragma once

// ============================================================================
//  collision/broadphase.hpp —— Morton LBVH + 确定性 pair 生成
// ----------------------------------------------------------------------------
//  【确定性契约】
//   1) 叶子按 (morton_code, body_index) 排序。同码时用 body_index 破平局，
//      因此结果不依赖 std::sort 的实现细节。
//   2) pair 统一存成 (min_body, max_body)，最后按 (a,b) 升序排序并去重。
//   3) 输出顺序只由 body 索引与 Morton 码决定，不依赖线程/哈希/时间。
//
//  【拥有关系】
//   Broadphase 拥有 BVH 节点和 pair 缓冲；AABB 由调用方以 span 传入，不复制。
//   没有任何 Rigid* / void* / 全局池。
// ============================================================================

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "collision/aabb.hpp"
#include "collision/bvh/node_storage.hpp"

namespace collision {

// Broadphase 本身不关心索引代表 body 还是 proxy；在 CollisionWorld 里现在喂的是
// proxy（body + slot）的 AABB，因此 pair 的两个下标也是 proxy id。
// 保留 BodyPair 这个名字作为历史别名，新增 ProxyPair 明确调用点语义。
struct BodyPair {
  std::size_t a = 0; // 恒有 a < b
  std::size_t b = 0;
};
using ProxyPair = BodyPair;

[[nodiscard]] inline bool operator<(const BodyPair &l, const BodyPair &r) noexcept {
  return l.a < r.a || (l.a == r.a && l.b < r.b);
}
[[nodiscard]] inline bool operator==(const BodyPair &l, const BodyPair &r) noexcept {
  return l.a == r.a && l.b == r.b;
}

struct BroadphaseConfig {
  // 增量 pair 是后续优化；第一阶段全量重建，但接口先固定下来。
  bool     incremental_pairs   = false;
  Scalar_t incremental_epsilon = Scalar_t(0);
};

class Broadphase {
public:
  Broadphase() = default;
  explicit Broadphase(BroadphaseConfig cfg) noexcept : cfg_(cfg) {}

  // aabbs[b] 是 body b 的世界 AABB；无效/NaN 的盒会被跳过。
  void update(std::span<const Aabb> aabbs);

  [[nodiscard]] std::span<const BodyPair> pairs() const noexcept { return pairs_; }
  [[nodiscard]] const bvh::NodeStorage   &nodes() const noexcept { return nodes_; }
  [[nodiscard]] std::size_t               node_count() const noexcept {
    return static_cast<std::size_t>(nodes_.size());
  }
  [[nodiscard]] const BroadphaseConfig &config() const noexcept { return cfg_; }

private:
  struct MortonEntry {
    std::uint64_t code = 0;
    std::size_t   body = 0;
  };

  static std::uint64_t expand_bits21(std::uint32_t v) noexcept;
  static std::uint64_t morton3d(Scalar_t x, Scalar_t y, Scalar_t z) noexcept;

  BroadphaseConfig         cfg_{};
  bvh::NodeStorage         nodes_{};
  std::vector<MortonEntry> coded_{}; // scratch，跨帧复用
  std::vector<BodyPair>    pairs_{};
  int                      leaf_count_ = 0;
  int                      first_leaf_ = 0;
};

} // namespace collision
