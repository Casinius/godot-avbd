#pragma once

// ============================================================================
//  collision/manifold.hpp —— 跨帧持久的接触流形 + 确定性存储
// ----------------------------------------------------------------------------
//  【两类状态必须分开（参考实现踩过的坑）】
//   1) cross-step：feature_key / rA / rB / lambda / penalty / beta_scale / stick
//      —— 它们要活到下一帧，用于 warm start 和摩擦锥；
//   2) step-scoped：当前帧的 basis/深度/见证点 —— 每帧由 narrowphase 重算。
//   本文件把 1) 放在 ManifoldContact，把 2) 放在 ContactPatch（narrowphase 输出），
//   两者只在 ManifoldStore::match 里合并。
//
//  【确定性契约】
//   ManifoldStore 只使用按 (proxy_a,proxy_b) 升序的 std::vector，双缓冲 current/next。
//   一个 body pair 可以对应多个 manifold（compound 的每个 shape pair 一个），
//   因为不同 shape pair 的法线/basis 不能共享。
//   **禁止 unordered_map**：哈希序会污染结果顺序，违反 state.hpp 不变量 10。
//   feature 匹配是固定 8 槽的线性扫描，冲突时取第一个匹配。
// ============================================================================

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "collision/broadphase.hpp"
#include "collision/narrowphase.hpp"

namespace collision {

inline constexpr Scalar_t kPenaltyFloor = Scalar_t(1);
inline constexpr Scalar_t kPenaltyCeil  = Scalar_t(1e10);

struct ManifoldContact {
  std::uint64_t feature_key = 0;
  Vec3          rA{0};
  Vec3          rB{0};
  Scalar_t      depth = 0;
  Vec3          lambda{0};      // normal / tangent0 / tangent1
  Vec3          penalty{1};     // 每行罚因子；第一阶段 AVBD 只读全局 k
  Scalar_t      beta_scale = 1; // 自适应 beta 预留
  bool          stick      = false;
};

struct Manifold {
  std::size_t                                          body_a  = 0;
  std::size_t                                          body_b  = 0;
  std::size_t                                          proxy_a = 0;
  std::size_t                                          proxy_b = 0;
  BasisFrame                                           basis{};
  Scalar_t                                             friction = Scalar_t(0.5);
  std::array<ManifoldContact, kMaxContactsPerManifold> contacts{};
  std::size_t                                          count = 0;

  // 下一次 motion gate / 增量 narrowphase 的位姿快照。第一阶段不消费，但字段先固定。
  Vec3 snapshot_a_lin{0};
  Vec3 snapshot_b_lin{0};
  Quat snapshot_a_ang = identity_quat();
  Quat snapshot_b_ang = identity_quat();

  [[nodiscard]] ProxyPair pair() const noexcept { return ProxyPair{proxy_a, proxy_b}; }
  [[nodiscard]] BodyPair  body_pair() const noexcept { return BodyPair{body_a, body_b}; }
};

class ManifoldStore {
public:
  // 每帧开始：清空 next 缓冲，游标归零。pairs 必须是 (a,b) 升序的 broadphase 输出。
  void begin(std::span<const BodyPair> pairs);

  // 按 ascending pair 顺序调用。patch.count == 0 时不应调用（该 manifold 会自然消失）。
  void match(const ContactPatch &patch);

  // 双缓冲交换。current 变成刚构建的、严格升序的 manifold 列表。
  void commit();

  [[nodiscard]] std::span<const Manifold> manifolds() const noexcept { return current_; }

  // 仅供 warm-start 回写：调用方只允许按 manifold_id/contact_index 改
  // lambda/penalty/beta_scale/stick，**不允许**增删、重排或改 pair/basis/count。
  // 不返回裸指针、不暴露内部 vector 的修改接口。
  [[nodiscard]] std::span<Manifold> manifolds_for_warm_start() noexcept { return current_; }

  [[nodiscard]] std::size_t size() const noexcept { return current_.size(); }

  void clear() {
    current_.clear();
    next_.clear();
    cursor_ = 0;
  }

private:
  std::vector<Manifold> current_{};
  std::vector<Manifold> next_{};
  std::size_t           cursor_ = 0;
};

} // namespace collision
