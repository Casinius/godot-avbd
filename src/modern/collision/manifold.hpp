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

// ==================== 跨帧热状态的释放 / 修复阈值 ====================
// 与 legacy Manifold::initialize 的 CONTACT_RELEASE_GAP / STICK_DRIFT_HEAL 同值。
//
// 释放：新鲜接触点在法线方向的真实间隙超过 margin + 2 mm 时，该对已经分开，双变量
// 里的 λ 是旧位姿的力，不能再跨过间隙施加。warm-start 衰减只有 ~0.99/步，一次冲击
// 的冲量会残留几百步，接触的罚项又每轮把它重新写进求解，表现为接触分开很久后仍被
// 弹飞。释放只丢 λ/stick，不动 penalty：刚度留着重合时一轮就能重建支撑。
//
// 粘连修复：stick 把 rA/rB 冻在静摩擦生效时的材料点上。两锚点互相粘着，它们的间距
// （STICK_THRESH 每轮量的量）永远≈0，所以哪怕整对已经滑了几十厘米，stick 也不会自己
// 解除，冻住的锚点却早已滑离真实接触面，C0 变成幻影穿透。漂移超过 2 mm 就重新采用
// 本帧的新锚点并清掉摩擦历史；否则保持冻结锚点（稳定静摩擦）。
inline constexpr Scalar_t kContactReleaseGap  = Scalar_t(2e-3);
inline constexpr Scalar_t kStickDriftHeal     = Scalar_t(2e-3);

// ==================== TGS Soft（Box2D v3）接触参数与柔度 ====================
// b2MakeSoft 的逐式移植；见 collision/tgs_soft.hpp 的注释。massScale +
// impulseScale == 1（hertz == 0 的刚体特例除外，那时两者都是 0）。
struct TgsSoftness {
  Scalar_t biasRate     = 0;
  Scalar_t massScale    = 0;
  Scalar_t impulseScale = 0;
};

struct TgsSoftParams {
  Scalar_t contact_hertz         = Scalar_t(30);
  Scalar_t contact_damping_ratio = Scalar_t(10);
  Scalar_t contact_speed         = Scalar_t(3);
};

[[nodiscard]] inline TgsSoftness make_soft(Scalar_t hertz, Scalar_t zeta,
                                           Scalar_t h) noexcept {
  if (!(hertz > Scalar_t(0)))
    return TgsSoftness{};
  const Scalar_t omega = Scalar_t(2) * Scalar_t(3.141592653589793) * hertz;
  const Scalar_t a1    = Scalar_t(2) * zeta + h * omega;
  const Scalar_t a2    = h * omega * a1;
  const Scalar_t a3    = Scalar_t(1) / (Scalar_t(1) + a2);
  return TgsSoftness{omega / a1, a2 * a3, a3};
}

struct ManifoldContact {
  std::uint64_t feature_key = 0;
  Vec3          rA{0};
  Vec3          rB{0};
  Scalar_t      depth = 0;
  Vec3          lambda{0};      // normal / tangent0 / tangent1
  Vec3          penalty{1};     // 每行罚因子；第一阶段 AVBD 只读全局 k
  Scalar_t      beta_scale = 1; // 自适应 beta 预留
  bool          stick      = false;

  // —— TGS Soft 的跨子步状态（prepare 每帧写一次，push/solve 只读/更新）——
  // lambda[0..2] 是持久化的法线/切向冲量（kg·m/s）；下面这些是同一组量的
  // per-step 工作副本。分开存是为了不污染 AVBD/B-PXPBD 的 warm-start 口径。
  Scalar_t base_separation     = 0; // 步首的 witness-point 分离（含 margin）
  Scalar_t normal_impulse      = 0;
  Scalar_t tangent_impulse[2]  = {0, 0};
  Scalar_t total_normal_impulse = 0;
  Vec3     anchorA_world{0};        // 步首 rA 的世界向量（相对质心）
  Vec3     anchorB_world{0};
  Scalar_t normal_mass         = 0; // 1 / (Σ gᵀ M⁻¹ g)，法线方向
  Scalar_t tangent_mass[2]     = {0, 0};
  TgsSoftness softness{};           // static/dynamic softness（prepare 时选择）
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

  // 清除跨帧热状态（λ / stick），保留 manifold 本身与罚因子。求解器切换时用：
  // AVBD 的 λ 是力（N），B-PXPBD 的 λ 是位置乘子（kg·m），量纲不同，不能互相
  // 暖启动。只在切换时调用，不能每帧调用（那会丢掉正常 warm start）。
  void reset_warm_start() noexcept {
    for (Manifold &m : current_) {
      for (std::size_t i = 0; i < m.count; ++i) {
        ManifoldContact &c  = m.contacts[i];
        c.lambda            = Vec3{0};
        c.stick             = false;
        c.normal_impulse    = 0;
        c.tangent_impulse[0] = 0;
        c.tangent_impulse[1] = 0;
        c.total_normal_impulse = 0;
      }
    }
  }

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
