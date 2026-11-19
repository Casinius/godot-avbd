#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

#include "base_concept.hpp"

// ==================== 4. 具体约束 + traits 特化 ====================
//
// 【决策记录 / decision —— 为什么这里没有 Manifold / ManifoldStore】
//   曾有 phys_core/manifold.hpp 描述"manifold 头 + feature key warm start 表"。
//   结论：删除（而不是补全）。理由按权重：
//     a) 零消费者：全仓没有任何 TU include 它（xmake 只编 src/*.cpp），同一时期的
//        DualUnit 也是零使用 ⇒ "3 个 contact 作为一组联合求解"从未实现；
//     b) 字段重复 / 错位：bodies 在 Manifold 与 ContactConstraint 各存一份；friction
//        只存在于已删除的 Manifold 里，ContactConstraint 侧没有对应字段，
//        "contact 状态在 ConstraintList 里"这条口径因此没有唯一权威；
//     c) API 不可用：ManifoldStore::clear() 只有声明没有定义；prev_by_feature 用
//        unordered_map，与本仓库"结果不得依赖哈希序、每步不分配"的确定性契约冲突；
//     d) 缺 narrowphase 契约：没有接触点、法线、局部锚点、life_time、feature 对
//        稳定性约定 ⇒ 它是"索引头"，不是 manifold。
//   ⇒ 将来引入 narrowphase 时**新建**一份完整契约（4 点上限、feature 对稳定性、
//     跨帧持久化、确定性 warm start），不要复活这个外壳。
//   ⇒ warm start 落地时的硬约束（不是建议）：不得用 std::unordered_map 之类的
//     哈希容器做"上一帧索引表"——本仓库的确定性契约（state.hpp 不变量 10）要求
//     结果与调度/哈希序无关；用按 feature key 排序的 vector + begin/commit 双缓冲
//     （对齐 BodyCache::commit_prev），并在帧末固定顺序 commit。
//
// 【接触分组的不变量 —— 现在由 ContactConstraint 自己承担】
//   1) 同一 manifold 的所有 ContactConstraint 共享同一组 basis（normal + 2 tangents）；
//   2) manifold_id 相同的约束在 ConstraintList 中必须**连续**，且 manifold_id 非降序
//      ⇒ DualUnit::first 可由 3 * contact_count 的前缀和给出，不需要任何 store；
//   3) manifold_id == kStandaloneManifold 表示独立约束（不参与分组）；
//   4) 约束列表的拥有者（scene / 未来的 narrowphase）是这两条的唯一维护者：
//      改动顺序或增删就必须重新分配 manifold_id。
//   5) 同一个接触点的行按 normal, tangent0, tangent1 顺序**紧邻**发射
//      （contact_constraints.hpp 保证）；AVBDStep::post_dual 按
//      (manifold_id, contact_index) 把它们聚成一组做联合夹紧。μ == 0 的接触只发
//      法线一行。

// —— 距离约束（两个 body，纯线性 DOF） ——
struct DistanceConstraint {
  std::array<std::size_t, 2> bodies_     = {0, 0};
  Scalar_t                   rest_length = 1.0;
  Scalar_t                   lambda      = 0.0; // XPBD / AVBD 的 dual 变量
  Scalar_t                   compliance  = 0.0;

  Scalar_t evaluate(const State &s) const {
    // 注意：glm 的 vec::length() 返回分量个数(3)，这里必须用 glm::length
    return glm::length(s.pose.position(bodies_[0]) - s.pose.position(bodies_[1])) - rest_length;
  }
  std::array<Vec3, 2> gradient(const State &s) const {
    Vec3           d   = s.pose.position(bodies_[0]) - s.pose.position(bodies_[1]);
    const Scalar_t len = glm::length(d);
    const Vec3     n   = (len > Scalar_t(1e-12)) ? d * (Scalar_t(1) / len) : Vec3{0, 0, 0};
    return {n, n * Scalar_t(-1)};
  }
  // 无角向自由度：显式给空梯度，让 Policy 的 if constexpr 走"纯线性"分支
  std::array<Vec3, 2> gradient_angular(const State &) const {
    return {Vec3{0, 0, 0}, Vec3{0, 0, 0}};
  }
  std::span<const std::size_t> bodies() const { return bodies_; }

  void update_dual(const State &s, Scalar_t k) { lambda += k * evaluate(s); }
};

template <> struct constraint_traits<DistanceConstraint> {
  using Scalar_t   = ::Scalar_t;
  using Gradient_t = std::array<Vec3, 2>;
};

// —— 接触约束（两个 body，单边：法线 + 两条切线共用一个 basis） ——
//
// manifold_id 的哨兵值：表示“这条约束不属于任何 manifold”（独立约束）。
// DualUnit::manifold 用同一个值，所以集中定义在这里，不要在别处再写 UINT32_MAX。
inline constexpr std::uint32_t kStandaloneManifold = 0xFFFFFFFFu;

enum class ContactKind : std::uint8_t { Normal, Tangent0, Tangent1 };

// ==================== 接触的罚函数 / 摩擦常数 ====================
// 全部对齐参考实现 savant117/avbd-demo3d（source/manifold.cpp + solver.h）：
//   alpha = 0.99 (stabilization), gamma = 0.999 (warmstart decay),
//   betaLin = 1e4 (penalty ramp), PENALTY_MIN/MAX = 1 / 1e10, STICK_THRESH = 1e-5。
//
// 约束值用一阶稳定化：
//     C_eff = C_current - alpha * C0        （C0 = 步首原始违反量 = initial_c）
// 它等价于参考实现的  C = C0*(1-alpha) + J*dq  —— 因为 basis/锚点在步内固定，
// C_current = C0 + J*dq 精确成立。C0 每步只在 capture 时取一次。
inline constexpr Scalar_t kContactStabilization = Scalar_t(0.99);
inline constexpr Scalar_t kWarmStartGamma       = Scalar_t(0.999);
// warm start 的 λ 衰减率 = alpha * gamma（参考实现 Eq.19）。
inline constexpr Scalar_t kWarmStartDecay = kContactStabilization * kWarmStartGamma;
inline constexpr Scalar_t kContactPenaltyFloor = Scalar_t(1);
inline constexpr Scalar_t kContactPenaltyCeil  = Scalar_t(1e10);
// 每迭代的罚因子 ramp 速率（单位 1/m：K += betaLin*|C|）。
inline constexpr Scalar_t kContactPenaltyBeta = Scalar_t(1e4);
// 静摩擦判定：切向 C 的模小于它 ⇒ 仍然粘着，下一帧保留旧锚点。
inline constexpr Scalar_t kContactStickThreshold = Scalar_t(1e-5);

// ==================== 关节行的罚函数 / 稳定化常数 ====================
// 与接触同一套语义（alpha / gamma / floor / ceil），但角的 beta 用参考实现的 betaAng，
// 且两者都受 penalty_limit() 的频率上限约束。单独命名是为了让"改接触会顺手改到关节"
// 不会发生；值本身必须与接触一致，下面的 static_assert 盯着。
inline constexpr Scalar_t kJointStabilization   = kContactStabilization;
inline constexpr Scalar_t kJointWarmStartGamma  = kWarmStartGamma;
inline constexpr Scalar_t kJointWarmStartDecay  = kWarmStartDecay;
inline constexpr Scalar_t kJointPenaltyFloor    = kContactPenaltyFloor;
inline constexpr Scalar_t kJointPenaltyCeil     = kContactPenaltyCeil;
// 新行的起始罚因子。接触从 floor 起步（它的 ramp 由 |C|*1e4 一步拉满），关节的
// C 被稳定化压到 1% 后 ramp 太慢，靠它跨帧爬会让锚点长期停在毫米级。直接用材料
// 刚度起步，再由 penalty_limit() 夹到该自由度能承受的上限。
inline constexpr Scalar_t kJointInitialPenalty = Scalar_t(1e4);
inline constexpr Scalar_t kJointPenaltyBetaLin  = Scalar_t(1e4);
inline constexpr Scalar_t kJointPenaltyBetaAng  = Scalar_t(1e2);
// 约束频率上限 = 该比例 * 步率。参考实现的 0.25 是给"λ 直接反馈"的显式迭代留的余量；
// 现代核心的 primal 是带 m/dt² 的隐式 6x6 块解、并做了 mass-split，能承受与接触同级
// 的刚度。4.0 给 1 kg @ 60 Hz 的线性自由度留 5.76e4 N/m 的上限（接触默认 k0=1e4），
// 关节在头几帧就能把重力位移压到亚毫米，而不是靠 λ 跨帧慢慢爬。
inline constexpr Scalar_t kJointMaxFrequencyRatio = Scalar_t(4);

struct ContactConstraint {
  // bodies_[0] = A（被推开的 body）；bodies_[1] = B（支撑面/另一侧）。
  // gradient 返回 {+n, -n}，所以 A 沿 +n 被推开——符号约定与 DistanceConstraint 一致。
  std::array<std::size_t, 2> bodies_       = {0, 0};
  std::uint32_t              manifold_id   = 0; // 分组用，dual 阶段按它聚合
  std::uint8_t               contact_index = 0; // 在 manifold 内的 contact 索引，决定 dual 顺序
  ContactKind                kind          = ContactKind::Normal;

  // 冗余存 basis（3 个 Vec3，36 B）：scatter 因此是纯函数，不查任何 store。
  // 分组规则见文件头的【接触分组的不变量】——由列表拥有者保证 constraints 连续。
  Vec3 basis[3]; // [0] = normal, [1..2] = tangents

  Vec3     rA{0}, rB{0};   // 局部接触点，角向梯度用
  Scalar_t target     = 0; // Normal: -margin（静止时的允许重叠量）；Tangent: 0
  Scalar_t compliance = 0;

  // —— dual / 罚因子状态（每行一份；三行是同一次 manifold 级 dual 的联合状态）——
  // lambda   ：上一次 dual 夹紧后的**力** F = K*C + λ（法线 ≤ 0；切线在锥内）。
  //            注意它存的是力，不是 λ 递增量——增量形式（λ += k*C）会发散。
  // penalty  ：该行的自适应罚因子 K。接触行只用它，非接触行用 Policy 的全局 k。
  // initial_c：步首的原始 C0，稳定化用（每步由 capture_step_start 重取一次）。
  // friction ：该接触点的 μ（同一 contact 的三行相同；行数也因此决定）。
  // force_min/max：primal 里该行力的界（参考实现 updatePrimal 也夹紧：法线 F ≤ 0，
  //            切向 |F| ≤ μ|F_n|）。必须夹在 RHS 上：否则一个"先粘住、后滑动"的
  //            接触会拿 K_t*C_t（= K_t·滑移量）当摩擦力，滑移量的量级是 cm ⇒ 上万 N，
  //            表现为 box 一落地就被焊住、再被甩飞。锥界用**上一轮**的 F_n（lag 一轮），
  //            因为 gather 是 per-row 的，没有"三行同时"的位置；差别只有一轮迭代。
  Scalar_t lambda     = 0;
  Scalar_t penalty    = 0;
  Scalar_t initial_c  = 0;
  Scalar_t friction   = 0;
  Scalar_t force_min  = -kForceUnclamped;
  Scalar_t force_max  = kForceUnclamped;
  Scalar_t beta_scale = 1;     // 预留：per-constraint ramp 缩放（未消费）
  bool     stick      = false; // 上一帧是否粘着：决定下一帧是否重锚定 rA/rB

  // 原始约束值（不掺稳定化）。C 的符号约定与 DistanceConstraint 相同：C < 0 ⇒ 需要沿 +g 推开。
  //   法线：C = n·(xA - xB) + margin   （穿透超过 margin ⇒ C<0）
  //   切线：C = t·(xA - xB)            （一旦相对滑动 ⇒ C≠0）
  Scalar_t value(const State &s) const {
    const Vec3 xA = s.pose.position(bodies_[0]) + rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 xB = s.pose.position(bodies_[1]) + rotate(s.pose.rotation(bodies_[1]), rB);
    return glm::dot(xA - xB, basis[static_cast<int>(kind)]) - target;
  }

  // 求解器看到的约束值：一阶稳定化后的 C_eff（见文件头的常数说明）。
  // initial_c == 0（非接触 / XPBD 路径）时与 value() 完全相同。
  Scalar_t evaluate(const State &s) const {
    return value(s) - kContactStabilization * initial_c;
  }

  std::array<Vec3, 2> gradient(const State &) const {
    const Vec3 n = basis[static_cast<int>(kind)];
    return {n, -n};
  }

  std::array<Vec3, 2> gradient_angular(const State &s) const {
    const Vec3 n   = basis[static_cast<int>(kind)];
    const Vec3 rAw = rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 rBw = rotate(s.pose.rotation(bodies_[1]), rB);
    // C = n·(xA - xB) - target ⇒
    //   ∂C/∂θ_A = +(rAw × n)
    //   ∂C/∂θ_B = -(rBw × n)      ← B 的贡献带负号（与线性梯度 {+n,-n} 一致）
    // 这里曾经漏掉 B 的负号：接触力作用在 B 上的**力矩方向**因此是反的，表现为
    // 平地上纯滑动的 box 会朝"翻车"的反方向持续自旋、再靠法线重分布去纠——而
    // 1~2 点 manifold 根本纠不回来，于是整摞/整块持续被泵能量。
    return {glm::cross(rAw, n), -glm::cross(rBw, n)};
  }

  std::span<const std::size_t> bodies() const { return bodies_; }

  // 为什么**没有** update_dual：
  // 接触的三行（normal/t0/t1）是一个联合问题——法线单边夹紧和切向摩擦锥夹紧必须
  // 同时做，逐行更新会把切线拉出锥外（也是上一版发散的原因之一）。
  // 联合更新在 AVBDStep::post_dual（串行，按 manifold 分组）里完成。
  // ConstraintRef::update_dual 因此对接触为 nullptr ⇒ 逐约束 dual 阶段跳过它们。
};
template <> struct constraint_traits<ContactConstraint> {
  using Scalar_t   = ::Scalar_t;
  using Gradient_t = std::array<Vec3, 2>;
};

// ==================== 接触行的 block-Jacobi mass-split 权重 ====================
//
// 【为什么要它】本仓库的 accumulate 是 scatter + per-body gather（block-Jacobi）：
// 两条约束行对同一对 body 的贡献在**同一个线性化点**上被折进两个独立的 6x6
// 局部系统。对一条刚度 K、两侧可动质量倒数 em_self / em_other 的行，
// 每侧各走一步 ⇒ 相对位移是
//     u_jacobi = -2F/(K + 1/em)        （K ≫ 1/em 时 ≈ -2F/K）
// 而耦合精确解是 u = -F/(K + 1/(em_self+em_other)) ≈ -F/K。
// 没有 λ 时这只是"每迭代乘 -1"的临界振荡（旧行为，有界）；λ 一旦进 RHS，
// 迭代变成 e → -3e（特征值 2+√3），几轮就把接触推开、下一帧再砸回来。
//
// 【修法】把每个 body 块的 H 乘上 (1 + em_other/em_self)：
//     Δ_self = -F / (K(1+em_other/em_self) + 1/em_self)
//            = em_self/(em_self+em_other) · u
// 正好是耦合精确解按可动性分配的那一份。静/动（另一侧 em_other == 0）时系数
// 为 1，保持"半步收敛"的行为；两侧都动（等质量）时系数为 2，相对位移一步到位。
// RHS 的力不受影响（仍是 K·C + λ），所以平衡点、λ 的语义与摩擦锥的界都不变。
[[nodiscard]] inline Scalar_t contact_hessian_scale(const State &s, std::size_t body_self,
                                                   std::size_t body_other, Vec3 g_self,
                                                   Vec3 ga_self, Vec3 g_other,
                                                   Vec3 ga_other) noexcept {
  const Scalar_t em_self  = s.inertial.effective_mass(body_self, g_self, ga_self);
  const Scalar_t em_other = s.inertial.effective_mass(body_other, g_other, ga_other);
  if (!(em_self > Scalar_t(0)) || !(em_other > Scalar_t(0)))
    return Scalar_t(1); // 任一侧静态：没有对侧过冲
  return Scalar_t(1) + em_other / em_self;
}

// （已删除：contact_sort_key / feature_id 排序草案。它引用了一个并不存在的
//  feature 概念，且排序 key 只在 narrowphase 侧有意义——等 narrowphase 落地时
//  连同 manifold 契约一起设计，见文件头的决策记录。）

// ==================== 关节轴约束（每个自由度一行） ====================
//
// 【为什么是一行一个约束，而不是一个 6-DOF 约束】
//   驱动器的 ConstraintList 每个 alternative 只有**一个标量 λ** 和一组 body 上的梯度。
//   一个带 Vec6 λ 的复合关节无法塞进这个形态：duel 阶段无法逐行并行，也无法复用既有的
//   scatter/gather。所以关节按"每轴一行"发射，和接触"每点三行"完全同构：
//     Locked  → 双边等式（λ 不夹紧）
//     Limited → 单边不等式（λ 夹到指向可行域的一侧；在范围内时 λ = 0）
//     Free    → 不发射行
//   每行独立，所以既有的并行 dual / 拓扑缓存 / 位精确性全部原样成立。
//
// 【线性 vs 角向】
//   * 线性行：C = axis·(anchorA - anchorB) - target，轴与锚点在 body A 的局部系里。
//     梯度同时有线性项和角向项（锚点偏离质心时，转动也会移动锚点）。
//   * 角向行：C 是绕 joint frame 某轴的转角（弧度），测量自关节建立时的相对姿态。
//     必须提供 gradient_angular；线性梯度为 0。
//   * body A / body B 都是真实 State 下标；世界侧由调用方用一个静态 body 表示
//     （现代服务器的 ModernSpace 常驻一个 world body），因此这里不需要 nullptr 分支。
//
// 【角向测量的 twist 分解】
//   锁定轴不能直接读相对旋转向量的分量：绕自由轴转一整圈会让它每半圈扫过整圈，铰链
//   会在几圈后把自己撕开。所以先把 dev 中绕每个**非锁定**轴的 twist 拆出来（那个轴的
//   读数就是 twist 角，限位/弹簧读它），剩下的 swing 才交给锁定轴。拆一个 twist 是精确的；
//   多个自由轴时逐个拆是通常近似。
//
// 【稳定化 / 罚因子 / 暖启动】
//   与接触同一套常数与语义：
//     Locked : C_eff = (C - target) - alpha * C0    （C0 = 步首原始值）
//              λ ← K*C_eff + λ（存力，不是增量），K ← K + beta|C_eff|
//     Limited: 只在越界时产生 C = value - bound；λ 强制指向可行域一侧；范围内 λ=0、无刚度。
//   暖启动在 capture_step_start 里做一次/步（λ *= alpha*gamma，K *= gamma）。
//   刚度的上限见 penalty_limit()：罚因子涨到求解器在一个步长内分辨不了时就只会发散。
//
// 【已知缺口】spring / softness 没有映射：一行只有一个 (K, λ, C)，装不下"硬约束 + 弹簧"
//   两套力。Limited/Locked/Free 与 target 已覆盖 6DOF / pin / hinge / slider 的服务端
//   映射；弹簧参数仍只记录在引擎记录里（与 legacy 服务器同一口径）。
enum class JointAxisMode : std::uint8_t {
  Free    = 0,
  Locked  = 1,
  Limited = 2,
};

namespace joint_axis_detail {

[[nodiscard]] inline Vec3 principal_axis(int p_axis) noexcept {
  static const std::array<Vec3, 3> axes = {Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
  return axes[static_cast<std::size_t>(p_axis)];
}

[[nodiscard]] inline Scalar_t quat_component(const Quat &p_q, int p_axis) noexcept {
  switch (p_axis) {
    case 0: return p_q.x;
    case 1: return p_q.y;
    default: return p_q.z;
  }
}

// 把 dev 中绕 axis 的 twist 拆掉，角度写进 r_angle。四元数按 glm 的 {w,x,y,z} 构造。
// 纯 twist = normalize((n*(v·n), w))；近半圈时 twist 退化（vector 与 w 同时趋零），
// 此时没有可拆的 twist。
inline void remove_twist(Quat &p_dev, int p_axis, Scalar_t &r_angle) noexcept {
  const Vec3 n = principal_axis(p_axis);
  const Vec3 v{p_dev.x, p_dev.y, p_dev.z};
  const Scalar_t proj = glm::dot(v, n);
  const Quat twist{p_dev.w, n * proj}; // glm: qua(w, xyz)
  if (glm::length(twist) <= Scalar_t(1.0e-12)) {
    r_angle = Scalar_t(0);
    return;
  }
  r_angle = Scalar_t(2) * std::atan2(proj, p_dev.w);
  p_dev  = glm::normalize(p_dev) * glm::conjugate(glm::normalize(twist));
}

// 三个角向自由度的读数（弧度，相对建立时的姿态，joint frame 下）。
// Locked 轴读 swing 分量；Free/Limited 轴读拆出来的 twist 角。
inline void angular_values(const State &s, const std::array<std::size_t, 2> &bodies,
                           const Quat &frame_a, const Quat &frame_b, const Quat &rest,
                           const std::array<JointAxisMode, 3> &modes,
                           Scalar_t p_out[3]) noexcept {
  const Quat q_a = glm::normalize(s.pose.rotation(bodies[0]) * frame_a);
  const Quat q_b = glm::normalize(s.pose.rotation(bodies[1]) * frame_b);
  const Quat relative_now = glm::normalize(glm::conjugate(q_a) * q_b);
  // dev 是 **B 相对 A** 自建立以来多转的量（先取 rest⁻¹ 再乘当前相对姿态）。
  // 这与 Godot 的 hinge/slider 角度同号，服务器可以直接把限位抄过来；角向梯度则
  // 相应取 {-n, +n}，与 dev 的符号一致。
  Quat dev = glm::normalize(glm::conjugate(rest) * relative_now);

  // 非锁定轴的 twist 先拆掉（顺序固定为 0,1,2），锁定轴只看剩下的 swing。
  for (int i = 0; i < 3; ++i) {
    if (modes[static_cast<std::size_t>(i)] != JointAxisMode::Locked)
      remove_twist(dev, i, p_out[i]);
  }
  for (int i = 0; i < 3; ++i) {
    if (modes[static_cast<std::size_t>(i)] == JointAxisMode::Locked)
      p_out[i] = quat_component(dev, i) * Scalar_t(2);
  }
}

} // namespace joint_axis_detail

// 一个自由度对应的一行。自由度是 body A 的 joint frame 的主轴。
struct JointAxisConstraint {
  std::array<std::size_t, 2> bodies_ = {0, 0};
  Vec3     rA{0}, rB{0};        // 局部锚点（线性行用；角向行也可携带，便于调试）
  Vec3     axis{1, 0, 0};       // 轴方向，在 body A 的**局部系**里（joint frame 已折进来）
  bool     angular = false;     // true: 角向自由度；false: 线性
  JointAxisMode mode = JointAxisMode::Free;
  Scalar_t target = 0;          // Locked 的目标值（线性：米；角向：弧度）
  Scalar_t lower  = 0;          // Limited 的可行区间（与测量同单位）
  Scalar_t upper  = 0;

  // 角向测量用：joint frame 在各自 body 局部系里的朝向，建立时的相对姿态，以及三个
  // 角向轴的模式（twist 分解需要知道哪些轴非锁定）。线性行不使用这些字段。
  Quat frame_a = identity_quat();
  Quat frame_b = identity_quat();
  Quat rest    = identity_quat();
  std::array<JointAxisMode, 3> angular_modes{{JointAxisMode::Free, JointAxisMode::Free,
                                              JointAxisMode::Free}};
  std::uint8_t axis_index = 0;

  Scalar_t lambda     = 0; // 上一次 dual 夹紧后的**力** F（不是增量）
  Scalar_t penalty    = 0; // 自适应罚因子；<= 0 表示未初始化
  Scalar_t initial    = 0; // 步首的原始 C（稳定化基准）
  Scalar_t compliance = 0; // 仅 XPBD 路径使用；AVBD 不读

  [[nodiscard]] Scalar_t value(const State &s) const {
    if (angular) {
      Scalar_t values[3];
      joint_axis_detail::angular_values(s, bodies_, frame_a, frame_b, rest, angular_modes,
                                        values);
      return values[axis_index];
    }
    const Vec3 x_a = s.pose.position(bodies_[0]) + rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 x_b = s.pose.position(bodies_[1]) + rotate(s.pose.rotation(bodies_[1]), rB);
    const Vec3 n   = rotate(s.pose.rotation(bodies_[0]), axis);
    return glm::dot(x_a - x_b, n);
  }

  // 求解器看到的约束值：Locked 做一阶稳定化，Limited 只在越界时返回违反量（范围内为 0）。
  [[nodiscard]] Scalar_t evaluate(const State &s) const {
    const Scalar_t v = value(s);
    switch (mode) {
      case JointAxisMode::Locked:
        return (v - target) - kJointStabilization * initial;
      case JointAxisMode::Limited:
        if (v > upper) return v - upper;
        if (v < lower) return v - lower;
        return Scalar_t(0);
      case JointAxisMode::Free:
      default:
        return Scalar_t(0);
    }
  }

  [[nodiscard]] std::array<Vec3, 2> gradient(const State &s) const {
    if (angular)
      return {Vec3{0}, Vec3{0}};
    const Vec3 n = rotate(s.pose.rotation(bodies_[0]), axis);
    return {n, -n};
  }

  [[nodiscard]] std::array<Vec3, 2> gradient_angular(const State &s) const {
    const Vec3 n = rotate(s.pose.rotation(bodies_[0]), axis);
    // 角向测量是 B 相对 A（见 angular_values），所以符号与线性行相反。
    if (angular)
      return {-n, n};
    const Vec3 r_a = rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 r_b = rotate(s.pose.rotation(bodies_[1]), rB);
    return {glm::cross(r_a, n), -glm::cross(r_b, n)};
  }

  [[nodiscard]] std::span<const std::size_t> bodies() const { return bodies_; }

  // 这一行的 RHS 力允许的区间：Limited 只许朝可行域内推；范围内 [0,0] 表示不施力。
  void force_bounds(const State &s, Scalar_t &r_min, Scalar_t &r_max) const {
    if (mode != JointAxisMode::Limited) {
      r_min = -kForceUnclamped;
      r_max = kForceUnclamped;
      return;
    }
    const Scalar_t v = value(s);
    if (v > upper) {
      r_min = Scalar_t(0);
      r_max = kForceUnclamped;
    } else if (v < lower) {
      r_min = -kForceUnclamped;
      r_max = Scalar_t(0);
    } else {
      r_min = Scalar_t(0);
      r_max = Scalar_t(0);
    }
  }

  // 该自由度能承受的罚因子上限。约束的频率是 sqrt(K * 有效可动性)，接近步率时迭代不再
  // 收敛：误差涨 ⇒ K 涨 ⇒ 更不收敛。上限取"频率仍是步率一小部分"的 K，用两侧 body 里
  // 更轻的那一侧（参考实现的 GenericJoint::penaltyLimit 同一公式）。
  [[nodiscard]] Scalar_t penalty_limit(const State &s) const {
    // 用参考实现的量纲：线性轴看质量，角向轴看绕轴的转动惯量。注意**不要**用整行的
    // gᵀM⁻¹g：离轴锚点的线性行还耦合着转动，那个有效质量会把上限压到 1/7，关节在
    // 头几帧会像没有约束一样下坠。
    Scalar_t cap = kForceUnclamped;
    const Vec3 n = rotate(s.pose.rotation(bodies_[0]), axis);
    for (std::size_t i = 0; i < 2; ++i) {
      Scalar_t inertia = 0;
      if (angular) {
        const Scalar_t inv_along =
                glm::dot(n, s.inertial.world_inv_inertia(bodies_[i]) * n);
        if (!(inv_along > Scalar_t(0)))
          continue;
        inertia = Scalar_t(1) / inv_along;
      } else {
        const Scalar_t inv_mass = s.inertial.inv_mass(bodies_[i]);
        if (!(inv_mass > Scalar_t(0)))
          continue;
        inertia = Scalar_t(1) / inv_mass;
      }
      const Scalar_t limit = inertia * kJointMaxFrequencyRatio * kJointMaxFrequencyRatio /
                             (s.dt * s.dt);
      cap = std::min(cap, limit);
    }
    return cap;
  }

  // 每步一次，在迭代之前：记录稳定化基准并暖启动对偶状态（Eq.19）。
  void capture_step_start(const State &s) {
    if (mode == JointAxisMode::Free)
      return;
    initial = value(s) - target;
    lambda *= kJointWarmStartDecay;
    const Scalar_t cap = std::max(penalty_limit(s), kJointPenaltyFloor);
    if (!(penalty > Scalar_t(0)))
      penalty = std::min(kJointInitialPenalty, cap);
    else
      penalty = std::clamp(penalty * kJointWarmStartGamma, kJointPenaltyFloor, cap);
    // 互补松弛：不在边界上的限位没有力；留着上一帧的 λ 会在回到边界时踢一脚。
    if (mode == JointAxisMode::Limited) {
      const Scalar_t v = value(s);
      if (!(v > upper) && !(v < lower))
        lambda = Scalar_t(0);
    }
  }

  void update_dual(const State &s, Scalar_t /*p_k*/) {
    if (mode == JointAxisMode::Free)
      return;
    const Scalar_t C = evaluate(s);
    Scalar_t F       = penalty * C + lambda;
    if (mode == JointAxisMode::Limited) {
      if (C > Scalar_t(0))
        F = std::max(F, Scalar_t(0));
      else if (C < Scalar_t(0))
        F = std::min(F, Scalar_t(0));
      else
        F = Scalar_t(0);
    }
    lambda = F;

    const bool active = mode == JointAxisMode::Locked || C != Scalar_t(0);
    if (active) {
      const Scalar_t beta = angular ? kJointPenaltyBetaAng : kJointPenaltyBetaLin;
      const Scalar_t cap  = std::max(penalty_limit(s), kJointPenaltyFloor);
      penalty = std::clamp(penalty + std::abs(C) * beta, Scalar_t(0), cap);
    }
  }
};

template <> struct constraint_traits<JointAxisConstraint> {
  using Scalar_t   = ::Scalar_t;
  using Gradient_t = std::array<Vec3, 2>;
};

// ==================== 5. 异质约束容器 ====================
// 注意：新增 alternative 只能在**末尾**追加。variant index 被 collision digest 哈希，
// 改变既有顺序会无谓地换掉已有场景的摘要。
using ConstraintVariant =
    std::variant<DistanceConstraint, ContactConstraint, JointAxisConstraint>;
using ConstraintList = std::vector<ConstraintVariant>;

// ==================== 8. Policy 相关的类型擦除视图 ====================
// variant 里某个类型的 alternative 下标（编译期）。
template <typename C, typename V> struct variant_index;
template <typename C, typename... Ts> struct variant_index<C, std::variant<Ts...>> {
  static constexpr std::size_t value = [] {
    std::size_t i = 0, out = 0;
    ((std::is_same_v<C, Ts> ? (out = i) : 0, ++i), ...);
    return out;
  }();
};

template <typename Policy> struct ConstraintRef {
  void       *obj        = nullptr;
  std::size_t type_index = 0; // variant alternative 下标，用于按类型分批

  Scalar_t (*evaluate)(const void *, const State &)              = nullptr;
  void (*gradient)(const void *, const State &, std::span<Vec3>) = nullptr;
  // 只有 6-DOF 约束才提供；constraint 无角向自由度时为 nullptr，scatter 会填
  // 0。
  void (*gradient_angular)(const void *, const State &, std::span<Vec3>) = nullptr;
  std::span<const std::size_t> (*bodies)(const void *)                   = nullptr;

  // 只有 AVBD 用到，nullable
  void (*update_dual)(void *, const State &, Scalar_t) = nullptr;

  // 单约束原子操作（转发到 Policy::step<C>）
  void (*step_fn)(void *, State &, Policy &) = nullptr;
  void step(State &s, Policy &p) { step_fn(obj, s, p); }

  template <typename C>
    requires Constraint<C, State>
  static ConstraintRef from(C &c) {
    ConstraintRef r;
    r.obj        = &c;
    r.type_index = variant_index<C, ConstraintVariant>::value;
    r.evaluate   = [](const void *p, const State &s) {
      return static_cast<const C *>(p)->evaluate(s);
    };
    r.gradient = [](const void *p, const State &s, std::span<Vec3> out) {
      auto g = static_cast<const C *>(p)->gradient(s);
      std::copy(g.begin(), g.end(), out.begin());
    };
    r.bodies = [](const void *p) { return static_cast<const C *>(p)->bodies(); };

    if constexpr (HasAngularGradient<C, State>) {
      r.gradient_angular = [](const void *p, const State &s, std::span<Vec3> out) {
        auto ga = static_cast<const C *>(p)->gradient_angular(s);
        std::copy(ga.begin(), ga.end(), out.begin());
      };
    }

    if constexpr (requires(C &cc, const State &ss, Scalar_t k) { cc.update_dual(ss, k); }) {
      r.update_dual = [](void *p, const State &s, Scalar_t k) {
        static_cast<C *>(p)->update_dual(s, k);
      };
    }
    r.step_fn = [](void *p, State &s, Policy &pol) { pol.step(s, *static_cast<C *>(p)); };
    return r;
  }
};
