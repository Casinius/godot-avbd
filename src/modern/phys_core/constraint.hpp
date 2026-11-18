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

// ==================== 5. 异质约束容器 ====================
using ConstraintVariant = std::variant<DistanceConstraint, ContactConstraint>;
using ConstraintList    = std::vector<ConstraintVariant>;

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
