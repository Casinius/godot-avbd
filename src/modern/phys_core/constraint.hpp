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

struct ContactConstraint {
  // bodies_[0] = A（被推开的 body）；bodies_[1] = B（支撑面/另一侧）。
  // gradient 返回 {+n, -n}，所以 A 沿 +n 被推开——符号约定与 DistanceConstraint 一致。
  std::array<std::size_t, 2> bodies_       = {0, 0};
  std::uint32_t              manifold_id   = 0; // 分组用，dual 阶段按它聚合
  std::uint8_t               contact_index = 0; // 在 manifold 内的 contact 索引，决定 dual 顺序
  ContactKind                kind          = ContactKind::Normal;

  // 冗余存 basis（3 个 Vec3，36 B）：scatter 因此是纯函数，不查任何 store。
  // 分组规则见文件头的【接触分组的不变量】——由列表拥有者保证 constraints 连续。
  // 整条约束 108 B（16 + 4 + 1 + 1 + 2 pad + 36 + 12 + 12 + 20 + 1 + 3 pad），
  // 放在 variant 里约 116 B/条；改字段顺序会改变这个数字，也改变 digest 的稳定性。
  Vec3 basis[3]; // [0] = normal, [1..2] = tangents

  Vec3     rA{0}, rB{0};   // 局部接触点，角向梯度用
  Scalar_t target     = 0; // Normal: depth（静止时的允许重叠量，>0）；Tangent: 0
  Scalar_t lambda     = 0; // dual 变量；单边 ⇒ λ ≤ 0（见 update_dual）
  Scalar_t compliance = 0;
  // 下面三个字段当前**没有任何消费者**：AVBD 用的是全局 k（Policy::k），摩擦锥也还没做。
  // 保留它们是因为 digest 已经把它们纳入哈希（见 digest.hpp），删掉会改变基线；
  // 它们的真实语义要等 manifold 级聚合与摩擦一起设计。
  Scalar_t penalty    = 0;     // 预留：per-constraint penalty
  Scalar_t beta_scale = 1;     // 预留：per-constraint ramp 缩放
  bool     stick      = false; // 预留：静摩擦标记

  // C 的符号约定与 DistanceConstraint 完全相同：C < 0 ⇒ 需要沿 +g 推开。
  //   法线：C = n·(xA - xB) - depth   （比允许重叠还深 ⇒ C<0）
  //   切线：C = t·(xA - xB) - target  （target 恒为 0，一旦相对滑动 ⇒ C≠0）
  // 早期版本的法线分支直接 `return target`，C 与位姿无关（梯度是常量）——那等于
  // 把接触写成了单向推进器，永远是同一个位移，现在改为位姿相关。
  Scalar_t evaluate(const State &s) const {
    const Vec3 xA = s.pose.position(bodies_[0]) + rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 xB = s.pose.position(bodies_[1]) + rotate(s.pose.rotation(bodies_[1]), rB);
    return glm::dot(xA - xB, basis[static_cast<int>(kind)]) - target;
  }

  std::array<Vec3, 2> gradient(const State &) const {
    const Vec3 n = basis[static_cast<int>(kind)];
    return {n, -n};
  }

  std::array<Vec3, 2> gradient_angular(const State &s) const {
    const Vec3 n   = basis[static_cast<int>(kind)];
    const Vec3 rAw = rotate(s.pose.rotation(bodies_[0]), rA);
    const Vec3 rBw = rotate(s.pose.rotation(bodies_[1]), rB);
    return {glm::cross(rAw, n), glm::cross(rBw, n)};
  }

  std::span<const std::size_t> bodies() const { return bodies_; }

  // 单边接触的 dual 更新：只推不拉 ⇒ λ ∈ (-∞, 0]。
  // （符号回推：g_A = +n，delta = g·(inv_mass·dλ)；C<0 要推出 ⇒ dλ>0 ⇒ λ 不能为正。）
  // 切向的摩擦锥 |λ_t| ≤ μ|λ_n| 需要 manifold 级聚合，在它落地前**不做伪更新**，
  // 所以切向直接返回（见文件头的决策记录）。AVBD 目前不消费 λ（只用全局 k 与 C），
  // 所以这里维护 λ 主要是为了 digest 的完整性与将来摩擦锥的容量。
  void update_dual(const State &s, Scalar_t k) {
    if (kind != ContactKind::Normal)
      return;
    lambda = std::min(Scalar_t(0), lambda + k * evaluate(s));
  }
};
template <> struct constraint_traits<ContactConstraint> {
  using Scalar_t   = ::Scalar_t;
  using Gradient_t = std::array<Vec3, 2>;
};

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
