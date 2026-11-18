#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>
#include <variant>

#include "base_concept.hpp"
#include "constraint.hpp"

// ==================== 三种 Policy ====================
//
// 【Policy 契约】—— 并发化之后唯一还允许存在的接口形态。
// 每个 Policy 只提供"可以被并行调用的纯计算"，**完全不知道线程/执行器是什么**：
//
//   static constexpr bool uses_solver_buffers;
//   template <Constraint<State> C> void step(State&, C&);              // 单约束原子操作
//   void accumulate_entry(State&, std::size_t body, const ConstraintEntry&) const noexcept;
//                                                                     // 单 entry → 局部系统（仅
//                                                                     block-Jacobi 类 Policy）
//   template <typename Exec, typename Sched>
//     void begin_iterations(State&, Exec&, const Sched&);            // 每步一次
//   template <typename Exec, typename Sched>
//     void solve_primal(State&, Exec&, const Sched&);                // 逐 body
//   template <typename Ref, typename Exec, typename Sched>
//     void solve_dual(State&, std::span<Ref>, Exec&, const Sched&);  // 逐约束
//   void ramp_stiffness() noexcept;                                  // 串行标量
//
// 分工原则：Policy 决定"某个阶段内部怎么算"，driver（solve）决定"谁在算、
// 切多大块、什么时候必须停下来对齐"。所以加线程只改 driver 和 Executor，
// 不动任何一行物理代码。

// —— AVBD：分四个阶段，primal 解局部系统，dual 更新 λ，最后 ramp 刚度 ——
class AVBDStep {
public:
  // 只有 AVBD 会用到 SolverBuffers；XPBD/PBD 不走局部系统。
  static constexpr bool uses_solver_buffers = true;

  // 刚度 ramp 的起始值与步内增长率。k 是**步内**状态：每步回到 k0（见 begin_iterations）。
  // k0 可配置（引擎按材质/接触刚度映射），默认值保持历史基线不变。
  static constexpr Scalar_t kStiffness0 = Scalar_t(1e4);
  Scalar_t                  k0          = kStiffness0;
  Scalar_t                  k           = kStiffness0;
  Scalar_t                  beta        = Scalar_t(1e-2);

  // 单约束：把 k·g gᵀ / k·C g 累加进 per-body 系统（只走 buffers 接口）
  // 接触行不一样：它有自己的 K（自适应罚因子）和 λ（上一轮 dual 夹紧后的力），
  // 所以走 accumulate_row()。其余约束仍用 policy 的全局 k 且 λ = 0 ⇒ 与旧行为逐位相同。
  template <Constraint<State> C> void step(State &s, C &c) {
    const Scalar_t Cv     = c.evaluate(s);
    auto           g      = c.gradient(s);
    auto           bodies = c.bodies();

    Scalar_t K_row = k;
    Scalar_t l_row = Scalar_t(0);
    Scalar_t h_row = k;
    if constexpr (std::is_same_v<C, ContactConstraint>) {
      auto           g_ang = angular_gradient_of(c, s);
      K_row                = c.penalty;
      l_row                = c.lambda;
      // 见 contact_hessian_scale()：H 的权重按 body 做 mass-split，RHS 仍是物理力。
      const Scalar_t row_scale = Scalar_t(1);
      for (std::size_t n = 0; n < g.size(); ++n) {
        const std::size_t j     = (n == 0) ? 1 : 0;
        const Scalar_t    scale = contact_hessian_scale(s, c.bodies_[n], c.bodies_[j], g[n],
                                                        g_ang[n], g[j], g_ang[j]);
        s.buffers.accumulate_row(bodies[n], g[n], g_ang[n], K_row, Cv, l_row,
                                 K_row * scale * row_scale, c.force_min, c.force_max);
      }
      return;
    }

    if constexpr (HasAngularGradient<C, State>) {
      auto g_ang = c.gradient_angular(s);
      for (std::size_t n = 0; n < g.size(); ++n)
        s.buffers.accumulate_row(bodies[n], g[n], g_ang[n], K_row, Cv, l_row, h_row);
    } else {
      for (std::size_t n = 0; n < g.size(); ++n)
        s.buffers.accumulate_row(bodies[n], g[n], Vec3{0}, K_row, Cv, l_row, h_row);
    }
  }

  // —— 阶段 1b：把单个 entry 折进局部系统 ——
  // gather 按 body 并行调用这里，每个调用只写 buffers[b]，且同 body 内的调用顺序
  // 由 body_entries（约束升序）决定 ⇒ 与串行逐条累加位精确。
  // 注意 k 只在 accumlate 阶段读，改写它的是 ramp_stiffness()。
  // entry.penalty > 0 ⇒ 这一行自带罚因子（接触），否则用 policy 的 k、λ = 0。
  void accumulate_entry(State &s, std::size_t b, const ConstraintEntry &e) const noexcept {
    const Scalar_t K_row = e.penalty > Scalar_t(0) ? e.penalty : k;
    // hessian < 0 是"这一行不提供刚度"的哨兵（滑动的摩擦行）；== 0 表示沿用 K_row。
    const Scalar_t h_row =
        e.hessian > Scalar_t(0) ? e.hessian : (e.hessian < Scalar_t(0) ? Scalar_t(0) : K_row);
    s.buffers.accumulate_row(b, e.g_lin, e.g_ang, K_row, e.C, e.lambda, h_row, e.force_min,
                             e.force_max);
  }

  // —— 阶段 2：逐 body 解局部系统 ——
  // 并行安全性：每个 body 只写自己的 pose/cache、只读自己的 buffers/inertial，
  // 与分区无关 ⇒ **位精确**地等于串行，跟线程数无关。这也是并发收益的主战场。
  // buffers.reset(i) 塞在同一个 body 的处理尾部：原实现是在整个循环之后统一 reset，
  // 中间没有任何人读 buffers，所以逐 body 清与一次性清最终状态完全相同。
  template <typename Exec, typename Sched>
  void solve_primal(State &s, Exec &exec, const Sched &sched) {
    run_phase(exec, sched.parallel_primal, s.size(), sched.grain_body, sched.min_items_body,
              [&s](std::size_t b, std::size_t e) {
                for (std::size_t i = b; i < e; ++i) {
                  Vec3                d_lin{0}, d_ang{0};
                  const BodyResiduals r = s.inertial_residual(i);
                  if (s.buffers.solve_body(i, s.inertial, r, s.dt, d_lin, d_ang))
                    s.pose.apply_delta(i, d_lin, d_ang);
                  s.cache.refresh_rotmat(i, s.pose);
                  s.cache.refresh_dq(i, s.pose);
                  s.buffers.reset(i);
                }
              });
  }

  // —— 阶段 3：dual（λ）更新 ——
  // 并行安全性：每条约束只写自己的 λ，对 pose 只读 ⇒ 可并行，且结果与执行顺序无关。
  template <typename Ref, typename Exec, typename Sched>
  void solve_dual(State &s, std::span<Ref> refs, Exec &exec, const Sched &sched) {
    // k 在整个 dual 阶段只读；改写它的是下个阶段的 ramp_stiffness()（串行）。
    const Scalar_t k_now = k;
    run_phase(exec, sched.parallel_dual, refs.size(), sched.grain_constraint,
              sched.min_items_constraint, [&s, refs, k_now](std::size_t b, std::size_t e) {
                for (std::size_t i = b; i < e; ++i) {
                  auto &r = refs[i];
                  if (r.update_dual)
                    r.update_dual(r.obj, s, k_now);
                }
              });
  }

  // —— 阶段 3b：manifold 级 dual（**串行**，接触专用）——
  //
  // 为什么必须是串行、且按"每点三行"成组：
  //   法线的单边夹紧（F_n ≤ 0）与切线的摩擦锥夹紧（|F_t| ≤ μ|F_n|）是同一个
  //   投影的两个部分。逐行做只能各夹各的，结果不属于锥内；参考实现
  //   （avbd-demo3d Manifold::updateDual）也是每点一次算完三行。
  // 为什么放在逐约束 parallel dual 之后：接触行没有 update_dual（见 constraint.hpp），
  //   所以并行 dual 阶段对它们什么都不做，这里不会重复更新。
  // 确定性：只写每条约束自己的字段，与调度无关；约定同(manifold_id, contact_index)
  //   的行连续（见 constraint.hpp 【接触分组的不变量】5）。
  void post_dual(State &s, ConstraintList &cs) const {
    const std::size_t n = cs.size();
    std::size_t       i = 0;
    while (i < n) {
      if (!std::holds_alternative<ContactConstraint>(cs[i])) {
        ++i;
        continue;
      }
      const ContactConstraint &first = std::get<ContactConstraint>(cs[i]);
      std::size_t              end   = i + 1;
      while (end < n && end < i + 3 && std::holds_alternative<ContactConstraint>(cs[end])) {
        const ContactConstraint &next = std::get<ContactConstraint>(cs[end]);
        if (next.manifold_id != first.manifold_id || next.contact_index != first.contact_index)
          break;
        ++end;
      }
      ContactConstraint *rows[3] = {nullptr, nullptr, nullptr};
      for (std::size_t j = i; j < end; ++j) {
        ContactConstraint &c     = std::get<ContactConstraint>(cs[j]);
        const int          row_k = static_cast<int>(c.kind);
        if (row_k >= 0 && row_k < 3 && rows[row_k] == nullptr)
          rows[row_k] = &c;
      }
      dual_contact(s, rows);
      i = end;
    }
  }

  // 接触罚因子的上限：材料刚度 k0（引擎映射 beta_linear），再夹到参考实现的
  // PENALTY_MAX / PENALTY_MIN 之间。k0 设成极大值时退化成参考实现的 1e10。
  [[nodiscard]] Scalar_t contact_penalty_max() const noexcept {
    return std::min(kContactPenaltyCeil, std::max(kContactPenaltyFloor, k0));
  }

  // —— 步首准备（每步一次，在迭代循环之前）——
  // 对齐参考实现 Manifold::initialize() 的 C0 部分：取 C0 = 步首的**原始**约束值，
  // 作为整步的稳定化基准（步内不变，见 ContactConstraint::evaluate）。
  // 罚因子不在这里初始化：它由 manifold 的 warm start 带来（emit 时已按 gamma 衰减），
  // 新接触就是 floor，随后由 dual_contact 的 ramp 自适应增长。
  void capture_step_start(State &s, ConstraintList &cs) {
    for (ConstraintVariant &v : cs) {
      if (!std::holds_alternative<ContactConstraint>(v))
        continue;
      ContactConstraint &c = std::get<ContactConstraint>(v);
      c.initial_c          = c.value(s);
      // 罚因子一律从 manifold 带来的值起步（新接触 = floor），K 由 dual_contact 的
      // ramp 自适应增长——这是参考实现的做法，也是唯一同时满足下面两条的起点：
      //   a) 深度穿透（box 生成在斜坡里 9 cm / 落地 11 m/s）不能产生上万 N 的
      //      瞬时法线力：F = K·(1-alpha)·C0 与 K 成正比，K 从 floor 起步就压得住；
      //   b) 滑动的接触不能带着 K_t = 1e5 的 Hessian：局部系统会把切向自由度钉死
      //      （Δ = -β·r/(K_t+β)），μ 小的接触根本滑不动。
      // 上限用材料刚度 k0（引擎的 beta_linear）：既给 ramp 一个有物理含义的天花板，
      // 也让"设成 1e10 就永远是参考实现的 PENALTY_MAX"这种映射仍然可用。
      if (!(c.penalty > Scalar_t(0)))
        c.penalty = kContactPenaltyFloor;
      c.penalty = std::min(c.penalty, contact_penalty_max());
    }
  }

  // —— 阶段 4：刚度 ramp —— 一个标量，串行，不该并行。
  void ramp_stiffness() noexcept { k *= (Scalar_t(1) + beta); }

  // 每步一次（solve() 调用一次 = 一帧一次）：
  //  1) 刚度 k 归位到 k0。ramp 是**步内**语义；
  //  2) pose 留在步首位置 —— AVBD 是隐式残差形式，靠 solve_body 的 -M/dt²·r 把它
  //     拉向惯性目标，所以这里不要"顺手"去调用 move_body_to_inertial_target()。
  //
  // 第 1 条是修的一个真 bug：早期版本从不在步首复位 k，跨步累乘 1.01^(iters·steps)，
  // 1000 帧 @60Hz 之后 k → inf，接着 solve_body 的 10 m clamp 会把 DOF 逐个冻死
  // （表现为"跑一会儿就僵住"，且完全不报错）。现在有断言盯着它。
  template <typename Exec, typename Sched>
  void begin_iterations(State &, Exec &, const Sched &) noexcept {
    k = k0;
  }

private:
  // 单个接触点（至多三行）的联合 dual。rows[k] == nullptr 表示该行不存在（μ == 0）。
  // 公式（与参考实现逐项对应）：
  //   F      = K*C_eff + λ                       （逐行）
  //   F_n    = min(F_n, 0)                       单边
  //   F_t    = F_t * min(1, μ|F_n| / |F_t|)       径向夹紧到摩擦锥
  //   λ      = F                                 （存力，不是增量）
  //   K_row += betaLin*|C_row|                    仅在 F_n<0 / 切向仍在锥内时 ramp
  void dual_contact(State &s, ContactConstraint *const rows[3]) const {
    if (rows[0] == nullptr)
      return; // 没有法线行 ⇒ 没有锥可夹，切向也就没有意义

    Scalar_t C[3] = {0, 0, 0};
    Scalar_t K[3] = {0, 0, 0};
    Scalar_t F[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) {
      if (rows[k] == nullptr)
        continue;
      C[k] = rows[k]->evaluate(s);
      K[k] = rows[k]->penalty;
      F[k] = K[k] * C[k] + rows[k]->lambda;
    }

    F[0] = std::min(F[0], Scalar_t(0)); // 只推不拉

    const bool     tangent = rows[1] != nullptr && rows[2] != nullptr;
    const Scalar_t bounds  = std::abs(F[0]) * rows[0]->friction;
    const Scalar_t scale   = std::sqrt(F[1] * F[1] + F[2] * F[2]);
    if (tangent && scale > bounds && scale > Scalar_t(0)) {
      const Scalar_t f = bounds / scale;
      F[1] *= f;
      F[2] *= f;
    }

    // 写回：λ 是夹紧后的力；同时把 primal 的力界刷成同一组（下一轮迭代生效）。
    // 法线行**不夹** RHS：夹成 F=0 但保留 H=K 会让局部系统的最小化目标与 RHS 不再
    // 一致，二体块 Jacobi 的平移/转动耦合项就变成凭空扭矩（实测平地上纯滑动的 box
    // 会自己仰头、越转越快）。单边性由 dual 的 F_n ≤ 0 保证；旧实现也一直是两侧
    // 弹簧（C>0 时轻微吸附），保持它既一致又稳定。
    // 切向行的 ±μ|F_n| 不一致只是"欠松弛"（H 与 RHS 同向、只是被缩小），不会注入
    // 能量；而它正是防止"先粘后滑"的接触把 K_t·滑移量当摩擦力的关键。
    rows[0]->lambda    = F[0];
    rows[0]->force_min = -kForceUnclamped;
    rows[0]->force_max = kForceUnclamped;
    for (int k = 1; k < 3; ++k) {
      if (rows[k] == nullptr)
        continue;
      rows[k]->force_min = -bounds;
      rows[k]->force_max = bounds;
      rows[k]->lambda    = F[k];
    }

    // 罚因子 ramp（Eq.16）：法线只在受压时涨；切向只在**夹紧前**还在锥内时涨。
    if (F[0] < Scalar_t(0))
      rows[0]->penalty =
          std::min(K[0] + kContactPenaltyBeta * std::abs(C[0]), contact_penalty_max());
    if (tangent && scale <= bounds) {
      rows[1]->penalty =
          std::min(K[1] + kContactPenaltyBeta * std::abs(C[1]), contact_penalty_max());
      rows[2]->penalty =
          std::min(K[2] + kContactPenaltyBeta * std::abs(C[2]), contact_penalty_max());
      // 静摩擦标记：切向基本没滑 ⇒ 下一帧保留旧锚点（ManifoldStore::match 消费）。
      const bool stick = std::sqrt(C[1] * C[1] + C[2] * C[2]) < kContactStickThreshold;
      rows[1]->stick   = stick;
      rows[2]->stick   = stick;
    }
  }
};

// —— XPBD：step 内完成投影 + dual 微更新（逐约束细粒度） ——
class XPBDStep {
public:
  // 不走 per-body 局部系统：进度完全在 step() 里的直接位置投影（Gauss-Seidel）。
  static constexpr bool uses_solver_buffers = false;

  template <Constraint<State> C> void step(State &s, C &c) {
    const Scalar_t Cv     = c.evaluate(s);
    auto           g      = c.gradient(s);
    auto           bodies = c.bodies();

    if constexpr (HasAngularGradient<C, State>) {
      auto     g_ang = c.gradient_angular(s);
      Scalar_t denom = 0;
      for (std::size_t n = 0; n < g.size(); ++n)
        denom += s.inertial.effective_mass(bodies[n], g[n], g_ang[n]);

      const Scalar_t alpha_tilde = c.compliance / (s.dt * s.dt);
      const Scalar_t dlambda     = (-Cv - alpha_tilde * c.lambda) / (denom + alpha_tilde);
      c.lambda += dlambda;
      for (std::size_t n = 0; n < g.size(); ++n) {
        const Vec3 d_ang = (s.inertial.world_inv_inertia(bodies[n]) * g_ang[n]) * dlambda;
        s.pose.apply_delta(bodies[n], g[n] * (s.inertial.inv_mass(bodies[n]) * dlambda), d_ang);
      }
    } else {
      Scalar_t denom = 0;
      for (std::size_t n = 0; n < g.size(); ++n)
        denom += s.inertial.effective_mass(bodies[n], g[n], Vec3{0});

      const Scalar_t alpha_tilde = c.compliance / (s.dt * s.dt);
      const Scalar_t dlambda     = (-Cv - alpha_tilde * c.lambda) / (denom + alpha_tilde);
      c.lambda += dlambda;
      for (std::size_t n = 0; n < g.size(); ++n)
        s.pose.apply_delta(bodies[n], g[n] * (s.inertial.inv_mass(bodies[n]) * dlambda), Vec3{0});
    }
  }

  template <typename Ref> void batch_update(State &, std::span<Ref>) {} // 无批次阶段

  // XPBD 是显式位置投影：必须先把 pose 搬到惯性目标（含本次重力/外力的位移），
  // 再在它上面逐约束投影。不搬就等于把惯性项丢了。
  // 逐 body 独立 ⇒ 可并行且位精确。
  template <typename Exec, typename Sched>
  void begin_iterations(State &s, Exec &exec, const Sched &sched) {
    run_phase(exec, true, s.size(), sched.grain_body, sched.min_items_body,
              [&s](std::size_t b, std::size_t e) {
                for (std::size_t i = b; i < e; ++i)
                  s.move_body_to_inertial_target(i);
              });
  }

  // XPBD 没有批次阶段：累积和投影都在 step() 里完成。
  template <typename Exec, typename Sched> void solve_primal(State &, Exec &, const Sched &) {}
  template <typename Ref, typename Exec, typename Sched>
  void solve_dual(State &, std::span<Ref>, Exec &, const Sched &) {}
  void ramp_stiffness() noexcept {}
};

// —— PBD：XPBD 去掉 compliance/λ 项之后的退化形式 ——
// 它存在的唯一理由是给 main.cpp 的 check_xpbd_equals_pbd_when_rigid 提供对照：
// “XPBD(compliance=0) 与 PBD 位精确相同”是公式级事实（dlambda = -C/denom），
// 任何把 XPBD 重写成 Jacobi/block 形式的改动都必须先过这条断言。
// 所以这里**故意**不写成 `using PBDStep = XPBDStep`：别名会把那条断言变成恒真。
//
// 与 XPBD 的差异（只有这两处，且顺序必须先 denom 后 dlambda，才能位精确）：
//   1) 没有 λ：PBD 的位置投影不累积对偶变量；
//   2) 没有 alpha_tilde：分母就是 denom，不做 compliance 除以 dt²。
class PBDStep {
public:
  static constexpr bool uses_solver_buffers = false;

  template <Constraint<State> C> void step(State &s, C &c) {
    const Scalar_t Cv     = c.evaluate(s);
    auto           g      = c.gradient(s);
    auto           bodies = c.bodies();

    if constexpr (HasAngularGradient<C, State>) {
      auto     g_ang = c.gradient_angular(s);
      Scalar_t denom = 0;
      for (std::size_t n = 0; n < g.size(); ++n)
        denom += s.inertial.effective_mass(bodies[n], g[n], g_ang[n]);

      // 与 XPBD 在 compliance==0 时逐位一致：denom + 0.0f 与 denom 相等，
      // (-Cv - 0.0f * lambda) 与 -Cv 相等，除法走同一条指令。
      const Scalar_t dlambda = -Cv / denom;
      for (std::size_t n = 0; n < g.size(); ++n) {
        const Vec3 d_ang = (s.inertial.world_inv_inertia(bodies[n]) * g_ang[n]) * dlambda;
        s.pose.apply_delta(bodies[n], g[n] * (s.inertial.inv_mass(bodies[n]) * dlambda), d_ang);
      }
    } else {
      Scalar_t denom = 0;
      for (std::size_t n = 0; n < g.size(); ++n)
        denom += s.inertial.effective_mass(bodies[n], g[n], Vec3{0});

      const Scalar_t dlambda = -Cv / denom;
      for (std::size_t n = 0; n < g.size(); ++n)
        s.pose.apply_delta(bodies[n], g[n] * (s.inertial.inv_mass(bodies[n]) * dlambda), Vec3{0});
    }
  }

  // PBD 同样是显式位置投影：先搬到惯性目标，再逐约束投影。
  template <typename Exec, typename Sched>
  void begin_iterations(State &s, Exec &exec, const Sched &sched) {
    run_phase(exec, true, s.size(), sched.grain_body, sched.min_items_body,
              [&s](std::size_t b, std::size_t e) {
                for (std::size_t i = b; i < e; ++i)
                  s.move_body_to_inertial_target(i);
              });
  }

  template <typename Exec, typename Sched> void solve_primal(State &, Exec &, const Sched &) {}
  template <typename Ref, typename Exec, typename Sched>
  void solve_dual(State &, std::span<Ref>, Exec &, const Sched &) {}
  void ramp_stiffness() noexcept {}
};
