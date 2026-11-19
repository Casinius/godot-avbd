#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

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
  // 拓扑布局（per-constraint entry 邻接 / block 分区）只有 AVBD/BPXPBD 需要。
  // 两者共用同一个 TopologyKey，但 block 只在 constraint_blocks 策略下构建。
  static constexpr bool needs_topology = true;
  // AVBD 的约束更新在 scatter/gather 里，不扫 block。
  static constexpr bool constraint_blocks = false;

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

    // 关节行与接触行同一形态：自有 K/λ/力界；Hessian 走 mass-split，限位内不提供刚度。
    // 走 solve() 时它们由 scatter_typed 的专用重载处理，这条分支是给直接调 step() 的
    // Policy 用法（以及 selftest 的 scatter/gather 对照）对齐同一组数字。
    if constexpr (std::is_same_v<C, JointAxisConstraint>) {
      auto     g_ang = c.gradient_angular(s);
      Scalar_t force_min = -kForceUnclamped;
      Scalar_t force_max = kForceUnclamped;
      c.force_bounds(s, force_min, force_max);
      const Scalar_t row_k = c.penalty > Scalar_t(0) ? c.penalty : k;
      const bool inactive_limit = c.mode == JointAxisMode::Limited && Cv == Scalar_t(0);
      for (std::size_t n = 0; n < g.size(); ++n) {
        const std::size_t j = (n == 0) ? 1 : 0;
        const Scalar_t h = inactive_limit
                ? Scalar_t(0)
                : row_k * contact_hessian_scale(s, c.bodies_[n], c.bodies_[j], g[n], g_ang[n],
                                                g[j], g_ang[j]);
        s.buffers.accumulate_row(bodies[n], g[n], g_ang[n], row_k, Cv, c.lambda, h, force_min,
                                 force_max);
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
      // 关节行：记录本步的稳定化基准并暖启动 λ/K，一次/步。
      if (auto *joint = std::get_if<JointAxisConstraint>(&v)) {
        joint->capture_step_start(s);
        continue;
      }
      if (!std::holds_alternative<ContactConstraint>(v))
        continue;
      ContactConstraint &c = std::get<ContactConstraint>(v);
      c.capture_contact_c0(s);
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
  static constexpr bool needs_topology = false;
  static constexpr bool constraint_blocks = false;

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

// —— B-PXPBD：Blended Primal XPBD（Chen et al., MIG'23） ——
//
// 在 XPBD 的位置更新上多加一个全局残差修正项：
//
//   x̃ = x + dt·(v + M⁻¹f_ext)           （begin_iterations 已把 pose 搬到 x̃）
//   for k in 0..K-1:
//     R_i = Σ_c λ_c ∇C_c(x)             （begin_sweep，串行、约束升序）
//     for each unit（着色 block 扫掠）:
//       C = C_d(x); g = ∇C_d(x); w = Σ gᵀM⁻¹g
//       α̃ = compliance/dt²
//       Δλ = -(C + α̃λ) / (w + α̃)
//       λ += Δλ
//       for i: B_i = (x̃_i - x_i) + M_i⁻¹R_i
//              x_i += M_i⁻¹g_iΔλ + (1-ζ)B_i
//              R_i += g_iΔλ
//   v = (x - x_prev)/dt                 （finish_step）
//
// ζ = 1 时 blend 项恒为 0 ⇒ 与 XPBDStep 的串行 Gauss-Seidel **逐位相同**（已有断言）。
// ζ = 0 是一次精确 fixed-point 步；默认 0.5。
//
// 与 AVBD 的语义差：λ 是**乘子**（力 = λ∇C），不是 AVBD 存的力 F = K·C + λ。
// 关节 Limited 行的可行方向因此与 AVBD 相反，夹紧区间取 (-force_max, -force_min)。
// 关节初始稳定化基准 `initial` 必须为 0（不跑 AVBD 的 capture_step_start）。
class BPXPBDStep : public XPBDStep {
public:
  static constexpr bool needs_topology = true;
  static constexpr bool uses_solver_buffers = false; // 不需要 per-body 局部系统
  static constexpr bool constraint_blocks = true;    // 走着色 block 扫掠

  Scalar_t zeta = Scalar_t(0.5);

  // 接触行的附加柔度（单位 m/N，即 1/K）。0 = 刚体接触，与发射路径一致
  // （ContactEmitOptions::compliance 默认 0）；>0 时用 α = 1/K 把 AVBD 的罚刚度
  // 映射成 XPBD 柔度，用来做"软/硬" A/B。它加在 ContactConstraint::compliance 上，
  // 不修改发射路径。
  Scalar_t contact_compliance = Scalar_t(0);

  // R = Σ λ∇C，每步/每迭代刷新。resize 只在 body 数变化时发生。
  std::vector<Vec3> R_lin;
  std::vector<Vec3> R_ang;

  // 本步是否已经清过接触的 λ/initial_c（每步第一次 begin_sweep 清一次）。
  bool contacts_initialized_ = false;

  // 继承 XPBD 的 begin_iterations（把 pose 搬到惯性目标），额外准备 R。
  template <typename Exec, typename Sched>
  void begin_iterations(State &s, Exec &exec, const Sched &sched) {
    XPBDStep::begin_iterations(s, exec, sched);
    resize_R(s.size());
    contacts_initialized_ = false;
  }

  // 每次迭代开始时调用一次：串行、按约束索引升序 gather R。
  // 这里同时把关节行的 XPBD 语义钉死：
  //   * initial 必须为 0（XPBD 不做 AVBD 的一阶稳定化）；
  //   * Limited 行在可行域内 λ = 0（互补松弛；λ 是乘子，可行方向一侧见 step）。
  void begin_sweep(State &s, ConstraintList &cs) {
    resize_R(s.size());
    for (auto &r : R_lin) r = Vec3{0};
    for (auto &r : R_ang) r = Vec3{0};

    // 接触 λ 的量纲是乘子（kg·m），AVBD 的 λ 是力（N）。两者不能互相暖启动，
    // 所以每步**第一次** begin_sweep 把接触的 λ/initial_c 清零；之后各迭代保留
    // 本步累积的乘子。manifold store 里的 AVBD λ 只能通过 reset_warm_start()
    // 在切换求解器时清掉（见 CollisionWorld）。
    if (!contacts_initialized_) {
      for (ConstraintVariant &v : cs) {
        if (auto *c = std::get_if<ContactConstraint>(&v)) {
          c->lambda        = Scalar_t(0);
          c->initial_c     = Scalar_t(0);
          c->initial_c_raw = Scalar_t(0);
        }
      }
      contacts_initialized_ = true;
    }

    for (ConstraintVariant &v : cs) {
      if (auto *j = std::get_if<JointAxisConstraint>(&v)) {
        j->initial = Scalar_t(0);
        if (j->mode == JointAxisMode::Limited) {
          const Scalar_t val = j->value(s);
          if (!(val > j->upper) && !(val < j->lower)) j->lambda = Scalar_t(0);
        }
      }
      std::visit([&](const auto &c) { gather_R(s, c); }, v);
    }
  }

  // 单行内核。接触点的 1~3 行必须由 step_unit() 的 unit 内核联合求解（法线夹紧
  // 决定摩擦锥上界），所以直接逐行调用这里时接触是 no-op；driver 的 block 扫掠
  // 总是走 step_unit()。
  template <Constraint<State> C> void step(State &s, C &c) {
    if constexpr (std::is_same_v<C, ContactConstraint>) {
      (void)s;
      (void)c;
      return;
    } else {
      const Scalar_t Cv     = c.evaluate(s);
      const auto     g      = c.gradient(s);
      const auto     bodies = c.bodies();

      if constexpr (HasAngularGradient<C, State>) {
        const auto ga = c.gradient_angular(s);

        Scalar_t denom = 0;
        for (std::size_t n = 0; n < g.size(); ++n)
          denom += s.inertial.effective_mass(bodies[n], g[n], ga[n]);

        const Scalar_t alpha_tilde = c.compliance / (s.dt * s.dt);
        const Scalar_t denom_all   = denom + alpha_tilde;
        if (!(denom_all > Scalar_t(0))) return; // 全静态且无 compliance：没有可动 DOF
        Scalar_t dlambda = (-Cv - alpha_tilde * c.lambda) / denom_all;
        apply_lambda(s, c, dlambda);

        const Scalar_t blend = Scalar_t(1) - zeta;
        for (std::size_t n = 0; n < g.size(); ++n) {
          const std::size_t b = bodies[n];
          if (s.inertial.is_static(b)) continue; // 无限质量 ⇒ 不写 pose，也不写 R

          Vec3 d_lin = g[n] * (s.inertial.inv_mass(b) * dlambda);
          Vec3 d_ang = (s.inertial.world_inv_inertia(b) * ga[n]) * dlambda;
          if (blend != Scalar_t(0)) {
            // B_i = (x̃_i - x_i) + M_i⁻¹R_i，角向用同一阶旋转向量表示。
            const Vec3 B_lin = (s.cache.inertial_position(b) - s.pose.position(b)) +
                               s.inertial.inv_mass(b) * R_lin[b];
            const Vec3 B_ang =
                    rotvec_between(s.pose.rotation(b), s.cache.inertial_rotation(b)) +
                    s.inertial.world_inv_inertia(b) * R_ang[b];
            d_lin += B_lin * blend;
            d_ang += B_ang * blend;
          }
          s.pose.apply_delta(b, d_lin, d_ang);
          R_lin[b] += g[n] * dlambda;
          R_ang[b] += ga[n] * dlambda;
        }
      } else {
        Scalar_t denom = 0;
        for (std::size_t n = 0; n < g.size(); ++n)
          denom += s.inertial.effective_mass(bodies[n], g[n], Vec3{0});

        const Scalar_t alpha_tilde = c.compliance / (s.dt * s.dt);
        const Scalar_t denom_all   = denom + alpha_tilde;
        if (!(denom_all > Scalar_t(0))) return;
        Scalar_t dlambda = (-Cv - alpha_tilde * c.lambda) / denom_all;
        apply_lambda(s, c, dlambda);

        const Scalar_t blend = Scalar_t(1) - zeta;
        for (std::size_t n = 0; n < g.size(); ++n) {
          const std::size_t b = bodies[n];
          if (s.inertial.is_static(b)) continue;

          Vec3 d_lin = g[n] * (s.inertial.inv_mass(b) * dlambda);
          if (blend != Scalar_t(0)) {
            const Vec3 B_lin = (s.cache.inertial_position(b) - s.pose.position(b)) +
                               s.inertial.inv_mass(b) * R_lin[b];
            d_lin += B_lin * blend;
          }
          s.pose.apply_delta(b, d_lin, Vec3{0});
          R_lin[b] += g[n] * dlambda;
        }
      }
    }
  }

  // —— unit 内核：driver 的 block 扫掠以一个 unit 为单位调用 ——
  //
  // 普通约束的 unit 只有一行，转给 step()（与 XPBD 同序同式）。接触点的 unit 是
  // 同一 (manifold_id, contact_index) 的 1~3 行，必须一起解：法线的单边夹紧决定
  // 切线的摩擦锥上界，逐行解会把切线拉出锥外。
  template <typename Ref>
  void step_unit(State &s, Ref *refs, std::uint32_t first, std::uint32_t rows) {
    if (refs[first].type_index == variant_index<ContactConstraint, ConstraintVariant>::value) {
      contact_unit(s, refs, first, rows);
      return;
    }
    for (std::uint32_t r = 0; r < rows; ++r) refs[first + r].step(s, *this);
  }

private:
  // 单个接触点（法线 + 两条切线）的联合 XPBD 行。
  //
  //   C        = ContactConstraint::value(s)      （**原始几何**，不掺 AVBD 稳定化）
  //   α̃        = compliance/dt²
  //   w        = Σ gᵀM⁻¹g（线性 + 角向，两侧 body）
  //   法线     Δλ = -(C_n + α̃λ_n)/(w_n + α̃)，λ_n = max(0, λ_n + Δλ)
  //   切线     Δλ = -(C_t + α̃λ_t)/(w_t + α̃)，然后 |λ_t| ≤ μλ_n 径向夹紧
  //   位置     x += M⁻¹g Δλ_clamped（夹紧后的增量），blend 项每个 body 每 unit 一次
  //
  // 单边性：C_n > 0 且 λ_n == 0 时整行跳过（分离、没有累积冲量）。这保证分离的
  // 接触绝不产生吸引力（AVBD 的 value 符号：C<0 需要沿 +n 推开 A）。
  //
  // R 的更新用 Δλ_clamped，且 blend 的 B_i = (x̃-x) + M⁻¹R 在**本 unit 开始前**的
  // R 上算一次（一个 unit 的三行是同一个整体；按行各加一次会把 blend 放大三倍）。
  template <typename Ref>
  void contact_unit(State &s, Ref *refs, std::uint32_t first, std::uint32_t rows) {
    ContactConstraint *row[3] = {nullptr, nullptr, nullptr};
    for (std::uint32_t r = 0; r < rows; ++r) {
      auto *c = static_cast<ContactConstraint *>(refs[first + r].obj);
      const auto k = static_cast<int>(c->kind);
      if (k >= 0 && k < 3 && row[k] == nullptr)
        row[k] = c;
    }
    if (row[0] == nullptr)
      return; // 没有法线行 ⇒ 没有锥可夹，切线没有意义

    const auto bodies = row[0]->bodies_;
    const Scalar_t inv_dt2 = Scalar_t(1) / (s.dt * s.dt);

    bool     have[3] = {false, false, false};
    Scalar_t C[3]    = {0, 0, 0};
    Scalar_t w[3]    = {0, 0, 0};
    Scalar_t alpha[3] = {0, 0, 0};
    Vec3     g[3][2];
    Vec3     ga[3][2];
    for (int k = 0; k < 3; ++k) {
      if (row[k] == nullptr)
        continue;
      have[k] = true;
      C[k]    = row[k]->value(s); // raw geometry: no alpha*initial_c stabilization
      const auto gl = row[k]->gradient(s);
      const auto glang = row[k]->gradient_angular(s);
      for (int n = 0; n < 2; ++n) {
        g[k][n]  = gl[n];
        ga[k][n] = glang[n];
      }
      w[k] = s.inertial.effective_mass(bodies[0], g[k][0], ga[k][0]) +
             s.inertial.effective_mass(bodies[1], g[k][1], ga[k][1]);
      alpha[k] = (row[k]->compliance + contact_compliance) * inv_dt2;
    }

    // —— 法线：单边 XPBD dual ——
    Scalar_t dlam[3] = {0, 0, 0};
    Scalar_t lam_n   = row[0]->lambda;
    if (!(C[0] > Scalar_t(0) && lam_n == Scalar_t(0))) {
      const Scalar_t denom = w[0] + alpha[0];
      if (denom > Scalar_t(0)) {
        const Scalar_t next = std::max(Scalar_t(0), lam_n - (C[0] + alpha[0] * lam_n) / denom);
        dlam[0]             = next - lam_n;
        row[0]->lambda      = next;
        lam_n               = next;
      }
    }

    // —— 切线：先各自做 XPBD dual，再按**新鲜**的 λ_n 做径向摩擦锥夹紧 ——
    Scalar_t cand[3] = {0, 0, 0};
    for (int k = 1; k < 3; ++k) {
      if (!have[k])
        continue;
      const Scalar_t old   = row[k]->lambda;
      const Scalar_t denom = w[k] + alpha[k];
      cand[k] = denom > Scalar_t(0) ? old - (C[k] + alpha[k] * old) / denom : old;
    }
    if (have[1] && have[2]) {
      const Scalar_t bound = std::max(Scalar_t(0), row[0]->friction) * lam_n;
      const Scalar_t tau   = std::sqrt(cand[1] * cand[1] + cand[2] * cand[2]);
      if (tau > bound && tau > Scalar_t(0)) {
        const Scalar_t scale = bound / tau;
        cand[1] *= scale;
        cand[2] *= scale;
      }
    } else {
      // μ == 0 只发一行法线；万一发射了单条切线，也要守锥。
      const Scalar_t bound = std::max(Scalar_t(0), row[0]->friction) * lam_n;
      for (int k = 1; k < 3; ++k)
        if (have[k]) cand[k] = std::clamp(cand[k], -bound, bound);
    }
    for (int k = 1; k < 3; ++k) {
      if (!have[k])
        continue;
      dlam[k]        = cand[k] - row[k]->lambda;
      row[k]->lambda = cand[k];
    }

    // —— 位置 + R：每个 body 在本 unit 内只写一次，blend 也只加一次 ——
    const Scalar_t blend = Scalar_t(1) - zeta;
    for (int n = 0; n < 2; ++n) {
      const std::size_t b = bodies[n];
      if (s.inertial.is_static(b))
        continue; // 静态 body 不写 pose，也不写 R（与 XPBD/距离内核一致）

      Vec3 d_lin{0}, d_ang{0}, r_lin{0}, r_ang{0};
      for (int k = 0; k < 3; ++k) {
        if (!have[k] || dlam[k] == Scalar_t(0))
          continue;
        d_lin += g[k][n] * (s.inertial.inv_mass(b) * dlam[k]);
        d_ang += (s.inertial.world_inv_inertia(b) * ga[k][n]) * dlam[k];
        r_lin += g[k][n] * dlam[k];
        r_ang += ga[k][n] * dlam[k];
      }
      if (blend != Scalar_t(0)) {
        const Vec3 B_lin = (s.cache.inertial_position(b) - s.pose.position(b)) +
                           s.inertial.inv_mass(b) * R_lin[b];
        const Vec3 B_ang =
                rotvec_between(s.pose.rotation(b), s.cache.inertial_rotation(b)) +
                s.inertial.world_inv_inertia(b) * R_ang[b];
        d_lin += B_lin * blend;
        d_ang += B_ang * blend;
      }
      s.pose.apply_delta(b, d_lin, d_ang);
      R_lin[b] += r_lin;
      R_ang[b] += r_ang;
    }
  }

  void resize_R(std::size_t n) {
    if (R_lin.size() != n) {
      R_lin.assign(n, Vec3{0});
      R_ang.assign(n, Vec3{0});
    }
  }

  // λ 是乘子。普通行直接累加；Limited 行夹到与 AVBD 力界相反的区间。
  // 注意：位置更新必须用**夹紧后**的 Δλ，否则约束会施加超过限位允许的力。
  template <typename C> void apply_lambda(State &s, C &c, Scalar_t &dlambda) {
    if constexpr (std::is_same_v<C, JointAxisConstraint>) {
      if (c.mode == JointAxisMode::Limited) {
        Scalar_t force_min = -kForceUnclamped;
        Scalar_t force_max = kForceUnclamped;
        c.force_bounds(s, force_min, force_max);
        const Scalar_t lo = -force_max;
        const Scalar_t hi = -force_min;
        const Scalar_t nl = std::clamp(c.lambda + dlambda, lo, hi);
        dlambda           = nl - c.lambda;
        c.lambda          = nl;
        return;
      }
    }
    c.lambda += dlambda;
  }

  // R_i += λ_c ∇C_c,i。串行、按 ConstraintList 顺序 ⇒ 与线程/调度无关。
  // 接触行的 λ 在步首已被 begin_sweep 清成乘子，这里和普通行同式累加；它把上一
  // 迭代的乘子重新纳入 R，使 R 始终等于"当前 λ 的 Σλ∇C"。
  template <typename C> void gather_R(const State &s, const C &c) {
    if constexpr (std::is_same_v<C, ContactConstraint>) {
      if (c.lambda == Scalar_t(0))
        return;
      const auto g      = c.gradient(s);
      const auto ga     = c.gradient_angular(s);
      const auto bodies = c.bodies();
      for (std::size_t n = 0; n < g.size(); ++n) {
        const std::size_t b = bodies[n];
        if (s.inertial.is_static(b))
          continue;
        R_lin[b] += g[n] * c.lambda;
        R_ang[b] += ga[n] * c.lambda;
      }
    } else {
      if constexpr (std::is_same_v<C, JointAxisConstraint>) {
        if (c.mode == JointAxisMode::Free) return;
      }
      if (c.lambda == Scalar_t(0)) return; // 0 乘任何梯度都是 0，跳过不改变结果

      const auto g      = c.gradient(s);
      const auto bodies = c.bodies();
      if constexpr (HasAngularGradient<C, State>) {
        const auto ga = c.gradient_angular(s);
        for (std::size_t n = 0; n < g.size(); ++n) {
          const std::size_t b = bodies[n];
          if (s.inertial.is_static(b)) continue;
          R_lin[b] += g[n] * c.lambda;
          R_ang[b] += ga[n] * c.lambda;
        }
      } else {
        for (std::size_t n = 0; n < g.size(); ++n) {
          const std::size_t b = bodies[n];
          if (s.inertial.is_static(b)) continue;
          R_lin[b] += g[n] * c.lambda;
        }
      }
    }
  }
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
  static constexpr bool needs_topology = false;
  static constexpr bool constraint_blocks = false;

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
