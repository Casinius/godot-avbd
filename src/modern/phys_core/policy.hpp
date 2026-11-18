#pragma once

#include <cstddef>
#include <span>

#include "base_concept.hpp"

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
  static constexpr Scalar_t kStiffness0 = Scalar_t(1e4);
  Scalar_t                  k           = kStiffness0;
  Scalar_t                  beta        = Scalar_t(1e-2);

  // 单约束：把 k·g gᵀ / k·C g 累加进 per-body 系统（只走 buffers 接口）
  template <Constraint<State> C> void step(State &s, C &c) {
    const Scalar_t Cv     = c.evaluate(s);
    auto           g      = c.gradient(s);
    auto           bodies = c.bodies();

    if constexpr (HasAngularGradient<C, State>) {
      auto g_ang = c.gradient_angular(s);
      for (std::size_t n = 0; n < g.size(); ++n)
        s.buffers.accumulate_penalty(bodies[n], g[n], g_ang[n], k, Cv);
    } else {
      for (std::size_t n = 0; n < g.size(); ++n)
        s.buffers.accumulate_penalty(bodies[n], g[n], Vec3{0}, k, Cv);
    }
  }

  // —— 阶段 1b：把单个 entry 折进局部系统 ——
  // gather 按 body 并行调用这里，每个调用只写 buffers[b]，且同 body 内的调用顺序
  // 由 body_entries（约束升序）决定 ⇒ 与串行逐条累加位精确。
  // 注意 k 只在 accumlate 阶段读，改写它的是 ramp_stiffness()。
  void accumulate_entry(State &s, std::size_t b, const ConstraintEntry &e) const noexcept {
    s.buffers.accumulate_penalty(b, e.g_lin, e.g_ang, k, e.C);
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
    k = kStiffness0;
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
