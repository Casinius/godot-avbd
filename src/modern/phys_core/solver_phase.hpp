#pragma once
#include "constraint.hpp"
#include "sched.hpp"
// ==================== 10. 统一求解循环（driver） ====================
//
// 这里是"调度"的唯一归宿。四件事全在这里，Policy 一概不管：
//   1) 相位顺序（谁必须在谁之后，barrier 在哪）；
//   2) 每个相位能不能并行（SolveSchedule）；
//   3) 切块粒度；
//   4) 观察点/计时点。
// 想换调度器（线程池/TBB/单线程/串行回归）只需换 Exec；想验证并发正确性只需
// 用 SerialExecutor 跑同一条路径——两者位级结果必须一致。

// —— 相位词汇表首——
// 观察点用它们命名；以后做 per-phase 计时、per-phase 采样都挂在这里。
enum class SolvePhase { IterationBegin, Scatter, Gather, Primal, Dual, Ramp };

// 默认观察者：四个空函数，全部内联掉，串行路径不付代价。
struct NullPhaseObserver {
  void on_phase(SolvePhase, int, const State &) noexcept {}
};

template <typename Policy, typename Exec, typename Observer = NullPhaseObserver>
void solve(State &s, ConstraintList &cs, Policy &p, int iters, Exec &exec,
           const SolveSchedule &sched = {}, Observer *obs = nullptr,
           std::uint64_t topology_version = 0) {
  std::vector<ConstraintRef<Policy>> refs;
  refs.reserve(cs.size());
  for (auto &v : cs)
    std::visit([&](auto &c) { refs.push_back(ConstraintRef<Policy>::from(c)); }, v);

  // 每步一次：把策略特有的“准备”（如 XPBD/PBD 的显式位姿积分）放在迭代循环外。
  p.begin_iterations(s, exec, sched);

  // 按类型分批视图。依赖 scr.begin（拓扑），所以在 layout 之后建。
  ConstraintBatches<Policy> batches;

  // 拓扑缓存的唯一消费点。topology_version == 0 时每步重建（安全默认）；
  // 约束列表的拥有者声明“拓扑没变”才会真的复用。
  if constexpr (Policy::uses_solver_buffers) {
    const LayoutResult layout =
        ensure_layout(s.scratch, std::span<const ConstraintRef<Policy>>{refs},
                      TopologyKey{topology_version, refs.size(), s.size()}, sched.verify_topology);
    // stale_detected 意味着有人改了约束却没 bump version —— 错得很安静，必须暴露：
    // 计一次数（测试断言 stale_count == 0），并主动作废缓存，避免一步错、步步错。
    if (layout.stale_detected) {
      ++s.scratch.stale_count;
      s.scratch.key = TopologyKey{};
    }
    batches.build(std::span<const ConstraintRef<Policy>>{refs}, s.scratch);
  }

  // 步首一次：约束的稳定化基准（C0）与接触罚因子初值。必须在**任何迭代之前**，
  // 且必须在碰撞 sync 之后（接触锚点/基已经就位）——调用方保证这个顺序。
  if constexpr (requires { p.capture_step_start(s, cs); })
    p.capture_step_start(s, cs);

  for (int it = 0; it < iters; ++it) {
    if (obs) obs->on_phase(SolvePhase::IterationBegin, it, s); //阶段hook 用来观察的

    // —— 相位 1：累积 ——
    if constexpr (Policy::uses_solver_buffers) {
      // 块 Jacobi：先 scatter 到私有 slot，再按 body 折叠。两步都可并行。
      // 按类型分批 scatter：每批都是单态循环（见 9b 的说明）。
      // 分批只改写入顺序，不改每条约束写到哪一段 slot ⇒ 位精确性不受影响。
      for (auto &batch : batches.batches) {
        if (batch.size() == 0) continue;
        run_phase(exec, sched.parallel_accumulate, batch.size(), sched.grain_constraint,
                  sched.min_items_constraint,
                  [&s, &batch](std::size_t b, std::size_t e) {
                    batch.scatter(s, batch.slots.data(), batch.objs.data(), s.scratch, b, e);
                  });
      }
      if (obs) obs->on_phase(SolvePhase::Scatter, it, s);

      gather_entries(s, p, s.scratch, exec, sched);
      if (obs) obs->on_phase(SolvePhase::Gather, it, s);
    } else {
      // XPBD / PBD：step() 直接写 pose（Gauss-Seidel），**并行化会错**。
      // 要并行必须改算法（约束着色或 Jacobi 化），那是行为变更，不在本轮范围。
      for (auto &r : refs) r.step(s, p);
      if (obs) obs->on_phase(SolvePhase::Scatter, it, s);
    }

    // —— 相位 2：primal（逐 body；buffers 的清空也在这里逐 body 完成）——
    p.solve_primal(s, exec, sched);
    if (obs) obs->on_phase(SolvePhase::Primal, it, s);

    // —— 相位 3：dual（逐约束）——
    p.solve_dual(s, std::span{refs}, exec, sched);
    // 相位 3b：manifold 级 dual（串行）。接触的三行必须联合夹紧（单边 + 摩擦锥），
    // 逐行做不出来；接触没有 update_dual，所以与上面的并行 dual 不重叠。
    if constexpr (requires { p.post_dual(s, cs); })
      p.post_dual(s, cs);
    if (obs) obs->on_phase(SolvePhase::Dual, it, s);

    // —— 相位 4：刚度 ramp（串行标量，并行它只会更慢）——
    p.ramp_stiffness();
    if (obs) obs->on_phase(SolvePhase::Ramp, it, s);
  }
}
