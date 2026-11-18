#pragma once

// ============================================================================
//  frame.hpp —— 一帧的流水线（task 层，建在求解器之上，不耦合进去）
// ----------------------------------------------------------------------------
//  【这一层做什么、不做什么】
//  做：声明一帧里有哪些 job、它们的读写集、谁必须等谁；把 job 派给 Executor。
//  不做：任何物理。求解只是其中一个 job，`solve()` 不知道也没必要知道自己被谁调用。
//  这就是"task 可以构建在其上而不耦合进去"的具体形态：求解器只认 Executor，
//  帧层只认 State 的命名接口，两边都不需要知道对方的调度细节。
//
//  【一帧的 DAG】
//
//     prepare ──▶ solve ──┬─▶ finish  ┐
//                         │           ├─▶ （同一一次 fork，写集不相交）
//                         └─▶ consume ┘
//
//  prepare：复位外力累加器（含"重力 = uniform_accel"）→ begin_step。
//  solve  ：iters 次迭代，内部四个相位（scatter/gather/primal/dual），已并行且位精确。
//  finish ：由位移反推速度 + 补齐 dq/rotMat —— **只写 velocity / cache**。
//  consume：任意只读消费者（渲染 AABB、蒙皮、网络同步）—— **只读 pose**。
//
//  finish 与 consume 的写集不相交、且互不读写对方的数据，所以可以塞进同一次
//  parallel_for：把两件事的参数空间拼成 [0, n_body + n_consume)，前半段跑 finish，
//  后半段跑 consume。这不是"看起来并行"，是真的在同一次 fork 里跑。
//
//  【为什么这个位置可以接更多东西】
//  真实引擎在这里还要挂 broadphase / narrowphase / CCD / sleep。它们的读写集是：
//    broadphase  ：读 pose          → 可与 finish/consume 同批
//    narrowphase ：读 pose          → 同上
//    CCD         ：读 pose+velocity → 必须排在 finish 之后
//  只要把"读集/写集"声明清楚，组合规则就只有一个：**写集相交 ⇒ 必须排序**。
//  求解器永远不需要知道这些。
// ============================================================================

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "phys_core/base_concept.hpp"
#include "phys_core/sched.hpp"
#include "phys_core/solver_phase.hpp"
#include "task.hpp"

namespace frame {

// 一帧的调度参数：帧层自己关心的东西（并行阈值/粒度）+ 求解参数。
struct FrameSchedule {
  SolveSchedule solve{};
  // finish/consume 这一批的粒度。它们每 body 的活很轻（几百 ns），
  // 所以阈值要比 solve 阶段更保守。
  std::size_t grain_tail     = kDefaultGrainBody;
  std::size_t min_items_tail = 128;
};

// 一帧的产出（供上层/诊断用；不参与物理）。
struct FrameStats {
  std::size_t bodies           = 0;
  std::size_t constraints      = 0;
  std::size_t tail_items       = 0; // finish + consume 合并后的参数空间大小
  bool        topology_rebuilt = false;
};

// 跑一帧。
//   consumer(b) 会在**求解完成之后**被调用，body 下标 b ∈ [0, s.size())。
//
//   【consumer 的硬契约 —— 违反它就是一个真数据竞争】
//     a) 只**读** State（不许写 pose/velocity/...）；
//     b) 只写**自己那一格**的输出（按 b 索引），不许累加进任何共享变量。
//   为什么：consumer 与 finish 挤在同一次 parallel_for 里，多个线程会同时跑
//   consumer(b)。我第一版就是让每个 b 去 merge 同一个 AABB，TSan 立刻报
//   data race（float std::min in SceneAabb::merge_point）。
//   正确地做聚合：先按 b 落到 points[b]，再在 fork **之外**做一次固定顺序归约
//   （见 SceneAabb::reduce）—— 这样既无竞争又保持确定性。
template <typename Policy, typename Exec, typename Observer, typename Consumer>
FrameStats step_frame(State &s, ConstraintList &cs, Policy &policy, Exec &exec,
                      const FrameSchedule &sched, int iters, std::uint64_t topology_version,
                      Observer *obs, Consumer *consumer) {
  // —— prepare ——
  s.reset_external(); // 外力累加器复位到 uniform_accel（重力只是它的默认值）
  s.begin_step();

  // —— solve ——（求解器眼中的世界到此为止）
  solve(s, cs, policy, iters, exec, sched.solve, obs, topology_version);

  FrameStats st;
  st.bodies           = s.size();
  st.constraints      = cs.size();
  st.topology_rebuilt = s.scratch.build_count != 0;

  // —— finish ∥ consume ——
  const std::size_t nb = s.size();
  const std::size_t nc = (consumer != nullptr) ? nb : 0;
  st.tail_items        = nb + nc;

  run_phase(exec, true, nb + nc, sched.grain_tail, sched.min_items_tail,
            [&s, consumer, nb](std::size_t b0, std::size_t b1) {
              for (std::size_t i = b0; i < b1; ++i) {
                if (i < nb) {
                  s.finish_body(i); // 只写 velocity / cache[b]
                } else {
                  (*consumer)(i - nb); // 只读 pose
                }
              }
            });

  s.cache.commit_prev(s.velocity);
  return st;
}

// 不带观察者的便捷重载（observer 无法从 nullptr 推导，所以单列一个）。
template <typename Policy, typename Exec, typename Consumer>
FrameStats step_frame(State &s, ConstraintList &cs, Policy &p, Exec &exec, const FrameSchedule &fs,
                      int iters, std::uint64_t topology_version, Consumer *consumer) {
  return step_frame<Policy, Exec, NullPhaseObserver, Consumer>(s, cs, p, exec, fs, iters,
                                                               topology_version, nullptr, consumer);
}

// ==================== 一个具体的只读消费者：场景 AABB ====================
// 渲染/裁剪/网络剔除都要它。它只读 pose，所以能和 finish 同批跑。
//
// 两段式，正是为了满足 consumer 契约：
//   store(b, p) —— 在融合 fork 里被调用，只写 points[b]（并行安全）；
//   reduce()    —— 在 fork **之外**按 b 升序归约（确定性；不并行）。
struct SceneAabb {
  std::vector<Vec3> points;
  Vec3              lo{0}, hi{0};
  bool              valid = false;

  void reset(std::size_t n) {
    points.assign(n, Vec3{0});
    lo    = Vec3{0};
    hi    = Vec3{0};
    valid = false;
  }

  // 融合 fork 内调用：只写自己那一格。
  void store(std::size_t b, Vec3 p) noexcept { points[b] = p; }

  // fork 外调用：固定顺序归约，结果与线程数无关。
  void reduce() {
    valid = false;
    for (Vec3 p : points) {
      if (!valid) {
        lo    = p;
        hi    = p;
        valid = true;
        continue;
      }
      lo.x = std::min(lo.x, p.x);
      lo.y = std::min(lo.y, p.y);
      lo.z = std::min(lo.z, p.z);
      hi.x = std::max(hi.x, p.x);
      hi.y = std::max(hi.y, p.y);
      hi.z = std::max(hi.z, p.z);
    }
  }
};

} // namespace frame
