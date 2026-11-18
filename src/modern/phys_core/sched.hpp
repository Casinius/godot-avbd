#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <tuple>
#include <type_traits>
#include <variant>

#include "constraint.hpp"
// 调度参数。**唯一**声明"什么能并行、切多大块"的地方。
// 默认值是"可证明位一致"的保守配置：一改这里的 flag，就必须用
// SerialExecutor 跑一遍基线验证结果没变。
struct SolveSchedule {
  // —— 相位开关（A/B 用；默认全开）——
  // parallel_accumulate 对 XPBD/PBD **无效**：那两个 Policy 的投影是 Gauss-Seidel，
  // 只能串行（见 solve() 的 else 分支）。
  bool parallel_accumulate = true;
  bool parallel_primal = true;
  bool parallel_dual = true;

  // —— 分区粒度 ——
  // grain_body 必须是 kDefaultGrainBody（16）的倍数，原因见 state.hpp 的伪共享推导：
  // 132 B/body × 16 = 33 条完整 cacheline，chunk 起点正好落在行边界上。
  std::size_t grain_body = kDefaultGrainBody;
  // 约束项更小（一个 entry 36 B），但逐约束的活也更轻，用 32 让每块有足够工作。
  std::size_t grain_constraint = 32;

  // —— 串行快速路径阈值 ——
  // 一次 fork-join 的固定成本约 1.5 us。item 数低于阈值时并行纯亏，
  // 所以直接串行跑。这两个数字由 `avbd-mo --sweep` 在真实场景上标定，
  // 不是拍脑袋：见 tests/README 与 step 7 的报告。
  std::size_t min_items_body = 64;
  std::size_t min_items_constraint = 96;

  // 拓扑缓存校验（测试用）：复用缓存时把布局重算一遍比对。
  // 开着会多一份 O(entries) 的临时分配，所以默认关。
  bool verify_topology = false;
};

// —— 唯一的“该不该并行”判定 ——
// driver 与 Policy 共用这一处规则，避免"有人忘了判阈值"导致小场景白白 fork。
template <typename Exec, typename F>
void run_phase(Exec &exec, bool allow_parallel, std::size_t n, std::size_t grain,
               std::size_t min_items, F &&f) noexcept {
  if (allow_parallel && n >= min_items) {
    exec.parallel_for(n, grain, f);
  } else {
    f(std::size_t(0), n); // 串行快速路径：与并行路径同代码、同顺序 ⇒ 位一致
  }
}

// ==================== 9. scatter / gather：消除 accumulate 的写冲突 ====================
//
// 问题：串行版把每条约束的贡献直接 `buffers.accumulate_penalty(b, ...)` 进去，
// 这是对同一个 body 的 read-modify-write。同一 body 被 k 条约束碰，并行就必错。
//
// 解法（方案 A：scatter-then-gather）：
//   1) scatter（可并行）：约束 c 把 (C, g_lin, g_ang) 写进 entries[begin[c]..) —— 私有 slot；
//   2) gather（可并行）：按 body 把 slot 折进 buffers，加总顺序 = 约束索引升序。
// 第 2 步的顺序是**数据布局**保证的，不是线程调度保证的 ⇒ 与串行位精确。
// 比“约束着色 + 锁/原子”好在：无原子、无锁、且不只正确而且确定。

// entries 的临时梯度缓冲上界由 ConstraintVariant 自动推导，不写死 ——
// 以后加一种 4-body 约束时，这里自动跟上，不会静默截断梯度。
template <typename V> struct max_gradient_size;
template <typename... Ts>
struct max_gradient_size<std::variant<Ts...>>
    : std::integral_constant<std::size_t,
                             std::max({std::tuple_size_v<constraint_gradient_t<Ts>>...})> {};

inline constexpr std::size_t kMaxConstraintBodies = max_gradient_size<ConstraintVariant>::value;

// —— 布局 + body 邻接 ——
// 遍历约束必须**升序**：body_entries 的排列就是 gather 的加总顺序。
// 全部用 resize/assign（不是 assign(n, 0) 新建）—— 缓存复用后容量不再增长，
// 每步不会有 malloc；这是 step 5 真正的收益（O(n) 重建成本本身只有 accumulate 的 1/iters）。
template <typename Ref>
void build_entry_layout(ConstraintScratch &scr, std::span<const Ref> refs, std::size_t body_count) {
  const std::size_t nc = refs.size();

  scr.begin.resize(nc + 1);
  scr.begin[0] = 0;
  for (std::size_t c = 0; c < nc; ++c)
    scr.begin[c + 1] = scr.begin[c] + static_cast<std::uint32_t>(refs[c].bodies(refs[c].obj).size());
  scr.entries.resize(scr.begin[nc]);

  scr.body_begin.resize(body_count + 1);
  for (std::size_t b = 0; b <= body_count; ++b) scr.body_begin[b] = 0;

  // 第一遍：统计每个 body 有多少 entry   counting sort
  for (std::size_t c = 0; c < nc; ++c)
    for (std::size_t b : refs[c].bodies(refs[c].obj)) ++scr.body_begin[b + 1];
  //
  //前缀和 → body_begin[b] 是 body b 的桶起点
  for (std::size_t b = 0; b < body_count; ++b) scr.body_begin[b + 1] += scr.body_begin[b];

  scr.body_entries.resize(scr.body_begin[body_count]);
  scr.cursor.resize(body_count);
  // 第二遍：按约束索引升序扫描，把 entry 塞进各自 body 的桶
  for (std::size_t b = 0; b < body_count; ++b) scr.cursor[b] = scr.body_begin[b];
  for (std::size_t c = 0; c < nc; ++c) {
    auto bodies = refs[c].bodies(refs[c].obj);
    for (std::size_t i = 0; i < bodies.size(); ++i)
      scr.body_entries[scr.cursor[bodies[i]]++] = scr.begin[c] + static_cast<std::uint32_t>(i);
  }
}

// —— 拓扑缓存：拓扑不变就不重建 ——
// verify = true 时（测试用）会把布局重算一遍逐字节比对，把"拓扑没变"从“相信调用方”
// 变成“可验证”。stale_detected 一旦为真，就是 version 维护错了。
struct LayoutResult {
  bool rebuilt = false;
  bool stale_detected = false;
};

template <typename Ref>
LayoutResult ensure_layout(ConstraintScratch &scr, std::span<const Ref> refs, TopologyKey want,
                           bool verify) {
  LayoutResult out;
  if (scr.key.same_as(want)) {
    if (verify) {
      ConstraintScratch fresh;
      build_entry_layout(fresh, refs, want.body_count);
      out.stale_detected = fresh.begin != scr.begin || fresh.body_begin != scr.body_begin ||
                           fresh.body_entries != scr.body_entries;
    }
    return out; // 复用
  }

  build_entry_layout(scr, refs, want.body_count);
  scr.key = want;
  ++scr.build_count;
  out.rebuilt = true;
  return out;
}

// ==================== 9b. 按类型分批：把间接调用从“逐约束”降到“逐 (类型, 分块)” ====================
// 为什么需要它：ConstraintRef 的 gradient/evaluate 是函数指针，**逐条约束**一次间接调用，
// 编译器永远看不到循环体内部，既不能内联也不能向量化。
// 把约束按类型分组之后，每组的循环体是单态的（`const DistanceConstraint&`），
// evaluate/gradient 全部内联 —— 这才是“typed pool”真正换来的东西。
//
// 位精确性不受影响：分批只改变**写入顺序**，而每条约束写的还是它自己那一段 slot。
// 真正决定浮点结果的是 gather 的加总顺序，那由 body_entries 决定，与分批无关。

// 没有角向自由度的约束类型，补一个同尺寸的全零梯度，让下面的循环保持单态。
template <typename C>
constexpr auto angular_gradient_of(const C &c, const State &s) noexcept {
  if constexpr (HasAngularGradient<C, State>) {
    return c.gradient_angular(s);
  } else {
    std::array<Vec3, std::tuple_size_v<constraint_gradient_t<C>>> z{};
    for (auto &v : z) v = Vec3{0};
    return z;
  }
}

template <typename C>
void scatter_typed(const State &s, const std::uint32_t *slots, const void *const *objs,
                   ConstraintScratch &scr, std::size_t b, std::size_t e) noexcept {
  for (std::size_t i = b; i < e; ++i) {
    const auto &c = *static_cast<const C *>(objs[i]);
    const Scalar_t Cv = c.evaluate(s);          // ← 单态，可内联
    const auto g = c.gradient(s);               // ← 单态，可内联/向量化
    const auto ga = angular_gradient_of(c, s);
    constexpr std::size_t n = std::tuple_size_v<decltype(g)>;
    const std::uint32_t slot = slots[i];
    for (std::size_t k = 0; k < n; ++k) {
      ConstraintEntry &ent = scr.entries[slot + k];
      ent.g_lin = g[k];
      ent.g_ang = ga[k];
      ent.C = Cv;
      // 非接触行没有自己的罚因子 / 对偶力：0 是"用 Policy 的全局 k、λ=0"的哨兵。
      ent.penalty = Scalar_t(0);
      ent.lambda  = Scalar_t(0);
      ent.hessian = Scalar_t(0);
      ent.force_min = -kForceUnclamped;
      ent.force_max = kForceUnclamped;
    }
  }
}

template <typename V> struct BatchDispatch;
template <typename... Ts> struct BatchDispatch<std::variant<Ts...>> {
  template <typename Batch> static void assign(std::vector<Batch> &batches) {
    std::size_t i = 0;
    ((batches[i++].scatter = &scatter_typed<Ts>), ...);
  }
};

// 每个 variant alternative 一个 Batch。成员顺序 = 约束的全局索引升序。
template <typename Policy> struct ConstraintBatches {
  struct Batch {
    void (*scatter)(const State &, const std::uint32_t *, const void *const *,
                    ConstraintScratch &, std::size_t, std::size_t) = nullptr;
    std::vector<std::uint32_t> slots; // 每个成员约束的 entry 起始槽位
    std::vector<const void *> objs;
    std::size_t size() const noexcept { return objs.size(); }
  };
  std::vector<Batch> batches;

  ConstraintBatches() {
    batches.resize(std::variant_size_v<ConstraintVariant>);
    BatchDispatch<ConstraintVariant>::assign(batches);
  }

  // 按 refs 重建（保持每类内部按全局索引升序）。
  // 只做 O(nc) 的 push_back，容量复用，不分配新内存。
  template <typename Ref> void build(std::span<const Ref> refs, const ConstraintScratch &scr) {
    for (auto &b : batches) {
      b.slots.clear();
      b.objs.clear();
    }
    for (std::size_t c = 0; c < refs.size(); ++c) {
      Batch &b = batches[refs[c].type_index];
      b.slots.push_back(scr.begin[c]);
      b.objs.push_back(refs[c].obj);
    }
  }
};

// —— scatter 的一个连续区间 ——
// 只调 evaluate/gradient，**不写 pose**，所以对任何 Policy 都安全可并行。
template <typename Ref>
void scatter_range(const State &s, std::span<const Ref> refs, ConstraintScratch &scr, std::size_t c0,
                   std::size_t c1) noexcept {
  std::array<Vec3, kMaxConstraintBodies> g_lin{};
  std::array<Vec3, kMaxConstraintBodies> g_ang{};
  for (std::size_t c = c0; c < c1; ++c) {
    const auto &r = refs[c];
    const std::size_t n = std::min(r.bodies(r.obj).size(), kMaxConstraintBodies);
    const Scalar_t Cv = r.evaluate(r.obj, s);
    r.gradient(r.obj, s, std::span<Vec3>{g_lin.data(), n});
    if (r.gradient_angular) {
      r.gradient_angular(r.obj, s, std::span<Vec3>{g_ang.data(), n});
    } else {
      for (std::size_t i = 0; i < n; ++i) g_ang[i] = Vec3{0};
    }
    for (std::size_t i = 0; i < n; ++i) {
      ConstraintEntry &e = scr.entries[scr.begin[c] + i];
      e.g_lin = g_lin[i];
      e.g_ang = g_ang[i];
      e.C = Cv;
      e.penalty = Scalar_t(0); // 通用路径：用 Policy 的全局 k、λ = 0
      e.lambda  = Scalar_t(0);
      e.hessian = Scalar_t(0);
      e.force_min = -kForceUnclamped;
      e.force_max = kForceUnclamped;
    }
  }
}

// —— gather 的一个连续 body 区间 ——
// 每个 body 只写自己的 buffers[b]，按 body_entries 的（约束升序）依次累加。
template <typename Policy, typename Exec>
void gather_entries(State &s, Policy &p, const ConstraintScratch &scr, Exec &exec,
                    const SolveSchedule &sched) noexcept {
  auto body_range = [&s, &p, &scr](std::size_t b0, std::size_t b1) {
    for (std::size_t b = b0; b < b1; ++b)
      for (std::uint32_t k = scr.body_begin[b]; k < scr.body_begin[b + 1]; ++k)
        p.accumulate_entry(s, b, scr.entries[scr.body_entries[k]]);
  };
  run_phase(exec, sched.parallel_accumulate, s.size(), sched.grain_body, sched.min_items_body,
            body_range);
}


//scatter_typed<ContactConstraint>
// 接触行的三个特例：
//   1) C 是**稳定化后**的 C_eff（ContactConstraint::evaluate 已减去 alpha*initial_c）；
//   2) 每行带自己的罚因子 K 与上一轮夹紧后的力 λ，一起进 entry（gather 直接用）；
//   3) Hessian 权重按 body 做 block-Jacobi mass-split（见 contact_hessian_scale），
//      与 AVBDStep::step 的直写路径必须给出同一组数字。
template <>
inline void scatter_typed<ContactConstraint>(const State &s, const std::uint32_t *slots,
                                      const void *const *objs, ConstraintScratch &scr,
                                      std::size_t b, std::size_t e) noexcept {
  for (std::size_t i = b; i < e; ++i) {
    const auto &c = *static_cast<const ContactConstraint *>(objs[i]);
    const Scalar_t C = c.evaluate(s);           // 单态，内联
    const auto g = c.gradient(s);
    const auto ga = c.gradient_angular(s);
    const std::uint32_t slot = slots[i];
    for (std::size_t k = 0; k < 2; ++k) {
      ConstraintEntry &ent = scr.entries[slot + k];
      ent.g_lin = g[k];
      ent.g_ang = ga[k];
      ent.C = C;
      ent.penalty = c.penalty;
      ent.lambda = c.lambda;
      ent.force_min = c.force_min;
      ent.force_max = c.force_max;
      const std::size_t j = (k == 0) ? 1 : 0;
      ent.hessian = c.penalty * contact_hessian_scale(s, c.bodies_[k], c.bodies_[j], g[k], ga[k],
                                                      g[j], ga[j]);
    }
  }
}
