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

  // 测试专用：强制"一 unit 一 block，block 升序"。这是 ζ=1 时与 XPBD 串行
  // Gauss-Seidel **逐位等价**的唯一分区形态（见 oracle 1），用来在固定的 pinned
  // 场景上验证 block 扫掠本身没有改变更新顺序/算术。生产路径保持 false。
  bool force_single_unit_blocks = false;
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

// —— 块分区的 body 可动性签名 ——
// block 的冲突判定会跳过静态 body，因此分区依赖"哪些 body 是静态的"。
// 静态性不在 TopologyKey 里（改质量不算拓扑变更），所以 block 缓存额外带一个
// 可动性签名：它变了就只重建 block 分区，entry 布局照常复用。
[[nodiscard]] inline std::uint64_t static_mask_hash(const State &s) noexcept {
  std::uint64_t h = 1469598103934665603ull;
  const auto *  p = s.inertial.mass_inv.data();
  for (std::size_t b = 0; b < s.inertial.mass_inv.size(); ++b) {
    const unsigned char v = p[b] == Scalar_t(0) ? 1u : 0u;
    h = (h ^ static_cast<std::uint64_t>(v)) * 1099511628211ull;
  }
  return h;
}

// —— 块分区（constraint_blocks 策略）——
// 与 entry 布局同属拓扑状态：只依赖约束→body 邻接与"哪两行是同一 unit"，
// 所以和 entry 布局一起挂在 TopologyKey 下缓存、一起在 verify 时重算比对。
//
// 规则（确定性）：
//   1) unit 按 ConstraintList 索引升序枚举；
//   2) 每个 unit first-fit 到**编号最小**的、所有动态 body 都空闲的 block；
//   3) 静态 body 不参与冲突判定 —— 否则常驻的 world body 会把所有约束塞进同一个
//      block，着色退化成完全串行。
// 结果只由输入顺序决定，没有任何哈希序/调度参与。
template <typename Ref>
void build_block_layout(ConstraintScratch &scr, std::span<const Ref> refs, const State &s,
                        bool force_single_unit = false) {
  scr.unit_begin.clear();
  scr.unit_rows.clear();
  scr.block_begin.clear();
  scr.block_units.clear();

  const std::size_t nc = refs.size();

  // 1) unit：相邻且属于同一 unit 的行（接触点的 1~3 行）合在一起。
  // 必须先比 type_index：joins_unit 是**按 refs[end] 的类型**强转两端的函数指针，
  // 跨类型调用会把 Distance/JointAxis 当成 ContactConstraint 读（UB）。
  std::size_t u = 0;
  while (u < nc) {
    std::size_t end = u + 1;
    while (end < nc && end - u < 3 && refs[end].type_index == refs[end - 1].type_index &&
           refs[end].joins_unit != nullptr &&
           refs[end].joins_unit(refs[end - 1].obj, refs[end].obj))
      ++end;
    scr.unit_begin.push_back(static_cast<std::uint32_t>(u));
    scr.unit_rows.push_back(static_cast<std::uint32_t>(end - u));
    u = end;
  }

  const std::size_t nu = scr.unit_rows.size();

  // 2) first-fit。block_units_tmp / block_bodies 都是拓扑重建期的临时量，
  //    不进缓存；块数通常很小（网格着色），复杂度不是热路径。
  std::vector<std::vector<std::uint32_t>> block_units_tmp;
  std::vector<std::vector<std::size_t>>   block_bodies;

  for (std::size_t ui = 0; ui < nu; ++ui) {
    if (force_single_unit) {
      block_units_tmp.emplace_back();
      block_units_tmp.back().push_back(static_cast<std::uint32_t>(ui));
      continue;
    }
    // 收集这个 unit 的动态 body（去重；一个 unit 最多 3 行 × 2 body）。
    std::array<std::size_t, 6> bodies{};
    std::size_t                nb = 0;
    const std::uint32_t        r0 = scr.unit_begin[ui];
    for (std::uint32_t r = 0; r < scr.unit_rows[ui]; ++r) {
      const Ref &ref = refs[r0 + r];
      for (std::size_t b : ref.bodies(ref.obj)) {
        if (s.inertial.is_static(b)) continue; // 静态 body 不制造冲突
        bool dup = false;
        for (std::size_t k = 0; k < nb; ++k) dup = dup || bodies[k] == b;
        if (!dup) bodies[nb++] = b;
      }
    }

    std::size_t chosen = block_units_tmp.size();
    for (std::size_t bi = 0; bi < block_units_tmp.size(); ++bi) {
      bool conflict = false;
      for (std::size_t k = 0; k < nb && !conflict; ++k)
        for (std::size_t x : block_bodies[bi]) conflict = conflict || x == bodies[k];
      if (!conflict) {
        chosen = bi;
        break;
      }
    }
    if (chosen == block_units_tmp.size()) {
      block_units_tmp.emplace_back();
      block_bodies.emplace_back();
    }
    block_units_tmp[chosen].push_back(static_cast<std::uint32_t>(ui));
    for (std::size_t k = 0; k < nb; ++k) block_bodies[chosen].push_back(bodies[k]);
  }

  // 3) 展平成 block_begin/block_units（块内 unit 索引保持升序）。
  scr.block_begin.resize(block_units_tmp.size() + 1);
  scr.block_begin[0] = 0;
  for (std::size_t bi = 0; bi < block_units_tmp.size(); ++bi)
    scr.block_begin[bi + 1] =
        scr.block_begin[bi] + static_cast<std::uint32_t>(block_units_tmp[bi].size());
  scr.block_units.resize(scr.block_begin.back());
  std::size_t w = 0;
  for (const auto &bl : block_units_tmp)
    for (std::uint32_t ui : bl) scr.block_units[w++] = ui;
}

// —— 拓扑缓存：拓扑不变就不重建 ——
// verify = true 时（测试用）会把布局重算一遍逐字节比对，把"拓扑没变"从“相信调用方”
// 变成“可验证”。stale_detected 一旦为真，就是 version 维护错了。
//
// entry 布局与 block 分区共用同一个 TopologyKey；block 是**按需**构建的
// （with_blocks = Policy::constraint_blocks），并用 block_valid 记录，避免
// "同一个 Scene 先跑 AVBD 再跑 BPXPBD"时命中 entry 缓存却没有 block 分区。
struct LayoutResult {
  bool rebuilt = false;
  bool stale_detected = false;
};

template <typename Ref>
LayoutResult ensure_layout(ConstraintScratch &scr, std::span<const Ref> refs, TopologyKey want,
                           bool verify, const State &s, bool with_blocks,
                           bool force_single_unit = false) {
  LayoutResult  out;
  std::uint64_t mass_hash = static_mask_hash(s);
  if (force_single_unit) mass_hash ^= 0x9E3779B97F4A7C15ull; // 分区模式进缓存签名

  if (scr.key.same_as(want)) {
    if (verify) {
      ConstraintScratch fresh;
      build_entry_layout(fresh, refs, want.body_count);
      out.stale_detected = fresh.begin != scr.begin || fresh.body_begin != scr.body_begin ||
                           fresh.body_entries != scr.body_entries;
      // 分区依赖静态性：只有签名一致时才能拿 fresh 和缓存比。签名变了会重建，
      // 不算"拓扑没 bump"的 stale。
      if (!out.stale_detected && with_blocks && scr.block_valid &&
          scr.block_mass_hash == mass_hash) {
        build_block_layout(fresh, refs, s, force_single_unit);
        out.stale_detected = fresh.unit_begin != scr.unit_begin ||
                             fresh.unit_rows != scr.unit_rows ||
                             fresh.block_begin != scr.block_begin ||
                             fresh.block_units != scr.block_units;
      }
    }
    if (with_blocks && (!scr.block_valid || scr.block_mass_hash != mass_hash)) {
      build_block_layout(scr, refs, s, force_single_unit);
      scr.block_mass_hash = mass_hash;
      scr.block_valid     = true;
    }
    return out; // 复用
  }

  build_entry_layout(scr, refs, want.body_count);
  if (with_blocks) {
    build_block_layout(scr, refs, s, force_single_unit);
    scr.block_mass_hash = mass_hash;
    scr.block_valid     = true;
  } else {
    scr.unit_begin.clear();
    scr.unit_rows.clear();
    scr.block_begin.clear();
    scr.block_units.clear();
    scr.block_mass_hash = 0;
    scr.block_valid     = false;
  }
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

// 关节行：与接触行相同的三条特例（自有 K / λ / 力界），外加两件事：
//   1) Hessian 同样做 block-Jacobi mass-split。二体行在 per-body 局部系统里会被两侧
//      各走一步，没有这项修正时相对位移是精确解的 2 倍，λ 在 RHS 里会把 2 倍放大成
//      每迭代 3 倍（与接触注释里推导的完全一样）。
//   2) 落在可行区间内的 Limited 行不提供任何刚度：hessian < 0 是"这一行没有 H"的
//      哨兵（见 accumulate_entry）。若照常给 K，局部系统会把轴钉在当前点，限位内的
//      自由度就动不了了。
template <>
inline void scatter_typed<JointAxisConstraint>(const State &s, const std::uint32_t *slots,
                                               const void *const *objs, ConstraintScratch &scr,
                                               std::size_t b, std::size_t e) noexcept {
  for (std::size_t i = b; i < e; ++i) {
    const auto &c = *static_cast<const JointAxisConstraint *>(objs[i]);
    const Scalar_t C = c.evaluate(s);
    const auto g = c.gradient(s);
    const auto ga = c.gradient_angular(s);
    Scalar_t force_min = -kForceUnclamped;
    Scalar_t force_max = kForceUnclamped;
    c.force_bounds(s, force_min, force_max);
    // Limited 且未越界 ⇒ C == 0 且没有力：既不提供刚度也不提供 RHS。
    const bool inactive_limit = c.mode == JointAxisMode::Limited && C == Scalar_t(0);
    const std::uint32_t slot = slots[i];
    for (std::size_t k = 0; k < 2; ++k) {
      ConstraintEntry &ent = scr.entries[slot + k];
      ent.g_lin = g[k];
      ent.g_ang = ga[k];
      ent.C = C;
      ent.penalty = c.penalty;
      ent.lambda = c.lambda;
      ent.force_min = force_min;
      ent.force_max = force_max;
      if (inactive_limit) {
        ent.hessian = Scalar_t(-1);
      } else {
        const std::size_t j = (k == 0) ? 1 : 0;
        ent.hessian = c.penalty * contact_hessian_scale(s, c.bodies_[k], c.bodies_[j], g[k], ga[k],
                                                        g[j], ga[j]);
      }
    }
  }
}
