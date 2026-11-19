// ============================================================================
//  main.cpp —— 回归基线 + 结构自检 + 基准驱动
// ----------------------------------------------------------------------------
//  用法：
//    avbd-mo                  跑全部用例，stdout 打印规范化摘要行，stderr 打印耗时
//    avbd-mo --check FILE     与录制好的基线逐行比对，任何一位不同即返回 1
//    avbd-mo --selftest       跑结构不变量断言（不比对数值），失败返回 1
//    avbd-mo --trace          额外把每次迭代的三个阶段摘要打到 stderr
//    avbd-mo --save FILE      把当前摘要写成基线（只在确认行为变更时用）
//    avbd-mo --collision-digest
//                             打印 collision 层可观测摘要（不跑求解器用例）
//    avbd-mo --check-collision FILE
//                             对比 collision 可观测基线；变化时返回 1
//    avbd-mo --save-collision FILE
//                             记录 collision 可观测基线（算法变更时重新录）
//
//  为什么要有 --selftest（摘要之外还要断言）：
//    摘要只能回答"和上次是不是一样"，不能回答"这样对不对"。并发重构最容易犯的
//    错，恰恰是让摘要稳定地错下去（例如 Acc 阶段漏刷新 dq、batch 后漏 reset）。
//    所以这里把 state.hpp 顶部那几条 invariants 逐条变成可执行的断言。
//
//  退出码：0 = 全部通过 / 无检查请求；1 = 检查失败；2 = 参数或 IO 错误。
// ============================================================================

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <limits>
#include <print>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

#include "digest.hpp"
#include "frame.hpp"
#include "phys_core/constraint.hpp"
#include "phys_core/joint.hpp"
#include "phys_core/policy.hpp"
#include "phys_core/sched.hpp"
#include "phys_core/solver_phase.hpp"
#include "scene.hpp"
#include "task.hpp"
#include "collision_selftest.hpp"

namespace {

// —— 阶段追踪器：把每次迭代的每个相位结束后的摘要记下来 ——
// 只用于定位漂移发生在哪个相位，不参与基线比对。
struct PhaseTrace {
  const ConstraintList        *cs = nullptr;
  Vec3                         consumer_lo{0}; // 只读消费者在融合 fork 里看到的第一个 AABB 角
  static constexpr std::size_t kPhases = 6;    // IterationBegin/Scatter/Gather/Primal/Dual/Ramp
  std::vector<std::uint64_t>   trace[kPhases];

  void on_phase(SolvePhase ph, int, const State &s) {
    trace[static_cast<std::size_t>(ph)].push_back(digest::parts(s, *cs).all);
  }

  void dump(std::string_view case_name) const {
    static const char *kNames[kPhases] = {"iter", "scatter", "gather", "primal", "dual", "ramp"};
    for (std::size_t p = 0; p < kPhases; ++p)
      for (std::size_t it = 0; it < trace[p].size(); ++it)
        std::println(stderr, "  trace {} it={:3d} {:>7}={}", case_name, it, kNames[p],
                     digest::to_hex(trace[p][it]));
  }
};

// 一次完整的"多步"执行：与真实驱动的调用顺序一致。
// 注意 reset_external() 必须在 begin_step() 之前：begin_step 只消费外力不复位。
// 想要风力/爆炸力，就在这两行之间 add_force()，与重力走同一条路。
template <typename Policy, typename Exec>
void drive(scene::Scene &sc, Policy policy, Exec &exec, const SolveSchedule &sched, int steps,
           int iters, PhaseTrace *tracer = nullptr) {
  State               &s = sc.state;
  frame::FrameSchedule fs{};
  fs.solve = sched;

  // 只读消费者：渲染侧的 AABB。它只读 pose，所以被 frame.hpp 融进同一次 fork。
  frame::SceneAabb aabb;

  // finish ∥ consume：consumer 只读 pose，且只写 points[b]（并行安全），
  // 归约放到 fork 之外做 —— 见 frame.hpp 里 consumer 契约那一段。
  auto consume = [&aabb, &s](std::size_t b) { aabb.store(b, s.pose.position(b)); };

  for (int step = 0; step < steps; ++step) {
    aabb.reset(s.size());
    frame::step_frame(s, sc.constraints, policy, exec, fs, iters, sc.topology_version,
                      step == 0 ? tracer : nullptr, &consume);
    aabb.reduce();
    if (tracer && step == 0 && aabb.valid)
      tracer->consumer_lo = aabb.lo;
  }
}

struct Result {
  std::string line;     // stdout 的规范化行
  double      ms = 0.0; // stderr 的耗时
};

template <typename Exec, typename Policy>
Result run_case(Exec &exec, const SolveSchedule &sched, bool trace, std::string_view case_name,
                std::string_view policy_name, scene::Scene sc, Policy policy, int steps,
                int iters) {
  PhaseTrace tracer;
  tracer.cs = &sc.constraints;

  const auto t0 = std::chrono::steady_clock::now();
  drive(sc, policy, exec, sched, steps, iters, trace ? &tracer : nullptr);
  const auto t1 = std::chrono::steady_clock::now();

  if (trace)
    tracer.dump(case_name);

  const digest::Parts d = digest::parts(sc.state, sc.constraints);
  Result              r;
  r.ms   = std::chrono::duration<double, std::milli>(t1 - t0).count();
  r.line = std::format(
      "{} {} steps={} iters={} all={} pose={} vel={} inert={} ext={} cache={} buf={} lam={}",
      case_name, policy_name, steps, iters, digest::to_hex(d.all), digest::to_hex(d.pose),
      digest::to_hex(d.velocity), digest::to_hex(d.inertial), digest::to_hex(d.external),
      digest::to_hex(d.cache), digest::to_hex(d.buffers), digest::to_hex(d.constraints));
  return r;
}

// —— 用例表：顺序即基线文件的行序 ——
// 覆盖目标（缺一不可）：
//   chain x3 policy  : 三个 Policy 的更新公式确实互不相同（除 XPBD≡PBD，见自检）
//   chain@rest       : 静止金丝雀（strain==1 ⇒ C≡0，任何位移都是 bug）
//   fan x2 policy    : 单个静态 body 被 n 条约束累加 —— accumulate 写冲突最点
//   lattice 1D/3D    : 共享密度更高，且是后续扩展性测试的载体
//   mixed x2 policy  : Distance + Contact 两种类型同在热路径上
//
// 整张表是一个**以 Exec 为模板参数**的函数：这样同一份用例既能在串行执行器上跑
// （参考语义），也能在线程池上跑，再逐位比对 —— 并发正确性的验证不需要第二套代码。
template <typename Exec>
std::vector<Result> run_all(Exec &exec, const SolveSchedule &sched, bool trace,
                            std::uint64_t topology_version_override = 1) {
  auto run = [&](std::string_view case_name, std::string_view policy_name, scene::Scene sc,
                 auto policy, int steps, int iters) {
    sc.topology_version = topology_version_override;
    return run_case(exec, sched, trace, case_name, policy_name, std::move(sc), policy, steps,
                    iters);
  };

  std::vector<Result> out;
  out.push_back(run("chain/32", "avbd", scene::chain(32), AVBDStep{}, 4, 50));
  out.push_back(run("chain/32", "xpbd", scene::chain(32), XPBDStep{}, 4, 50));
  out.push_back(run("chain/32", "pbd", scene::chain(32), PBDStep{}, 4, 50));
  {
    // strain==1 且关掉重力：每个约束的 C 恒为 0，全文应该一动不动。
    scene::Scene sc = scene::chain(32, Scalar_t(0.5), Scalar_t(1.0), false);
    scene::make_gravity_free(sc);
    out.push_back(run("chain/32@rest+nograv", "avbd", std::move(sc), AVBDStep{}, 4, 50));
  }
  out.push_back(run("freefall/4", "avbd", scene::freefall(4), AVBDStep{}, 8, 10));
  out.push_back(run("freefall/4", "xpbd", scene::freefall(4), XPBDStep{}, 8, 10));
  out.push_back(run("fan/64", "avbd", scene::fan(64), AVBDStep{}, 3, 30));
  out.push_back(run("fan/64", "xpbd", scene::fan(64), XPBDStep{}, 3, 30));
  out.push_back(run("lattice/5x4x1", "avbd", scene::lattice(5, 4, 1), AVBDStep{}, 3, 30));
  out.push_back(run("lattice/4x3x2", "avbd", scene::lattice(4, 3, 2), AVBDStep{}, 3, 30));
  out.push_back(run("mixed/20", "avbd", scene::mixed(20), AVBDStep{}, 3, 30));
  out.push_back(run("mixed/20", "xpbd", scene::mixed(20), XPBDStep{}, 3, 30));
  // 缩放样本：480 body / ~2100 约束。也是后续 step 6/7/11 的测量载体。
  out.push_back(run("lattice/12x10x4", "avbd", scene::lattice(12, 10, 4), AVBDStep{}, 3, 30));
  return out;
}

// ==================== 结构自检 ====================
// 自检关心的是语义（是不是位一致、有没有凭空位移），不是吞吐，所以统一用串行
// 执行器 —— 它同时就是"参考实现"。并发的正确性在 check_parallel_* 里单独验。
inline task::SerialExecutor g_serial;
inline const SolveSchedule  g_sched{};

struct Check {
  std::string name;
  bool        ok = true;
  std::string detail;
};

inline bool bitwise_eq(Vec3 a, Vec3 b) noexcept { return std::memcmp(&a, &b, sizeof(Vec3)) == 0; }
inline bool is_zero(Vec3 v) noexcept { return v.x == 0.0f && v.y == 0.0f && v.z == 0.0f; }
inline bool is_zero(const Mat3 &m) noexcept {
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r)
      if (m[c][r] != 0.0f)
        return false;
  return true;
}

// invariant 4：buffers 是"每步累积、每步清空"；batch_update 消费完就必须清空。
// 并发化之后这条最容易破（重复累加 / 漏清 / 竞态留下的残值），所以单独一条断言。
Check check_buffers_cleared_after_solve() {
  Check        c{"buffers reset after batch_update", true, ""};
  scene::Scene sc = scene::lattice(4, 3, 2);
  drive(sc, AVBDStep{}, g_serial, g_sched, 1, 10);
  const auto &buf = sc.state.buffers;
  for (std::size_t b = 0; b < buf.size(); ++b) {
    if (!is_zero(buf.hessian_lin[b]) || !is_zero(buf.hessian_lin_ang[b]) ||
        !is_zero(buf.hessian_ang[b]) || !is_zero(buf.grad_lin[b]) || !is_zero(buf.grad_ang[b])) {
      c.ok     = false;
      c.detail = std::format("body {} not cleared", b);
      break;
    }
  }
  return c;
}

// invariant 2：dq 必须等于 pose 相对 initial 的差。任何写了 pose 却忘记
// refresh_dq 的地方都会在这里被抓到（并发重构的头号回归来源）。
Check check_dq_consistent() {
  Check        c{"dq == pose - initial (linear & angular)", true, ""};
  scene::Scene sc = scene::chain(32);
  drive(sc, AVBDStep{}, g_serial, g_sched, 2, 20);
  const auto &s = sc.state;
  for (std::size_t b = 0; b < s.size(); ++b) {
    const Vec3 expect_lin = s.pose.position(b) - s.cache.initial_position(b);
    const Vec3 expect_ang = rotvec_between(s.cache.initial_rotation(b), s.pose.rotation(b));
    if (!bitwise_eq(expect_lin, s.cache.dq_linear(b)) ||
        !bitwise_eq(expect_ang, s.cache.dq_angular(b))) {
      c.ok     = false;
      c.detail = std::format("body {} dq mismatch", b);
      break;
    }
  }
  return c;
}

// invariant 5：静态 body 永不移动、速度恒为 0、逆惯量为 0。
Check check_static_bodies_frozen() {
  Check        c{"static bodies never move", true, ""};
  scene::Scene sc = scene::fan(64);
  // fan 的 body 0 是静态枢纽，且被全部约束累加。
  std::vector<Vec3> before(sc.state.size());
  for (std::size_t b = 0; b < sc.state.size(); ++b)
    before[b] = sc.state.pose.position(b);

  drive(sc, AVBDStep{}, g_serial, g_sched, 2, 30);

  for (std::size_t b = 0; b < sc.state.size(); ++b) {
    if (!sc.state.inertial.is_static(b))
      continue;
    if (!bitwise_eq(before[b], sc.state.pose.position(b))) {
      c.ok     = false;
      c.detail = std::format("static body {} moved", b);
      break;
    }
    if (!is_zero(sc.state.velocity.velocity(b)) ||
        !is_zero(sc.state.velocity.angular_velocity(b))) {
      c.ok     = false;
      c.detail = std::format("static body {} has nonzero velocity", b);
      break;
    }
  }
  return c;
}

// 静止金丝雀（零重力）：C 初始恒为 0 且无外力 ⇒ 全文不动、λ 全 0。
// 它存在的价值是防止并发化引入"凭空产生位移"的 bug。
Check check_gravity_free_rest_scene_is_frozen() {
  Check        c{"gravity-free strain==1 scene is a true no-op", true, ""};
  scene::Scene sc = scene::chain(32, Scalar_t(0.5), Scalar_t(1.0), false);
  scene::make_quiescent(sc);
  const auto pose_before = digest::pose(sc.state.pose);

  drive(sc, AVBDStep{}, g_serial, g_sched, 2, 50);

  if (digest::pose(sc.state.pose) != pose_before) {
    c.ok     = false;
    c.detail = "pose changed with C==0 and no gravity";
    return c;
  }
  for (std::size_t b = 0; b < sc.state.size(); ++b) {
    if (!is_zero(sc.state.velocity.velocity(b)) ||
        !is_zero(sc.state.velocity.angular_velocity(b))) {
      c.ok     = false;
      c.detail = std::format("body {} gained velocity from nothing", b);
      return c;
    }
  }
  for (const auto &v : sc.constraints) {
    const bool nonzero = std::visit([](const auto &cc) { return cc.lambda != 0.0f; }, v);
    if (nonzero) {
      c.ok     = false;
      c.detail = "lambda became nonzero with C==0";
      return c;
    }
  }
  return c;
}

// 惯性与重力的端到端验收：半隐式积分下 x_n = x0 + g·dt²·n(n+1)/2，v_n = n·g·dt。
// 这条断言一旦通过，就说明"外力→速度→惯性目标→solve_body 残差→位姿→速度"
// 整条链路都是通的；任何一环断掉这里都会立刻翻车。
//
// 关于容差（不是拍脑袋）：模拟是逐步 float 累加，与一次性解析式必然有尾差。
// 误差源头是坐标量级 |y|≈10 处的 ulp，而速度 = Δx/dt 会把位置误差放大 1/dt 倍。
// 所以容差写成"32 ulp 的位置误差，再除以 dt"，而不是写一个魔法常数。
Check check_free_fall_matches_analytic() {
  Check          c{"free fall matches semi-implicit analytic solution", true, ""};
  scene::Scene   sc = scene::freefall(4, Scalar_t(10));
  const int      n  = 8;
  const Vec3     x0 = sc.state.pose.position(0);
  const Scalar_t g  = sc.state.external.uniform_accel.y;
  const Scalar_t dt = sc.state.dt;

  drive(sc, AVBDStep{}, g_serial, g_sched, n, 10); // 无约束，迭代数无关

  const Scalar_t expect_v = g * dt * Scalar_t(n);
  const Scalar_t expect_y = x0.y + g * dt * dt * Scalar_t(n * (n + 1)) * Scalar_t(0.5);

  constexpr Scalar_t kEps  = std::numeric_limits<Scalar_t>::epsilon();
  const Scalar_t     ulp_y = kEps * std::max(std::abs(expect_y), Scalar_t(1));
  const Scalar_t     tol_y = Scalar_t(32) * ulp_y;
  const Scalar_t     tol_v = tol_y / dt + Scalar_t(32) * kEps * std::abs(expect_v);

  for (std::size_t b = 0; b < sc.state.size(); ++b) {
    const Vec3 p = sc.state.pose.position(b);
    const Vec3 v = sc.state.velocity.velocity(b);
    if (std::abs(p.y - expect_y) > tol_y || std::abs(v.y - expect_v) > tol_v) {
      c.ok     = false;
      c.detail = std::format("body {}: y={} (want {}+-{}), vy={} (want {}+-{})", b, p.y, expect_y,
                             tol_y, v.y, expect_v, tol_v);
      return c;
    }
    if (std::abs(p.x - x0.x - Scalar_t(b) * Scalar_t(0.3)) > 1e-6f || p.z != 0.0f) {
      c.ok     = false;
      c.detail = std::format("body {} drifted sideways", b);
      return c;
    }
  }
  return c;
}

// 直接验收用户提出的要求："重力不应该是特殊内部力"。
// uniform_accel = g 与"对每个动态 body 调 add_force(m·g)"必须产生**逐位相同**的结果。
// 两个 mass_inv==1 的相等让 f·m⁻¹ 与 g 恰好同位，因此这里可以要求严格相等。
Check check_gravity_is_not_special() {
  Check c{"gravity == add_force(m*g), no special path", true, ""};

  auto run = [](bool express_as_force) {
    scene::Scene sc = scene::fan(16);
    State       &s  = sc.state;
    const Vec3   g  = s.external.uniform_accel;
    if (express_as_force)
      s.external.uniform_accel = Vec3{0};

    AVBDStep policy;
    for (int step = 0; step < 3; ++step) {
      s.reset_external();
      if (express_as_force)
        for (std::size_t b = 0; b < s.size(); ++b) {
          if (s.inertial.is_static(b))
            continue;
          const Scalar_t mass = Scalar_t(1) / s.inertial.inv_mass(b);
          s.add_force(b, g * mass);
        }
      s.begin_step();
      solve(s, sc.constraints, policy, 20, g_serial, g_sched);
      s.finish_step();
      s.cache.commit_prev(s.velocity);
    }
    return std::tuple{digest::pose(s.pose), digest::velocity(s.velocity), digest::cache(s.cache)};
  };

  const auto [pos_u, vel_u, cache_u] = run(false);
  const auto [pos_f, vel_f, cache_f] = run(true);
  // 只比物理结果，不比 external 字段本身 —— 两种写法本来就是不同的 uniform_accel，
  // 把它们算进比较只会让这条断言恒假，钓不到任何真问题。
  if (pos_u != pos_f || vel_u != vel_f || cache_u != cache_f) {
    c.ok = false;
    c.detail =
        std::format("pose {} vs {}, vel {} vs {} — gravity IS special", digest::to_hex(pos_u),
                    digest::to_hex(pos_f), digest::to_hex(vel_u), digest::to_hex(vel_f));
  }
  return c;
}

// XPBD 在 compliance==0 时与 PBD 位置更新代数等价（dlambda = -C/denom 完全相同）。
// 这不是巧合，而是并发重构时必须保住的"公式级"事实：任何把 XPBD 改成
// Jacobi/block 形式的重写，都会先在这条断言上现形。
Check check_xpbd_equals_pbd_when_rigid() {
  Check        c{"XPBD(compliance=0) pose == PBD pose", true, ""};
  scene::Scene a = scene::chain(32);
  scene::Scene b = scene::chain(32);
  drive(a, XPBDStep{}, g_serial, g_sched, 2, 30);
  drive(b, PBDStep{}, g_serial, g_sched, 2, 30);
  const auto pa = digest::pose(a.state.pose);
  const auto pb = digest::pose(b.state.pose);
  if (pa != pb) {
    c.ok     = false;
    c.detail = std::format("{} != {}", digest::to_hex(pa), digest::to_hex(pb));
  }
  return c;
}

// ==================== B-PXPBD（blended primal XPBD）====================
// 这些断言覆盖 Phase 0–2 的三条硬契约：
//   1) ζ = 1 必须逐位退化为 XPBD 的串行 Gauss-Seidel；
//   2) ζ = 0 的 fixed-point 修正必须让真实残差单调不增，且优于 ζ = 1；
//   3) 着色 block 扫掠必须对线程数/粒度/串行阈值不敏感；
//   4) 静止链条不许出现由 stale λ/initial 引起的"幻影拉力"。

// 单步内的残差探针：每次迭代结束（Ramp 之后）记一次
//     ‖ M(x - x̃) - Σ λ∇C(x) ‖       （对所有动态 body 取平方和再开方）
// 这是 XPBD/B-PXPBD 的乘子平衡式：M(x-x̃) = Σ λ∇C（λ 的量纲是 kg·m，力 = λ∇C/dt²）。
// 注意不要再除 dt²：那会把乘子项的量级弄错 3600 倍，残差看起来"越迭代越大"。
struct ResidualProbe {
  const ConstraintList *cs = nullptr;
  std::vector<Scalar_t> residual;

  void on_phase(SolvePhase ph, int, const State &s) {
    if (ph != SolvePhase::Ramp || cs == nullptr) return;
    std::vector<Vec3> R(s.size(), Vec3{0});
    for (const auto &v : *cs) {
      std::visit(
          [&](const auto &c) {
            if constexpr (std::is_same_v<std::decay_t<decltype(c)>, ContactConstraint>) {
              return; // Phase 3
            }
            if (c.lambda == Scalar_t(0)) return;
            const auto g      = c.gradient(s);
            const auto bodies = c.bodies();
            for (std::size_t n = 0; n < g.size(); ++n)
              if (!s.inertial.is_static(bodies[n])) R[bodies[n]] += g[n] * c.lambda;
          },
          v);
    }
    Scalar_t       sum     = 0;
    for (std::size_t b = 0; b < s.size(); ++b) {
      if (s.inertial.is_static(b)) continue;
      const Vec3 d = (s.pose.position(b) - s.cache.inertial_position(b)) *
                         s.inertial.inv_mass(b) - R[b];
      sum += glm::dot(d, d);
    }
    residual.push_back(std::sqrt(sum));
  }
};

// 所有动态 body 共享一个可动 hub 的星形场景。每条约束都碰 hub ⇒ first-fit 必然
// 给每个 unit 一个独立 block，且 block 顺序 = 约束顺序（正是 oracle 1 的前提）。
scene::Scene hub_star(std::size_t n) {
  scene::Scene sc;
  sc.name     = "hub_star";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(n + 1);
  scene::add_body(sc.state, 0, Vec3{0, 0, 0}, 1.0f, Vec3{0.1f, 0, 0}, Vec3{0, 0, 0.05f});
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3 p{Scalar_t(1) + Scalar_t(0.1) * Scalar_t(i), Scalar_t(0.2) * Scalar_t(i % 3),
                 Scalar_t(0.05) * Scalar_t(i % 5)};
    scene::add_body(sc.state, i + 1, p, 0.0f); // 静态叶子：只有 hub 是动态 body
  }
  for (std::size_t i = 0; i < n; ++i)
    sc.constraints.emplace_back(
        DistanceConstraint{{0, i + 1}, Scalar_t(0.9) * scene::jitter(i)});
  return sc;
}

// 一行 JointAxis 的最小场景：body 0 是静态世界，body 1 是动态滑块。
scene::Scene joint_slider(JointAxisMode mode, Scalar_t x, Scalar_t lower, Scalar_t upper) {
  scene::Scene sc;
  sc.name     = "joint_slider";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(2);
  scene::add_body(sc.state, 0, Vec3{0, 0, 0}, 0.0f);
  scene::add_body(sc.state, 1, Vec3{x, 0, 0}, 1.0f);
  scene::make_gravity_free(sc); // 只让关节轴上的自由度参与

  JointAxisConstraint row;
  row.bodies_ = {1, 0};
  row.axis    = Vec3{1, 0, 0};
  row.angular = false;
  row.mode    = mode;
  row.target  = Scalar_t(0.5);
  row.lower   = lower;
  row.upper   = upper;
  sc.constraints.emplace_back(row);
  return sc;
}

// 角向 JointAxis 的最小场景：body 0 静态世界 + body 1 动态，绕 X 轴 Locked。
// 关节建立后再把 body 1 转 0.2 rad ⇒ C=0.2。用 Joint::rebuild 生成 frame/rest/
// angular_modes，避免手写四元数把测试写错。
scene::Scene angular_joint() {
  scene::Scene sc;
  sc.name     = "angular_joint";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(2);
  scene::add_body(sc.state, 0, Vec3{0, 0, 0}, 0.0f);
  scene::add_body(sc.state, 1, Vec3{0, 0, 0}, 1.0f);
  scene::make_gravity_free(sc);

  JointDesc desc;
  desc.body_a          = 0;
  desc.body_b          = 1;
  desc.angular[0].mode = JointAxisMode::Locked;
  Joint joint(desc, sc.state);
  for (const auto &row : joint.rows())
    if (row.angular) sc.constraints.emplace_back(row);

  // rest 已在 rebuild 时捕获为"当前相对姿态"；现在制造一个角向违反。
  sc.state.pose.set_rotation(1, quat_from_rotvec(Vec3{0.2f, 0, 0}));
  return sc;
}

// oracle 1a：在"一 unit 一 block、block 升序"的分区下，ζ = 1 必须与 XPBDStep 的
// 串行 GS 逐位相同。这里用 SolveSchedule 的测试开关强制该分区，覆盖 pinned 的
// chain/fan（Distance-only 且 XPBD 有有限解）：
// ζ = 1 时 blend 项被显式跳过，行内核的 denom/Δλ/apply_delta 与 XPBD 同序同式；
// 强制分区只是把全局约束顺序原样搬进 block 顺序，不改变任何算术。
Check check_bpxpbd_zeta_one_equals_xpbd() {
  Check c{"BPXPBD(zeta=1) == XPBD (single-unit blocks, pinned scenes)", true, ""};

  SolveSchedule sched  = g_sched;
  sched.force_single_unit_blocks = true;

  auto run = [&](auto policy, scene::Scene sc) {
    drive(sc, policy, g_serial, sched, 4, 50);
    return digest::parts(sc.state, sc.constraints).all;
  };

  BPXPBDStep bpxpbd;
  bpxpbd.zeta = Scalar_t(1);

  const std::pair<const char *, scene::Scene> cases[] = {
      {"chain/32", scene::chain(32)},
      {"fan/64", scene::fan(64)},
  };
  // 注：lattice 不在等价场景里 —— 它含有**两端都是静态**的 Distance 行
  // （底面层的 x/y 约束），XPBDStep 在 w==0 且 compliance==0 时除零得到 NaN，
  // 参考实现本身没有有限值。BPXPBD 对全静态行的处理是显式跳过（见 row kernel 的
  // denom 守卫），这是有意的健壮性差异，不属于"block 扫掠改变顺序"。
  for (const auto &[name, sc] : cases) {
    const auto a = run(XPBDStep{}, sc);
    const auto b = run(bpxpbd, sc);
    if (a != b) {
      c.ok     = false;
      c.detail = std::format("{}: xpbd={} bpxpbd={}", name, digest::to_hex(a), digest::to_hex(b));
      return c;
    }
  }
  return c;
}

// oracle 1b：显式构造"一 unit 一 block、block 升序"的分区，再比一次。
// 这条同时检查 first-fit 的分配结果本身（block_units[k] == k）。
Check check_bpxpbd_one_unit_per_block_equals_xpbd() {
  Check c{"BPXPBD(zeta=1) == XPBD with one unit per block", true, ""};

  scene::Scene sc_bp = hub_star(24);
  BPXPBDStep   bpxpbd;
  bpxpbd.zeta = Scalar_t(1);
  drive(sc_bp, bpxpbd, g_serial, g_sched, 3, 30);

  const ConstraintScratch &scr = sc_bp.state.scratch;
  const std::size_t        nb  = scr.block_begin.empty() ? 0 : scr.block_begin.size() - 1;
  if (nb != sc_bp.constraints.size()) {
    c.ok     = false;
    c.detail = std::format("{} blocks for {} units (want one per unit)", nb,
                           sc_bp.constraints.size());
    return c;
  }
  for (std::size_t b = 0; b < nb; ++b) {
    if (scr.block_begin[b] != b || scr.block_units[b] != b) {
      c.ok     = false;
      c.detail = std::format("block {} is not unit {} (partition not ascending)", b,
                             scr.block_units.empty() ? 0 : scr.block_units[b]);
      return c;
    }
  }

  scene::Scene sc_x = hub_star(24);
  drive(sc_x, XPBDStep{}, g_serial, g_sched, 3, 30);

  const auto a = digest::parts(sc_x.state, sc_x.constraints).all;
  const auto b = digest::parts(sc_bp.state, sc_bp.constraints).all;
  if (a != b) {
    c.ok     = false;
    c.detail = std::format("xpbd={} bpxpbd={}", digest::to_hex(a), digest::to_hex(b));
  }
  return c;
}

// oracle 2：ζ = 0 的 blend 必须让 fixed-point 残差整体单调下降（允许 float 尾差），
// 且明显优于 ζ = 1。
Check check_bpxpbd_blend_reduces_residual() {
  Check c{"BPXPBD(zeta=0) residual decreases and beats zeta=1", true, ""};

  // 用有 compliance 的链条：B-PXPBD 的 blend 项是为柔性约束设计的。刚度取 3e-4
  // 让 ζ=0 的 fixed-point 映射保持收缩，同时 residual 仍远高于 float 噪声底。
  auto probes = [](Scalar_t zeta) {
    scene::Scene sc = scene::chain(12, Scalar_t(0.5), Scalar_t(0.85), true);
    for (auto &v : sc.constraints) std::get<DistanceConstraint>(v).compliance = Scalar_t(3e-4);

    BPXPBDStep    p;
    p.zeta = zeta;
    ResidualProbe probe;
    probe.cs = &sc.constraints;

    sc.state.reset_external();
    sc.state.begin_step();
    solve(sc.state, sc.constraints, p, 10, g_serial, g_sched, &probe, sc.topology_version);
    sc.state.finish_step();
    return probe.residual;
  };

  const auto r_blend = probes(Scalar_t(0));
  const auto r_xpbd  = probes(Scalar_t(1));
  if (r_blend.size() != 10 || r_xpbd.size() != 10) {
    c.ok     = false;
    c.detail = std::format("expected 10 samples, got {} / {}", r_blend.size(), r_xpbd.size());
    return c;
  }

  // ζ=0：整体单调下降。允许一处 float 尾差回升（实测 ~1.35×），阈值放到 1.5×。
  for (std::size_t i = 1; i < r_blend.size(); ++i) {
    if (r_blend[i] > r_blend[i - 1] * Scalar_t(1.5)) {
      c.ok     = false;
      c.detail = std::format("residual grew at iter {}: {:.3e} -> {:.3e}", i, r_blend[i - 1],
                             r_blend[i]);
      return c;
    }
  }
  if (!(r_blend.back() < r_blend.front() * Scalar_t(0.1))) {
    c.ok     = false;
    c.detail = std::format("residual did not drop by 10x: {:.3e} -> {:.3e}", r_blend.front(),
                           r_blend.back());
    return c;
  }
  // ζ=0 必须优于 ζ=1（ζ=1 在这条链条上停滞在 ~1e-4）。
  if (!(r_blend.back() <= r_xpbd.back())) {
    c.ok     = false;
    c.detail = std::format("zeta=0 residual {:.3e} > zeta=1 residual {:.3e}", r_blend.back(),
                           r_xpbd.back());
    return c;
  }
  return c;
}

// 拓扑缓存校验：BPXPBD 走的是 block 分区，verify_topology 必须同时覆盖
// entry 邻接与 block 分区，且缓存要真的被复用（build_count == 1）。
Check check_bpxpbd_verify_topology() {
  Check c{"BPXPBD topology cache verifies clean and is reused", true, ""};

  SolveSchedule strict{};
  strict.verify_topology = true;

  scene::Scene sc = scene::lattice(4, 3, 2);
  drive(sc, BPXPBDStep{}, g_serial, strict, 3, 20);

  if (sc.state.scratch.stale_count != 0) {
    c.ok     = false;
    c.detail = std::format("stale_count={} (block verifier disagreed with the cache)",
                           sc.state.scratch.stale_count);
    return c;
  }
  if (sc.state.scratch.build_count != 1) {
    c.ok     = false;
    c.detail = std::format("build_count={} after 3 steps (want 1)", sc.state.scratch.build_count);
    return c;
  }
  if (!sc.state.scratch.block_valid) {
    c.ok     = false;
    c.detail = "block partition is not marked valid after a BPXPBD step";
    return c;
  }
  return c;
}

// 定义在接触用例一节；这里提前声明，供确定性检查复用。
scene::Scene contact_fan_scene(std::size_t n, Scalar_t friction);

// oracle 3：BPXPBD 的着色扫掠对线程数、grain、min_items 逐位不敏感。
Check check_bpxpbd_determinism() {  Check c{"BPXPBD digests invariant across threads x grain x thresholds", true, ""};

  auto run = [](auto &exec, const SolveSchedule &sched, scene::Scene sc) {
    BPXPBDStep p; // 默认 zeta = 0.5，blend 项真的启用
    drive(sc, p, exec, sched, 3, 30);
    return digest::parts(sc.state, sc.constraints).all;
  };

  const auto ref_chain   = run(g_serial, SolveSchedule{}, scene::chain(32));
  const auto ref_lattice = run(g_serial, SolveSchedule{}, scene::lattice(4, 3, 2));
  // 接触 unit 的并行路径：一个 block 内的 unit 只写自己的 body ⇒ 与线程数无关。
  const auto ref_contact = run(g_serial, SolveSchedule{}, contact_fan_scene(12, Scalar_t(0.5)));

  for (unsigned threads : {1u, 8u}) {
    for (std::size_t grain : {std::size_t(1), std::size_t(7), std::size_t(64)}) {
      SolveSchedule cfg{};
      cfg.grain_body           = grain;
      cfg.grain_constraint     = grain;
      cfg.min_items_body       = 1;
      cfg.min_items_constraint = 1;

      task::ThreadPoolExecutor pool(threads);
      if (run(pool, cfg, scene::chain(32)) != ref_chain) {
        c.ok     = false;
        c.detail = std::format("chain diverged: threads={} grain={}", threads, grain);
        return c;
      }
      if (run(pool, cfg, scene::lattice(4, 3, 2)) != ref_lattice) {
        c.ok     = false;
        c.detail = std::format("lattice diverged: threads={} grain={}", threads, grain);
        return c;
      }
      if (run(pool, cfg, contact_fan_scene(12, Scalar_t(0.5))) != ref_contact) {
        c.ok     = false;
        c.detail = std::format("contact_fan diverged: threads={} grain={}", threads, grain);
        return c;
      }
    }
  }
  return c;
}

// oracle 4：静止 + 无力 ⇒ BPXPBD 必须一动不动（不允许 stale λ / initial 造成幻影拉力）。
Check check_bpxpbd_rest_chain_is_frozen() {
  Check        c{"BPXPBD gravity-free rest chain is a true no-op", true, ""};
  scene::Scene sc = scene::chain(32, Scalar_t(0.5), Scalar_t(1.0), false);
  scene::make_quiescent(sc);
  const auto pose_before = digest::pose(sc.state.pose);

  BPXPBDStep p;
  drive(sc, p, g_serial, g_sched, 2, 50);

  if (digest::pose(sc.state.pose) != pose_before) {
    c.ok     = false;
    c.detail = "pose changed with C==0 and no gravity";
    return c;
  }
  for (std::size_t b = 0; b < sc.state.size(); ++b) {
    if (!is_zero(sc.state.velocity.velocity(b)) ||
        !is_zero(sc.state.velocity.angular_velocity(b))) {
      c.ok     = false;
      c.detail = std::format("body {} gained velocity from nothing", b);
      return c;
    }
  }
  for (const auto &v : sc.constraints) {
    if (std::visit([](const auto &cc) { return cc.lambda != 0.0f; }, v)) {
      c.ok     = false;
      c.detail = "lambda became nonzero with C==0";
      return c;
    }
  }
  return c;
}

// Phase 2：JointAxis 的 XPBD 语义。
//   * Locked：λ 是乘子，更新式与 XPBDStep 完全一致 ⇒ ζ=1 时逐位相同；
//   * Limited：λ 夹到可行域一侧；在限位内 λ=0、无位移。
Check check_bpxpbd_joint_axis_locked_equals_xpbd() {
  Check c{"BPXPBD(zeta=1) == XPBD on a Locked joint axis", true, ""};

  scene::Scene a = joint_slider(JointAxisMode::Locked, Scalar_t(1), Scalar_t(0), Scalar_t(0));
  scene::Scene b = joint_slider(JointAxisMode::Locked, Scalar_t(1), Scalar_t(0), Scalar_t(0));

  drive(a, XPBDStep{}, g_serial, g_sched, 3, 20);
  BPXPBDStep p;
  p.zeta = Scalar_t(1);
  drive(b, p, g_serial, g_sched, 3, 20);

  const auto da = digest::parts(a.state, a.constraints).all;
  const auto db = digest::parts(b.state, b.constraints).all;
  if (da != db) {
    c.ok     = false;
    c.detail = std::format("xpbd={} bpxpbd={}", digest::to_hex(da), digest::to_hex(db));
    return c;
  }
  const Scalar_t x = b.state.pose.position(1).x;
  if (std::abs(x - Scalar_t(0.5)) > Scalar_t(1e-4)) {
    c.ok     = false;
    c.detail = std::format("locked joint settled at x={} (want 0.5)", x);
  }
  return c;
}

// 角向行与线性行走同一条内核；这里额外覆盖 gradient_angular / 角向 apply_delta。
Check check_bpxpbd_angular_joint_locked_equals_xpbd() {
  Check c{"BPXPBD(zeta=1) == XPBD on an angular Locked joint", true, ""};

  scene::Scene a = angular_joint();
  scene::Scene b = angular_joint();

  drive(a, XPBDStep{}, g_serial, g_sched, 3, 20);
  BPXPBDStep p;
  p.zeta = Scalar_t(1);
  drive(b, p, g_serial, g_sched, 3, 20);

  const auto da = digest::parts(a.state, a.constraints).all;
  const auto db = digest::parts(b.state, b.constraints).all;
  if (da != db) {
    c.ok     = false;
    c.detail = std::format("xpbd={} bpxpbd={}", digest::to_hex(da), digest::to_hex(db));
    return c;
  }
  const auto &row = std::get<JointAxisConstraint>(b.constraints.front());
  if (std::abs(row.value(b.state)) > Scalar_t(1e-4)) {
    c.ok     = false;
    c.detail = std::format("angular joint settled at {} rad (want 0)", row.value(b.state));
  }
  return c;
}

Check check_bpxpbd_joint_axis_limited() {  Check c{"BPXPBD Limited joint stays in range and has no phantom lambda", true, ""};

  // 限位内：C==0 ⇒ λ 必须保持 0，且 pose 完全不动。
  scene::Scene inside = joint_slider(JointAxisMode::Limited, Scalar_t(0.1), Scalar_t(0),
                                     Scalar_t(0.2));
  const Vec3   inside_before = inside.state.pose.position(1);
  BPXPBDStep   p;
  drive(inside, p, g_serial, g_sched, 2, 20);
  if (!bitwise_eq(inside.state.pose.position(1), inside_before)) {
    c.ok     = false;
    c.detail = "Limited joint moved while inside its limits";
    return c;
  }
  for (const auto &v : inside.constraints) {
    const auto &j = std::get<JointAxisConstraint>(v);
    if (j.lambda != Scalar_t(0) || j.initial != Scalar_t(0)) {
      c.ok     = false;
      c.detail = std::format("inside limits: lambda={} initial={}", j.lambda, j.initial);
      return c;
    }
  }

  // 越上限：一步位置投影把它压回 x=upper，λ 必须是**负**的乘子（AVBD 的力是 F=K·C+λ，
  // 与 XPBD 乘子反号）。用 1 次迭代直接观察夹紧后的 λ；多步会因为 v=Δx/dt 的惯性
  // 把 body 甩过下界，反而测不到夹紧符号。
  scene::Scene outside = joint_slider(JointAxisMode::Limited, Scalar_t(0.25), Scalar_t(0),
                                      Scalar_t(0.2));
  drive(outside, p, g_serial, g_sched, 1, 1);
  const Scalar_t x = outside.state.pose.position(1).x;
  if (std::abs(x - Scalar_t(0.2)) > Scalar_t(1e-5)) {
    c.ok     = false;
    c.detail = std::format("Limited joint ended at x={} (want upper=0.2)", x);
    return c;
  }
  for (const auto &v : outside.constraints) {
    const auto &j = std::get<JointAxisConstraint>(v);
    if (!(j.lambda < Scalar_t(0))) {
      c.ok     = false;
      c.detail = std::format("Limited lambda={} (want < 0 for an upper-limit violation)",
                             j.lambda);
      return c;
    }
  }
  return c;
}

// ==================== B-PXPBD 接触行（Phase 3） ====================
// 这些用例直接合成 ContactConstraint（不跑窄相），把锚点放在两 body 的界面上。
// penetration > 0 是几何穿透（C_n = -penetration + margin < 0），< 0 是间隙。
void add_synthetic_contact(scene::Scene &sc, std::size_t body_a, std::size_t body_b,
                           std::uint32_t manifold_id, Vec3 normal, Scalar_t penetration,
                           Scalar_t friction, Vec3 slip = Vec3{0}) {
  normal = glm::normalize(normal);
  const Vec3 helper = std::abs(normal.y) < Scalar_t(0.9) ? Vec3{0, 1, 0} : Vec3{1, 0, 0};
  const Vec3 t0     = glm::normalize(glm::cross(normal, helper));
  const Vec3 t1     = glm::normalize(glm::cross(normal, t0));

  const Vec3 pA  = sc.state.pose.position(body_a);
  const Vec3 pB  = sc.state.pose.position(body_b);
  const Vec3 mid = (pA + pB) * Scalar_t(0.5);
  const Vec3 xA  = mid - normal * (penetration * Scalar_t(0.5)) + slip;
  const Vec3 xB  = mid + normal * (penetration * Scalar_t(0.5));
  const Vec3 rA  = transpose(sc.state.pose.rotation_matrix(body_a)) * (xA - pA);
  const Vec3 rB  = transpose(sc.state.pose.rotation_matrix(body_b)) * (xB - pB);

  for (int k = 0; k < 3; ++k) {
    ContactConstraint c;
    c.bodies_       = {body_a, body_b};
    c.manifold_id   = manifold_id;
    c.contact_index = 0;
    c.kind          = static_cast<ContactKind>(k);
    c.basis[0]      = normal;
    c.basis[1]      = t0;
    c.basis[2]      = t1;
    c.rA            = rA;
    c.rB            = rB;
    c.target        = (k == 0) ? -Scalar_t(0.005) : Scalar_t(0);
    c.friction      = friction;
    c.compliance    = Scalar_t(0);
    c.lambda        = Scalar_t(0);
    c.penalty       = Scalar_t(0);
    c.initial_c     = Scalar_t(0);
    c.initial_c_raw = Scalar_t(0);
    sc.constraints.emplace_back(c);
  }
}

// A 在 B 上方 1 m，锚点落在界面 y=0.5；法线 +y 指向 A。重力可选。
scene::Scene contact_pair_scene(Scalar_t penetration, Scalar_t friction, Vec3 slip,
                                bool gravity) {
  scene::Scene sc;
  sc.name     = "contact_pair";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(2);
  scene::add_body(sc.state, 0, Vec3{0, 0, 0}, 0.0f); // B: static support
  scene::add_body(sc.state, 1, Vec3{0, 1, 0}, 1.0f); // A: dynamic, above B
  add_synthetic_contact(sc, 1, 0, 0, Vec3{0, 1, 0}, penetration, friction, slip);
  if (!gravity)
    scene::make_quiescent(sc);
  return sc;
}

// n 个独立的动态 body 各自与同一个静态支撑接触：所有 unit 只共享静态 body，
// 所以自然分区把它们放进同一个 block，且块内顺序 = 约束顺序。这正是一 unit
// 一 block 的"串行顺序"，两者必须逐位相同。
scene::Scene contact_fan_scene(std::size_t n, Scalar_t friction) {
  scene::Scene sc;
  sc.name     = "contact_fan";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(n + 1);
  scene::add_body(sc.state, 0, Vec3{0, 0, 0}, 0.0f); // static support
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3 p{Scalar_t(0.5) * Scalar_t(i), 1, Scalar_t(0.25) * Scalar_t(i % 5)};
    scene::add_body(sc.state, i + 1, p, 1.0f);
  }
  for (std::size_t i = 0; i < n; ++i) {
    const Vec3 slip{Scalar_t(0.01) * Scalar_t(i % 3), 0, 0};
    add_synthetic_contact(sc, i + 1, 0, static_cast<std::uint32_t>(i), Vec3{0, 1, 0},
                          Scalar_t(0.01), friction, slip);
  }
  return sc;
}

// Phase 3a：分离的接触（C_n > 0，λ_n == 0）不得产生任何吸引力，法线行必须整行跳过。
Check check_bpxpbd_contact_unilateral_no_pull() {
  Check c{"BPXPBD separated contact exerts no pull", true, ""};

  scene::Scene sc = contact_pair_scene(Scalar_t(-0.1), Scalar_t(0.5), Vec3{0}, false);
  const Vec3   before = sc.state.pose.position(1);
  BPXPBDStep   p;
  drive(sc, p, g_serial, g_sched, 1, 10);
  const Vec3 after = sc.state.pose.position(1);
  if (!bitwise_eq(before, after)) {
    c.ok     = false;
    c.detail = std::format("A moved from ({},{},{}) to ({},{},{})", before.x, before.y, before.z,
                           after.x, after.y, after.z);
    return c;
  }
  for (const ConstraintVariant &v : sc.constraints) {
    const ContactConstraint &cc = std::get<ContactConstraint>(v);
    if (cc.lambda != Scalar_t(0)) {
      c.ok     = false;
      c.detail = std::format("row {} got lambda={} on a separated contact",
                             static_cast<int>(cc.kind), cc.lambda);
      return c;
    }
  }
  return c;
}

// Phase 3b：穿透接触只推不拉（λ_n > 0 ⇒ A 沿 +n 移动），切向乘子必须落在
// μλ_n 的锥内；slip 足够大时锥必须是紧的（不是把切线白白放过）。
Check check_bpxpbd_contact_friction_cone() {
  Check c{"BPXPBD contact friction cone bounds tangent multipliers", true, ""};

  scene::Scene sc = contact_pair_scene(Scalar_t(0.01), Scalar_t(0.5), Vec3{0.1f, 0, 0}, false);
  const Vec3   before = sc.state.pose.position(1);
  BPXPBDStep   p;
  drive(sc, p, g_serial, g_sched, 1, 1);
  const Vec3 after = sc.state.pose.position(1);
  if (!(after.y > before.y)) {
    c.ok     = false;
    c.detail = std::format("penetrating contact pushed A down: y {} -> {}", before.y, after.y);
    return c;
  }

  Scalar_t lam[3] = {0, 0, 0};
  for (const ConstraintVariant &v : sc.constraints) {
    const ContactConstraint &cc = std::get<ContactConstraint>(v);
    lam[static_cast<int>(cc.kind)] = cc.lambda;
  }
  if (!(lam[0] > Scalar_t(0))) {
    c.ok     = false;
    c.detail = std::format("normal multiplier is {} (want > 0)", lam[0]);
    return c;
  }
  const Scalar_t tau   = std::sqrt(lam[1] * lam[1] + lam[2] * lam[2]);
  const Scalar_t bound = Scalar_t(0.5) * lam[0];
  if (!(tau <= bound * (Scalar_t(1) + Scalar_t(1e-4)) + Scalar_t(1e-7))) {
    c.ok     = false;
    c.detail = std::format("tangent |lambda|={} exceeds mu*lambda_n={}", tau, bound);
    return c;
  }
  if (!(tau >= bound * Scalar_t(0.99))) {
    c.ok     = false;
    c.detail = std::format("cone is not active: |lambda_t|={} < 0.99*mu*lambda_n={}", tau, bound);
    return c;
  }
  return c;
}

// Phase 3c：接触 unit（法线 + 两切线）在 block 扫掠与"一 unit 一 block"的串行
// 顺序下必须逐位相同；一个接触点必须是一个 unit（3 行同块），不能被拆开。
Check check_bpxpbd_contact_block_oracle() {
  Check c{"BPXPBD contact block sweep == one unit per block", true, ""};

  constexpr std::size_t kContacts = 12;

  scene::Scene natural = contact_fan_scene(kContacts, Scalar_t(0.5));
  BPXPBDStep   p_nat;
  drive(natural, p_nat, g_serial, g_sched, 3, 20);

  SolveSchedule forced = g_sched;
  forced.force_single_unit_blocks = true;
  scene::Scene serial = contact_fan_scene(kContacts, Scalar_t(0.5));
  BPXPBDStep   p_ser;
  drive(serial, p_ser, g_serial, forced, 3, 20);

  const auto d_nat = digest::parts(natural.state, natural.constraints).all;
  const auto d_ser = digest::parts(serial.state, serial.constraints).all;
  if (d_nat != d_ser) {
    c.ok     = false;
    c.detail = std::format("block sweep={} serial={}", digest::to_hex(d_nat),
                           digest::to_hex(d_ser));
    return c;
  }

  // 自然分区：这些 unit 只共享静态 body ⇒ 全部落在同一个 block，顺序 = 约束顺序。
  const auto &scr_nat = natural.state.scratch;
  const std::size_t nb_nat = scr_nat.block_begin.empty() ? 0 : scr_nat.block_begin.size() - 1;
  if (nb_nat != 1) {
    c.ok     = false;
    c.detail = std::format("independent contacts used {} blocks (want 1)", nb_nat);
    return c;
  }

  // 强制一 unit 一 block：每个接触点必须仍是 3 行一个 unit，且块数 = 接触点数。
  const auto &scr_ser = serial.state.scratch;
  const std::size_t nb_ser = scr_ser.block_begin.empty() ? 0 : scr_ser.block_begin.size() - 1;
  if (scr_ser.unit_rows.size() != kContacts) {
    c.ok     = false;
    c.detail = std::format("forced partition made {} units (want {})", scr_ser.unit_rows.size(),
                           kContacts);
    return c;
  }
  for (std::size_t u = 0; u < scr_ser.unit_rows.size(); ++u) {
    if (scr_ser.unit_rows[u] != 3) {
      c.ok     = false;
      c.detail = std::format("contact unit {} has {} rows (want 3: normal + 2 tangents)", u,
                             scr_ser.unit_rows[u]);
      return c;
    }
  }
  if (nb_ser != kContacts) {
    c.ok     = false;
    c.detail = std::format("forced partition made {} blocks (want {})", nb_ser, kContacts);
    return c;
  }

  // 别让这个 oracle 空转：接触必须真的在承力。3 步之后物体可能已经被位置投影
  // 的"修正变速度"顶得短暂分离（C>0 ⇒ 单边投影把 λ 放回 0），所以这里用一个
  // 从静止开始的单步探针检查内核确实产生了正的法线乘子。
  scene::Scene probe = contact_fan_scene(4, Scalar_t(0.5));
  BPXPBDStep   p_probe;
  drive(probe, p_probe, g_serial, g_sched, 1, 1);
  Scalar_t lam_sum = 0;
  for (const ConstraintVariant &v : probe.constraints) {
    const ContactConstraint &cc = std::get<ContactConstraint>(v);
    if (cc.kind == ContactKind::Normal)
      lam_sum += cc.lambda;
  }
  if (!(lam_sum > Scalar_t(0))) {
    c.ok     = false;
    c.detail = "all contact normal multipliers are zero (oracle is vacuous)";
  }
  return c;
}

// Phase 3d：深度穿透 cap。投影看到的法线 C 被夹在 -kContactMaxPenetration，
// 所以**一次约束求解**最多修正 2 cm。5 cm / 9 cm 的初始重叠必须给出相同的位移，
// 而不是各自一步弹出去（9 cm 一步会把 0.085/dt ≈ 5 m/s 直接写进速度，正是斜坡
// box 被"发射"的机制）。切向不受影响（这条只测法线）。
Check check_bpxpbd_contact_deep_penetration_cap() {
  Check c{"BPXPBD deep-overlap contact caps one solve at 2 cm", true, ""};

  auto displacement = [](Scalar_t penetration, int iters) {
    scene::Scene sc = contact_pair_scene(penetration, Scalar_t(0), Vec3{0}, false);
    const Vec3   before = sc.state.pose.position(1);
    BPXPBDStep   p;
    drive(sc, p, g_serial, g_sched, 1, iters);
    return sc.state.pose.position(1) - before;
  };

  // 单轮迭代：cap 生效时两个深度给出同一个 2 cm 修正，未超 cap 的给出原始深度。
  const Vec3 shallow = displacement(Scalar_t(0.015), 1); // C=-0.01：cap 不生效
  const Vec3 mid     = displacement(Scalar_t(0.05), 1);  // C=-0.045：cap 生效
  const Vec3 deep    = displacement(Scalar_t(0.09), 1);  // C=-0.085：cap 生效

  constexpr Scalar_t kTol = Scalar_t(1e-6);
  if (std::abs(shallow.y - Scalar_t(0.01)) > kTol) {
    c.ok     = false;
    c.detail = std::format("uncapped 1 cm overlap moved {} (want 0.01)", shallow.y);
    return c;
  }
  if (std::abs(deep.y - kContactMaxPenetration) > kTol) {
    c.ok     = false;
    c.detail = std::format("9 cm overlap moved {} (want cap {})", deep.y,
                           kContactMaxPenetration);
    return c;
  }
  // cap 内的两个深度必须给出**逐位相同**的修正量。
  if (!bitwise_eq(mid, deep)) {
    c.ok     = false;
    c.detail = std::format("5 cm moved {} but 9 cm moved {} under the same cap", mid.y, deep.y);
    return c;
  }
  // 没有 cap 的旧行为是 0.085 m；这里必须明显更小（发射速度 0.02/dt vs 0.085/dt）。
  if (!(deep.y < Scalar_t(0.03))) {
    c.ok     = false;
    c.detail = std::format("9 cm correction {} is not capped", deep.y);
    return c;
  }
  if (std::abs(deep.x) > kTol || std::abs(deep.z) > kTol) {
    c.ok     = false;
    c.detail = std::format("cap leaked into tangents: dx={} dz={}", deep.x, deep.z);
    return c;
  }

  // 这是真正的"不发射"断言：默认的 8 轮迭代也必须只走 2 cm。单纯把实时 C 夹到
  // -2 cm 的实现会在这里露馅 —— 每轮都读到 -2 cm，8 轮下来仍然修完 8.5 cm。
  const Vec3 deep8 = displacement(Scalar_t(0.09), 8);
  if (!(deep8.y > Scalar_t(0)) || !(deep8.y <= kContactMaxPenetration + kTol)) {
    c.ok     = false;
    c.detail = std::format("9 cm overlap moved {} over 8 iterations (want ≤ cap)", deep8.y);
    return c;
  }

  // 撞击步**不能**夹：静息深重叠才夹（见 begin_sweep 的速度门）。用 2 m/s 的
  // 接近速度喂给同一个 9 cm 重叠，位移必须是完整的去穿透（≈ 8.5 cm），否则
  // 冲击只能以 2 cm/步退穿透，整摞会被压着穿过地板（双金字塔实测的第一处
  // 背离）。这条同时钉住"cap 没有偷偷变成全局的软化"。
  {
    scene::Scene sc = contact_pair_scene(Scalar_t(0.09), Scalar_t(0), Vec3{0}, false);
    sc.state.velocity.set_velocity(1, Vec3{0, -2, 0});
    const Vec3 before = sc.state.pose.position(1);
    BPXPBDStep p;
    drive(sc, p, g_serial, g_sched, 1, 1);
    const Vec3 d = sc.state.pose.position(1) - before;
    if (!(d.y > Scalar_t(0.05))) {
      c.ok     = false;
      c.detail = std::format("impact (2 m/s) was depenetrated only {} m (cap must be off)", d.y);
      return c;
    }
  }
  return c;
}

// Phase 3e：步末速度级摩擦必须是**纯耗散**的。
// 手工给一个真实重叠的接触一点法向乘子 λ_n > 0 和切向滑块速度，然后只跑
// post_step：切向相对速度必须下降，系统动能不许上升。
Check check_bpxpbd_friction_pass_is_dissipative() {
  Check c{"BPXPBD friction pass removes tangential speed without adding energy", true, ""};

  auto kinetic = [](const State &s) {
    Scalar_t e = 0;
    for (std::size_t b = 0; b < s.size(); ++b) {
      if (s.inertial.is_static(b))
        continue;
      const Vec3 v = s.velocity.velocity(b);
      const Vec3 w = s.velocity.angular_velocity(b);
      const Mat3 I = glm::inverse(s.inertial.world_inv_inertia(b));
      e += Scalar_t(0.5) * (glm::dot(v, v) / s.inertial.inv_mass(b) + glm::dot(w, I * w));
    }
    return e;
  };

  // 接触点相对切向速度（与 policy 的 post_step 同一几何）。
  auto tangential_speed = [](const State &s, const ContactConstraint &row) {
    const Vec3 rAw = rotate(s.pose.rotation(row.bodies_[0]), row.rA);
    const Vec3 rBw = rotate(s.pose.rotation(row.bodies_[1]), row.rB);
    const Vec3 vA =
        s.velocity.velocity(row.bodies_[0]) +
        glm::cross(s.velocity.angular_velocity(row.bodies_[0]), rAw);
    const Vec3 vB =
        s.velocity.velocity(row.bodies_[1]) +
        glm::cross(s.velocity.angular_velocity(row.bodies_[1]), rBw);
    const Vec3     v_rel = vA - vB;
    const Scalar_t x     = glm::dot(v_rel, row.basis[1]);
    const Scalar_t y     = glm::dot(v_rel, row.basis[2]);
    return std::sqrt(x * x + y * y);
  };

  auto probe = [&](Scalar_t lambda_n) {
    scene::Scene sc = contact_pair_scene(Scalar_t(0.01), Scalar_t(0.5), Vec3{0}, false);
    for (ConstraintVariant &v : sc.constraints) {
      ContactConstraint &cc = std::get<ContactConstraint>(v);
      if (cc.kind == ContactKind::Normal)
        cc.lambda = lambda_n;
    }
    sc.state.velocity.set_velocity(1, Vec3{1, 0, 0}); // 1 m/s across the contact plane
    // post_step 读世界逆惯量；真实帧里 begin_step 已经刷新，这里直接调钩子要先补上。
    sc.state.inertial.refresh_all(sc.state.pose);
    const ContactConstraint &normal = std::get<ContactConstraint>(sc.constraints.front());
    const Scalar_t           before_vt = tangential_speed(sc.state, normal);
    const Scalar_t           before_e  = kinetic(sc.state);

    BPXPBDStep p;
    p.post_step(sc.state, sc.constraints);

    const Scalar_t after_vt = tangential_speed(sc.state, normal);
    const Scalar_t after_e  = kinetic(sc.state);
    return std::tuple{before_vt, after_vt, before_e, after_e};
  };

  // λ_n = 0.1 kg·m ⇒ μλ_n/dt 远超 1 m/s：静摩擦必须把切向速度压到 ~0。
  const auto [b0, a0, e0, f0] = probe(Scalar_t(0.1));
  if (!(b0 > Scalar_t(0.5))) {
    c.ok     = false;
    c.detail = std::format("probe is vacuous: initial tangential speed {}", b0);
    return c;
  }
  if (!(f0 <= e0 + Scalar_t(1e-6))) {
    c.ok     = false;
    c.detail = std::format("friction added energy: {} -> {}", e0, f0);
    return c;
  }
  if (!(a0 < b0 * Scalar_t(0.01)) || !(a0 < Scalar_t(1e-4))) {
    c.ok     = false;
    c.detail = std::format("static friction left tangential speed {} -> {}", b0, a0);
    return c;
  }

  // λ_n = 1e-4 kg·m ⇒ 锥界只允许削掉一点点：速度必须下降但仍有剩余。
  const auto [b1, a1, e1, f1] = probe(Scalar_t(1e-4));
  if (!(a1 < b1) || !(a1 > Scalar_t(0.5) * b1)) {
    c.ok     = false;
    c.detail = std::format("bounded impulse is wrong: {} -> {}", b1, a1);
    return c;
  }
  if (!(f1 <= e1 + Scalar_t(1e-6))) {
    c.ok     = false;
    c.detail = std::format("bounded friction added energy: {} -> {}", e1, f1);
    return c;
  }
  return c;
}

// 数值卫生：所有解算结果必须有限，且不许超出物理量级（solve_body 的冻结阈值是 10m）。
Check check_all_finite_and_bounded() {
  Check c{"all finite and within physical bound", true, ""};
  for (scene::Scene sc :
       {scene::chain(32), scene::fan(64), scene::lattice(4, 3, 2), scene::mixed(20)}) {
    drive(sc, AVBDStep{}, g_serial, g_sched, 2, 30);
    for (std::size_t b = 0; b < sc.state.size(); ++b) {
      const Vec3 p = sc.state.pose.position(b);
      if (!all_finite(p) || !all_finite(sc.state.velocity.velocity(b)) ||
          !all_finite(sc.state.velocity.angular_velocity(b))) {
        c.ok     = false;
        c.detail = std::format("non-finite at body {}", b);
        return c;
      }
      if (std::abs(p.x) > 1.0e4f || std::abs(p.y) > 1.0e4f || std::abs(p.z) > 1.0e4f) {
        c.ok     = false;
        c.detail = std::format("body {} ran away to ({}, {}, {})", b, p.x, p.y, p.z);
        return c;
      }
    }
  }
  return c;
}

// 场景必须真的有初始约束违反，否则整张用例表测不出任何东西（见 scene.hpp 的坑）。
Check check_scenes_are_violated() {
  Check c{"scenes start with violated constraints", true, ""};
  auto  max_abs_C = [](const scene::Scene &sc) {
    Scalar_t worst = 0;
    for (const auto &v : sc.constraints)
      worst = std::max(
          worst, std::abs(std::visit([&](const auto &cc) { return cc.evaluate(sc.state); }, v)));
    return worst;
  };
  const Scalar_t w_chain = max_abs_C(scene::chain(32));
  const Scalar_t w_fan   = max_abs_C(scene::fan(64));
  const Scalar_t w_mixed = max_abs_C(scene::mixed(20));
  if (!(w_chain > 1e-3f) || !(w_fan > 1e-3f) || !(w_mixed > 1e-3f)) {
    c.ok     = false;
    c.detail = std::format("chain={} fan={} mixed={}", w_chain, w_fan, w_mixed);
  }
  return c;
}

// 无约束时，AVBD 的隐式残差与 XPBD 的显式位姿积分必须给出**逐位相同**的结果：
// 前者解出 Δ = -r（即 x ← x_pred），后者直接把 pose 搬到 x_pred。两条完全不同的
// 代码路径在"纯惯性"极限下重合，说明残差项的符号/系数/参考位姿都没搞错。
Check check_inertial_limit_policies_agree() {
  Check        c{"no constraints: AVBD implicit == XPBD explicit", true, ""};
  scene::Scene a = scene::freefall(4);
  scene::Scene b = scene::freefall(4);
  drive(a, AVBDStep{}, g_serial, g_sched, 6, 10);
  drive(b, XPBDStep{}, g_serial, g_sched, 6, 10);
  const auto pa = digest::pose(a.state.pose), pb = digest::pose(b.state.pose);
  const auto va = digest::velocity(a.state.velocity), vb = digest::velocity(b.state.velocity);
  if (pa != pb || va != vb) {
    c.ok     = false;
    c.detail = std::format("pose {} vs {}, vel {} vs {}", digest::to_hex(pa), digest::to_hex(pb),
                           digest::to_hex(va), digest::to_hex(vb));
  }
  return c;
}

// 重力方向/符号的端到端验收。
//
// 教训：第一版写的是" strained chain 的 tip 必须往下沉 "，它在 k 从不复位的旧代码上
// 恰好成立，纯属巧合 —— chain() 的 strain=0.85 让约束把链条强烈地往锚点收（x 方向），
// 竖直方向只是收缩+转角耦合出来的瞬态，方向和重力没有必然关系（实测 +0.028）。
// 断言写错了场景，会把"物理正常"误判成"物理坏了"。
//
// 现在的做法是两条互补的断言：
//   1) 用**静止**链条（strain==1 ⇒ C≡0）单独隔离重力：此时唯一的驱动就是重力，
//      tip 必须严格下沉；没有重力时逐位不动（金丝雀已在别处覆盖）。
//   2) 对 strained chain 只做**差分**断言：有重力和无重力必须不同，
//      不假定方向 —— 这样收缩瞬态怎么演化都不会误报。
Check check_gravity_pulls_an_at_rest_chain_down() {
  Check        c{"gravity pulls an at-rest chain down (sign)", true, ""};
  scene::Scene sc            = scene::chain(8, Scalar_t(0.5), Scalar_t(1.0), false);
  const Vec3   anchor_before = sc.state.pose.position(0);

  drive(sc, AVBDStep{}, g_serial, g_sched, 10, 30);

  if (!bitwise_eq(anchor_before, sc.state.pose.position(0))) {
    c.ok     = false;
    c.detail = "static anchor moved under gravity";
    return c;
  }
  const Scalar_t tip_y = sc.state.pose.position(7).y;
  if (!(tip_y < -1e-3f)) {
    c.ok     = false;
    c.detail = std::format("tip did not sag: y={} (expected < -1e-3)", tip_y);
    return c;
  }
  return c;
}

Check check_gravity_changes_a_strained_chain() {
  Check c{"gravity measurably affects a strained chain", true, ""};
  auto  run = [](bool gravity) {
    scene::Scene sc = scene::chain(8);
    if (!gravity)
      scene::make_gravity_free(sc);
    drive(sc, AVBDStep{}, g_serial, g_sched, 10, 30);
    return sc.state.pose.position(7);
  };
  const Vec3 with_g  = run(true);
  const Vec3 without = run(false);
  if (std::abs(with_g.y - without.y) < 1e-3f && std::abs(with_g.x - without.x) < 1e-3f) {
    c.ok     = false;
    c.detail = std::format("gravity had no effect: ({}, {}) vs ({}, {})", with_g.x, with_g.y,
                           without.x, without.y);
  }
  return c;
}

// ==================== 执行器契约自检 ====================
// parallel_for 的契约是"每个 index 恰好被执行一次"。这不是小事：后面所有阶段都
// 建在它上面，多跑一次少跑一次在物理里表现为"略有点不对"，几乎不可能用肉眼看出来。
// 所以用"每个 index 被命中次数"的计数器直接验证，而不是间接看结果。
Check check_executor_covers_exactly_once(task::Executor auto &exec) {
  Check c{std::format("executor: exact-once coverage ({} threads, serial={})", exec.thread_count(),
                      std::remove_reference_t<decltype(exec)>::is_serial()),
          true, ""};

  // (n, grain) 组合要覆盖：n==0、n<grain、n%grain!=0、chunks>threads。
  const std::size_t cases[][2] = {{0, 1},  {1, 1},     {1, 8},    {7, 3},    {64, 1},
                                  {64, 7}, {1000, 13}, {999, 64}, {4096, 32}};
  for (const auto &n_grain : cases) {
    const std::size_t n     = n_grain[0];
    const std::size_t grain = n_grain[1];
    std::vector<int>  hits(n, 0);
    exec.parallel_for(n, grain, [&hits](std::size_t b, std::size_t e) {
      for (std::size_t i = b; i < e; ++i)
        hits[i] += 1;
    });
    for (std::size_t i = 0; i < n; ++i) {
      if (hits[i] != 1) {
        c.ok     = false;
        c.detail = std::format("n={} grain={}: index {} hit {} times", n, grain, i, hits[i]);
        return c;
      }
    }
  }

  // 区间非空且不越界（契约 3）。
  bool bad_range = false;
  exec.parallel_for(37, 5, [&](std::size_t b, std::size_t e) {
    if (b >= e || e > 37)
      bad_range = true;
  });
  if (bad_range) {
    c.ok     = false;
    c.detail = "empty or out-of-range chunk";
  }
  return c;
}

// 反复 fork-join 不得泄漏/死锁/重复（generation 竞态的守门测试）。
Check check_executor_reuse_is_stable(task::Executor auto &exec) {
  Check c{std::format("executor: repeated fork-join stable ({} threads)", exec.thread_count()),
          true, ""};
  for (int round = 0; round < 200; ++round) {
    std::atomic<long> sum{0};
    const std::size_t n = 1000;
    exec.parallel_for(n, 1, [&sum](std::size_t b, std::size_t e) {
      long local = 0;
      for (std::size_t i = b; i < e; ++i)
        local += static_cast<long>(i);
      sum.fetch_add(local, std::memory_order_relaxed);
    });
    const long want = static_cast<long>(n) * (static_cast<long>(n) - 1) / 2;
    if (sum.load() != want) {
      c.ok     = false;
      c.detail = std::format("round {}: sum {} != {}", round, sum.load(), want);
      return c;
    }
  }
  return c;
}

// 参考实现：把 AVBD 的累积阶段强制回"直接写 buffers"的老路子。
// 靠把 uses_solver_buffers 打成 false 走 driver 的 else 分支 —— 那条分支正是
// step 4 之前的行为：逐约束 read-modify-write 到 per-body 局部系统。
//
// 这就是 AVBDStep::step 的存活理由：它是 scatter/gather 的**语义基准**。
// 两条路径直接对比，比"和录下来的数字比"更能说明问题，也不依赖基线文件。
struct AVBDDirectAccumulate : AVBDStep {
  static constexpr bool uses_solver_buffers = false; // 走直接累积分支
};

// scatter/gather 必须与直接累积**位精确**相等。这是 step 4 的核心断言：
// 它同时证明了"分配顺序 == 串行顺序"，以及 accumulate_entry 的权重没错。
Check check_scatter_gather_equals_direct_accumulation() {
  Check c{"scatter/gather == direct per-body accumulation", true, ""};

  auto run_with = [](auto policy, scene::Scene sc, int steps, int iters) {
    drive(sc, policy, g_serial, g_sched, steps, iters);
    return digest::parts(sc.state, sc.constraints).all;
  };

  const std::pair<scene::Scene, std::pair<int, int>> cases[] = {
      {scene::chain(32), {4, 50}},
      {scene::fan(64), {3, 30}},
      {scene::lattice(4, 3, 2), {3, 30}},
      {scene::mixed(20), {3, 30}},
  };
  for (const auto &[sc, si] : cases) {
    const auto via_entries = run_with(AVBDStep{}, sc, si.first, si.second);
    const auto via_direct  = run_with(AVBDDirectAccumulate{}, sc, si.first, si.second);
    if (via_entries != via_direct) {
      c.ok     = false;
      c.detail = std::format("{} diverged: entries={} direct={}", sc.name,
                             digest::to_hex(via_entries), digest::to_hex(via_direct));
      return c;
    }
  }
  return c;
}

// gather 的分块边界不得影响结果。这是"顺序来自数据布局、不来自调度"的直接检验：
// 只要有人把 gather 改成"按线程收集再合并"，块边界一变结果就会变。
Check check_gather_partitioning_is_bit_exact(task::Executor auto &pool) {
  Check c{std::format("gather chunking invariant ({} threads)", pool.thread_count()), true, ""};

  SolveSchedule whole{};
  whole.parallel_accumulate = false; // 一次跑完整个 body 区间

  std::uint64_t reference = 0;
  {
    auto results = run_all(pool, whole, false);
    reference    = std::hash<std::string>{}(results[0].line + results[6].line + results[9].line);
  }

  for (std::size_t grain : {1u, 3u, 7u, 64u}) {
    SolveSchedule chunked{};
    chunked.parallel_accumulate = true;
    chunked.grain_body          = grain;
    chunked.grain_constraint    = grain;
    const auto results          = run_all(pool, chunked, false);
    const auto h = std::hash<std::string>{}(results[0].line + results[6].line + results[9].line);
    if (h != reference) {
      c.ok     = false;
      c.detail = std::format("grain={} changed results", grain);
      return c;
    }
    if (grain == 7u) { // 顺带把整个表也比一遍，避免只看三个用例
      const auto whole_results = run_all(pool, whole, false);
      for (std::size_t i = 0; i < results.size(); ++i) {
        if (results[i].line != whole_results[i].line) {
          c.ok     = false;
          c.detail = std::format("grain=7 case {} diverged", i + 1);
          return c;
        }
      }
    }
  }
  return c;
}

// 刚度 ramp 必须是**步内**语义：跨步不许累积。
// 没有这条，k 会跨步累乘 (1+beta)^(iters·steps)，几十万帧后 → inf，
// 然后 10 m 的修正上限会把每个 DOF 逐个冻结 —— 物理静默僵死，不报任何错。
Check check_stiffness_ramp_is_bounded_across_steps() {
  Check          c{"stiffness ramp resets every step", true, ""};
  scene::Scene   sc = scene::chain(8);
  AVBDStep       p;
  const Scalar_t k0    = AVBDStep::kStiffness0;
  const int      iters = 30;
  const int      steps = 40;

  for (int step = 0; step < steps; ++step) {
    sc.state.reset_external();
    sc.state.begin_step();
    solve(sc.state, sc.constraints, p, iters, g_serial, g_sched,
          static_cast<NullPhaseObserver *>(nullptr), sc.topology_version);
    sc.state.finish_step();
    sc.state.cache.commit_prev(sc.state.velocity);
  }

  // 步内 ramp iters 次 ⇒ k <= k0·(1+beta)^iters。跨步不许多涨一次。
  const Scalar_t bound = k0 * std::pow(Scalar_t(1) + p.beta, Scalar_t(iters + 1));
  if (!(p.k > 0) || !(p.k <= bound)) {
    c.ok     = false;
    c.detail = std::format("k={} after {} steps (want 0 < k <= {}); ramp leaked across steps?", p.k,
                           steps, bound);
  }
  return c;
}

// ==================== 帧流水线（step 10）====================

// 帧层必须是纯重组：走 frame::step_frame 与手工按序调用（reset_external →
// begin_step → solve → finish_step → commit_prev）必须逐位相同。
// 这条把"帧层不偷偷改变物理"钉死。
Check check_frame_pipeline_equals_manual_drive() {
  Check c{"frame::step_frame == manual lifecycle", true, ""};

  auto manual = [] {
    scene::Scene sc = scene::lattice(4, 3, 2);
    State       &s  = sc.state;
    AVBDStep     p;
    for (int step = 0; step < 3; ++step) {
      s.reset_external();
      s.begin_step();
      solve(s, sc.constraints, p, 25, g_serial, g_sched, static_cast<NullPhaseObserver *>(nullptr),
            sc.topology_version);
      s.finish_step();
      s.cache.commit_prev(s.velocity);
    }
    return digest::parts(s, sc.constraints).all;
  };

  auto via_frame = [] {
    scene::Scene         sc = scene::lattice(4, 3, 2);
    AVBDStep             p;
    frame::FrameSchedule fs{};
    fs.solve  = g_sched;
    auto noop = [](std::size_t) {};
    for (int step = 0; step < 3; ++step)
      frame::step_frame(sc.state, sc.constraints, p, g_serial, fs, 25, sc.topology_version, &noop);
    return digest::parts(sc.state, sc.constraints).all;
  };

  const auto a = manual();
  const auto b = via_frame();
  if (a != b) {
    c.ok     = false;
    c.detail = std::format("manual={} frame={}", digest::to_hex(a), digest::to_hex(b));
  }
  return c;
}

// 融合进来的 consumer 必须在**求解完成后**看到最终位姿，且与"求解后串行遍历"一致。
// 这条防的是"融合把 consumer 排到了 solve 前面"这种顺序错误 —— 那种错误下
// 物理摘要完全正常，只有消费者拿到的是上一帧的数据。
Check check_frame_consumer_sees_final_pose() {
  Check c{"fused consumer observes the post-solve pose", true, ""};

  scene::Scene         sc = scene::lattice(4, 3, 2);
  AVBDStep             p;
  frame::FrameSchedule fs{};
  fs.solve          = g_sched;
  fs.min_items_tail = 1; // 强制融合分叉真的并行

  frame::SceneAabb fused;
  fused.reset(sc.state.size());
  auto consume = [&fused, &sc](std::size_t b) { fused.store(b, sc.state.pose.position(b)); };
  task::ThreadPoolExecutor pool(4);
  frame::step_frame(sc.state, sc.constraints, p, pool, fs, 20, sc.topology_version, &consume);
  fused.reduce();

  frame::SceneAabb reference;
  reference.reset(sc.state.size());
  for (std::size_t b = 0; b < sc.state.size(); ++b)
    reference.store(b, sc.state.pose.position(b));
  reference.reduce();

  if (!fused.valid || !reference.valid) {
    c.ok     = false;
    c.detail = "consumer produced no AABB";
    return c;
  }
  if (!bitwise_eq(fused.lo, reference.lo) || !bitwise_eq(fused.hi, reference.hi)) {
    c.ok = false;
    c.detail =
        std::format("fused=({},{},{})..({},{},{}) vs serial=({},{},{})..({},{},{})", fused.lo.x,
                    fused.lo.y, fused.lo.z, fused.hi.x, fused.hi.y, fused.hi.z, reference.lo.x,
                    reference.lo.y, reference.lo.z, reference.hi.x, reference.hi.y, reference.hi.z);
  }
  return c;
}

// ==================== 确定性矩阵（step 8）====================
// 不变量 10 的可执行版本：把"线程数 × 分区粒度 × 串行阈值"的笛卡尔积跑一遍，
// 每一位都必须等于串行参考。刻意包含 grain=1 / min_items=1 这种"最碎、最并行"的
// 极端配置 —— 那是最容易暴露"求和顺序依赖调度"的地方。
Check check_determinism_matrix() {
  Check                     c{"digests invariant across threads x grain x thresholds", true, ""};
  const std::vector<Result> reference = run_all(g_serial, SolveSchedule{}, false);

  for (unsigned threads : {2u, 3u, 8u}) {
    for (std::size_t grain : {std::size_t(1), std::size_t(5), std::size_t(16), std::size_t(64)}) {
      for (std::size_t min_items : {std::size_t(1), std::size_t(4096)}) {
        SolveSchedule cfg{};
        cfg.grain_body           = grain;
        cfg.grain_constraint     = grain;
        cfg.min_items_body       = min_items;
        cfg.min_items_constraint = min_items;

        task::ThreadPoolExecutor  pool(threads);
        const std::vector<Result> got = run_all(pool, cfg, false);
        for (std::size_t i = 0; i < reference.size(); ++i) {
          if (reference[i].line != got[i].line) {
            c.ok     = false;
            c.detail = std::format("threads={} grain={} min_items={} case {} diverged", threads,
                                   grain, min_items, i + 1);
            return c;
          }
        }
      }
    }
  }
  return c;
}

// 同一配置重复跑必须完全相同 —— 抓任何一处"读了未初始化/上一轮的残留"。
Check check_runs_are_reproducible(task::Executor auto &pool) {
  Check                     c{"repeated runs are bit-identical", true, ""};
  const std::vector<Result> a = run_all(pool, SolveSchedule{}, false);
  const std::vector<Result> b = run_all(pool, SolveSchedule{}, false);
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (a[i].line != b[i].line) {
      c.ok     = false;
      c.detail = std::format("case {} differs between runs: {} vs {}", i + 1, a[i].line, b[i].line);
      return c;
    }
  }
  return c;
}

// ==================== 拓扑缓存（step 5）====================

// 缓存必须真的被复用：跑多步只允许重建一次。
// 这条断言防的是"每次 solve 都换了个 version / 签名永远失配"的退化——那种情况下
// 功能完全正确，缓存一行都没生效，性能却悄悄退回到每步重建。
Check check_topology_cache_is_reused() {
  Check        c{"topology layout rebuilt exactly once for N steps", true, ""};
  scene::Scene sc = scene::lattice(5, 4, 1);
  drive(sc, AVBDStep{}, g_serial, g_sched, 6, 20); // 6 步
  if (sc.state.scratch.build_count != 1) {
    c.ok     = false;
    c.detail = std::format("build_count={} after 6 steps (want 1)", sc.state.scratch.build_count);
  }
  return c;
}

// 复用出来的布局必须与"每步重建"逐字节相同。
// 做法：同一场景跑两遍，一遍 version=1（全程复用），一遍 version=0（每步重建），
// 摘要与 build_count 双双比对——既证明缓存内容没错，也证明缓存确实在起作用。
Check check_cached_layout_equals_rebuilt() {
  Check        c{"cached layout == freshly rebuilt layout", true, ""};
  scene::Scene reused      = scene::lattice(4, 3, 2);
  scene::Scene rebuilt     = scene::lattice(4, 3, 2);
  rebuilt.topology_version = 0; // 0 = 永远重建

  drive(reused, AVBDStep{}, g_serial, g_sched, 4, 25);
  drive(rebuilt, AVBDStep{}, g_serial, g_sched, 4, 25);

  const auto a = digest::parts(reused.state, reused.constraints).all;
  const auto b = digest::parts(rebuilt.state, rebuilt.constraints).all;
  if (a != b) {
    c.ok     = false;
    c.detail = std::format("reused={} rebuilt={}", digest::to_hex(a), digest::to_hex(b));
    return c;
  }
  if (reused.state.scratch.build_count != 1 || rebuilt.state.scratch.build_count != 4) {
    c.ok     = false;
    c.detail = std::format("build_count: reused={} (want 1), rebuilt={} (want 4)",
                           reused.state.scratch.build_count, rebuilt.state.scratch.build_count);
    return c;
  }
  return c;
}

// 缓存失效必须走对路：body 数不变时 resize_bodies() 必须是个完全的空操作。
// （这里踩过一个坑：早期版本无条件清 body_entries，却不动 key ⇒ 下一轮
// "缓存命中"拿到的是空邻接，gather 会越界读。这条断言就是为了钉死这个坑。）
Check check_resize_with_same_body_count_keeps_cache() {
  Check         c{"resize to same body count is a no-op for the topology cache", true, ""};
  SolveSchedule strict{};
  strict.verify_topology = true;

  scene::Scene sc = scene::lattice(4, 3, 2);
  drive(sc, AVBDStep{}, g_serial, strict, 1, 5);
  const std::uint64_t builds_before = sc.state.scratch.build_count;
  const std::size_t   adj_before    = sc.state.scratch.body_entries.size();

  sc.state.resize(sc.state.size()); // 同尺寸重分配
  drive(sc, AVBDStep{}, g_serial, strict, 1, 5);

  if (sc.state.scratch.build_count != builds_before) {
    c.ok     = false;
    c.detail = std::format("build_count {} -> {} (cache was needlessly invalidated)", builds_before,
                           sc.state.scratch.build_count);
  } else if (sc.state.scratch.body_entries.size() != adj_before || adj_before == 0) {
    c.ok     = false;
    c.detail = std::format("adjacency size {} -> {} (cleared under a valid key!)", adj_before,
                           sc.state.scratch.body_entries.size());
  } else if (sc.state.scratch.stale_count != 0) {
    c.ok     = false;
    c.detail = "verifier flagged a stale layout";
  }
  return c;
}

// 校验模式必须报 0 次 stale：只要有人改了约束却没 bump version，这里就会非零。
// 这条同时也自检了 verifier 本身不是"永远说 OK"的摆设。
Check check_topology_version_bookkeeping_strict() {
  Check         c{"verify_topology reports 0 stale cache hits", true, ""};
  SolveSchedule strict{};
  strict.verify_topology = true;

  scene::Scene sc = scene::lattice(4, 3, 2);
  drive(sc, AVBDStep{}, g_serial, strict, 5, 20);
  if (sc.state.scratch.stale_count != 0) {
    c.ok     = false;
    c.detail = std::format("stale_count={} (someone mutated constraints without bumping "
                           "topology_version)",
                           sc.state.scratch.stale_count);
  }
  return c;
}

// verifier 必须真的能抓到错误版本号 —— 否则上面那条断言是永真的摆设。
// 故意在改完拓扑之后不复位缓存、也不换 version，看 verifier 会不会红。
Check check_strict_topology_verifier_has_teeth() {
  Check         c{"verify_topology detects a deliberately stale cache", true, ""};
  SolveSchedule strict{};
  strict.verify_topology = true;

  scene::Scene sc = scene::lattice(4, 3, 2);
  // 先建立缓存（version=1）
  drive(sc, AVBDStep{}, g_serial, strict, 1, 5);
  if (sc.state.scratch.build_count != 1) {
    c.ok     = false;
    c.detail = "cache was not built on the first step";
    return c;
  }

  // 在**不改 version、不清缓存**的前提下改拓扑：把最后一条约束挪到别的 body 上。
  // 这正是"忘记 bump"的真实形态。
  std::get<DistanceConstraint>(sc.constraints.back()).bodies_ = {0, 2};
  drive(sc, AVBDStep{}, g_serial, strict, 1, 5);

  if (sc.state.scratch.stale_count == 0) {
    c.ok     = false;
    c.detail = "verifier accepted a stale layout — it is a no-op";
  }
  return c;
}

// 并发正确性的主断言：同一张用例表在线程池上跑出来的每一位都必须等于串行结果。
//
// 为什么这条现在就能过：step 3 只把 dispatch 抽出来了，primal（逐 body）与 dual
// （逐约束）本来就是每项独立写的，分区不影响算术序列。唯一还串行的 accumulate
// 由 SolveSchedule::parallel_accumulate 守着（step 4 用 scatter/gather 解开）。
// 将来任何"顺手并行化"只要碰上真实写冲突，这条立刻翻车。
Check check_parallel_suite_is_bit_identical(task::Executor auto &pool) {
  Check c{std::format("parallel suite bit-identical to serial ({} threads)", pool.thread_count()),
          true, ""};
  const std::vector<Result> serial   = run_all(g_serial, g_sched, false);
  const std::vector<Result> parallel = run_all(pool, g_sched, false);
  if (serial.size() != parallel.size()) {
    c.ok     = false;
    c.detail = "case count mismatch";
    return c;
  }
  for (std::size_t i = 0; i < serial.size(); ++i) {
    if (serial[i].line != parallel[i].line) {
      c.ok     = false;
      c.detail = std::format("case {} diverged\n     serial: {}\n   parallel: {}", i + 1,
                             serial[i].line, parallel[i].line);
      return c;
    }
  }
  return c;
}

// 换线程数也必须位一致。逐 body / 逐约束的相位天然满足；一旦有人引入
// "按线程数决定求和顺序"的写法（例如固定 range 分片 + 有序 reduce），这条就会抓到。
Check check_parallel_is_thread_count_invariant() {
  Check                     c{"digests invariant across 1/2/4/8 threads", true, ""};
  const std::vector<Result> reference = run_all(g_serial, g_sched, false);
  for (unsigned threads : {1u, 2u, 4u, 8u}) {
    task::ThreadPoolExecutor  pool(threads);
    const std::vector<Result> got = run_all(pool, g_sched, false);
    for (std::size_t i = 0; i < reference.size(); ++i) {
      if (reference[i].line != got[i].line) {
        c.ok     = false;
        c.detail = std::format("threads={} case {} diverged", threads, i + 1);
        return c;
      }
    }
  }
  return c;
}

std::vector<Check> run_collision_selftest() {
  std::vector<Check> out;
  for (auto &rc : collision_selftest::run())
    out.push_back(Check{rc.name, rc.ok, rc.detail});
  return out;
}

std::vector<Check> run_selftest() {
  std::vector<Check> checks{
      check_scenes_are_violated(),
      check_buffers_cleared_after_solve(),
      check_dq_consistent(),
      check_static_bodies_frozen(),
      check_gravity_free_rest_scene_is_frozen(),
      check_free_fall_matches_analytic(),
      check_gravity_is_not_special(),
      check_xpbd_equals_pbd_when_rigid(),
      check_bpxpbd_zeta_one_equals_xpbd(),
      check_bpxpbd_one_unit_per_block_equals_xpbd(),
      check_bpxpbd_blend_reduces_residual(),
      check_bpxpbd_verify_topology(),
      check_bpxpbd_determinism(),
      check_bpxpbd_rest_chain_is_frozen(),
      check_bpxpbd_joint_axis_locked_equals_xpbd(),
      check_bpxpbd_angular_joint_locked_equals_xpbd(),
      check_bpxpbd_joint_axis_limited(),
      check_bpxpbd_contact_unilateral_no_pull(),
      check_bpxpbd_contact_friction_cone(),
      check_bpxpbd_contact_block_oracle(),
      check_bpxpbd_contact_deep_penetration_cap(),
      check_bpxpbd_friction_pass_is_dissipative(),
      check_scatter_gather_equals_direct_accumulation(),
      check_topology_cache_is_reused(),
      check_cached_layout_equals_rebuilt(),
      check_resize_with_same_body_count_keeps_cache(),
      check_topology_version_bookkeeping_strict(),
      check_strict_topology_verifier_has_teeth(),
      check_inertial_limit_policies_agree(),
      check_gravity_pulls_an_at_rest_chain_down(),
      check_gravity_changes_a_strained_chain(),
      check_all_finite_and_bounded(),
  };
  for (auto &c : run_collision_selftest())
    checks.push_back(std::move(c));
  return checks;
}

// 执行器自检要同时跑串行与多线程两个实现。
std::vector<Check> run_executor_selftest(task::Executor auto &serial, task::Executor auto &pool) {
  return {
      check_executor_covers_exactly_once(serial),
      check_executor_reuse_is_stable(serial),
      check_executor_covers_exactly_once(pool),
      check_executor_reuse_is_stable(pool),
      check_parallel_suite_is_bit_identical(pool),
      check_parallel_is_thread_count_invariant(),
      check_gather_partitioning_is_bit_exact(pool),
      check_determinism_matrix(),
      check_runs_are_reproducible(pool),
      check_frame_pipeline_equals_manual_drive(),
      check_frame_consumer_sees_final_pose(),
      check_stiffness_ramp_is_bounded_across_steps(),
  };
}

// ==================== 布局报告：伪共享推导的依据 ====================
void print_layout() {
  std::println(stderr, "\n[layout] sizeof Vec3={} Quat={} Mat3={} ConstraintEntry={}", sizeof(Vec3),
               sizeof(Quat), sizeof(Mat3), sizeof(ConstraintEntry));
  std::println(stderr, "  SolverBuffers/body = {}", kBodyBufferBytes);
  std::println(stderr, "  bodies/cacheline : buffers={:.2f} pos={:.2f} rot={:.2f} entry={:.2f}",
               double(kCacheLineBytes) / double(kBodyBufferBytes),
               double(kCacheLineBytes) / double(sizeof(Vec3)),
               double(kCacheLineBytes) / double(sizeof(Quat)),
               double(kCacheLineBytes) / double(kEntryBytes));
  std::println(stderr, "  grain_body={} ⇒ chunk={} B = {} cacheline{}", kDefaultGrainBody,
               kDefaultGrainBody * kBodyBufferBytes,
               (kDefaultGrainBody * kBodyBufferBytes) / kCacheLineBytes,
               (kDefaultGrainBody * kBodyBufferBytes) % kCacheLineBytes == 0 ? " (aligned)"
                                                                             : " (MISALIGNED)");
}

// ==================== 阈值/粒度标定 ====================
// 一次 fork-join 约 1.5 us（--bench 可验）。所以阈值不是"拍一个数"：
// 这里在真实场景上把串行 / 各线程数的耗时打出来，阈值取"并行开始稳定不亏"的点。
void run_sweep(const SolveSchedule &base, std::uint64_t topology_version) {
  struct Sample {
    const char *name;
    scene::Scene (*make)();
  };
  const Sample samples[] = {
      {"lattice/4x3x2  (24b)", [] { return scene::lattice(4, 3, 2); }},
      {"lattice/12x10x4(480b)", [] { return scene::lattice(12, 10, 4); }},
  };

  std::println(stderr, "\n[sweep] grain_body × threads (best of 5, ms)");
  for (const auto &s : samples) {
    std::println(stderr, "  {}", s.name);
    for (std::size_t grain : {4u, 8u, 16u, 32u, 64u}) {
      std::println(stderr, "    grain_body={:<3}", grain);
      for (unsigned threads : {1u, 2u, 4u, 8u}) {
        SolveSchedule cfg = base;
        cfg.grain_body    = grain;
        double best       = 1e9;
        for (int rep = 0; rep < 5; ++rep) {
          scene::Scene sc     = s.make();
          sc.topology_version = topology_version;
          task::SerialExecutor     serial;
          task::ThreadPoolExecutor pool(threads);
          double                   ms;
          if (threads == 1) {
            ms = run_case(serial, cfg, false, s.name, "avbd", std::move(sc), AVBDStep{}, 3, 30).ms;
          } else {
            ms = run_case(pool, cfg, false, s.name, "avbd", std::move(sc), AVBDStep{}, 3, 30).ms;
          }
          best = std::min(best, ms);
        }
        std::println(stderr, "      threads={:<2} {:8.3f} ms", threads, best);
      }
    }
  }
}

// ==================== 裸 fork-join 延迟 ====================
// 阈值不能拍脑袋：一次 fork-join 的固定成本（唤醒 worker + 原子取块 + 主线程等待）
// 决定了"多小的活不值得并行"。step 7 的串行快速路径阈值直接用这里的数字定。
void bench_fork_join(unsigned threads) {
  task::ThreadPoolExecutor pool(threads);
  constexpr int            kIters = 2000;
  std::println(stderr, "\n[fork-join latency] threads={} empty body, {} iters", pool.thread_count(),
               kIters);
  for (std::size_t n : {64u, 256u, 1024u, 4096u, 16384u}) {
    for (std::size_t grain : {8u, 32u, 128u}) {
      volatile std::size_t sink = 0;
      const auto           t0   = std::chrono::steady_clock::now();
      for (int i = 0; i < kIters; ++i)
        pool.parallel_for(n, grain, [&sink](std::size_t b, std::size_t e) { sink = b + e; });
      const auto        t1 = std::chrono::steady_clock::now();
      const double      us = std::chrono::duration<double, std::micro>(t1 - t0).count() / kIters;
      const std::size_t chunks = (n + grain - 1) / grain;
      std::println(stderr, "  n={:6} grain={:4} chunks={:5}  {:8.2f} us/call  {:6.2f} us/chunk", n,
                   grain, chunks, us, us / static_cast<double>(chunks));
    }
  }
}

} // namespace

int main(int argc, char **argv) {
  bool        trace            = false;
  bool        selftest         = false;
  bool        bench            = false;
  bool        sweep            = false;
  bool        rebuild_topology = false;
  bool        collision_digest = false;
  std::size_t grain_body = 0, grain_constraint = 0, min_items_body = 0, min_items_constraint = 0;
  unsigned    threads = 1;
  std::string check_path;
  std::string save_path;
  std::string check_collision_path;
  std::string save_collision_path;

  for (int i = 1; i < argc; ++i) {
    const std::string_view a = argv[i];
    if (a == "--trace") {
      trace = true;
    } else if (a == "--selftest") {
      selftest = true;
    } else if (a == "--bench") {
      bench = true;
    } else if (a == "--rebuild-topology") {
      rebuild_topology = true;
    } else if (a == "--sweep") {
      sweep = true;
    } else if (a == "--grain-body" && i + 1 < argc) {
      grain_body = static_cast<std::size_t>(std::strtoul(argv[++i], nullptr, 10));
    } else if (a == "--grain-constraint" && i + 1 < argc) {
      grain_constraint = static_cast<std::size_t>(std::strtoul(argv[++i], nullptr, 10));
    } else if (a == "--min-items-body" && i + 1 < argc) {
      min_items_body = static_cast<std::size_t>(std::strtoul(argv[++i], nullptr, 10));
    } else if (a == "--min-items-constraint" && i + 1 < argc) {
      min_items_constraint = static_cast<std::size_t>(std::strtoul(argv[++i], nullptr, 10));
    } else if (a == "--threads" && i + 1 < argc) {
      threads = static_cast<unsigned>(std::max(1L, std::strtol(argv[++i], nullptr, 10)));
    } else if (a == "--check" && i + 1 < argc) {
      check_path = argv[++i];
    } else if (a == "--save" && i + 1 < argc) {
      save_path = argv[++i];
    } else if (a == "--collision-digest") {
      collision_digest = true;
    } else if (a == "--check-collision" && i + 1 < argc) {
      check_collision_path = argv[++i];
    } else if (a == "--save-collision" && i + 1 < argc) {
      save_collision_path = argv[++i];
    } else {
      std::println(stderr,
                   "unknown argument: {}\nusage: {} [--trace] [--selftest] [--threads N] "
                   "[--check FILE] [--save FILE] [--collision-digest] "
                   "[--check-collision FILE] [--save-collision FILE]",
                   a, argv[0]);
      return 2;
    }
  }

  int exit_code = 0;

  // collision 层可观测摘要：与主基线解耦。它不参与 tests/run.sh 的强制门禁，
  // 只用于“统一 polytope 管线前后到底变了哪些接触”的对比。
  if (collision_digest || !save_collision_path.empty() || !check_collision_path.empty()) {
    const std::vector<std::string> lines = collision_selftest::observable_lines();

    if (collision_digest)
      for (const std::string &line : lines)
        std::println("{}", line);

    if (!check_collision_path.empty()) {
      std::ifstream in(check_collision_path, std::ios::binary);
      if (!in) {
        std::println(stderr, "cannot open {} for reading", check_collision_path);
        return 2;
      }
      std::vector<std::string> expected;
      for (std::string l; std::getline(in, l);)
        if (!l.empty())
          expected.push_back(l);

      int failures = 0;
      for (std::size_t i = 0; i < lines.size(); ++i) {
        if (i >= expected.size()) {
          std::println(stderr, "[FAIL] {} : baseline has only {} lines but {} cases ran",
                       check_collision_path, expected.size(), lines.size());
          ++failures;
          break;
        }
        if (expected[i] != lines[i]) {
          std::println(stderr, "[FAIL-COLLISION] line {}\n   expected: {}\n   actual  : {}",
                       i + 1, expected[i], lines[i]);
          ++failures;
        }
      }
      if (expected.size() > lines.size()) {
        std::println(stderr, "[FAIL-COLLISION] baseline has {} lines, run produced {}",
                     expected.size(), lines.size());
        ++failures;
      }
      if (failures == 0)
        std::println(stderr, "[OK] {} collision cases match {}", lines.size(),
                     check_collision_path);
      else
        exit_code = 1;
    }

    if (!save_collision_path.empty()) {
      std::ofstream out(save_collision_path, std::ios::binary | std::ios::trunc);
      if (!out) {
        std::println(stderr, "cannot open {} for writing", save_collision_path);
        return 2;
      }
      for (const std::string &line : lines)
        out << line << '\n';
      std::println(stderr, "[collision-baseline] wrote {} lines to {}", lines.size(),
                   save_collision_path);
    }

    return exit_code;
  }

  if (bench) {
    print_layout();
    bench_fork_join(threads > 1 ? threads : 8);
  }

  if (selftest) {
    task::SerialExecutor     serial;
    task::ThreadPoolExecutor pool(4);

    std::vector<Check> checks = run_selftest();
    for (auto &extra : run_executor_selftest(serial, pool))
      checks.push_back(std::move(extra));

    int failures = 0;
    for (const auto &c : checks) {
      std::println(stderr, "[{}] {}{}", c.ok ? "PASS" : "FAIL", c.name,
                   c.ok ? std::string{} : std::format("  <-- {}", c.detail));
      if (!c.ok)
        ++failures;
    }
    std::println(stderr, "[selftest] {}/{} passed",
                 checks.size() - static_cast<std::size_t>(failures), checks.size());
    if (failures != 0)
      exit_code = 1;
  }

  // 同一张用例表，只换执行器。--threads 1 走串行参考实现，>1 走线程池。
  // 两条路径的 stdout 必须一模一样（run.sh 会用 --check 在多个线程数下各跑一遍）。
  // --rebuild-topology 把拓扑缓存关掉（等价于 version=0），用来度量缓存的真实收益。
  const std::uint64_t topology_version_override = rebuild_topology ? 0 : 1;

  SolveSchedule sched{};
  if (grain_body)
    sched.grain_body = grain_body;
  if (grain_constraint)
    sched.grain_constraint = grain_constraint;
  if (min_items_body)
    sched.min_items_body = min_items_body;
  if (min_items_constraint)
    sched.min_items_constraint = min_items_constraint;

  if (sweep) {
    run_sweep(sched, topology_version_override);
    return 0;
  }

  std::vector<Result>  results;
  task::SerialExecutor serial_exec;
  if (threads <= 1) {
    results = run_all(serial_exec, sched, trace, topology_version_override);
  } else {
    task::ThreadPoolExecutor pool_exec(threads);
    std::println(stderr, "[executor] thread pool with {} threads", pool_exec.thread_count());
    results = run_all(pool_exec, sched, trace, topology_version_override);
  }

  for (const auto &r : results) {
    std::println("{}", r.line);
    const std::string case_name   = r.line.substr(0, r.line.find(' '));
    const std::string policy_name = r.line.substr(
        case_name.size() + 1, r.line.find(' ', case_name.size() + 1) - case_name.size() - 1);
    std::println(stderr, "[timing] {:<24} {:8.3f} ms", case_name + " " + policy_name, r.ms);
  }

  if (!save_path.empty()) {
    std::ofstream out(save_path, std::ios::binary | std::ios::trunc);
    if (!out) {
      std::println(stderr, "cannot open {} for writing", save_path);
      return 2;
    }
    for (const auto &r : results)
      out << r.line << '\n';
    std::println(stderr, "[baseline] wrote {} lines to {}", results.size(), save_path);
  }

  if (!check_path.empty()) {
    std::ifstream in(check_path, std::ios::binary);
    if (!in) {
      std::println(stderr, "cannot open {} for reading", check_path);
      return 2;
    }
    std::vector<std::string> expected;
    for (std::string l; std::getline(in, l);)
      if (!l.empty())
        expected.push_back(l);

    int failures = 0;
    for (std::size_t i = 0; i < results.size(); ++i) {
      if (i >= expected.size()) {
        std::println(stderr, "[FAIL] {} : baseline has only {} lines but {} cases ran", check_path,
                     expected.size(), results.size());
        ++failures;
        break;
      }
      if (expected[i] != results[i].line) {
        std::println(stderr, "[FAIL] line {}\n   expected: {}\n   actual  : {}", i + 1, expected[i],
                     results[i].line);
        ++failures;
      }
    }
    if (expected.size() > results.size()) {
      std::println(stderr, "[FAIL] baseline has {} lines, run produced {}", expected.size(),
                   results.size());
      ++failures;
    }
    if (failures == 0) {
      std::println(stderr, "[OK] {} cases bit-identical to {}", results.size(), check_path);
    } else {
      exit_code = 1;
    }
  }

  return exit_code;
}
