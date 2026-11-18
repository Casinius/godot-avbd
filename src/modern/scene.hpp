#pragma once

// ============================================================================
//  scene.hpp —— 确定性测试场景
// ----------------------------------------------------------------------------
//  规则：
//   1) 禁止随机数、禁止时间/哈希序 —— 同样的参数必须产出逐位同样的 State。
//      基线的全部意义就在于此。所有"扰动"都由 body/约束索引的整数表达式算出。
//   2) 场景只通过 State 的命名接口写（set_position/set_inv_mass/…），不允许直接
//      摸 pose.pos[]。
//   3) **场景必须有初始约束违反**（strain != 1）。这条是踩过的坑：本仓库的
//      solve_body 的 RHS 只有 -grad，没有 (x - x_inertial) 惯性项，也没有重力
//      （cache.inertial_position() 目前无人消费）。所以 C≡0 的场景里三个 Policy
//      全都零位移、λ 全 0，摘要会完全相同 —— 那样的"基线"什么都测不出来。
//      保留一个 strain==1 的静止金丝雀用例，专门抓"不该动却动了"的漂移。
// ============================================================================

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "phys_core/constraint.hpp"

namespace scene {

struct Scene {
  const char    *name = "";
  State          state;
  ConstraintList constraints;

  // 拓扑版本号：约束列表在场景构建期一次性建好、之后不再增删 ⇒ 默认 1。
  // **如果**你在拿到 Scene 之后又往 constraints 里 emplace/push，必须 bump 这个值，
  // 否则 gather 会继续用旧邻接（开 SolveSchedule::verify_topology 能抓到这个错）。
  std::uint64_t topology_version = 1;
};

// 统一 body 初始化：位置/朝向/质量/惯量/速度一次写完，避免漏字段导致
// 摘要里出现"某些 body 恰好全零"的假通过。
inline void add_body(State &s, std::size_t b, Vec3 p, Scalar_t m_inv, Vec3 v = Vec3{0},
                     Vec3 omega = Vec3{0}, Quat q = identity_quat()) noexcept {
  s.pose.set_position(b, p);
  s.pose.set_rotation(b, q);
  s.velocity.set_velocity(b, v);
  s.velocity.set_omega(b, omega);
  s.inertial.set_inv_mass(b, m_inv);
  if (m_inv > 0.0f)
    s.inertial.set_local_inv_inertia(b, Mat3(1));
}

// 确定性"抖动"：把约束静止长度按索引打散成 0.8x~1.2x。
// 作用有两个：让每个约束的违反量互不相同（否则所有 C 相等时，累加顺序对浮点
// 结果没有区分度），以及让不同位置的约束给出不同的 λ。
inline constexpr Scalar_t jitter(std::size_t i) noexcept {
  return Scalar_t(0.8) + Scalar_t(0.1) * Scalar_t(i % 5); // 0.8, 0.9, 1.0, 1.1, 1.2
}

// —— 链条：n 个 body 沿 x 排开，body 0 是静态锚点 ——
// strain < 1 → 链条被"压紧"，约束初始违反；
// strain == 1 且 vary == false → 位置严格等距且 rest 就是间距，C 恒为 0（金丝雀）。
// vary 必须显式关掉：jitter 会乘在 rest 上，光把 strain 设成 1 并不会让 C 归零。
inline Scene chain(std::size_t n, Scalar_t spacing = Scalar_t(0.5),
                   Scalar_t strain = Scalar_t(0.85), bool vary = true) {
  Scene sc;
  sc.name     = "chain";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(n);
  for (std::size_t b = 0; b < n; ++b) {
    const Scalar_t x      = Scalar_t(b) * spacing;
    const bool     anchor = (b == 0);
    add_body(sc.state, b, Vec3{x, 0, 0}, anchor ? 0.0f : 1.0f, anchor ? Vec3{0} : Vec3{0, 0, 0},
             anchor ? Vec3{0} : Vec3{0, 0, Scalar_t(0.1) * Scalar_t(b + 1)});
  }
  for (std::size_t b = 0; b + 1 < n; ++b) {
    const Scalar_t vary_k = vary ? jitter(b) : Scalar_t(1);
    sc.constraints.emplace_back(DistanceConstraint{{b, b + 1}, spacing * strain * vary_k});
  }
  return sc;
}

// —— 扇面：1 个静态中心 + n 个动态叶子，一圈排开 ——
// 这是 accumulate 阶段的写冲突热点：中心 body 被 n 条约束同时累加。
// 多线程化后只要这里不竞态，别处基本不会竞态。
inline Scene fan(std::size_t n, Scalar_t radius = Scalar_t(1)) {
  Scene sc;
  sc.name                 = "fan";
  const Scalar_t dt_angle = Scalar_t(6.2831853) / Scalar_t(n);
  sc.state.dt             = Scalar_t(1) / 60;
  sc.state.resize(n + 1);

  add_body(sc.state, 0, Vec3{0, 0, 0}, 0.0f); // 静态枢纽
  for (std::size_t i = 0; i < n; ++i) {
    const Scalar_t a = dt_angle * Scalar_t(i);
    // std::cos/std::sin 是确定性的（同 libm 同参数同结果），但为了让不同 libm 的
    // 构建之间也不互相污染基线，这里只在 x-y 上加一点由整数索引算出的偏移。
    const Vec3 p{std::cos(a) * radius, std::sin(a) * radius, Scalar_t(0.05) * Scalar_t(i % 3)};
    add_body(sc.state, i + 1, p, 1.0f, Vec3{0, 0, 0}, Vec3{0, 0, Scalar_t(0.02) * Scalar_t(i % 4)});
  }
  for (std::size_t i = 0; i < n; ++i)
    sc.constraints.emplace_back(DistanceConstraint{{0, i + 1}, radius * Scalar_t(0.9) * jitter(i)});
  return sc;
}

// —— 关节格：nx*ny*nz 网格，沿三个轴向各连一条距离约束 ——
// 每个格点都被 3~6 条约束共享，是 accumulate 冲突密度的上限场景，也是后面
// 扩展性测试（1k / 10k body）的载体。底面一层静当作地面。
inline Scene lattice(std::size_t nx, std::size_t ny, std::size_t nz,
                     Scalar_t spacing = Scalar_t(0.5), Scalar_t strain = Scalar_t(0.9)) {
  Scene sc;
  sc.name             = "lattice";
  sc.state.dt         = Scalar_t(1) / 60;
  const std::size_t n = nx * ny * nz;
  sc.state.resize(n);

  auto index = [&](std::size_t i, std::size_t j, std::size_t k) { return (k * ny + j) * nx + i; };

  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t b = index(i, j, k);
        const Vec3        p{Scalar_t(i) * spacing, Scalar_t(j) * spacing, Scalar_t(k) * spacing};
        const Scalar_t    m_inv = (j == 0) ? 0.0f : 1.0f; // 底面一层静态
        const Vec3        v{0, 0, Scalar_t(0.05) * Scalar_t((b % 7) - 3)};
        add_body(sc.state, b, p, m_inv, v, Vec3{0, 0, Scalar_t(0.02) * Scalar_t(b % 5)});
      }

  for (std::size_t k = 0; k < nz; ++k)
    for (std::size_t j = 0; j < ny; ++j)
      for (std::size_t i = 0; i < nx; ++i) {
        const std::size_t b = index(i, j, k);
        if (i + 1 < nx)
          sc.constraints.emplace_back(
              DistanceConstraint{{b, index(i + 1, j, k)}, spacing * strain * jitter(i + j + k)});
        if (j + 1 < ny)
          sc.constraints.emplace_back(
              DistanceConstraint{{b, index(i, j + 1, k)}, spacing * strain * jitter(i + 2 * j)});
        if (k + 1 < nz)
          sc.constraints.emplace_back(
              DistanceConstraint{{b, index(i, j, k + 1)}, spacing * strain * jitter(j + 3 * k)});
      }
  return sc;
}

// —— 混合：Distance 与 Contact 交错 ——
// ContactConstraint 在此仓库里是"预给法线+深度的单边约束"（无碰撞检测），
// 所以这里手动造出法线与深度。重点是让两种约束类型同时活跃在同一条热路径上，
// 且 Distance 侧也带 strain（=1 的话距离约束会被接触推着走，反而不是初始违反）。
inline Scene mixed(std::size_t n, Scalar_t spacing = Scalar_t(0.4),
                   Scalar_t strain = Scalar_t(0.8)) {
  Scene sc;
  sc.name     = "mixed";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(n);
  for (std::size_t b = 0; b < n; ++b) {
    const Vec3 p{Scalar_t(b) * spacing, Scalar_t(0), Scalar_t(0)};
    add_body(sc.state, b, p, (b == 0) ? 0.0f : 1.0f, Vec3{0, Scalar_t(-0.1) * Scalar_t(b % 3), 0},
             Vec3{0, 0, Scalar_t(0.03) * Scalar_t(b % 4)});
  }
  for (std::size_t b = 0; b + 1 < n; ++b) {
    sc.constraints.emplace_back(DistanceConstraint{{b, b + 1}, spacing * strain * jitter(b)});
    // 每两个 body 加一条接触：body b 相对静态锚点 body 0 的 +y 支撑面。
    // 深度随索引变化，避免“全同深度”掩盖梯度 bug。
    // basis 必须完整给出（normal + 两条正交切线），rA/rB 是局部锚点，这里取质心。
    if (b % 2 == 1) {
      ContactConstraint cc;
      cc.bodies_       = {b, 0};
      cc.manifold_id   = kStandaloneManifold; // 不参与 manifold 分组
      cc.contact_index = 0;
      cc.kind          = ContactKind::Normal;
      cc.basis[0]      = Vec3{0, 1, 0}; // normal
      cc.basis[1]      = Vec3{1, 0, 0}; // tangent0
      cc.basis[2]      = Vec3{0, 0, 1}; // tangent1
      cc.rA            = Vec3{0};
      cc.rB            = Vec3{0};
      cc.target        = Scalar_t(0.01) * Scalar_t(1 + b % 5); // depth
      sc.constraints.emplace_back(cc);
    }
  }
  return sc;
}

// —— 自由落体：无约束，只测重力（=默认均匀外力）是否真的进了惯性目标 ——
// 这是唯一能用解析解验证的场景：半隐式积分下 x_n = x0 + g·dt²·n(n+1)/2。
inline Scene freefall(std::size_t n, Scalar_t height = Scalar_t(10)) {
  Scene sc;
  sc.name     = "freefall";
  sc.state.dt = Scalar_t(1) / 60;
  sc.state.resize(n);
  for (std::size_t b = 0; b < n; ++b)
    add_body(sc.state, b, Vec3{Scalar_t(b) * Scalar_t(0.3), height, 0}, 1.0f);
  return sc;
}

// 零重力：把 uniform_accel 置零。仅用于"不该动却动了"的金丝雀。
inline void make_gravity_free(Scene &sc) { sc.state.external.uniform_accel = Vec3{0}; }

// 完全静默：零重力 + 零速度。
// 为什么必须也清速度：一旦求解器有了惯性残差项，非零 ω 会通过
// rotvec(θ_pred → θ) 合法地把 body 转过去。那不是 bug，是角速度积分。
// 所以“C==0 且不动”这个金丝雀只有在速度也归零时才能成立。
inline void make_quiescent(Scene &sc) {
  make_gravity_free(sc);
  for (std::size_t b = 0; b < sc.state.size(); ++b) {
    sc.state.velocity.set_velocity(b, Vec3{0});
    sc.state.velocity.set_omega(b, Vec3{0});
  }
}

} // namespace scene
