#pragma once

// ============================================================================
//  collision/tgs_soft.hpp —— Box2D v3 TGS Soft 接触核（速度级）
// ----------------------------------------------------------------------------
//  与 AVBD / B-PXPBD 的根本区别：接触不产生位置乘子。每帧 prepare 一次（拷贝
//  锚点、算 base separation、选 softness、暖启动冲量），然后每子步按
//  integrate velocities → warm start → push (bias) → integrate positions →
//  solve (bias-free + friction) 的顺序推进。本文件只负责接触核；位置积分由
//  TgsSoftStep::begin_iterations 完成（见 phys_core/policy.hpp）。
//
//  【法线约定】本仓库的 manifold normal 指向 B→A（basis[0]），与 Box2D 的
//  normal（A→B）恰好相反。因此 vn = dot(vA - vB, n) 与 Box2D 的
//  dot(vB - vA, n_b2) 数值相同（负 = 接近），正冲量把 A 沿 +n 推开、B 沿 -n
//  推开。分离量 s 的增量项也相应取 (dpA-dpB) 与 (rot(dqA,rA)-rot(dqB,rB))。
//
//  【确定性】所有循环按 manifold / contact 索引升序串行执行；不引入任何线程
//  原语，结果与调度无关。多点的 k<j 耦合通过"逐点立即施加冲量"实现（见 push）。
// ============================================================================

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>

#include "collision/manifold.hpp"
#include "collision/narrowphase.hpp" // kMaxContactsPerManifold
#include "phys_core/state.hpp"

namespace collision {
namespace tgs_soft {

namespace detail {

inline Vec3 contact_velocity(const State &s, std::size_t b, Vec3 r_world) noexcept {
  return s.velocity.velocity(b) + glm::cross(s.velocity.angular_velocity(b), r_world);
}

// ds = (pA - pA0) - (pB - pB0) + (dqA·rA - dqB·rB)，即 Box2D 的 dp/dq 项在
// 本仓库法线约定下的形式。dqLin/dqAng 由 State::cache 相对步首位姿给出。
inline Vec3 separation_delta(const State &s, std::size_t a, std::size_t b,
                             Vec3 rA_world, Vec3 rB_world) noexcept {
  const Quat dqA = quat_from_rotvec(s.cache.dq_angular(a));
  const Quat dqB = quat_from_rotvec(s.cache.dq_angular(b));
  return (s.cache.dq_linear(a) - s.cache.dq_linear(b)) +
         (rotate(dqA, rA_world) - rotate(dqB, rB_world));
}

inline Scalar_t effective_mass(const State &s, std::size_t b, Vec3 g_lin,
                               Vec3 g_ang) noexcept {
  return s.inertial.effective_mass(b, g_lin, g_ang);
}

// 两体有效质量 gᵀM⁻¹g 之和的倒数（= Box2D 的 normalMass / tangentMass）。
// 注意**不能**写成 1/wA + 1/wB：那会把二体接触质量放大 4 倍，堆叠在几步内
// 冲到 4000 m/s。
inline Scalar_t contact_mass(const State &s, std::size_t a, std::size_t b, Vec3 g_lin,
                             Vec3 g_ang_a, Vec3 g_ang_b) noexcept {
  const Scalar_t w = effective_mass(s, a, g_lin, g_ang_a) +
                     effective_mass(s, b, g_lin, g_ang_b);
  return w > Scalar_t(0) ? Scalar_t(1) / w : Scalar_t(0);
}

// 累加一个冲量：P_ours 是把 A 沿 +n 推的向量；角冲量是 r × P。
inline void apply_impulse(State &s, std::size_t a, std::size_t b, Vec3 P,
                          Vec3 ang_a, Vec3 ang_b) noexcept {
  const Scalar_t inv_ma = s.inertial.inv_mass(a);
  const Scalar_t inv_mb = s.inertial.inv_mass(b);
  if (inv_ma > Scalar_t(0)) {
    s.velocity.add_velocity(a, P * inv_ma);
    s.velocity.add_omega(a, s.inertial.world_inv_inertia(a) * ang_a);
  }
  if (inv_mb > Scalar_t(0)) {
    s.velocity.add_velocity(b, P * (-inv_mb));
    s.velocity.add_omega(b, s.inertial.world_inv_inertia(b) * ang_b);
  }
}

} // namespace detail

// 每帧一次：准备锚点/分离量/有效质量/柔度，并把持久冲量施加到速度上。
inline void prepare(State &s, std::span<Manifold> manifolds, Scalar_t margin,
                    const TgsSoftParams &params) {
  if (manifolds.empty())
    return;

  // Box2D：hertz 先夹到 0.125/dt，静态接触再用 2×hertz（"avoid bodies pushed
  // through the ground"）。make_soft 保证 massScale + impulseScale == 1。
  const Scalar_t h     = s.dt;
  const Scalar_t hertz = std::min(params.contact_hertz, Scalar_t(0.125) / h);
  const TgsSoftness dynamic_softness = make_soft(hertz, params.contact_damping_ratio, h);
  const TgsSoftness static_softness  = make_soft(Scalar_t(2) * hertz,
                                                 params.contact_damping_ratio, h);

  for (Manifold &m : manifolds) {
    const std::size_t a = m.body_a;
    const std::size_t b = m.body_b;
    const bool is_static = s.inertial.is_static(a) || s.inertial.is_static(b);
    const TgsSoftness soft = is_static ? static_softness : dynamic_softness;
    const Vec3 n  = m.basis.normal;
    const Vec3 t0 = m.basis.tangent0;
    const Vec3 t1 = m.basis.tangent1;

    for (std::size_t i = 0; i < m.count; ++i) {
      ManifoldContact &c = m.contacts[i];
      c.softness         = soft;
      c.anchorA_world    = s.cache.initial_rot_matrix(a) * c.rA;
      c.anchorB_world    = s.cache.initial_rot_matrix(b) * c.rB;
      // Box2D: baseSeparation = separation - dot(rB - rA, normal). 本仓库的
      // normal 是 B→A，代入后就是下面这个形式；少了锚点项，push 在 dq=0 时
      // 会把 ~1 m 的质心-锚点距离误当成穿透，直接把堆叠发射出去。
      c.base_separation  = (margin - c.depth) -
                           glm::dot(c.anchorA_world - c.anchorB_world, n);

      const Vec3 rna = glm::cross(c.anchorA_world, n);
      const Vec3 rnb = glm::cross(c.anchorB_world, n);
      c.normal_mass = detail::contact_mass(s, a, b, n, rna, rnb);
      for (int k = 0; k < 2; ++k) {
        const Vec3 t   = (k == 0) ? t0 : t1;
        const Vec3 rta = glm::cross(c.anchorA_world, t);
        const Vec3 rtb = glm::cross(c.anchorB_world, t);
        c.tangent_mass[k] = detail::contact_mass(s, a, b, t, rta, rtb);
      }

      // warm start（scale = 1）：持久的 normal/tangent 冲量在本帧起点施加一次。
      c.normal_impulse       = c.lambda[0];
      c.tangent_impulse[0]   = c.lambda[1];
      c.tangent_impulse[1]   = c.lambda[2];
      c.total_normal_impulse = c.normal_impulse;
      // 静态接触用 static softness，但**暖启动冲量仍然要施加到动态侧**：
      // apply_impulse 逐 body 按 inv_mass == 0 跳过静态 body。若这里整体跳过，
      // 累积冲量会逐帧线性增长（静息盒子 5 s 内从 0.1 涨到 12 N·s），支撑
      // 却只靠 push/solve 当帧补，盒子缓慢下沉。
      const Vec3 P = n * c.normal_impulse + t0 * c.tangent_impulse[0] +
                     t1 * c.tangent_impulse[1];
      detail::apply_impulse(s, a, b, P, glm::cross(c.anchorA_world, P),
                            glm::cross(c.anchorB_world, -P));
    }
  }
}

// push（bias = true）：软性去穿透，不带摩擦。每子步一次，在位置积分之前。
//
// 【逐点立即施加】Box2D 对同一 manifold 的多个点用 k<j 耦合：后面的点在算 vn 时
// 要包含前面点的冲量。我们没有实现显式耦合项，而是把每个点的冲量**立即**写回
// body 速度，于是下一点读到的 vn 自然包含前面点的影响——代数上等价，且顺序固定
// （contact 索引升序），仍是确定性的。若先把所有点的冲量加起来最后一次性施加，
// 四点接触会拿到约 4 倍冲量，堆叠立刻爆炸。
inline void push(State &s, std::span<Manifold> manifolds,
                 const TgsSoftParams &params) {
  if (manifolds.empty())
    return;
  const Scalar_t inv_h = Scalar_t(1) / s.dt;

  for (Manifold &m : manifolds) {
    const std::size_t a = m.body_a;
    const std::size_t b = m.body_b;
    const Vec3 n = m.basis.normal;

    for (std::size_t i = 0; i < m.count; ++i) {
      ManifoldContact &c = m.contacts[i];
      const Scalar_t sep = c.base_separation +
              glm::dot(detail::separation_delta(s, a, b, c.anchorA_world,
                                                c.anchorB_world), n);

      Scalar_t velocity_bias, mass_scale, impulse_scale;
      if (sep > Scalar_t(0)) {
        // speculative bias：只阻止本步穿过间隙。no-pull 由下面对累积冲量的
        // 单边夹紧保证（vn >= -bias 时 impulse <= 0，被夹到 -normal_impulse）。
        velocity_bias = sep * inv_h;
        mass_scale    = Scalar_t(1);
        impulse_scale = Scalar_t(0);
      } else {
        velocity_bias = std::max(c.softness.massScale * c.softness.biasRate * sep,
                                 -params.contact_speed);
        mass_scale    = c.softness.massScale;
        impulse_scale = c.softness.impulseScale;
      }

      // 固定锚点上的相对法向速度（本仓库约定：负 = 接近）。
      const Vec3 vA = detail::contact_velocity(s, a, c.anchorA_world);
      const Vec3 vB = detail::contact_velocity(s, b, c.anchorB_world);
      const Scalar_t vn = glm::dot(vA - vB, n);

      const Scalar_t j =
              std::max(-c.normal_mass * (mass_scale * vn + velocity_bias) -
                               impulse_scale * c.normal_impulse,
                       -c.normal_impulse);
      c.normal_impulse += j;
      c.total_normal_impulse += j;
      const Vec3 P = n * j;
      detail::apply_impulse(s, a, b, P, glm::cross(c.anchorA_world, P),
                            glm::cross(c.anchorB_world, -P));
    }
  }
}

// solve（bias-free）：法线 + 双切线摩擦。每子步一次，在位置积分之后。
inline void solve(State &s, std::span<Manifold> manifolds,
                  const TgsSoftParams &) {
  if (manifolds.empty())
    return;
  const Scalar_t inv_h = Scalar_t(1) / s.dt;

  for (Manifold &m : manifolds) {
    const std::size_t a = m.body_a;
    const std::size_t b = m.body_b;
    const Vec3 n  = m.basis.normal;
    const Vec3 t0 = m.basis.tangent0;
    const Vec3 t1 = m.basis.tangent1;
    const Scalar_t mu = std::max(m.friction, Scalar_t(0)) * Scalar_t(2);

    for (std::size_t i = 0; i < m.count; ++i) {
      ManifoldContact &c = m.contacts[i];
      const Scalar_t sep = c.base_separation +
              glm::dot(detail::separation_delta(s, a, b, c.anchorA_world,
                                                c.anchorB_world), n);
      const Scalar_t velocity_bias = sep > Scalar_t(0) ? sep * inv_h : Scalar_t(0);

      const Vec3 vA = detail::contact_velocity(s, a, c.anchorA_world);
      const Vec3 vB = detail::contact_velocity(s, b, c.anchorB_world);
      const Scalar_t vn = glm::dot(vA - vB, n);
      const Scalar_t j = std::max(-c.normal_mass * (vn + velocity_bias),
                                  -c.normal_impulse);
      c.normal_impulse += j;
      c.total_normal_impulse += j;
      {
        const Vec3 P = n * j;
        detail::apply_impulse(s, a, b, P, glm::cross(c.anchorA_world, P),
                              glm::cross(c.anchorB_world, -P));
      }

      // 摩擦：切向冲量夹在 ±μ·normalImpulse 内，累积量不能越界。
      const Scalar_t max_friction = mu * c.normal_impulse;
      for (int k = 0; k < 2; ++k) {
        const Vec3 t = (k == 0) ? t0 : t1;
        const Vec3 va = detail::contact_velocity(s, a, c.anchorA_world);
        const Vec3 vb = detail::contact_velocity(s, b, c.anchorB_world);
        const Scalar_t vt = glm::dot(va - vb, t);
        const Scalar_t jt = std::clamp(-c.tangent_mass[k] * vt,
                                       -max_friction - c.tangent_impulse[k],
                                       max_friction - c.tangent_impulse[k]);
        c.tangent_impulse[k] += jt;
        const Vec3 P = t * jt;
        detail::apply_impulse(s, a, b, P, glm::cross(c.anchorA_world, P),
                              glm::cross(c.anchorB_world, -P));
      }
    }
  }
}

// 把工作副本写回持久冲量，供下一次 narrowphase 的 feature match 暖启动。
inline void store(std::span<Manifold> manifolds) noexcept {
  for (Manifold &m : manifolds) {
    const Scalar_t mu = m.friction;
    for (std::size_t i = 0; i < m.count; ++i) {
      ManifoldContact &c = m.contacts[i];
      c.lambda[0] = c.normal_impulse;
      c.lambda[1] = c.tangent_impulse[0];
      c.lambda[2] = c.tangent_impulse[1];
      const Scalar_t t2 = c.tangent_impulse[0] * c.tangent_impulse[0] +
                          c.tangent_impulse[1] * c.tangent_impulse[1];
      c.stick = c.normal_impulse > Scalar_t(0) &&
                t2 <= (Scalar_t(0.5) * mu * c.normal_impulse) *
                              (Scalar_t(0.5) * mu * c.normal_impulse);
    }
  }
}

} // namespace tgs_soft
} // namespace collision
