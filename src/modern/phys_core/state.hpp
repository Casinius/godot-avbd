#pragma once

// ============================================================================
//  state.hpp —— 求解器状态（SoA，按语义分组成 7 个明确接口）
// ----------------------------------------------------------------------------
//  【不变量 / invariants】—— 改动本文件前先读这段
//   1) pose.pos/rot 是唯一"真值"；cache 里的 dq* 都是它相对 initial* 的差，
//      任何写 pose 的地方都必须紧跟 cache.refresh_dq(b, pose)。
//   2) cache.dq_linear(b)  == pose.position(b) - cache.initial_position(b)
//      cache.dq_angular(b) == rotvec(cache.initial_rotation(b), pose.rotation(b))
//   3) inertial.inertia_world_inv 只在 begin_step() / 显式 refresh 时刷新，
//      即 R * I_local⁻¹ * Rᵀ（R 来自当前 pose，不是转置侧乘）。
//   4) buffers 是"每步累积、每步清空"的累加器：reset() 必须清全部 5 个字段，
//      且只能由 batch_update 在消费完之后调用。
//   5) mass_inv == 0 表示静态 body：inv inertia 也必须为 0，速度恒为 0，
//      solve_body 对它返回"无位移"。
//   6) 求解器只调用命名接口（position()/dq_linear()/accumulate_penalty()…），
//      不直接对 vector 下标 —— 这样换 SoA/AoS 后端时调用点不用改。
//   7) external 是"每步复位、每步消费"的外力累加器。**重力不是求解器里的一个
//      分支**，它只是 external.uniform_accel 的默认值 —— 也就是说重力跟风、爆炸
//      力走完全相同的 add_force/add_torque 路径，求解器永远不知道哪一份来自重力。
//      顺序硬约束：driver 必须先 reset_external()（+ 可选累加外力）再 begin_step()；
//      begin_step() 只消费、不复位。
//   8) cache.inertial_* 是 begin_step() 算出的惯性目标 x_pred。它的**唯一**消费者
//      是 solve_body 里的残差项 -M/dt²·(x - x_pred)。约束累积阶段禁止读它 ——
//      一读，求解器就"看见"了重力，重力就重新变成特殊力了。
//   9) 静态 body 的 external 加速度恒为 0、速度恒为 0。这是累加器与速度积分的
//      通用规则（无限质量 ⇒ 任何力都不产生加速度），不是对重力的特判。
//  10) 【确定性契约】同一份二进制、同一份场景、同一串调用顺序下，结果对以下三者
//      **完全不变**（逐位相同）：线程数、分区粒度 grain_*、串行/并行阈值 min_items_*。
//      保证来自三件事，缺一不可：
//        a) scatter 只写自己的 slot（无冲突，写入顺序无关）；
//        b) gather 的加总顺序由 body_entries（约束升序）决定，不由调度决定；
//        c) 串行快速路径与并行路径执行**同一段代码**（run_phase 只改分组，不改算术）。
//      **不**保证：跨编译器/跨 -O 级别/跨 FMA 设置逐位一致（浮点收缩会变 rounding），
//      也不保证 XPBD/PBD 能被并行化后仍与今天一致 —— 它们是 Gauss-Seidel，
//      并行化必须换算法（着色/Jacobi），那属于行为变更。
// ============================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <span>
#include <vector>

// ==================== 0. 标量与几何别名 ====================
using Scalar_t = float;
using Vec3 = glm::vec<3, Scalar_t>;
using Mat3 = glm::mat<3, 3, Scalar_t>;
using Quat = glm::qua<Scalar_t>;

// ==================== 1. 数学小工具（glm compat 层） ====================
inline Quat identity_quat() noexcept { return glm::identity<Quat>(); }

inline Mat3 outer(Vec3 a, Vec3 b) noexcept {
  return Mat3{a.x * b.x, a.y * b.x, a.z * b.x,
              a.x * b.y, a.y * b.y, a.z * b.y,
              a.x * b.z, a.y * b.z, a.z * b.z};
}

inline Mat3 rot_matrix(const Quat &q) noexcept { return glm::mat3_cast(q); }
inline Vec3 rotate(const Quat &q, Vec3 v) noexcept { return glm::mat3_cast(q) * v; }
inline Mat3 transpose(const Mat3 &m) noexcept { return glm::transpose(m); }

// 叉乘矩阵：skew(w) * v == cross(w, v)
inline Mat3 skew(Vec3 w) noexcept {
  return Mat3{0.0f, w.z, -w.y,
              -w.z, 0.0f, w.x,
              w.y, -w.x, 0.0f};
}

// 旋转向量 -> 四元数（指数映射，小角度走一阶展开避免 0/0）
inline Quat quat_from_rotvec(Vec3 rv) noexcept {
  const Scalar_t angle = glm::length(rv);
  if (angle < Scalar_t(1e-8))
    return glm::normalize(Quat{1.0f, rv.x * 0.5f, rv.y * 0.5f, rv.z * 0.5f});
  const Vec3 axis = rv / angle;
  const Scalar_t half = angle * Scalar_t(0.5);
  return Quat{std::cos(half), axis * std::sin(half)};
}

// 从 from 转到 to 的最短弧旋转向量
inline Vec3 rotvec_between(const Quat &from, const Quat &to) noexcept {
  Quat d = glm::normalize(to * glm::conjugate(from));
  if (d.w < 0.0f) d = -d;
  const Vec3 v{d.x, d.y, d.z};
  const Scalar_t s = glm::length(v);
  if (s < Scalar_t(1e-8)) return v * Scalar_t(2);
  return v * (Scalar_t(2) * std::atan2(s, d.w) / s);
}

// 世界系逆惯量：R * I_local⁻¹ * Rᵀ
inline Mat3 inertia_world_inv_from_local(const Mat3 &I_local_inv, const Mat3 &R) noexcept {
  return R * I_local_inv * transpose(R);
}

// 定长 6x6 对称正定/半正定系统的 LDLᵀ 求解（不依赖 GLM 的 mat6，GLM 只到 4x4）。
// 返回 false 表示遇到非正主元或解不合理，调用方应冻结该 DOF。
inline bool ldlt6_solve(const std::array<std::array<Scalar_t, 6>, 6> &a,
                        const std::array<Scalar_t, 6> &b,
                        std::array<Scalar_t, 6> &x) noexcept {
  std::array<std::array<Scalar_t, 6>, 6> l{};
  std::array<Scalar_t, 6> d{};
  for (int j = 0; j < 6; ++j) {
    Scalar_t dj = a[j][j];
    for (int k = 0; k < j; ++k) dj -= l[j][k] * l[j][k] * d[k];
    if (!(dj > Scalar_t(1e-12))) return false;
    d[j] = dj;
    l[j][j] = 1.0f;
    for (int i = j + 1; i < 6; ++i) {
      Scalar_t lij = a[i][j];
      for (int k = 0; k < j; ++k) lij -= l[i][k] * l[j][k] * d[k];
      l[i][j] = lij / dj;
    }
  }
  std::array<Scalar_t, 6> z{};
  for (int i = 0; i < 6; ++i) {
    Scalar_t v = b[i];
    for (int k = 0; k < i; ++k) v -= l[i][k] * z[k];
    z[i] = v;
  }
  for (int i = 0; i < 6; ++i) x[i] = z[i] / d[i];
  for (int i = 5; i >= 0; --i) {
    Scalar_t v = x[i];
    for (int k = i + 1; k < 6; ++k) v -= l[k][i] * x[k];
    x[i] = v;
  }
  return true;
}

inline bool all_finite(Vec3 v) noexcept {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

// ==================== 2. PoseState：位置 / 朝向 ====================
struct PoseState {
  std::vector<Vec3> pos;
  std::vector<Quat> rot;

  void resize(std::size_t n) {
    pos.assign(n, Vec3{0});
    rot.assign(n, identity_quat());
  }
  std::size_t size() const noexcept { return pos.size(); }
  void clear() {
    pos.clear();
    rot.clear();
  }

  Vec3 position(std::size_t b) const noexcept { return pos[b]; }
  Quat rotation(std::size_t b) const noexcept { return rot[b]; }
  Mat3 rotation_matrix(std::size_t b) const noexcept { return rot_matrix(rot[b]); }

  void set_position(std::size_t b, Vec3 p) noexcept { pos[b] = p; }
  void set_rotation(std::size_t b, Quat q) noexcept { rot[b] = glm::normalize(q); }
  void add_position(std::size_t b, Vec3 d) noexcept { pos[b] += d; }

  // 线位移 + 旋转向量增量（求解器唯一的写入口）
  void apply_delta(std::size_t b, Vec3 d_lin, Vec3 d_ang) noexcept {
    pos[b] += d_lin;
    if (glm::dot(d_ang, d_ang) > 0.0f)
      rot[b] = glm::normalize(quat_from_rotvec(d_ang) * rot[b]);
  }
};

// ==================== 3. InertialState：质量 / 惯量 ====================
struct InertialState {
  std::vector<Scalar_t> mass_inv;
  std::vector<Mat3> inertia_local_inv;
  std::vector<Mat3> inertia_world_inv;

  void resize(std::size_t n) {
    mass_inv.assign(n, Scalar_t(0));
    inertia_local_inv.assign(n, Mat3(0));
    inertia_world_inv.assign(n, Mat3(0));
  }
  std::size_t size() const noexcept { return mass_inv.size(); }
  void clear() {
    mass_inv.clear();
    inertia_local_inv.clear();
    inertia_world_inv.clear();
  }

  Scalar_t inv_mass(std::size_t b) const noexcept { return mass_inv[b]; }
  bool is_static(std::size_t b) const noexcept { return mass_inv[b] == 0.0f; }

  void set_inv_mass(std::size_t b, Scalar_t m_inv) noexcept {
    mass_inv[b] = m_inv;
    if (m_inv == 0.0f) inertia_local_inv[b] = Mat3(0); // 静态：逆惯量也为 0
  }
  void set_local_inv_inertia(std::size_t b, const Mat3 &I_inv) noexcept {
    inertia_local_inv[b] = I_inv;
  }
  Mat3 local_inv_inertia(std::size_t b) const noexcept { return inertia_local_inv[b]; }
  Mat3 world_inv_inertia(std::size_t b) const noexcept { return inertia_world_inv[b]; }

  void refresh_world_inertia(std::size_t b, const Mat3 &R) noexcept {
    inertia_world_inv[b] = inertia_world_inv_from_local(inertia_local_inv[b], R);
  }
  void refresh_all(const PoseState &pose) {
    for (std::size_t b = 0; b < mass_inv.size(); ++b)
      refresh_world_inertia(b, pose.rotation_matrix(b));
  }

  // gᵀ M⁻¹ g —— XPBD / PBD 的分母
  Scalar_t effective_mass(std::size_t b, Vec3 g_lin, Vec3 g_ang) const noexcept {
    return glm::dot(g_lin, g_lin) * mass_inv[b] +
           glm::dot(g_ang, inertia_world_inv[b] * g_ang);
  }

  // 惯性预测位姿：p + v·dt，q ⊗ exp(ω·dt)；静态 body 原地不动
  void predict_pose(std::size_t b, const PoseState &pose, Vec3 v, Vec3 omega,
                    Scalar_t dt, Vec3 &out_lin, Quat &out_ang) const noexcept {
    if (is_static(b)) {
      out_lin = pose.position(b);
      out_ang = pose.rotation(b);
      return;
    }
    out_lin = pose.position(b) + v * dt;
    out_ang = glm::normalize(quat_from_rotvec(omega * dt) * pose.rotation(b));
  }
};

// ==================== 4. ExternalState：外力累加器（重力只是它的默认值）====================
// 设计要点（不要退化）：
//   * 求解器只读 accel_lin/accel_ang，**从不读 uniform_accel**。
//   * uniform_accel 只是每步复位时的默认值。“加重力”和“加风力”在求解器眼里
//     没有任何区别，因为到了那个阶段两者都只是 accel 里的几个 float。
//   * 存的是加速度（已乘过 m⁻¹），不是力：这样复位时不需要拿质量做除法，
//     静态 body 也不需要特判（无限质量的加速度本来就是 0）。
struct ExternalState {
  std::vector<Vec3> accel_lin; // 已乘过 m⁻¹ 的线加速度 (m/s²)
  std::vector<Vec3> accel_ang; // 已乘过 Iw⁻¹ 的角加速度 (rad/s²)

  // “重力就是默认的均匀加速度”。想关掉重力就置零；想加风力就另调 add_force。
  Vec3 uniform_accel = Vec3{0.0f, -9.81f, 0.0f};

  void resize(std::size_t n) {
    accel_lin.assign(n, Vec3{0});
    accel_ang.assign(n, Vec3{0});
  }
  std::size_t size() const noexcept { return accel_lin.size(); }
  void clear() {
    accel_lin.clear();
    accel_ang.clear();
  }

  // 每步开头的复位：动态 body 置为 uniform_accel，静态 body 保持 0。
  // 静态清零是累加器的通用规则，不是对重力的特判。
  void reset(const InertialState &inertial) noexcept {
    for (std::size_t b = 0; b < accel_lin.size(); ++b) {
      accel_lin[b] = inertial.is_static(b) ? Vec3{0} : uniform_accel;
      accel_ang[b] = Vec3{0};
    }
  }

  Vec3 linear(std::size_t b) const noexcept { return accel_lin[b]; }
  Vec3 angular(std::size_t b) const noexcept { return accel_ang[b]; }
  void add_linear(std::size_t b, Vec3 a) noexcept { accel_lin[b] += a; }
  void add_angular(std::size_t b, Vec3 a) noexcept { accel_ang[b] += a; }
};

// ==================== 5. VelocityState：线速度 / 角速度 ====================
struct BodyCache; // 前置声明：from_pose_delta 需要它

struct VelocityState {
  std::vector<Vec3> v;
  std::vector<Vec3> omega;

  void resize(std::size_t n) {
    v.assign(n, Vec3{0});
    omega.assign(n, Vec3{0});
  }
  std::size_t size() const noexcept { return v.size(); }
  void clear() {
    v.clear();
    omega.clear();
  }

  Vec3 velocity(std::size_t b) const noexcept { return v[b]; }
  Vec3 angular_velocity(std::size_t b) const noexcept { return omega[b]; }
  void set_velocity(std::size_t b, Vec3 vel) noexcept { v[b] = vel; }
  void set_omega(std::size_t b, Vec3 w) noexcept { omega[b] = w; }
  void add_velocity(std::size_t b, Vec3 dv) noexcept { v[b] += dv; }
  void add_omega(std::size_t b, Vec3 dw) noexcept { omega[b] += dw; }
  void damp(std::size_t b, Scalar_t factor) noexcept {
    v[b] *= factor;
    omega[b] *= factor;
  }
  void zero(std::size_t b) noexcept {
    v[b] = Vec3{0};
    omega[b] = Vec3{0};
  }

  // 外力积分：v += dt·a，ω += dt·α。这是"惯性"进入求解器的唯一入口 ——
  // 它发生在预测之前，所以重力最终被烧进了 cache 的惯性目标，而不是变成
  // 求解循环里的一个额外项。静态 body 强制归零（invariant 5/9）。
  void apply_acceleration(const ExternalState &ext, const InertialState &inertial,
                          Scalar_t dt) noexcept {
    for (std::size_t b = 0; b < v.size(); ++b) {
      if (inertial.is_static(b)) {
        v[b] = Vec3{0};
        omega[b] = Vec3{0};
        continue;
      }
      v[b] += ext.linear(b) * dt;
      omega[b] += ext.angular(b) * dt;
    }
  }

  // 由位置增量反推速度（步末调用）。逐 body 版供并行/融合调用，算术相同。
  void from_pose_delta(std::size_t b, const PoseState &pose, const BodyCache &cache,
                       Scalar_t dt) noexcept;
  void from_pose_delta(const PoseState &pose, const BodyCache &cache, Scalar_t dt);
};

// ==================== 6. BodyCache：步内派生量（只读缓存） ====================
struct BodyCache {
  std::vector<Vec3> initialLin; // step 起始位置
  std::vector<Quat> initialAng; // step 起始朝向
  std::vector<Mat3> initialRotMat;
  std::vector<Vec3> inertialLin; // 惯性/预测位姿
  std::vector<Quat> inertialAng;
  std::vector<Mat3> rotMat;      // 当前朝向的矩阵形式
  std::vector<Vec3> dqLin;       // pos      - initialLin
  std::vector<Vec3> dqAng;       // 旋转向量(initialAng -> rot)
  std::vector<Vec3> prevVelocityLin; // extrapolation / CCD 用
  std::vector<Vec3> prevOmega;

  void resize(std::size_t n) {
    initialLin.assign(n, Vec3{0});
    initialAng.assign(n, identity_quat());
    initialRotMat.assign(n, Mat3(1));
    inertialLin.assign(n, Vec3{0});
    inertialAng.assign(n, identity_quat());
    rotMat.assign(n, Mat3(1));
    dqLin.assign(n, Vec3{0});
    dqAng.assign(n, Vec3{0});
    prevVelocityLin.assign(n, Vec3{0});
    prevOmega.assign(n, Vec3{0});
  }
  std::size_t size() const noexcept { return initialLin.size(); }
  void clear() {
    initialLin.clear();
    initialAng.clear();
    initialRotMat.clear();
    inertialLin.clear();
    inertialAng.clear();
    rotMat.clear();
    dqLin.clear();
    dqAng.clear();
    prevVelocityLin.clear();
    prevOmega.clear();
  }

  // —— 写：只在明确的生命周期点调用 ——
  void capture_initial(std::size_t b, const PoseState &pose) noexcept {
    initialLin[b] = pose.position(b);
    initialAng[b] = pose.rotation(b);
    initialRotMat[b] = pose.rotation_matrix(b);
  }

  void set_inertial(std::size_t b, Vec3 lin, Quat ang) noexcept {
    inertialLin[b] = lin;
    inertialAng[b] = ang;
  }

  void refresh_rotmat(std::size_t b, const PoseState &pose) noexcept {
    rotMat[b] = pose.rotation_matrix(b);
  }

  void refresh_dq(std::size_t b, const PoseState &pose) noexcept {
    dqLin[b] = pose.position(b) - initialLin[b];
    dqAng[b] = rotvec_between(initialAng[b], pose.rotation(b));
  }

  // 逐 body 版：并行时按 body 分块调用。算术序列与整循环版完全相同 ⇒ 位一致。
  void refresh_derived(std::size_t b, const PoseState &pose, Scalar_t dt) noexcept {
    refresh_rotmat(b, pose);
    refresh_dq(b, pose);
    prevVelocityLin[b] = dqLin[b] / dt;
    prevOmega[b] = dqAng[b] / dt;
  }
  void begin_step(const PoseState &pose, const VelocityState &vel,
                  const InertialState &inertial, Scalar_t dt) {
    if (size() != pose.size()) resize(pose.size());
    for (std::size_t b = 0; b < size(); ++b) {
      capture_initial(b, pose);
      rotMat[b] = initialRotMat[b];
      dqLin[b] = Vec3{0};
      dqAng[b] = Vec3{0};
      Vec3 lin;
      Quat ang;
      inertial.predict_pose(b, pose, vel.velocity(b), vel.angular_velocity(b), dt, lin, ang);
      set_inertial(b, lin, ang);
    }
  }
  void commit_prev(const VelocityState &vel) noexcept {
    for (std::size_t b = 0; b < size(); ++b) {
      prevVelocityLin[b] = vel.velocity(b);
      prevOmega[b] = vel.angular_velocity(b);
    }
  }

  // —— 读 ——
  Vec3 initial_position(std::size_t b) const noexcept { return initialLin[b]; }
  Quat initial_rotation(std::size_t b) const noexcept { return initialAng[b]; }
  Mat3 initial_rot_matrix(std::size_t b) const noexcept { return initialRotMat[b]; }
  Vec3 inertial_position(std::size_t b) const noexcept { return inertialLin[b]; }
  Quat inertial_rotation(std::size_t b) const noexcept { return inertialAng[b]; }
  Mat3 rot_matrix(std::size_t b) const noexcept { return rotMat[b]; }
  Vec3 dq_linear(std::size_t b) const noexcept { return dqLin[b]; }
  Vec3 dq_angular(std::size_t b) const noexcept { return dqAng[b]; }
  Vec3 prev_velocity(std::size_t b) const noexcept { return prevVelocityLin[b]; }
  Vec3 prev_omega(std::size_t b) const noexcept { return prevOmega[b]; }
};

inline void VelocityState::from_pose_delta(std::size_t b, const PoseState &pose,
                                           const BodyCache &cache, Scalar_t dt) noexcept {
  const Scalar_t inv_dt = Scalar_t(1) / dt;
  v[b] = (pose.position(b) - cache.initial_position(b)) * inv_dt;
  omega[b] = rotvec_between(cache.initial_rotation(b), pose.rotation(b)) * inv_dt;
}

inline void VelocityState::from_pose_delta(const PoseState &pose, const BodyCache &cache,
                                           Scalar_t dt) {
  for (std::size_t b = 0; b < v.size(); ++b) from_pose_delta(b, pose, cache, dt);
}

// 惯性目标残差：线向 x - x_pred，角向 rotvec(θ_pred -> θ(x))。
// 存在的意义就是让 solve_body 不去碰 pose / cache，只吃这个小结构体 ——
// 并发化时"算残差"和"解 6x6"可以被拆成两个阶段。
struct BodyResiduals {
  Vec3 lin{0};
  Vec3 ang{0};
};

// ==================== 7. SolverBuffers：per-body 局部线性系统 ====================
struct SolverBuffers {
  std::vector<Mat3> hessian_lin;     // ∂²E/∂p∂p
  std::vector<Mat3> hessian_lin_ang; // ∂²E/∂p∂θ
  std::vector<Mat3> hessian_ang;     // ∂²E/∂θ∂θ
  std::vector<Vec3> grad_lin;
  std::vector<Vec3> grad_ang;

  void resize(std::size_t n) {
    hessian_lin.assign(n, Mat3(0));
    hessian_lin_ang.assign(n, Mat3(0));
    hessian_ang.assign(n, Mat3(0));
    grad_lin.assign(n, Vec3{0});
    grad_ang.assign(n, Vec3{0});
  }
  std::size_t size() const noexcept { return hessian_lin.size(); }
  void clear() {
    hessian_lin.clear();
    hessian_lin_ang.clear();
    hessian_ang.clear();
    grad_lin.clear();
    grad_ang.clear();
  }

  // 唯一清空入口：5 个字段一起清，杜绝"漏一行"
  void reset() noexcept {
    for (auto &h : hessian_lin) h = Mat3(0);
    for (auto &h : hessian_lin_ang) h = Mat3(0);
    for (auto &h : hessian_ang) h = Mat3(0);
    for (auto &g : grad_lin) g = Vec3{0};
    for (auto &g : grad_ang) g = Vec3{0};
  }
  void reset(std::size_t b) noexcept {
    hessian_lin[b] = Mat3(0);
    hessian_lin_ang[b] = Mat3(0);
    hessian_ang[b] = Mat3(0);
    grad_lin[b] = Vec3{0};
    grad_ang[b] = Vec3{0};
  }

  Mat3 hessian(std::size_t b) const noexcept { return hessian_lin[b]; }
  Mat3 cross_hessian(std::size_t b) const noexcept { return hessian_lin_ang[b]; }
  Mat3 angular_hessian(std::size_t b) const noexcept { return hessian_ang[b]; }
  Vec3 gradient(std::size_t b) const noexcept { return grad_lin[b]; }
  Vec3 angular_gradient(std::size_t b) const noexcept { return grad_ang[b]; }

  // 罚项累积：H += k·g gᵀ，grad += k·C·g
  void accumulate_penalty(std::size_t b, Vec3 g_lin, Vec3 g_ang, Scalar_t k,
                          Scalar_t C) noexcept {
    hessian_lin[b] += outer(g_lin, g_lin) * k;
    hessian_lin_ang[b] += outer(g_lin, g_ang) * k;
    hessian_ang[b] += outer(g_ang, g_ang) * k;
    grad_lin[b] += g_lin * (k * C);
    grad_ang[b] += g_ang * (k * C);
  }

  // 块 Newton 步：
  //     (H + diag(m⁻¹, Iʷ⁻¹)/dt²) Δ = -(grad + diag(m⁻¹, Iʷ⁻¹)/dt² · r)
  // 第二项 r 就是惯性。没有约束（grad==0）时 Δ = -r，即 x ← x_pred，
  // 自由落体正是靠它成立的。注意这里**没有重力**：重力早已在 begin_step()
  // 之前被烧进了 x_pred，在这里出现的只会是"离 x_pred 有多远"。
  bool solve_body(std::size_t b, const InertialState &inertial, const BodyResiduals &r, Scalar_t dt,
                  Vec3 &d_lin, Vec3 &d_ang) const noexcept {
    d_lin = Vec3{0};
    d_ang = Vec3{0};
    if (inertial.is_static(b)) return false;

    const Scalar_t inv_dt2 = Scalar_t(1) / (dt * dt);
    const Scalar_t m_dt2 = inertial.inv_mass(b) * inv_dt2;
    const Mat3 a_ll = hessian_lin[b] + Mat3(m_dt2);
    const Mat3 a_lt = hessian_lin_ang[b];
    const Mat3 a_tt = hessian_ang[b] + inertial.world_inv_inertia(b) * inv_dt2;

    // 按块装配 6x6 行主序矩阵：A(lin,ang) = a_lt，A(ang,lin) = a_ltᵀ。
    // glm 矩阵是列主序（m[col][row]），这里一次性转成 row-major 免得再绕。
    std::array<std::array<Scalar_t, 6>, 6> a{};
    for (int r = 0; r < 3; ++r)
      for (int c = 0; c < 3; ++c) {
        a[r][c] = a_ll[c][r];         // element(row=r, col=c) = m[c][r]
        a[r][3 + c] = a_lt[c][r];     // A(lin, ang) = a_lt
        a[3 + r][c] = a_lt[c][r];     // A(ang, lin) = a_ltᵀ 的 element(r, c)
        a[3 + r][3 + c] = a_tt[c][r];
      }

    const Vec3 rlin = -grad_lin[b] - r.lin * m_dt2;
    const Vec3 rang = -grad_ang[b] - inertial.world_inv_inertia(b) * (r.ang * inv_dt2);
    const std::array<Scalar_t, 6> rhs{rlin.x, rlin.y, rlin.z, rang.x, rang.y, rang.z};
    std::array<Scalar_t, 6> x{};
    if (!ldlt6_solve(a, rhs, x)) return false;
    d_lin = Vec3{x[0], x[1], x[2]};
    d_ang = Vec3{x[3], x[4], x[5]};

    // sanity 策略：非有限或超出物理量级就冻结这个 DOF（物理上 10 m 已远超一个步长）
    constexpr Scalar_t kMaxCorrection = Scalar_t(10);
    const Scalar_t worst = std::max({std::abs(d_lin.x), std::abs(d_lin.y), std::abs(d_lin.z),
                                     std::abs(d_ang.x), std::abs(d_ang.y), std::abs(d_ang.z)});
    if (!all_finite(d_lin) || !all_finite(d_ang) || worst > kMaxCorrection) {
      d_lin = Vec3{0};
      d_ang = Vec3{0};
      return false;
    }
    return true;
  }
};

// ==================== 8. ConstraintScratch：scatter/gather 的中间缓冲 ====================
// 存在的**唯一**理由：让"多条约束累加到同一个 body"不再需要 read-modify-write。
//
//   scatter：每条约束把自己的 (C, g_lin, g_ang) 写进**自己的** slot（无冲突，可并行）；
//   gather ：逐个 body 按**约束索引升序**把 slot 折进 SolverBuffers。
//
// 两个关键性质（别改坏）：
//   1) entries[begin[c] .. begin[c+1]) 只属于约束 c ⇒ scatter 无数据冲突，不需要原子/锁。
//   2) body_entries 按**约束索引升序**排列 ⇒ 对某个 body 而言，gather 的加总顺序与
//      "串行逐约束累加"逐条相同 ⇒ 浮点结果**位精确**，与线程数、分区、调度全都无关。
//      顺序来自这张邻接表，而不是来自线程调度 —— 这是位一致的根源。
//
// 【缓存边界】（step 5）
//   entry_begin / body_begin / body_entries 只依赖"哪条约束碰哪些 body"，与位置、
//   rest_length、depth 全部无关 ⇒ 可作为**拓扑缓存**跨 step 复用。
//   entries[] 则是每迭代被 scatter 全量重写的派生量 ⇒ **不可**缓存。
//   两者放在同一个结构体里，但生命周期不同，改动时别搞混。
struct ConstraintEntry {
  Vec3 g_lin{0};
  Vec3 g_ang{0};
  Scalar_t C = 0;
};

// 拓扑签名。version 由**约束列表的拥有者**维护：增删/替换约束时必须自增。
// 只改约束参数（rest_length / depth / compliance / λ）**不算**拓扑变更，不用 bump。
// version == 0 是安全默认值：永远重建（“我不知道拓扑变没变”）。
struct TopologyKey {
  std::uint64_t version = 0;
  std::size_t constraint_count = 0;
  std::size_t body_count = 0;

  bool same_as(const TopologyKey &o) const noexcept {
    return version != 0 && version == o.version && constraint_count == o.constraint_count &&
           body_count == o.body_count;
  }
};



struct DualUnit {
  std::uint32_t first;   // ConstraintList 索引
  std::uint32_t count;   // 3 * contact_count（manifold）或 1（standalone）
  std::uint32_t manifold; // UINT32_MAX 表示 standalone
};


struct ConstraintScratch {
  TopologyKey key{};             // 当前缓存的拓扑签名
  std::uint64_t build_count = 0;  // 观测用：布局重建了几次（测试断言靠它）
  std::uint64_t stale_count = 0;  // 观测用：缓存校验发现“拓扑变了却没 bump”的次数

  std::vector<std::uint32_t> begin;        // n_constraint + 1
  std::vector<ConstraintEntry> entries;    // Σ bodies().size()
  std::vector<std::uint32_t> body_begin;   // n_body + 1
  std::vector<std::uint32_t> body_entries; // 每个 body 的 entry 索引，按约束索引升序
  std::vector<std::uint32_t> cursor;       // 建邻接用的临时写指针（复用，免每步分配）

  void clear() {
    key = TopologyKey{};
    build_count = 0;
    stale_count = 0;
    begin.clear();
    entries.clear();
    body_begin.clear();
    body_entries.clear();
    cursor.clear();
  }
  void resize_bodies(std::size_t n) {
    // 【缓存失效的唯一入口】
    // body 数没变 ⇒ 邻接仍然有效，**什么都不许动**（尤其是不能清 body_entries，
    // 那会让 key 仍然失配不了、而邻接已经是空的 —— gather 会越界读）。
    // body 数变了 ⇒ 尺寸类数据全作废，并主动把 key 打成无效强制重建。
    if (body_begin.size() == n + 1) return;
    body_begin.assign(n + 1, 0);
    body_entries.clear();
    begin.clear();
    entries.clear();
    cursor.clear();
    key = TopologyKey{};
  }
  std::size_t constraint_count() const noexcept {
    return begin.empty() ? 0 : begin.size() - 1;
  }
};

// ==================== 8b. 分区粒度与伪共享 ====================
// 为什么需要这组常数：SoA 布局下同一"字段"的相邻 body 是紧挨着的，所以**分区边界**
// 如果落在 cacheline 中间，两个线程就会互相写对方正在读的那条 64 B 行（false sharing），
// 表现为"并行反而更慢"，而且 profiler 上看不出来（没有共享变量、没有锁）。
//
// primal / gather 阶段每个 body 最热的数组是 SolverBuffers：3 个 Mat3 + 2 个 Vec3。
// 让 chunk 的起点落在 cacheline 边界上，需要 grain · kBodyBufferBytes ≡ 0 (mod 64)。
// glm 的 Mat3 是 36 B（列主序 3×float3，无填充），Vec3 是 12 B ⇒ 132 B/body；
// 132 = 4·33，33 是奇数 ⇒ grain 必须是 16 的倍数。16 个 body = 2112 B = 33 条整行。
inline constexpr std::size_t kCacheLineBytes = 64;
inline constexpr std::size_t kBodyBufferBytes = 3 * sizeof(Mat3) + 2 * sizeof(Vec3);
inline constexpr std::size_t kEntryBytes = sizeof(ConstraintEntry);
inline constexpr std::size_t kDefaultGrainBody = 16;

static_assert((kDefaultGrainBody * kBodyBufferBytes) % kCacheLineBytes == 0,
              "grain_body must keep chunk starts cacheline-aligned, otherwise the "
              "primal/gather phases self-inflict false sharing");

// ==================== 9. State：组装 + 一步生命周期 ====================
struct State {
  PoseState pose;
  VelocityState velocity;
  InertialState inertial;
  ExternalState external;
  BodyCache cache;
  SolverBuffers buffers;
  ConstraintScratch scratch;

  Scalar_t dt = Scalar_t(1) / 60;

  void resize(std::size_t n) {
    pose.resize(n);
    velocity.resize(n);
    inertial.resize(n);
    external.resize(n);
    cache.resize(n);
    buffers.resize(n);
    scratch.resize_bodies(n);
  }
  std::size_t size() const noexcept { return pose.size(); }
  void clear() {
    pose.clear();
    velocity.clear();
    inertial.clear();
    external.clear();
    cache.clear();
    buffers.clear();
    scratch.clear();
  }

  // —— 外力入口 ——
  // 全部走 external 累加器；重力与用户力在求解器里没有任何区别。
  // 顺序硬约束：reset_external() → （可选）add_force/add_torque → begin_step()。
  void reset_external() { external.reset(inertial); } // 复位到 uniform_accel（默认即重力）
  void add_force(std::size_t b, Vec3 f) {
    if (inertial.is_static(b)) return; // 通用规则：无限质量不做功，不是重力特判
    external.add_linear(b, f * inertial.inv_mass(b));
  }
  void add_torque(std::size_t b, Vec3 t) {
    if (inertial.is_static(b)) return;
    // 力矩 → 角加速度需要世界惯量；这里懒刷新，避免引出"必须先 refresh_all"的顺序坑。
    inertial.refresh_world_inertia(b, pose.rotation_matrix(b));
    external.add_angular(b, inertial.world_inv_inertia(b) * t);
  }

  // 惯性目标残差。唯一消费者是 solve_body 的 -M/dt²·r 项。
  BodyResiduals inertial_residual(std::size_t b) const noexcept {
    return BodyResiduals{pose.position(b) - cache.inertial_position(b),
                         rotvec_between(cache.inertial_rotation(b), pose.rotation(b))};
  }

  // step 开始：刷新惯量 → 把外力（含重力）烧进速度 → 预测惯性目标 → 清局部系统
  void begin_step() {
    inertial.refresh_all(pose);
    velocity.apply_acceleration(external, inertial, dt);
    cache.begin_step(pose, velocity, inertial, dt);
    buffers.reset();
  }

  // step 结束：先补齐派生量（dq/rotMat 必须与 pose 一致），再由位移反推速度。
  // 注意 from_pose_delta 是从 **initial**（步首位姿）算的，所以重力带来的位移会
  // 正确地留在速度里；这也是为什么不能改用 inertial 作基准。
  void finish_step() {
    for (std::size_t b = 0; b < size(); ++b) finish_body(b);
  }

  // 逐 body 的步末收尾。两个写动作都是 per-body 的，所以可以被外面融合进
  // 任意一次 fork（见 frame.hpp：finish 与只读消费者共用一次并行区间）。
  void finish_body(std::size_t b) {
    cache.refresh_derived(b, pose, dt);
    velocity.from_pose_delta(b, pose, cache, dt);
  }

  // 把 pose 直接搬到惯性目标。XPBD/PBD 的显式位置投影需要它；
  // AVBD 走隐式残差（solve_body 的 -M/dt²·r），**不调**这个。
  void move_pose_to_inertial_target() {
    for (std::size_t b = 0; b < size(); ++b) move_body_to_inertial_target(b);
  }
  // 逐 body 版（并行安全：只碰自己的 pos/rot/cache）。
  void move_body_to_inertial_target(std::size_t b) {
    pose.set_position(b, cache.inertial_position(b));
    pose.set_rotation(b, cache.inertial_rotation(b));
    cache.refresh_derived(b, pose, dt);
  }
};
