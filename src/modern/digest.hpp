#pragma once

// ============================================================================
//  digest.hpp —— 位级状态摘要（回归基线）
// ----------------------------------------------------------------------------
//  为什么存在：并发重构唯一的验收标准是"串行执行器下与重构前逐位相同"。
//  肉眼比对 pose 是不够的（dg/λ/prevVelocity/buffers 全都能悄悄漂），所以这里
//  把 State + ConstraintList 的每一个 float 按**固定顺序**喂进 FNV-1a。
//
//  规则（改本文件前先读）：
//   1) 只对 float 的字节做 hash，不对结构体整体 memcpy —— Mat3/glm 类型有对齐
//      填充，memcpy 会把未初始化字节混进来，摘要就不可复现了。
//   2) 遍历顺序即 ABI。任何"顺手改成按 body 分组"的优化都会换掉摘要，那不是
//      重构就是行为变更，必须重新录基线。
//   3) 摘要不含指针、不含容量、不含 size（size 由 scene 固定）。
// ============================================================================

#include <cstddef>
#include <cstdint>
#include <string>

#include "phys_core/base_concept.hpp"
#include "phys_core/constraint.hpp"

namespace digest {

inline constexpr std::uint64_t kFnvOffset = 1469598103934665603ull;
inline constexpr std::uint64_t kFnvPrime  = 1099511628211ull;

inline std::uint64_t mix_bytes(std::uint64_t h, const void *p, std::size_t n) noexcept {
  const auto *b = static_cast<const unsigned char *>(p);
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<std::uint64_t>(b[i]);
    h *= kFnvPrime;
  }
  return h;
}

inline std::uint64_t mix_scalar(std::uint64_t h, Scalar_t v) noexcept {
  return mix_bytes(h, &v, sizeof(v));
}

// glm::vec3 = 3 float，无填充；逐分量喂。
inline std::uint64_t mix_vec3(std::uint64_t h, Vec3 v) noexcept {
  h = mix_scalar(h, v.x);
  h = mix_scalar(h, v.y);
  return mix_scalar(h, v.z);
}

// glm::qua 的成员布局随 GLM_FORCE_QUAT_DATA_* 变，所以只走 .x/.y/.z/.w 访问器，
// 顺序固定为 (x, y, z, w)，与内存布局解耦。
inline std::uint64_t mix_quat(std::uint64_t h, Quat q) noexcept {
  h = mix_scalar(h, q.x);
  h = mix_scalar(h, q.y);
  h = mix_scalar(h, q.z);
  return mix_scalar(h, q.w);
}

// glm 矩阵列主序：m[col][row]，按 (col, row) 字典序喂，与内存顺序一致。
inline std::uint64_t mix_mat3(std::uint64_t h, const Mat3 &m) noexcept {
  for (int c = 0; c < 3; ++c)
    for (int r = 0; r < 3; ++r)
      h = mix_scalar(h, m[c][r]);
  return h;
}

inline std::string to_hex(std::uint64_t v) {
  static const char *kDigits = "0123456789abcdef";
  std::string        out(16, '0');
  for (int i = 15; i >= 0; --i) {
    out[static_cast<std::size_t>(i)] = kDigits[v & 0xFull];
    v >>= 4;
  }
  return out;
}

// ==================== 分部位摘要（用于定位漂移来自哪一段） ====================

inline std::uint64_t pose(const PoseState &s) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < s.size(); ++b) {
    h = mix_vec3(h, s.pos[b]);
    h = mix_quat(h, s.rot[b]);
  }
  return h;
}

inline std::uint64_t velocity(const VelocityState &s) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < s.size(); ++b) {
    h = mix_vec3(h, s.v[b]);
    h = mix_vec3(h, s.omega[b]);
  }
  return h;
}

inline std::uint64_t inertial(const InertialState &s) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < s.size(); ++b) {
    h = mix_scalar(h, s.mass_inv[b]);
    h = mix_mat3(h, s.inertia_local_inv[b]);
    h = mix_mat3(h, s.inertia_world_inv[b]);
  }
  return h;
}

// 外力累加器也必须进摘要：它决定了速度更新，而"每步该复位却没复位"是并发化
// 之后很难用肉眼看出来的类错——处于累加状态的 accel 会默默漂移。
inline std::uint64_t external(const ExternalState &s) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < s.size(); ++b) {
    h = mix_vec3(h, s.accel_lin[b]);
    h = mix_vec3(h, s.accel_ang[b]);
  }
  return mix_vec3(h, s.uniform_accel);
}

inline std::uint64_t cache(const BodyCache &c) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < c.size(); ++b) {
    h = mix_vec3(h, c.initialLin[b]);
    h = mix_quat(h, c.initialAng[b]);
    h = mix_mat3(h, c.initialRotMat[b]);
    h = mix_vec3(h, c.inertialLin[b]);
    h = mix_quat(h, c.inertialAng[b]);
    h = mix_mat3(h, c.rotMat[b]);
    h = mix_vec3(h, c.dqLin[b]);
    h = mix_vec3(h, c.dqAng[b]);
    h = mix_vec3(h, c.prevVelocityLin[b]);
    h = mix_vec3(h, c.prevOmega[b]);
  }
  return h;
}

// 迭代末尾 buffers 应当已被 batch_update 清空；把它纳入摘要，就是为了让
// "reset 漏了一行" 这种 bug 立刻让基线炸掉。
inline std::uint64_t buffers(const SolverBuffers &s) noexcept {
  std::uint64_t h = kFnvOffset;
  for (std::size_t b = 0; b < s.size(); ++b) {
    h = mix_mat3(h, s.hessian_lin[b]);
    h = mix_mat3(h, s.hessian_lin_ang[b]);
    h = mix_mat3(h, s.hessian_ang[b]);
    h = mix_vec3(h, s.grad_lin[b]);
    h = mix_vec3(h, s.grad_ang[b]);
  }
  return h;
}

// dual 变量（λ）不在 State 里，在约束对象里；它是 AVBD/XPBD 的隐式状态，
// 漏掉它就等于漏掉一半的求解器状态。
inline std::uint64_t constraints(const ConstraintList &cs) noexcept {
  std::uint64_t h = kFnvOffset;
  for (const auto &v : cs) {
    std::visit([&h](const auto &c) {
      h = mix_scalar(h, c.lambda);
      h = mix_scalar(h, c.compliance);
      for (std::size_t b : c.bodies()) h = mix_scalar(h, Scalar_t(b));
      if constexpr (std::is_same_v<std::decay_t<decltype(c)>, ContactConstraint>) {
        // 规则：约束的每一个标量都要进哈希。曾经漏过 target（depth）与 manifold_id，
        // 结果是“把接触深度改错”这类输入错误对基线完全不可见——而 target 恰恰是
        // 接触里唯一直接决定物理量级的标量。新增字段时同步这里。
        h = mix_scalar(h, c.target);
        h = mix_scalar(h, Scalar_t(c.manifold_id));
        h = mix_scalar(h, c.penalty);
        h = mix_scalar(h, c.initial_c);
        h = mix_scalar(h, c.friction);
        h = mix_scalar(h, c.beta_scale);
        h = mix_scalar(h, Scalar_t(c.stick));
        h = mix_scalar(h, Scalar_t(c.contact_index));
        h = mix_scalar(h, Scalar_t(static_cast<int>(c.kind)));
        // basis 和 rA/rB 也哈希，因为它们是约束状态的一部分
        for (int j = 0; j < 3; ++j) h = mix_vec3(h, c.basis[j]);
        h = mix_vec3(h, c.rA);
        h = mix_vec3(h, c.rB);
      } else if constexpr (std::is_same_v<std::decay_t<decltype(c)>, JointAxisConstraint>) {
        // 关节行的每一个标量都要进摘要：模式/限位/暖启动状态任一改错，基线才看得见。
        h = mix_vec3(h, c.axis);
        h = mix_scalar(h, Scalar_t(static_cast<int>(c.mode)));
        h = mix_scalar(h, Scalar_t(c.angular ? 1 : 0));
        h = mix_scalar(h, c.target);
        h = mix_scalar(h, c.lower);
        h = mix_scalar(h, c.upper);
        h = mix_scalar(h, c.penalty);
        h = mix_scalar(h, c.initial);
        h = mix_scalar(h, Scalar_t(c.axis_index));
        for (JointAxisMode m : c.angular_modes)
          h = mix_scalar(h, Scalar_t(static_cast<int>(m)));
        h = mix_vec3(h, c.rA);
        h = mix_vec3(h, c.rB);
        h = mix_quat(h, c.frame_a);
        h = mix_quat(h, c.frame_b);
        h = mix_quat(h, c.rest);
      }
    }, v);
  }
  return h;
}

struct Parts {
  std::uint64_t pose, velocity, inertial, external, cache, buffers, constraints, all;
};

inline Parts parts(const State &s, const ConstraintList &cs) noexcept {
  Parts p{};
  p.pose        = digest::pose(s.pose);
  p.velocity    = digest::velocity(s.velocity);
  p.inertial    = digest::inertial(s.inertial);
  p.external    = digest::external(s.external);
  p.cache       = digest::cache(s.cache);
  p.buffers     = digest::buffers(s.buffers);
  p.constraints = digest::constraints(cs);

  // 组合时必须按 8 字节喂，不能先收敛成 float —— 那会把 64 位摘要压成 24 位尾数，
  // 不同部位恰好尾数相同就互相抵消，回归测试会静默失效。
  std::uint64_t h = kFnvOffset;
  h               = mix_bytes(h, &p.pose, sizeof(p.pose));
  h               = mix_bytes(h, &p.velocity, sizeof(p.velocity));
  h               = mix_bytes(h, &p.inertial, sizeof(p.inertial));
  h               = mix_bytes(h, &p.external, sizeof(p.external));
  h               = mix_bytes(h, &p.cache, sizeof(p.cache));
  h               = mix_bytes(h, &p.buffers, sizeof(p.buffers));
  h               = mix_bytes(h, &p.constraints, sizeof(p.constraints));
  h               = mix_scalar(h, s.dt);
  p.all           = h;
  return p;
}

} // namespace digest
