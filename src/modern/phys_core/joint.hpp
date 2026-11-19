#pragma once

// ============================================================================
//  joint.hpp —— 关节聚合体：一份紧凑描述 → N 条 JointAxisConstraint 行
// ----------------------------------------------------------------------------
//  现代核心的约束容器每个 alternative 只有一行、一个标量 λ（见 joint_constraint
//  的决策记录）。所以一个 6-DOF 关节在这里被拆成 0..6 行：
//
//    * 线性轴 0..2：joint frame 的主轴，锚点差异投影到轴上（米）；
//    * 角向轴 0..2：joint frame 主轴上的转角（弧度），测量自建立时的相对姿态；
//    * Free 轴不发射行；
//    * Locked 行是双边等式，Limited 行是单边不等式，其余行为与接触行同构。
//
//  Joint 只拥有行，不参与求解：调用方把 rows() 拷进 ConstraintList（或直接持有
//  Joint 并按行更新）。rebuild() 可以反复调用，用于关节参数变化后重建行。
//
//  【坐标系】
//    frame_a / frame_b 是关节坐标系在各自 body 局部系里的朝向。锚点 rA / rB 与
//    轴都按 body 局部系表达；Joint 在发射行时把 frame_a 折进轴方向，把
//    frame_a / frame_b / 建立时的相对姿态折进角向测量，所以求解器只需 body pose。
//    服务器侧一般 frame 是单位四元数（关节节点的 ref basis 无旋转），但这里不假设。
//
//  【与 legacy 的差异】
//    legacy 的 GenericJoint 忽略 ref basis、并且只认 body 局部主轴（离轴铰链会被
//    snap 到最近的轴）。这里轴与角向读数都在 joint frame 里算，所以离轴关节不需要
//    snap，也不需要镜像限位符号。
// ============================================================================

#include <array>
#include <cstddef>
#include <vector>

#include "phys_core/constraint.hpp"

// 一个自由度的描述。
struct JointAxisSpec {
  JointAxisMode mode = JointAxisMode::Free;
  Scalar_t lower     = 0; // Limited 的可行区间（线性：米；角向：弧度）
  Scalar_t upper     = 0;
  Scalar_t target    = 0; // Locked 的目标值
};

// 一个关节的紧凑描述。body_a / body_b 都是 State 下标；世界侧用一个静态 body 表示。
struct JointDesc {
  std::size_t body_a = 0;
  std::size_t body_b = 0;
  Vec3 rA{0}, rB{0}; // 各自 body 局部系里的锚点
  Quat frame_a = identity_quat();
  Quat frame_b = identity_quat();
  std::array<JointAxisSpec, 3> linear{};  // 轴 0..2 = joint frame X..Z
  std::array<JointAxisSpec, 3> angular{};
};

class Joint {
public:
  Joint() = default;
  Joint(const JointDesc &p_desc, const State &p_state) { rebuild(p_desc, p_state); }

  // 从描述重建行。rest 由当前 pose 取得：关节建立时的相对姿态就是它的零点。
  // 罚因子 / λ 初始为暖启动的起点（λ=0，K=floor），下一次 capture_step_start 再衰减。
  void rebuild(const JointDesc &p_desc, const State &p_state) {
    desc_ = p_desc;

    const Quat q_a = glm::normalize(p_state.pose.rotation(p_desc.body_a) * p_desc.frame_a);
    const Quat q_b = glm::normalize(p_state.pose.rotation(p_desc.body_b) * p_desc.frame_b);
    const Quat rest = glm::normalize(glm::conjugate(q_a) * q_b);

    std::array<JointAxisMode, 3> angular_modes{};
    for (std::size_t i = 0; i < 3; ++i)
      angular_modes[i] = p_desc.angular[i].mode;

    rows_.clear();
    rows_.reserve(6);

    for (int i = 0; i < 3; ++i) {
      const JointAxisSpec &spec = p_desc.linear[static_cast<std::size_t>(i)];
      if (spec.mode == JointAxisMode::Free)
        continue;
      JointAxisConstraint row;
      row.bodies_ = {p_desc.body_a, p_desc.body_b};
      row.rA      = p_desc.rA;
      row.rB      = p_desc.rB;
      row.axis    = rotate(p_desc.frame_a, joint_axis_detail::principal_axis(i));
      row.angular = false;
      row.mode    = spec.mode;
      row.target  = spec.target;
      row.lower   = spec.lower;
      row.upper   = spec.upper;
      row.penalty = Scalar_t(0); // capture_step_start seeds the initial stiffness
      rows_.push_back(row);
    }

    for (int i = 0; i < 3; ++i) {
      const JointAxisSpec &spec = p_desc.angular[static_cast<std::size_t>(i)];
      if (spec.mode == JointAxisMode::Free)
        continue;
      JointAxisConstraint row;
      row.bodies_        = {p_desc.body_a, p_desc.body_b};
      row.rA             = p_desc.rA;
      row.rB             = p_desc.rB;
      row.axis           = rotate(p_desc.frame_a, joint_axis_detail::principal_axis(i));
      row.angular        = true;
      row.mode           = spec.mode;
      row.target         = spec.target;
      row.lower          = spec.lower;
      row.upper          = spec.upper;
      row.frame_a        = p_desc.frame_a;
      row.frame_b        = p_desc.frame_b;
      row.rest           = rest;
      row.angular_modes  = angular_modes;
      row.axis_index     = static_cast<std::uint8_t>(i);
      row.penalty        = Scalar_t(0); // capture_step_start seeds the initial stiffness
      rows_.push_back(row);
    }
  }

  [[nodiscard]] const std::vector<JointAxisConstraint> &rows() const noexcept { return rows_; }
  [[nodiscard]] const JointDesc &desc() const noexcept { return desc_; }

private:
  JointDesc desc_{};
  std::vector<JointAxisConstraint> rows_{};
};
