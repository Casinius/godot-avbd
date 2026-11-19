#pragma once

// ============================================================================
//  collision/contact_constraints.hpp —— Manifold → ContactConstraint
// ----------------------------------------------------------------------------
//  【发射契约】
//   1) 每个 manifold 发 3 × count 条 ContactConstraint：normal、tangent0、tangent1；
//   2) 同一 manifold 的约束在 ConstraintList 中**连续**，manifold_id 从 0 递增；
//   3) 接触约束只追加在列表尾部；调用方每帧先 erase_contact_constraints() 再 emit；
//   4) basis[0] = normal(B->A)，basis[1..2] = 两条切线；rA/rB 是局部锚点；
//   5) normal 的 target = -margin：
//        C = n·(xA-xB) - target = n·(xA-xB) + margin
//      穿透时 n·(xA-xB) < 0，C<0 触发沿 +n 推开 A；静息允许 margin 的重叠。
//   6) 切向 target = 0：摩擦是 penalty 型；摩擦锥与法线单边夹紧在同一处按
//      (manifold_id, contact_index) 成组做（AVBDStep::post_dual）。
//   7) 每个接触的 λ/K 在发射时做一次 warm start 衰减（alpha*gamma / gamma），
//      位置与参考实现的 Manifold::initialize 相同。
//
//  【为什么不新增约束类型】
//   现有 ContactConstraint 已经有 basis/rA/rB/target/lambda/penalty/beta_scale/stick，
//   它们正好覆盖一个接触点的三个方向。新增类型会让 batch/scatter/digest 全部改一遍，
//   而第一阶段并不需要新的数据形态。
// ============================================================================

#include <algorithm>
#include <cstddef>
#include <variant>

#include "collision/manifold.hpp"
#include "phys_core/constraint.hpp"

namespace collision {

struct ContactEmitOptions {
  // 静息接触允许的重叠量（法线 target = -margin）。与生成距离无关。
  Scalar_t margin     = Scalar_t(0.005);
  Scalar_t compliance = Scalar_t(0);
  // 预接触生成距离（米）。>0 时宽相把 AABB 外扩这么多，窄相为间隙不超过它的
  // 形状对生成 depth<0 的 speculative 接触点。默认 0：核心的钉死基线/摘要
  // 逐位不变，只有引擎显式打开时才产生新行为。
  Scalar_t speculative_distance = Scalar_t(0);
};

// 删除所有 ContactConstraint，保留 DistanceConstraint 等用户约束。
// 用 remove_if + erase 的稳定写法：接触约束总是尾部追加，所以前面的索引不变。
inline void erase_contact_constraints(ConstraintList &list) {
  list.erase(std::remove_if(list.begin(), list.end(),
                            [](const ConstraintVariant &v) {
                              return std::holds_alternative<ContactConstraint>(v);
                            }),
             list.end());
}

// 追加接触约束。返回发射的约束条数。manifold_id 从 0 开始按 manifolds 顺序分配。
//
// 【warm start 的初始化（与参考实现 Manifold::initialize 同一位置）】
//   λ ← λ * alpha*gamma          （Eq.19）
//   K ← clamp(K * gamma, floor, ceil)
// 放在这里而不是 absorb_warm_start：absorb 的契约是"把 ContactConstraint 的 dual
// 状态原样写回 manifold"（collision_selftest 直接断言逐位相等），参考实现也是在
// 每个接触被**初始化**成这一帧的行时衰减，而不是在写回时。两者都只衰减一次/帧。
inline std::size_t emit_contact_constraints(std::span<const Manifold> manifolds,
                                            const ContactEmitOptions &opt, ConstraintList &list) {
  std::size_t emitted = 0;
  for (std::size_t mi = 0; mi < manifolds.size(); ++mi) {
    const Manifold &m = manifolds[mi];
    for (std::size_t ci = 0; ci < m.count; ++ci) {
      const ManifoldContact &mc = m.contacts[ci];
      // speculative 接触（depth<0，几何上还隔着间隙）不发射求解行。
      // 理由：法线行的 primal RHS 是 K*C+λ，C = gap+margin > 0 且 λ=0 时
      // K*C 是吸引力；只靠 dual 的 F_n≤0 夹紧来不及，因为 primal 已经先动了。
      // 预接触只负责在 solve 之前做速度夹紧；真正的支撑等形状重叠后由普通
      // 接触行（depth≥0）建立。这样也与"不得加能量"硬约束一致。
      if (mc.depth < Scalar_t(0))
        continue;
      // μ == 0 时没有切向约束：不要发射两条永远不该有 λ 的切线。
      const int rows = m.friction > Scalar_t(0) ? 3 : 1;
      // primal 的力界（lagged 锥）：法线 F ≤ 0；切向 |F| ≤ μ|λ_n|（用发射时的 λ_n）。
      const Scalar_t lambda_n = mc.lambda[0] * kWarmStartDecay;
      const Scalar_t cone_cap = m.friction * std::abs(lambda_n);
      for (int k = 0; k < rows; ++k) {
        ContactConstraint cc;
        cc.bodies_       = {m.body_a, m.body_b};
        cc.manifold_id   = static_cast<std::uint32_t>(mi);
        cc.contact_index = static_cast<std::uint8_t>(ci);
        cc.kind          = static_cast<ContactKind>(k);
        cc.basis[0]      = m.basis.normal;
        cc.basis[1]      = m.basis.tangent0;
        cc.basis[2]      = m.basis.tangent1;
        cc.rA            = mc.rA;
        cc.rB            = mc.rB;
        cc.target        = (k == static_cast<int>(ContactKind::Normal)) ? -opt.margin : Scalar_t(0);
        cc.lambda        = mc.lambda[k] * kWarmStartDecay;
        cc.friction      = m.friction;
        cc.compliance    = opt.compliance;
        cc.penalty       = std::clamp(mc.penalty[k] * kWarmStartGamma, kContactPenaltyFloor,
                                      kContactPenaltyCeil);
        cc.initial_c     = 0; // 步首由 AVBDStep::capture_step_start 取（步作用域状态）
        cc.force_min     = (k == 0) ? -kForceUnclamped : -cone_cap;
        cc.force_max     = (k == 0) ? kForceUnclamped : cone_cap;
        cc.beta_scale    = mc.beta_scale;
        cc.stick         = mc.stick;
        list.emplace_back(cc);
        ++emitted;
      }
    }
  }
  return emitted;
}

// ContactSink 的具体实现。保持无状态：只是把自由函数包成策略对象。
struct ConstraintListSink {
  std::size_t emit(std::span<const Manifold> manifolds, const ContactEmitOptions &opt,
                   ConstraintList &list) const {
    return emit_contact_constraints(manifolds, opt, list);
  }
};

} // namespace collision
