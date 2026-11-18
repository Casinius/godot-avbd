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
//   6) 切向 target = 0：第一阶段用 penalty 型静摩擦，摩擦锥在 dual 阶段按
//      manifold_id 聚合（见 plan 的第二阶段）。
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
  Scalar_t margin     = Scalar_t(0.005);
  Scalar_t compliance = Scalar_t(0);
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
inline std::size_t emit_contact_constraints(std::span<const Manifold> manifolds,
                                            const ContactEmitOptions &opt, ConstraintList &list) {
  std::size_t emitted = 0;
  for (std::size_t mi = 0; mi < manifolds.size(); ++mi) {
    const Manifold &m = manifolds[mi];
    for (std::size_t ci = 0; ci < m.count; ++ci) {
      const ManifoldContact &mc = m.contacts[ci];
      // μ == 0 时没有切向约束：不要发射两条永远不该有 λ 的切线。
      const int rows = m.friction > Scalar_t(0) ? 3 : 1;
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
        cc.lambda        = mc.lambda[k];
        cc.compliance    = opt.compliance;
        cc.penalty       = mc.penalty[k];
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
