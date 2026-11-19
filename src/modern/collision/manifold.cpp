#include "collision/manifold.hpp"

#include <algorithm>

namespace collision {

namespace {
[[nodiscard]] inline bool pair_less(const Manifold &m, ProxyPair p) noexcept {
  return m.proxy_a < p.a || (m.proxy_a == p.a && m.proxy_b < p.b);
}
[[nodiscard]] inline bool pair_equal(const Manifold &m, ProxyPair p) noexcept {
  return m.proxy_a == p.a && m.proxy_b == p.b;
}
[[nodiscard]] inline bool pair_less_pair(BodyPair a, BodyPair b) noexcept {
  return a.a < b.a || (a.a == b.a && a.b < b.b);
}
} // namespace

void ManifoldStore::begin(std::span<const BodyPair> pairs) {
  (void)pairs; // pair 列表只用于 API 契约；match 按 patch 的 pair 自己推进游标。
  next_.clear();
  cursor_ = 0;
}

void ManifoldStore::match(const ContactPatch &patch) {
  if (patch.count == 0)
    return;
  const ProxyPair key{patch.proxy_a, patch.proxy_b};

  // 契约：match 必须按 pair 升序调用。乱序输入直接丢弃，避免破坏 current 的排序。
  if (!next_.empty() && !pair_less_pair(next_.back().pair(), key))
    return;

  const Manifold *old = nullptr;
  while (cursor_ < current_.size() && pair_less(current_[cursor_], key))
    ++cursor_;
  if (cursor_ < current_.size() && pair_equal(current_[cursor_], key)) {
    old = &current_[cursor_];
    ++cursor_;
  }

  Manifold m;
  m.body_a   = patch.body_a;
  m.body_b   = patch.body_b;
  m.proxy_a  = patch.proxy_a;
  m.proxy_b  = patch.proxy_b;
  m.basis    = patch.basis;
  m.friction = patch.friction;
  m.count    = patch.count;

  for (std::size_t i = 0; i < patch.count; ++i) {
    const ContactPoint &src = patch.points[i];
    ManifoldContact    &dst = m.contacts[i];
    dst.feature_key         = src.feature_key;
    dst.depth               = src.depth;
    dst.rA                  = src.rA;
    dst.rB                  = src.rB;
    dst.lambda              = Vec3{0};
    dst.penalty             = Vec3{kPenaltyFloor};
    dst.beta_scale          = 1;
    dst.stick               = false;

    if (old != nullptr) {
      for (std::size_t j = 0; j < old->count; ++j) {
        if (old->contacts[j].feature_key != src.feature_key)
          continue;
        const ManifoldContact &prev = old->contacts[j];
        dst.lambda                  = prev.lambda;
        dst.penalty                 = prev.penalty;
        dst.beta_scale              = prev.beta_scale;
        dst.stick                   = prev.stick;
        // 新鲜接触的几何间隙（法线方向，含 margin，与 ContactConstraint::value 同一口径）。
        // 超过 margin + kContactReleaseGap ⇒ 真的分开了：旧 λ 不能跨过间隙继续施力。
        // 只丢 λ/stick，penalty（刚度）保留：重合时一轮就能重建支撑。
        const Scalar_t fresh_c =
                glm::dot(patch.basis.normal, src.xA - src.xB) + patch.margin;
        if (fresh_c > patch.margin + kContactReleaseGap) {
          // 锚点保持本帧的新值（默认已写入），只清掉跨帧力与静摩擦标记。
          dst.lambda = Vec3{0};
          dst.stick  = false;
          break;
        }
        if (!prev.stick) {
          // 无静摩擦：直接用本帧的新锚点（默认值）。
          break;
        }
        // 静摩擦把锚点钉在上一帧的材料点上，避免接触点滑动导致“记忆”被重置。
        // 用局部锚点算漂移即可：同一 body 的旋转对两个局部向量等长保持。
        const Scalar_t drift = std::max(glm::length(prev.rA - src.rA),
                                        glm::length(prev.rB - src.rB));
        if (drift > kStickDriftHeal) {
          // 冻结锚点已滑离真实接触面：重新采用新锚点，并清掉绑在旧锚点上的摩擦
          // 历史（法线支撑是“压得多紧”的函数，不是锚点位置的函数，保留）。
          dst.rA     = src.rA;
          dst.rB     = src.rB;
          dst.stick  = false;
          dst.lambda = Vec3{dst.lambda[0], Scalar_t(0), Scalar_t(0)};
        } else {
          dst.rA = prev.rA;
          dst.rB = prev.rB;
        }
        break;
      }
    }
  }

  next_.push_back(m);
}

void ManifoldStore::commit() {
  current_.swap(next_);
  next_.clear();
  cursor_ = 0;
}

} // namespace collision
