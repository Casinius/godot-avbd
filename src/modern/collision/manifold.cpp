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
        if (prev.stick) {
          // 静摩擦把锚点钉在上一帧的材料点上，避免接触点滑动导致“记忆”被重置。
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
