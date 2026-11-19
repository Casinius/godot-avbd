#include "collision/collision_world.hpp"

#include <algorithm>
#include <cmath>

namespace collision {

namespace {
constexpr std::uint64_t kFnvOffset = 1469598103934665603ull;
constexpr std::uint64_t kFnvPrime  = 1099511628211ull;

inline std::uint64_t mix_u64(std::uint64_t h, std::uint64_t v) noexcept {
  for (int i = 0; i < 8; ++i) {
    h ^= (v >> (i * 8)) & 0xffull;
    h *= kFnvPrime;
  }
  return h;
}

// State + ShapeTable 的只读适配器。碰撞层不直接依赖 State 的内部布局：只要满足
// BodyShapeSource，就能被 step 消费。
struct StateShapeSource {
  const State      &state;
  const ShapeTable &shapes;

  [[nodiscard]] std::size_t size() const noexcept { return state.size(); }
  [[nodiscard]] std::size_t shape_count(std::size_t b) const noexcept {
    return shapes.shapes(b).size();
  }
  [[nodiscard]] ShapeInstance shape_instance(std::size_t b, std::size_t slot) const {
    return shapes.shapes(b)[slot];
  }
  [[nodiscard]] Vec3 position(std::size_t b) const noexcept { return state.pose.position(b); }
  [[nodiscard]] Quat rotation(std::size_t b) const noexcept { return state.pose.rotation(b); }
  [[nodiscard]] bool is_static(std::size_t b) const noexcept { return state.inertial.is_static(b); }
};

static_assert(BodyShapeSource<StateShapeSource>);
static_assert(BroadphaseModel<Broadphase>);
static_assert(NarrowphaseModel<DiscreteNarrowphase>);
static_assert(ContactSink<ConstraintListSink>);
} // namespace

void CollisionWorld::set_pair_filter(PairFilter filter, void *user) {
  // 即使 (filter, user) 与上次相同也作废：user 指向的规则表可能已经变了，而这里
  // 无从知道。调用方应在规则真正变化时才调用（每帧调用会清掉 warm start）。
  pair_filter_      = filter;
  pair_filter_user_ = user;
  invalidate_contacts();
}

void CollisionWorld::step(const State &state) {
  if (shapes_.body_count() != state.size())
    resize(state.size());

  stats_        = {};
  stats_.bodies = state.size();

  // —— 1) body pose + instance local transform → per-proxy ShapeView / AABB ——
  // proxy 顺序 = body-major + slot 顺序，因此 proxy id 的字典序就是
  // (body, slot) 的字典序。broadphase 的 pair 排序不需要再做二次排序。
  views_.clear();
  aabbs_.clear();
  proxy_body_.clear();
  proxy_slot_.clear();
  views_.reserve(shapes_.instance_count());
  aabbs_.reserve(shapes_.instance_count());
  proxy_body_.reserve(shapes_.instance_count());
  proxy_slot_.reserve(shapes_.instance_count());

  const StateShapeSource src{state, shapes_};
  for (std::size_t b = 0; b < src.size(); ++b) {
    const std::size_t count = src.shape_count(b);
    for (std::size_t slot = 0; slot < count; ++slot) {
      const ShapeInstance inst = src.shape_instance(b, slot);
      const ConvexMesh   *mesh = nullptr;
      if (inst.shape.type == ShapeType::Convex) {
        mesh = meshes_.get(inst.shape.mesh);
        if (mesh == nullptr)
          continue; // 未注册/已失效的 convex shape 不参与碰撞
      }
      proxy_body_.push_back(b);
      proxy_slot_.push_back(static_cast<std::uint32_t>(slot));
      views_.push_back(make_shape_view(inst.shape, src.position(b), src.rotation(b),
                                       inst.local_position, inst.local_rotation, mesh));
      aabbs_.push_back(shape_aabb(views_.back()));
    }
  }
  stats_.proxies = views_.size();

  // —— 2) proxy 级宽相 ——
  broadphase_.update(aabbs_);
  stats_.broad_pairs = broadphase_.pairs().size();

  // —— 3) proxy pair → body pair 过滤 → narrowphase → manifold ——
  store_.begin(broadphase_.pairs());
  for (const ProxyPair &pair : broadphase_.pairs()) {
    if (pair.a >= views_.size() || pair.b >= views_.size())
      continue;

    const std::size_t body_a = proxy_body_[pair.a];
    const std::size_t body_b = proxy_body_[pair.b];

    // 同一 body 的 shape 之间不做自碰撞；两个静态 body 也不产生约束。
    if (body_a == body_b)
      continue;
    if (src.is_static(body_a) && src.is_static(body_b))
      continue;

    // 引擎侧规则（layer/mask/例外表）。放在内建规则之后：内建规则永远不能被放开。
    if (pair_filter_ != nullptr && !pair_filter_(pair_filter_user_, body_a, body_b))
      continue;

    ContactPatch patch;
    patch.body_a   = body_a;
    patch.body_b   = body_b;
    patch.proxy_a  = pair.a;
    patch.proxy_b  = pair.b;
    patch.friction = std::sqrt(std::max(Scalar_t(0), views_[pair.a].friction) *
                               std::max(Scalar_t(0), views_[pair.b].friction));
    // 释放/粘连修复判定用同一个 margin（target = -margin），与 emit 保持一致。
    patch.margin = emit_options_.margin;

    const std::size_t n = narrowphase_.collide(views_[pair.a], views_[pair.b],
                                               std::span<ContactPoint>{patch.points}, patch.basis);
    if (n == 0)
      continue;
    patch.count = n;
    ++stats_.narrow_pairs;
    stats_.contact_points += n;
    store_.match(patch);
  }
  store_.commit();
  stats_.manifolds = store_.size();
  recompute_signature();
}

void CollisionWorld::emit(ConstraintList &list) const {
  erase_contact_constraints(list);
  sink_.emit(store_.manifolds(), emit_options_, list);
}

void CollisionWorld::absorb_warm_start(const ConstraintList &list) {
  std::span<Manifold> manifolds = store_.manifolds_for_warm_start();
  for (const ConstraintVariant &v : list) {
    if (!std::holds_alternative<ContactConstraint>(v))
      continue;
    const ContactConstraint &cc = std::get<ContactConstraint>(v);
    if (cc.manifold_id >= manifolds.size())
      continue;
    Manifold &m = manifolds[cc.manifold_id];
    if (cc.contact_index >= m.count)
      continue;
    const int k = static_cast<int>(cc.kind);
    if (k < 0 || k > 2)
      continue;
    ManifoldContact &mc = m.contacts[cc.contact_index];
    mc.lambda[k]        = cc.lambda;
    mc.penalty[k]       = cc.penalty;
    mc.beta_scale       = cc.beta_scale;
    mc.stick            = cc.stick;
  }
}

bool CollisionWorld::sync(const State &state, ConstraintList &list,
                          std::uint64_t &topology_version) {
  // list 里是上一帧求解后的约束；先把它承载的 dual 状态写回 manifold，
  // 再让 step() 的 feature 匹配去消费。这就是 warm start 的闭环。
  absorb_warm_start(list);
  step(state);
  emit(list);
  const bool changed = layout_signature_ != last_signature_;
  if (changed) {
    ++topology_version;
    last_signature_ = layout_signature_;
  }
  return changed;
}

void CollisionWorld::recompute_signature() {
  std::uint64_t h = kFnvOffset;
  h               = mix_u64(h, static_cast<std::uint64_t>(store_.size()));
  for (const Manifold &m : store_.manifolds()) {
    h = mix_u64(h, static_cast<std::uint64_t>(m.body_a));
    h = mix_u64(h, static_cast<std::uint64_t>(m.body_b));
    h = mix_u64(h, static_cast<std::uint64_t>(m.proxy_a));
    h = mix_u64(h, static_cast<std::uint64_t>(m.proxy_b));
    h = mix_u64(h, static_cast<std::uint64_t>(m.count));
  }
  layout_signature_ = h;
}

} // namespace collision
