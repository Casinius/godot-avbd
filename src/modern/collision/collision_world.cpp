#include "collision/collision_world.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

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
      Aabb box = shape_aabb(views_.back());
      // 预接触生成距离：宽相必须比接触距离更宽，否则窄相永远看不到间隙对。
      // 对"本步位移远超预接触距离"的快速 body，把 proxy 沿运动方向扫一遍：
      // 这不会生成假的接触（窄相仍按当前位姿的 spec 距离判定），只是保证
      // 快速 body 附近的 shape pair 不会被宽相提前剪掉。speculative_distance==0
      // 时这里逐位保持旧行为。
      const Scalar_t spec = emit_options_.speculative_distance;
      if (spec > Scalar_t(0)) {
        if (!src.is_static(b)) {
          const Vec3     motion = state.velocity.velocity(b) * state.dt;
          const Scalar_t travel = glm::length(motion);
          if (travel > Scalar_t(2) * spec) {
            box.expand(box.lo + motion);
            box.expand(box.hi + motion);
          }
        }
        box = box.expanded(spec);
      }
      aabbs_.push_back(box);
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
                                               std::span<ContactPoint>{patch.points}, patch.basis,
                                               emit_options_.speculative_distance);
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

namespace {

// 修改速度后必须把惯性预测位姿重算一遍：solve_body / XPBD 的 begin_iterations
// 只读 cache.inertial_*，它是在 State::begin_step() 里由旧速度算好的。CCD/预接触
// 夹紧发生在 begin_step 之后，不刷新就等于没夹。
//
// 基准必须是**步首位姿** cache.initial_*，不是当前 pose：B-PXPBD 在窄相之前就把
// pose 搬到了 x~，若仍从 pose 预测，x~ + dt*v_new 会把旧速度的位移重复积一次
// （预接触夹紧因此失效，finish_step 还会把夹紧掉的速度又算回来）。用初值预测
// 得到 x_initial + dt*v_new，正是速度被改写后本步应该到达的位置。AVBD 在调用点
// 的 pose 就是步首（其 begin_iterations 不搬 pose），所以这一步对 AVBD 逐位不变。
void refresh_inertial_target(State &s, std::size_t b, Scalar_t dt) noexcept {
  if (s.inertial.is_static(b))
    return;
  Vec3 lin;
  Quat ang;
  s.inertial.predict_pose(b, s.cache.initial_position(b), s.cache.initial_rotation(b),
                          s.velocity.velocity(b), s.velocity.angular_velocity(b), dt, lin, ang);
  s.cache.set_inertial(b, lin, ang);
}

} // namespace

void CollisionWorld::clamp_fast_bodies(State &s, Scalar_t dt) {
  const Scalar_t spec = emit_options_.speculative_distance;
  if (!(spec > Scalar_t(0)) || !(dt > Scalar_t(0)))
    return;

  const std::size_t nb = std::min(shapes_.body_count(), s.size());
  constexpr std::size_t kNone = std::numeric_limits<std::size_t>::max();

  auto body_active = [this](std::size_t b) noexcept {
    return body_active_ == nullptr || body_active_(body_active_user_, b);
  };

  auto build_views = [this, &s](std::size_t b, std::vector<ShapeView> &out, Scalar_t &r_min_half) {
    out.clear();
    r_min_half = std::numeric_limits<Scalar_t>::infinity();
    const std::span<const ShapeInstance> instances = shapes_.shapes(b);
    out.reserve(instances.size());
    for (const ShapeInstance &inst : instances) {
      const ConvexMesh *mesh = nullptr;
      if (inst.shape.type == ShapeType::Convex) {
        mesh = meshes_.get(inst.shape.mesh);
        if (mesh == nullptr)
          continue;
      }
      out.push_back(make_shape_view(inst.shape, s.pose.position(b), s.pose.rotation(b),
                                    inst.local_position, inst.local_rotation, mesh));
      const Vec3 half = out.back().local_half;
      r_min_half = std::min({r_min_half, std::abs(half.x), std::abs(half.y), std::abs(half.z)});
    }
  };

  std::vector<ShapeView> moving_views;
  std::vector<ShapeView> other_views;

  // 全局多遍：一次冲量可能把"对方"（尤其是编号更小、已经处理过的 body）推向
  // 另一个障碍。若只按 body 升序扫一遍，被推下去的下层 body 本步不会再检查
  // 地面，就会穿过去。每遍都重新按当前速度扫掠，通常 1~2 遍收敛，上限 4 遍。
  for (int pass = 0; pass < 4; ++pass) {
    bool clamped_any = false;
    for (std::size_t b = 0; b < nb; ++b) {
      if (s.inertial.is_static(b) || !body_active(b))
        continue;

      Scalar_t min_half = std::numeric_limits<Scalar_t>::infinity();
      build_views(b, moving_views, min_half);
      if (moving_views.empty() || !std::isfinite(min_half) || !(min_half > Scalar_t(0)))
        continue;

      const Vec3     vb     = s.velocity.velocity(b);
      const Vec3     motion = vb * dt;
      const Scalar_t travel = glm::length(motion);
      const Scalar_t threshold = std::max(Scalar_t(2) * spec, Scalar_t(0.5) * min_half);
      if (!(travel > threshold))
        continue;

      Scalar_t    best_t = Scalar_t(1);
      std::size_t best_c = kNone;
      for (std::size_t c = 0; c < nb; ++c) {
        if (c == b)
          continue;
        if (shapes_.shapes(c).empty())
          continue;
        if (pair_filter_ != nullptr && !pair_filter_(pair_filter_user_, b, c))
          continue;

        const bool c_movable = !s.inertial.is_static(c) && body_active(c);
        const Vec3 vc        = c_movable ? s.velocity.velocity(c) : Vec3{0};
        const Vec3 rel_motion = (vb - vc) * dt;
        if (!(glm::length(rel_motion) > Scalar_t(1e-9)))
          continue;

        Scalar_t ignored_half = 0;
        build_views(c, other_views, ignored_half);
        for (const ShapeView &mv : moving_views) {
          for (const ShapeView &ov : other_views) {
            // 形状级分离轴 TOI。这里没有先做 AABB 扫掠预筛：AABB 只能证明
            // "某轴不同时重叠"，而 SAT 的任一正间隙轴必须先闭合是更保守的
            // 必要条件；预筛掉这些近失（near-miss）会让高速体在堆叠边缘少掉
            // 一层保守减速。CCD 本来就是兜底，宁可早停，不允许穿过去。
            Scalar_t t = 0;
            if (sweep_shapes_toi(mv, ov, rel_motion, t) && t < best_t) {
              best_t = t;
              best_c = c;
            }
          }
        }
      }

      if (best_c == kNone || !(best_t < Scalar_t(1)))
        continue;

      // 留一点 skin（步长的 0.2%），不要正好停在接触面上。否则下一步的
      // 起始间隙为 0：exact 窄相用严格 `< radius`，speculative 又要求 gap>eps，
      // 两个路径都不生成接触；body 会带着被 CCD 允许的残余速度直接穿过去。
      // 留 skin 后，下一步由 speculative/普通接触接管。
      const Scalar_t safe_t = std::max(Scalar_t(0), best_t - Scalar_t(0.002));

      const bool c_movable = !s.inertial.is_static(best_c) && body_active(best_c);
      const Vec3 vc        = c_movable ? s.velocity.velocity(best_c) : Vec3{0};
      const Vec3 vrel      = vb - vc;
      const Vec3 dv        = vrel * (safe_t - Scalar_t(1)); // 相对速度减少量
      const Scalar_t inv_a = s.inertial.inv_mass(b);
      const Scalar_t inv_c = c_movable ? s.inertial.inv_mass(best_c) : Scalar_t(0);
      const Scalar_t inv_sum = inv_a + inv_c;
      if (!(inv_sum > Scalar_t(0)))
        continue;

      // 动量守恒的线性冲量：J 施加在 b，-J 施加在 c。旋转被忽略（见头文件）。
      const Vec3 J = dv * (Scalar_t(1) / inv_sum);
      s.velocity.add_velocity(b, J * inv_a);
      if (c_movable)
        s.velocity.add_velocity(best_c, -J * inv_c);
      refresh_inertial_target(s, b, dt);
      if (c_movable)
        refresh_inertial_target(s, best_c, dt);
      clamped_any = true;
    }
    if (!clamped_any)
      break;
  }
}

std::size_t CollisionWorld::apply_speculative_velocity_clamp(State &s, Scalar_t dt) {
  const Scalar_t spec = emit_options_.speculative_distance;
  if (!(spec > Scalar_t(0)) || !(dt > Scalar_t(0)))
    return 0;
  const Scalar_t margin = std::max(Scalar_t(0), emit_options_.margin);

  std::size_t clamped = 0;
  for (const Manifold &m : store_.manifolds()) {
    if (s.inertial.is_static(m.body_a) && s.inertial.is_static(m.body_b))
      continue;
    const Vec3 n = m.basis.normal;
    if (!(glm::dot(n, n) > Scalar_t(1e-12)))
      continue;

    for (std::size_t i = 0; i < m.count; ++i) {
      const ManifoldContact &mc = m.contacts[i];
      if (!(mc.depth < Scalar_t(0)))
        continue; // 只有预接触（间隙>0）需要速度夹紧

      const Vec3 rAw = rotate(s.pose.rotation(m.body_a), mc.rA);
      const Vec3 rBw = rotate(s.pose.rotation(m.body_b), mc.rB);
      const Vec3 xA  = s.pose.position(m.body_a) + rAw;
      const Vec3 xB  = s.pose.position(m.body_b) + rBw;
      const Scalar_t gap = glm::dot(n, xA - xB);
      if (!(gap > Scalar_t(0)) || gap > spec)
        continue;

      // 接触点相对法向速度：包含两侧角速度在锚点上的贡献，这样"旋转正在闭合"
      // 的预接触也能被识别。有效质量因此也含世界逆惯量项（与接触行同口径）。
      const Vec3 vA = s.velocity.velocity(m.body_a) +
                      glm::cross(s.velocity.angular_velocity(m.body_a), rAw);
      const Vec3 vB = s.velocity.velocity(m.body_b) +
                      glm::cross(s.velocity.angular_velocity(m.body_b), rBw);
      const Scalar_t vn = glm::dot(vA - vB, n);
      const Scalar_t target_vn = -(gap + margin) / dt;
      // 只在闭合过快时减速度；vn >= target 时不动（绝不加能量）。
      if (!(vn < target_vn))
        continue;

      const Vec3 gA = glm::cross(rAw, n);
      const Vec3 gB = -glm::cross(rBw, n);
      Scalar_t   inv_eff = 0;
      const Scalar_t inv_a = s.inertial.inv_mass(m.body_a);
      const Scalar_t inv_b = s.inertial.inv_mass(m.body_b);
      if (inv_a > Scalar_t(0))
        inv_eff += inv_a + glm::dot(gA, s.inertial.world_inv_inertia(m.body_a) * gA);
      if (inv_b > Scalar_t(0))
        inv_eff += inv_b + glm::dot(gB, s.inertial.world_inv_inertia(m.body_b) * gB);
      if (!(inv_eff > Scalar_t(0)))
        continue;

      const Scalar_t J = (target_vn - vn) / inv_eff;
      if (!(J > Scalar_t(0)))
        continue;

      // 预接触只施加**中心冲量**：用完整的含惯量有效质量算出冲量大小，但不在
      // 这个尚未真正接触的投影锚点上加力矩。原因是 speculative 锚点来自分离
      // 状态下的最近特征，一个 manifold 又可能有多个角点；对每个角点顺序施加
      // 力矩会把静止堆叠"拧"出横向速度（420 box 金字塔实测会滑出平板）。
      // 旋转响应留给真正重叠后的普通接触行——那时锚点已经是真实接触点，并且
      // 走求解器的联合 dual/primal，而不是这里的一次性预接触。
      if (inv_a > Scalar_t(0)) {
        s.velocity.add_velocity(m.body_a, n * (J * inv_a));
        refresh_inertial_target(s, m.body_a, dt);
      }
      if (inv_b > Scalar_t(0)) {
        s.velocity.add_velocity(m.body_b, n * (-J * inv_b));
        refresh_inertial_target(s, m.body_b, dt);
      }
      ++clamped;
    }
  }
  return clamped;
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
