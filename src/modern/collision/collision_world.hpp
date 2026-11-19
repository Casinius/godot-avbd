#pragma once

// ============================================================================
//  collision/collision_world.hpp —— 碰撞门面：shape → BVH → narrowphase → manifold
// ----------------------------------------------------------------------------
//  【职责边界】
//   拥有：ShapeTable（body → ShapeInstance）、Broadphase、ManifoldStore、统计信息。
//   不拥有：State / ConstraintList —— 它们以 span 或 const& 传入。
//   不写：pose、velocity、任何求解器状态。碰撞阶段是纯只读。
//
//  【每步顺序（硬契约）】
//     step(state)   : 从 pose 派生 ShapeView → broadphase → narrowphase → manifold.commit
//     emit(list)    : 先 erase 旧 ContactConstraint，再按 manifold 顺序追加
//     sync(...)     : step + emit，并在**拓扑签名变化**时 bump topology_version
//
//  【为什么 topology_version 由这里管】
//     接触约束每帧重建；只有 pair 集合/每 manifold 点数变化才需要重建 scatter 邻接。
//     签名只包含拓扑（pair + count），不包含深度/位置，因此静息接触不会每帧失效。
//
//  【body pair 过滤】
//     自碰撞与静态-静态是内建规则；引擎侧的 layer/mask/例外表通过 set_pair_filter()
//     注入。回调是普通函数指针 + 不透明 user 指针：不分配、不做类型擦除、不依赖
//     全局状态，因此结果仍只由 (body 索引, 过滤规则) 决定。过滤器在 broadphase
//     pair 转成 body pair 之后、narrowphase 之前咨询；返回 false 的 pair 直接丢弃。
//     修改过滤器会作废已缓存的 manifold（接触集合变了，warm start 不再有效）。
// ============================================================================

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "collision/broadphase.hpp"
#include "collision/concepts.hpp"
#include "collision/contact_constraints.hpp"
#include "collision/convex_mesh.hpp"
#include "collision/manifold.hpp"
#include "collision/narrowphase.hpp"
#include "collision/shape.hpp"
#include "collision/shape_instance.hpp"
#include "collision/shape_table.hpp"

namespace collision {

// 引擎侧过滤规则：返回 true 表示 body_a 与 body_b 这一对允许产生接触。
// 参数顺序不保证 a < b；实现方必须对称处理，且只读 (user, a, b)。
using PairFilter = bool (*)(void *user, std::size_t body_a, std::size_t body_b);

// CCD 的"醒着"查询。核心没有 sleep 概念；引擎通过它告诉 CCD 哪些 body 在
// 本步参与运动/受力。返回 false 的 body 在 CCD 里按静态处理（不会被推走），
// 这样显式 sleep 的 body 不会被穿过去。nullptr = 全部醒着。
using BodyActiveFilter = bool (*)(void *user, std::size_t body);

struct CollisionWorldStats {
  std::size_t bodies         = 0;
  std::size_t proxies        = 0;
  std::size_t broad_pairs    = 0;
  std::size_t narrow_pairs   = 0;
  std::size_t contact_points = 0;
  std::size_t manifolds      = 0;
};

class CollisionWorld {
public:
  void resize(std::size_t n) {
    // body 数量变化会让旧 pair/manifold 的 body 索引失去意义，warm start 必须作废。
    const bool body_count_changed = shapes_.body_count() != n;
    if (body_count_changed) {
      store_.clear();
      last_signature_ = 0;
    }
    shapes_.resize_bodies(n);
    clear_proxy_scratch();
  }
  [[nodiscard]] std::size_t size() const noexcept { return shapes_.body_count(); }

  // —— compound API：一个 body 可以挂 0..N 个 ShapeInstance ——
  // 单 shape API 已删除：调用方必须显式给出 ShapeInstance（local transform 可留零）。
  // 任何 shape 表变更都会作废 manifold warm start，并让下一次 sync() bump version。
  void set_shapes(std::size_t b, std::span<const ShapeInstance> instances) {
    if (b >= shapes_.body_count())
      resize(b + 1);
    shapes_.set_shapes(b, instances);
    invalidate_contacts();
  }

  void add_shapes(std::size_t b, std::span<const ShapeInstance> instances) {
    if (b >= shapes_.body_count())
      resize(b + 1);
    shapes_.add_shapes(b, instances);
    invalidate_contacts();
  }

  void clear_shapes(std::size_t b) {
    shapes_.clear_shapes(b);
    invalidate_contacts();
  }

  [[nodiscard]] std::span<const ShapeInstance> shapes(std::size_t b) const noexcept {
    return shapes_.shapes(b);
  }
  [[nodiscard]] std::size_t shape_count(std::size_t b) const noexcept {
    return shapes_.shapes(b).size();
  }
  [[nodiscard]] std::uint64_t geometry_epoch() const noexcept { return geometry_epoch_; }

  // —— convex mesh 注册表 ——
  // 客户端提供预构建 hull；Shape 只存 id。add/replace 都会失效当前 manifold。
  ConvexMeshId add_convex_mesh(ConvexMesh mesh) {
    const ConvexMeshId id = meshes_.add(std::move(mesh));
    if (id != kInvalidConvexMesh)
      invalidate_contacts();
    return id;
  }
  bool replace_convex_mesh(ConvexMeshId id, ConvexMesh mesh) {
    if (!meshes_.replace(id, std::move(mesh)))
      return false;
    invalidate_contacts();
    return true;
  }
  [[nodiscard]] const ConvexMesh *convex_mesh(ConvexMeshId id) const noexcept {
    return meshes_.get(id);
  }
  [[nodiscard]] const ConvexMeshStore &convex_meshes() const noexcept { return meshes_; }
  [[nodiscard]] std::size_t convex_mesh_count() const noexcept { return meshes_.size(); }

  void set_emit_options(ContactEmitOptions opt) noexcept { emit_options_ = opt; }
  [[nodiscard]] const ContactEmitOptions &emit_options() const noexcept { return emit_options_; }
  // 预接触生成距离的便捷入口（等价于改 emit_options().speculative_distance）。
  void set_speculative_distance(Scalar_t distance) noexcept {
    emit_options_.speculative_distance = std::max(Scalar_t(0), distance);
  }
  [[nodiscard]] Scalar_t speculative_distance() const noexcept {
    return emit_options_.speculative_distance;
  }

  void set_body_active_filter(BodyActiveFilter filter, void *user = nullptr) noexcept {
    body_active_      = filter;
    body_active_user_ = user;
  }

  // 引擎侧 pair 过滤（layer/mask/例外表）。nullptr 表示"除了内建规则外全部允许"。
  // 会作废 manifold：过滤结果变了，旧接触就不能继续 warm start。
  void set_pair_filter(PairFilter filter, void *user = nullptr);
  [[nodiscard]] PairFilter pair_filter() const noexcept { return pair_filter_; }

  // 从 State 的当前 pose 重建所有 ShapeView/AABB，跑宽相 + 窄相 + manifold 合并。
  void step(const State &state);

  // 先清空旧接触约束，再追加当前 manifold。非接触约束保持不变。
  void emit(ConstraintList &list) const;

  // 求解结束后，把 ContactConstraint 里的 dual 状态（lambda/penalty/stick）写回
  // manifold。下一次 step 的 feature 匹配才会真正 warm start。
  // sync() 在 step 之前自动调用它；手动用 step()/emit() 的调用方必须自己调。
  void absorb_warm_start(const ConstraintList &list);

  // 清空 manifold store 的跨帧对偶状态（λ/stick），保留 manifold 与罚因子。
  // 求解器切换时调用：AVBD 的 λ 是力，B-PXPBD 的 λ 是乘子，不能互相暖启动。
  // 不要每帧调用——那会关掉正常 warm start。
  void reset_warm_start() noexcept { store_.reset_warm_start(); }

  // absorb_warm_start + step + emit，并只在拓扑签名变化时 bump version。返回是否 bump。
  bool sync(const State &state, ConstraintList &list, std::uint64_t &topology_version);

  // —— Layer 2：保守扫掠 CCD 兜底 ——
  // 在 sync() 之后、求解之前调用。只处理"本步线位移超过
  // max(2*speculative_distance, 0.5*min_half_extent)"且醒着的动态 body；
  // 用平移 AABB 扫掠求最早 t∈(0,1]，把相对速度缩到 t（动量按逆质量分配）。
  // 不改 pose；下一步的接触求解负责剩余部分。**忽略旋转**（保守兜底，不是精确 TOI）。
  // speculative_distance <= 0 时是空操作。
  void clamp_fast_bodies(State &state, Scalar_t dt);

  // —— Layer 1：speculative 接触的速度夹紧 ——
  // 对 manifold 里 depth<0（几何间隙 gap>0）的接触点，如果本步闭合位移
  // -v_n*dt 超过 gap+margin，就沿法线施加冲量把闭合速度降到 (gap+margin)/dt。
  // 只用与接触行相同的逆质量/世界逆惯量；切向速度不动，且只减不增。
  // speculative_distance <= 0 时是空操作。返回被夹紧的接触点数。
  std::size_t apply_speculative_velocity_clamp(State &state, Scalar_t dt);

  [[nodiscard]] std::uint64_t layout_signature() const noexcept { return layout_signature_; }
  [[nodiscard]] std::span<const Manifold>  manifolds() const noexcept { return store_.manifolds(); }
  [[nodiscard]] const Broadphase          &broadphase() const noexcept { return broadphase_; }
  [[nodiscard]] const CollisionWorldStats &stats() const noexcept { return stats_; }

  void clear() {
    shapes_.clear();
    clear_proxy_scratch();
    broadphase_ = Broadphase{};
    meshes_.clear();
    store_.clear();
    layout_signature_ = 0;
    last_signature_   = 0;
    ++geometry_epoch_;
    pair_filter_      = nullptr;
    pair_filter_user_ = nullptr;
    body_active_      = nullptr;
    body_active_user_ = nullptr;
    stats_            = {};
  }

private:
  void recompute_signature();
  void clear_proxy_scratch() noexcept {
    views_.clear();
    aabbs_.clear();
    proxy_body_.clear();
    proxy_slot_.clear();
  }
  void invalidate_contacts() noexcept {
    store_.clear();
    layout_signature_ = 0;
    last_signature_   = 0;
    ++geometry_epoch_;
  }

  ShapeTable                  shapes_{};
  ConvexMeshStore             meshes_{};
  std::vector<ShapeView>      views_{};      // per-proxy scratch，跨帧复用
  std::vector<Aabb>           aabbs_{};      // per-proxy scratch，跨帧复用
  std::vector<std::size_t>    proxy_body_{}; // proxy -> body
  std::vector<std::uint32_t>  proxy_slot_{}; // proxy -> slot
  Broadphase             broadphase_{};
  DiscreteNarrowphase    narrowphase_{};
  ConstraintListSink     sink_{};
  ManifoldStore          store_{};
  ContactEmitOptions     emit_options_{};
  CollisionWorldStats    stats_{};
  std::uint64_t          layout_signature_ = 0;
  std::uint64_t          last_signature_   = 0;
  std::uint64_t          geometry_epoch_    = 0;
  PairFilter             pair_filter_      = nullptr;
  void                  *pair_filter_user_ = nullptr;
  BodyActiveFilter       body_active_      = nullptr;
  void                  *body_active_user_ = nullptr;
};

} // namespace collision
