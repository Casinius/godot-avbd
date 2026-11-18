#pragma once

// ============================================================================
//  collision/concepts.hpp —— 碰撞层的可替换接口约束
// ----------------------------------------------------------------------------
//  模仿 phys_core 的做法：接口只描述“能读什么/能调什么”，不暴露拥有关系，
//  也不把线程/求解器拉进来。四个 concept 对应四层：
//
//    BodyShapeSource  : pose + body→ShapeInstance 的只读视图（适配器在 world 里）
//    BroadphaseModel  : update(aabbs) + pairs() + nodes()
//    NarrowphaseModel : collide(a, b, out, basis) -> count
//    ContactSink      : emit(manifolds, options, list) -> count
//
//  它们只被 static_assert / requires 使用，不引入虚函数、不引入 type-erasure
//  的堆分配，也不允许实现方偷偷拥有 State。
// ============================================================================

#include <concepts>
#include <cstddef>
#include <span>

#include "collision/broadphase.hpp"
#include "collision/contact_constraints.hpp"
#include "collision/manifold.hpp"
#include "collision/narrowphase.hpp"
#include "collision/shape.hpp"
#include "collision/shape_instance.hpp"

namespace collision {

template <typename T>
concept BodyShapeSource = requires(const T &src, std::size_t b, std::size_t slot) {
  { src.size() } -> std::convertible_to<std::size_t>;
  { src.shape_count(b) } -> std::convertible_to<std::size_t>;
  { src.shape_instance(b, slot) } -> std::convertible_to<ShapeInstance>;
  { src.position(b) } -> std::convertible_to<Vec3>;
  { src.rotation(b) } -> std::convertible_to<Quat>;
  { src.is_static(b) } -> std::convertible_to<bool>;
};

template <typename T>
concept BroadphaseModel = requires(T &bp, std::span<const Aabb> aabbs) {
  { bp.update(aabbs) } -> std::same_as<void>;
  { bp.pairs() } -> std::convertible_to<std::span<const BodyPair>>;
  { bp.nodes() } -> std::convertible_to<const bvh::NodeStorage &>;
};

template <typename T>
concept NarrowphaseModel = requires(const T &np, const ShapeView &a, const ShapeView &b,
                                    std::span<ContactPoint> out, BasisFrame &basis) {
  { np.collide(a, b, out, basis) } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept ContactSink = requires(const T &sink, std::span<const Manifold> manifolds,
                               const ContactEmitOptions &opt, ConstraintList &list) {
  { sink.emit(manifolds, opt, list) } -> std::convertible_to<std::size_t>;
};

} // namespace collision
