#pragma once

// ============================================================================
//  base_concept.hpp —— 约束 / 类型擦除 / 三种求解 Policy
// ----------------------------------------------------------------------------
//  State 的位域与接口全部在 state.hpp；这里只允许通过命名接口访问状态，
//  不允许对 State 的 vector 直接下标（pre_compute/step 都一样）。
// ============================================================================

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <variant>
#include <vector>

#include "state.hpp"

// ==================== 1. Traits 主模板 ====================
template <typename C> struct constraint_traits; // 必须特化

template <typename C>
using constraint_scalar_t = typename constraint_traits<C>::Scalar_t;
template <typename C>
using constraint_gradient_t = typename constraint_traits<C>::Gradient_t;

// ==================== 2. 辅助 concept ====================
template <typename V, typename Elem>
concept Indexable = requires(V v, std::size_t i) {
  { v[i] } -> std::convertible_to<Elem>;
};

// 可选：角向梯度（6-DOF 约束才提供）
template <typename C, typename S>
concept HasAngularGradient = requires(const C &c, const S &s) {
  { c.gradient_angular(s) } -> Indexable<Vec3>;
};

// ==================== 3. Constraint concept ====================
template <typename C, typename S>
concept Constraint =
    requires {
      typename constraint_traits<C>::Scalar_t;
      typename constraint_traits<C>::Gradient_t;
    } && Indexable<constraint_gradient_t<C>, Vec3> &&
    requires(const C &c, const S &s) {
      { c.evaluate(s) } -> std::convertible_to<constraint_scalar_t<C>>;
      { c.gradient(s) } -> std::convertible_to<constraint_gradient_t<C>>;
      { c.bodies() } -> Indexable<std::size_t>;
    };

// dual / compliance 成员（XPBD、AVBD 用）
template <typename C>
concept DualConstraint = requires(C &c, Scalar_t v) {
  { c.lambda = v } -> std::convertible_to<Scalar_t &>;
  { c.compliance } -> std::convertible_to<Scalar_t>;
};
