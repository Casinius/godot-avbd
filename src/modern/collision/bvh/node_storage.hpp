#pragma once

// ============================================================================
//  collision/bvh/node_storage.hpp —— 隐式堆布局的 LBVH 节点存储
// ----------------------------------------------------------------------------
//  【布局契约】
//   节点 i 的孩子恒为 2i+1 / 2i+2，叶子位于 [padded-1, padded-1+n)。因此没有
//   left/right/parent 数组：布局本身就是树。
//   空子树用 lo=+inf / hi=-inf 编码（min > max），valid() 只比较 x 一对分量，
//   所以节点记录只有两个 Vec3，没有第四个 liveness 字段。
//   叶子 payload 是 body index（std::int32_t），不是 Rigid*：碰撞层不认识求解器。
// ============================================================================

#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "phys_core/state.hpp"

namespace collision::bvh {

struct Node {
  Vec3 lo{empty_lo()};
  Vec3 hi{empty_hi()};

  [[nodiscard]] bool valid() const noexcept { return lo.x <= hi.x; }

  [[nodiscard]] static constexpr Vec3 empty_lo() noexcept {
    constexpr Scalar_t inf = std::numeric_limits<Scalar_t>::infinity();
    return Vec3{inf, inf, inf};
  }
  [[nodiscard]] static constexpr Vec3 empty_hi() noexcept {
    constexpr Scalar_t inf = std::numeric_limits<Scalar_t>::infinity();
    return Vec3{-inf, -inf, -inf};
  }
};

class NodeStorage {
public:
  void clear() noexcept {
    nodes_.clear();
    body_index_.clear();
  }

  [[nodiscard]] int  size() const noexcept { return static_cast<int>(nodes_.size()); }
  [[nodiscard]] bool empty() const noexcept { return nodes_.empty(); }

  [[nodiscard]] std::span<const Node> nodes() const noexcept { return nodes_; }
  [[nodiscard]] const Node           &node(int index) const noexcept {
    return nodes_[static_cast<std::size_t>(index)];
  }
  [[nodiscard]] bool                valid(int index) const noexcept { return node(index).valid(); }
  [[nodiscard]] const std::int32_t &body_index(int index) const noexcept {
    return body_index_[static_cast<std::size_t>(index)];
  }

  // 隐式布局支持：按精确节点数 resize，所有节点先标成空。
  void resize_implicit(int node_count) noexcept {
    nodes_.resize(static_cast<std::size_t>(node_count)); // value-init：空节点
    body_index_.assign(static_cast<std::size_t>(node_count), std::int32_t(-1));
  }

  void set_leaf(int slot, Vec3 lo, Vec3 hi, std::int32_t body_index) noexcept {
    const auto i   = static_cast<std::size_t>(slot);
    nodes_[i].lo   = lo;
    nodes_[i].hi   = hi;
    body_index_[i] = body_index;
  }

  void set_internal(int slot, Vec3 lo, Vec3 hi, bool subtree_valid) noexcept {
    const auto i = static_cast<std::size_t>(slot);
    nodes_[i].lo = subtree_valid ? lo : Node::empty_lo();
    nodes_[i].hi = subtree_valid ? hi : Node::empty_hi();
  }

private:
  std::vector<Node>         nodes_;
  std::vector<std::int32_t> body_index_;
};

} // namespace collision::bvh
