#pragma once

// ============================================================================
//  collision/shape_table.hpp —— body → ShapeInstance 的扁平表
// ----------------------------------------------------------------------------
//  【布局】
//    instances_ 按 body-major、slot 顺序连续存放：
//      body 0 的 slot 0..n0-1，body 1 的 slot 0..n1-1，……
//    bodies_[b] = {begin, count} 指向 instances_ 的一段。
//
//  这个顺序是后面 proxy id 的稳定来源：proxy id 直接等于 instances_ 下标，
//  所以 broadphase 的 (proxy_a, proxy_b) 排序天然等价于
//  (body_a, slot_a, body_b, slot_b) 的字典序。
//
//  【为什么增删是 O(total instances)】
//  这是 setup-time 操作，不是每帧热路径。保持扁平存储 + 重建 range，换来
//  step() 里零分配、零指针追逐，比在 body 上挂 std::vector<ShapeInstance> 更值。
//
//  【确定性】
//  不依赖哈希容器；body 顺序与 slot 顺序都是插入顺序。revision() 每次结构变更
//  自增，CollisionWorld 用它来失效 manifold warm start / topology 缓存。
// ============================================================================

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "collision/shape_instance.hpp"

namespace collision {

class ShapeTable {
public:
  struct Range {
    std::size_t begin = 0;
    std::size_t count = 0;
  };

  void resize_bodies(std::size_t n) {
    if (bodies_.size() == n)
      return;

    std::vector<ShapeInstance> old_instances = std::move(instances_);
    std::vector<Range>         old_bodies    = std::move(bodies_);

    instances_.clear();
    instances_.reserve(old_instances.size());
    bodies_.assign(n, Range{});

    const std::size_t keep = std::min(n, old_bodies.size());
    for (std::size_t b = 0; b < keep; ++b) {
      const Range r    = old_bodies[b];
      bodies_[b].begin = instances_.size();
      bodies_[b].count = r.count;
      for (std::size_t i = 0; i < r.count; ++i)
        instances_.push_back(normalized(old_instances[r.begin + i]));
    }
    for (std::size_t b = keep; b < n; ++b)
      bodies_[b] = Range{instances_.size(), 0};

    ++revision_;
  }

  void set_shapes(std::size_t body, std::span<const ShapeInstance> src) { replace_body(body, src); }

  void add_shapes(std::size_t body, std::span<const ShapeInstance> src) {
    if (body >= bodies_.size())
      resize_bodies(body + 1);
    const std::span<const ShapeInstance> old = shapes(body);
    if (old.empty()) {
      replace_body(body, src);
      return;
    }

    std::vector<ShapeInstance> merged;
    merged.reserve(old.size() + src.size());
    merged.insert(merged.end(), old.begin(), old.end());
    merged.insert(merged.end(), src.begin(), src.end());
    replace_body(body, merged);
  }

  void clear_shapes(std::size_t body) {
    if (body >= bodies_.size())
      return;
    replace_body(body, std::span<const ShapeInstance>{});
  }

  [[nodiscard]] std::span<const ShapeInstance> shapes(std::size_t body) const noexcept {
    if (body >= bodies_.size())
      return {};
    const Range r = bodies_[body];
    return std::span<const ShapeInstance>{instances_.data() + r.begin, r.count};
  }

  [[nodiscard]] std::size_t   body_count() const noexcept { return bodies_.size(); }
  [[nodiscard]] std::size_t   instance_count() const noexcept { return instances_.size(); }
  [[nodiscard]] std::uint64_t revision() const noexcept { return revision_; }

  void clear() {
    instances_.clear();
    bodies_.clear();
    ++revision_;
  }

private:
  [[nodiscard]] static ShapeInstance normalized(ShapeInstance inst) {
    inst.local_rotation = glm::normalize(inst.local_rotation);
    return inst;
  }

  void replace_body(std::size_t body, std::span<const ShapeInstance> src) {
    if (body >= bodies_.size())
      resize_bodies(body + 1);

    std::vector<ShapeInstance> next;
    next.reserve(instances_.size() - bodies_[body].count + src.size());

    std::vector<Range> next_bodies(bodies_.size());
    for (std::size_t b = 0; b < bodies_.size(); ++b) {
      next_bodies[b].begin = next.size();
      if (b == body) {
        for (const ShapeInstance &inst : src)
          next.push_back(normalized(inst));
      } else {
        const Range r = bodies_[b];
        for (std::size_t i = 0; i < r.count; ++i)
          next.push_back(normalized(instances_[r.begin + i]));
      }
      next_bodies[b].count = next.size() - next_bodies[b].begin;
    }

    instances_.swap(next);
    bodies_.swap(next_bodies);
    ++revision_;
  }

  std::vector<ShapeInstance> instances_{};
  std::vector<Range>         bodies_{};
  std::uint64_t              revision_ = 0;
};

} // namespace collision
