#pragma once

// ============================================================================
//  collision/shape_instance.hpp —— 刚体上的一个形状实例
// ----------------------------------------------------------------------------
//  【为什么需要这一层】
//  求解器仍然以 body 为单位；但客户端（尤其是从引擎导入的模型）经常把多个
//  collider 挂在同一个刚体下。Shape 本身是纯局部几何，不能携带 local transform；
//  所以这里把“形状 + 相对 body COM 的局部位姿”打包成一个实例。
//
//  【不变量】
//   1) local_position/local_rotation 相对 **body COM**，不是相对 mesh 原点；
//   2) 不允许嵌套 compound：客户端/资产导入工具必须先把嵌套摊平；
//   3) 非均匀 scale 不在这一层支持，要求客户端烘进 Shape 的 half/vertices；
//   4) 该类型是纯值类型，不拥有 mesh；convex mesh 由 CollisionWorld 的注册表拥有。
// ============================================================================

#include <cstdint>

#include "collision/shape.hpp"

namespace collision {

struct ShapeInstance {
  Shape shape{};
  Vec3  local_position{0};
  Quat  local_rotation = identity_quat();
  // 可选：客户端自己的 collider/part id，碰撞层不解释它，只原样保存。
  std::uint32_t user_id = 0;
};

} // namespace collision
