/*
 * Space-query geometry for the AVBD solver core.
 *
 * This is an extension layer, not part of the solver. It only reads the public
 * avbd::Shape description and calls the existing collideShapes() narrow phase, so the
 * teaching-oriented core (solver.h/solver.cpp/collide.cpp) stays untouched.
 *
 * Shapes use exactly the proxy the solver simulates: a capsule query is a cylinder, a
 * convex-hull query is the hull's local AABB box. Keeping the query and the contact
 * model on the same geometry is deliberate - a raycast may not report a surface the
 * solver does not collide against.
 */

#pragma once

#include <cstddef>
#include <span>

#include "avbd/solver.h"

namespace avbd::query {

// A world-space (solver frame) ray. `direction` need not be unit length; `t` is in
// units of `direction`, so the usual callers pass a normalised direction.
struct Ray {
    float3 origin{0, 0, 0};
    float3 direction{0, 0, 1};
    float maxT = 0.0f; // <= 0 means unbounded
};

struct RayHit {
    float t = 0.0f;
    float3 point{0, 0, 0};
    float3 normal{0, 0, 1}; // outward surface normal at `point`
    bool inside = false;    // ray origin was inside the shape
};

// Exact ray/shape test for box, sphere and the solver's Z-axis cylinder. Returns false
// when the ray misses or the hit lies beyond `ray.maxT`.
[[nodiscard]] bool intersectRay(const Shape &shape, const Ray &ray, RayHit &r_hit) noexcept;

// Closest point on the shape's surface (or volume boundary when the point is inside).
[[nodiscard]] float3 closestPoint(const Shape &shape, const float3 &point) noexcept;

// Point containment, optionally inflated by `margin` (used by the point query).
[[nodiscard]] bool containsPoint(const Shape &shape, const float3 &point, float margin = 0.0f) noexcept;

// Support radius of a shape along `axis` (the box/cylinder extent projected on the axis).
[[nodiscard]] float supportRadius(const Shape &shape, const float3 &axis) noexcept;
// Penetration along a contact normal: support(a)+support(b)-|centre separation|. This is
// the SAT overlap on that axis and gives rest queries a real depth instead of a distance
// between two arbitrary anchor points.
[[nodiscard]] float penetrationDepth(const Shape &a, const Shape &b, const float3 &normal) noexcept;

// Tight world-space AABB of the shape. The cylinder bound is the bound of its local-axis
// bounding box, which is what the solver's broad phase uses.
void shapeAabb(const Shape &shape, float3 &r_min, float3 &r_max) noexcept;

// The existing narrow phase, with a name that reads naturally from query code.
[[nodiscard]] int shapeContacts(const Shape &a, const Shape &b,
        std::span<Manifold::Contact> contacts, float3x3 &r_basis,
        float margin = 0.0f);

// Translation sweep of `base` along `motion` against a static `target`.
//
// `safe` is the last tested fraction without a contact, `unsafe` the first with one;
// `margin` has the same meaning as collideShapes()' speculative margin. The sweep is
// bounded by the swept AABB overlap interval, so a fast motion cannot tunnel through a
// target just because its final pose is past it.
struct SweepHit {
    bool hit = false;
    float safe = 1.0f;
    float unsafe = 1.0f;
    float depth = 0.0f;
    float3 point{0, 0, 0};
    float3 normal{0, 0, 1}; // points from the target towards the moving shape
};

[[nodiscard]] SweepHit sweepShape(const Shape &base, const float3 &motion,
        const Shape &target, float margin = 0.0f);

} // namespace avbd::query
