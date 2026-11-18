/*
 * AVBDDirectSpaceState3D implementation. See the header.
 *
 * Ray and point tests are exact for box/sphere/cylinder in the shape's local frame;
 * shape-shape queries use conservative world AABBs. Both are the MVP query contract on
 * the modernize branch and are documented as such in the header.
 *
 * All coordinates are Godot world space, which is also the modern core's frame.
 */

#include "server/avbd_direct_space_state3d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "server/avbd_physics_server3d.hpp"
#include "server/modern_convert.hpp"

using namespace godot;

namespace {

struct WorldAabb {
    Vec3 lo{0};
    Vec3 hi{0};
};

inline bool overlaps(const WorldAabb &a, const WorldAabb &b) {
    return a.lo.x <= b.hi.x && a.hi.x >= b.lo.x && a.lo.y <= b.hi.y && a.hi.y >= b.lo.y &&
           a.lo.z <= b.hi.z && a.hi.z >= b.lo.z;
}

WorldAabb shape_aabb(const collision::Shape &p_shape, const Transform3D &p_xform) {
    const Vec3 half = p_shape.local_half_extents();
    const Mat3 R = rot_matrix(to_modern(p_xform.basis.get_rotation_quaternion()));
    const Vec3 world_half{
            std::abs(R[0][0]) * half.x + std::abs(R[1][0]) * half.y + std::abs(R[2][0]) * half.z,
            std::abs(R[0][1]) * half.x + std::abs(R[1][1]) * half.y + std::abs(R[2][1]) * half.z,
            std::abs(R[0][2]) * half.x + std::abs(R[1][2]) * half.y + std::abs(R[2][2]) * half.z};
    const Vec3 center = to_modern(p_xform.origin);
    return WorldAabb{center - world_half, center + world_half};
}

// Ray in the shape's local frame. Fills the nearest non-negative hit.
bool ray_local(const collision::Shape &s, Vec3 o, Vec3 d, Scalar_t max_t, Scalar_t &r_t, Vec3 &r_point,
        Vec3 &r_normal, bool &r_inside) {
    Scalar_t best = max_t;
    Vec3 best_point{0};
    Vec3 best_normal{0};
    bool hit = false;
    bool inside = false;

    const auto accept = [&](Scalar_t t, Vec3 point, Vec3 normal) {
        if (t < 0 || t > max_t)
            return;
        if (!hit || t < best) {
            hit = true;
            best = t;
            best_point = point;
            best_normal = normal;
        }
    };

    switch (s.type) {
    case collision::ShapeType::Sphere: {
        const Scalar_t r = s.radius;
        const Scalar_t b = glm::dot(o, d);
        const Scalar_t c = glm::dot(o, o) - r * r;
        if (c < 0) {
            inside = true; // origin inside
            hit = true;
            best = 0;
            best_point = o;
            best_normal = -d;
            break;
        }
        const Scalar_t disc = b * b - c; // |d| == 1
        if (disc < 0)
            break;
        const Scalar_t t = -b - std::sqrt(disc);
        if (t >= 0) {
            const Vec3 p = o + d * t;
            accept(t, p, glm::normalize(p));
        }
        break;
    }
    case collision::ShapeType::Cylinder: {
        const Scalar_t r = s.radius;
        const Scalar_t h = s.half_height;
        // Side: x^2 + z^2 = r^2.
        const Scalar_t a = d.x * d.x + d.z * d.z;
        if (a > 1e-12f) {
            const Scalar_t b = o.x * d.x + o.z * d.z;
            const Scalar_t c = o.x * o.x + o.z * o.z - r * r;
            const Scalar_t disc = b * b - a * c;
            if (disc >= 0) {
                const Scalar_t t = (-b - std::sqrt(disc)) / a;
                const Vec3 p = o + d * t;
                if (t >= 0 && std::abs(p.y) <= h) {
                    accept(t, p, Vec3{p.x / r, 0, p.z / r});
                }
            }
        }
        // Caps: y = +-h, point inside the disk.
        if (std::abs(d.y) > 1e-12f) {
            for (int side = -1; side <= 1; side += 2) {
                const Scalar_t t = (static_cast<Scalar_t>(side) * h - o.y) / d.y;
                if (t < 0)
                    continue;
                const Vec3 p = o + d * t;
                if (p.x * p.x + p.z * p.z <= r * r) {
                    accept(t, p, Vec3{0, static_cast<Scalar_t>(side), 0});
                }
            }
        }
        if (o.x * o.x + o.z * o.z < r * r && std::abs(o.y) < h) {
            inside = true;
            hit = true;
            best = 0;
            best_point = o;
            best_normal = -d;
        }
        break;
    }
    case collision::ShapeType::Box:
    default: {
        const Vec3 half = s.half;
        Scalar_t tmin = 0;
        Scalar_t tmax = max_t;
        int axis = -1;
        Scalar_t sign = 0;
        bool separated = false;
        for (int i = 0; i < 3; ++i) {
            if (std::abs(d[i]) < 1e-12f) {
                if (o[i] < -half[i] || o[i] > half[i])
                    separated = true;
                continue;
            }
            Scalar_t t1 = (-half[i] - o[i]) / d[i];
            Scalar_t t2 = (half[i] - o[i]) / d[i];
            if (t1 > t2)
                std::swap(t1, t2);
            if (t1 > tmin) {
                tmin = t1;
                axis = i;
                sign = d[i] > 0 ? Scalar_t(-1) : Scalar_t(1);
            }
            tmax = std::min(tmax, t2);
        }
        if (separated || tmin > tmax) {
            break;
        }
        if (axis < 0) {
            // Origin inside the box.
            inside = true;
            hit = true;
            best = 0;
            best_point = o;
            best_normal = -d;
        } else {
            Vec3 normal{0};
            normal[axis] = sign;
            accept(tmin, o + d * tmin, normal);
        }
        break;
    }
    }
    if (!hit)
        return false;
    r_t = best;
    r_point = best_point;
    r_normal = best_normal;
    r_inside = inside;
    return true;
}

bool ray_shape(const collision::Shape &p_shape, const Transform3D &p_xform, Vec3 p_origin, Vec3 p_dir,
        Scalar_t p_max_t, Scalar_t &r_t, Vec3 &r_point, Vec3 &r_normal, bool &r_inside) {
    const Quat q = to_modern(p_xform.basis.get_rotation_quaternion());
    const Mat3 R = rot_matrix(q);
    const Vec3 center = to_modern(p_xform.origin);
    const Vec3 local_o = glm::transpose(R) * (p_origin - center);
    const Vec3 local_d = glm::transpose(R) * p_dir;
    Scalar_t t = 0;
    Vec3 point, normal;
    bool inside = false;
    if (!ray_local(p_shape, local_o, local_d, p_max_t, t, point, normal, inside))
        return false;
    r_t = t;
    r_point = center + R * point;
    r_normal = R * normal;
    r_inside = inside;
    return true;
}

bool point_in_shape(const collision::Shape &s, const Transform3D &p_xform, Vec3 p_point) {
    const Mat3 R = rot_matrix(to_modern(p_xform.basis.get_rotation_quaternion()));
    const Vec3 local = glm::transpose(R) * (p_point - to_modern(p_xform.origin));
    switch (s.type) {
    case collision::ShapeType::Sphere:
        return glm::dot(local, local) <= s.radius * s.radius;
    case collision::ShapeType::Cylinder:
        return local.x * local.x + local.z * local.z <= s.radius * s.radius &&
               std::abs(local.y) <= s.half_height;
    case collision::ShapeType::Box:
    default:
        return std::abs(local.x) <= s.half.x && std::abs(local.y) <= s.half.y &&
               std::abs(local.z) <= s.half.z;
    }
}

Vec3 closest_point(const collision::Shape &s, const Transform3D &p_xform, Vec3 p_point) {
    const Mat3 R = rot_matrix(to_modern(p_xform.basis.get_rotation_quaternion()));
    const Vec3 center = to_modern(p_xform.origin);
    const Vec3 local = glm::transpose(R) * (p_point - center);
    Vec3 out = local;
    switch (s.type) {
    case collision::ShapeType::Sphere: {
        const Scalar_t len = glm::length(local);
        if (len > s.radius && len > 1e-9f)
            out = local * (s.radius / len);
        break;
    }
    case collision::ShapeType::Cylinder: {
        if (std::abs(out.y) > s.half_height)
            out.y = std::copysign(s.half_height, out.y);
        const Scalar_t radial = std::sqrt(out.x * out.x + out.z * out.z);
        if (radial > s.radius && radial > 1e-9f) {
            const Scalar_t k = s.radius / radial;
            out.x *= k;
            out.z *= k;
        }
        break;
    }
    case collision::ShapeType::Box:
    default:
        out.x = std::clamp(out.x, -s.half.x, s.half.x);
        out.y = std::clamp(out.y, -s.half.y, s.half.y);
        out.z = std::clamp(out.z, -s.half.z, s.half.z);
        break;
    }
    return center + R * out;
}

// AABB-vs-AABB sweep of the query box along `motion`. Returns true on a hit in [0, 1].
bool sweep_aabb(const WorldAabb &p_base, Vec3 p_motion, const WorldAabb &p_obstacle, Scalar_t &r_safe,
        Scalar_t &r_depth, Vec3 &r_normal) {
    Scalar_t enter = 0;
    Scalar_t exit = 1;
    int axis = -1;
    Scalar_t sign = 0;
    for (int i = 0; i < 3; ++i) {
        const Scalar_t m = p_motion[i];
        if (std::abs(m) < 1e-9f) {
            if (p_base.hi[i] < p_obstacle.lo[i] || p_base.lo[i] > p_obstacle.hi[i])
                return false;
            continue;
        }
        Scalar_t t1, t2;
        if (m > 0) {
            t1 = (p_obstacle.lo[i] - p_base.hi[i]) / m;
            t2 = (p_obstacle.hi[i] - p_base.lo[i]) / m;
        } else {
            t1 = (p_obstacle.hi[i] - p_base.lo[i]) / m;
            t2 = (p_obstacle.lo[i] - p_base.hi[i]) / m;
        }
        if (t1 > enter) {
            enter = t1;
            axis = i;
            sign = m > 0 ? Scalar_t(-1) : Scalar_t(1);
        }
        exit = std::min(exit, t2);
    }
    if (axis < 0 || enter > exit || enter > 1 || exit < 0)
        return false;
    r_safe = std::max<Scalar_t>(enter, 0);
    r_depth = std::max<Scalar_t>(0, std::min(p_base.hi[axis], p_obstacle.hi[axis]) -
                                        std::max(p_base.lo[axis], p_obstacle.lo[axis]));
    r_normal = Vec3{0};
    r_normal[axis] = sign;
    return true;
}

} // namespace

void AVBDDirectSpaceState3D::bind(AVBDPhysicsServer3D *p_server, const RID &p_space) {
    server = p_server;
    space = p_space;
}

std::vector<AVBDDirectSpaceState3D::QueryTarget> AVBDDirectSpaceState3D::collect_targets(
        const Vec3 &p_min, const Vec3 &p_max, uint32_t p_mask, bool p_bodies, bool p_areas,
        bool p_pick_ray) const {
    std::vector<QueryTarget> out;
    if (server == nullptr) {
        return out;
    }
    const WorldAabb query{p_min, p_max};

    if (p_bodies) {
        for (const AVBDPhysicsServer3D::QueryShapeHit &hit : server->query_body_shapes(space, p_mask)) {
            if (is_body_excluded_from_query(hit.rid)) {
                continue;
            }
            if (p_pick_ray && !hit.ray_pickable) {
                continue;
            }
            if (!overlaps(query, shape_aabb(hit.shape, hit.transform))) {
                continue;
            }
            out.push_back(QueryTarget{hit.shape, hit.transform, hit.rid, hit.collider_id,
                    hit.shape_index, hit.linear_velocity, false});
        }
    }
    if (p_areas) {
        for (const AVBDPhysicsServer3D::QueryShapeHit &hit : server->query_area_shapes(space, p_mask)) {
            if (is_body_excluded_from_query(hit.rid)) {
                continue;
            }
            if (p_pick_ray && !hit.ray_pickable) {
                continue;
            }
            if (!overlaps(query, shape_aabb(hit.shape, hit.transform))) {
                continue;
            }
            out.push_back(QueryTarget{hit.shape, hit.transform, hit.rid, hit.collider_id,
                    hit.shape_index, Vector3(), true});
        }
    }
    return out;
}

bool AVBDDirectSpaceState3D::_intersect_ray(const Vector3 &p_from, const Vector3 &p_to, uint32_t p_collision_mask,
        bool p_collide_with_bodies, bool p_collide_with_areas, bool p_hit_from_inside,
        bool p_hit_back_faces, bool p_pick_ray, PhysicsServer3DExtensionRayResult *r_result) {
    (void)p_hit_back_faces;
    if (server == nullptr || r_result == nullptr) {
        return false;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return false;
    }

    const Vec3 origin = to_modern(p_from);
    const Vec3 end = to_modern(p_to);
    Vec3 dir = end - origin;
    const Scalar_t max_distance = glm::length(dir);
    if (max_distance <= 1e-6f)
        return false;
    dir /= max_distance;
    const WorldAabb segment{glm::min(origin, end) - Vec3{1e-4f}, glm::max(origin, end) + Vec3{1e-4f}};

    struct Best {
        Scalar_t t = 0;
        Vec3 point{0};
        Vec3 normal{0};
        RID rid;
        uint64_t collider_id = 0;
        int32_t shape = 0;
        bool inside = false;
        bool valid = false;
    } best;

    const auto consider = [&](const QueryTarget &target) {
        const WorldAabb aabb = shape_aabb(target.shape, target.transform);
        if (!overlaps(segment, aabb))
            return;
        Scalar_t t = 0;
        Vec3 point, normal;
        bool inside = false;
        if (!ray_shape(target.shape, target.transform, origin, dir, max_distance, t, point, normal, inside))
            return;
        if (inside && !p_hit_from_inside)
            return;
        if (inside) {
            best = Best{0, origin, -dir, target.rid, target.collider_id, target.shape_index, true, true};
            return;
        }
        if (best.valid && best.inside)
            return;
        if (best.valid && t >= best.t)
            return;
        best = Best{t, point, normal, target.rid, target.collider_id, target.shape_index, false, true};
    };

    if (p_collide_with_bodies) {
        for (const QueryTarget &target :
                collect_targets(segment.lo, segment.hi, p_collision_mask, true, false, p_pick_ray)) {
            consider(target);
            if (best.valid && best.inside)
                break;
        }
    }
    if (p_collide_with_areas && !(best.valid && best.inside)) {
        for (const QueryTarget &target :
                collect_targets(segment.lo, segment.hi, p_collision_mask, false, true, p_pick_ray)) {
            consider(target);
            if (best.valid && best.inside)
                break;
        }
    }
    if (!best.valid)
        return false;

    r_result->position = to_godot(best.point);
    Vector3 normal = to_godot(best.normal);
    if (!best.inside && normal.dot(p_to - p_from) > 0.0) {
        normal = -normal;
    }
    r_result->normal = normal;
    r_result->rid = best.rid;
    r_result->collider_id = best.collider_id;
    r_result->collider = nullptr;
    r_result->shape = best.shape;
    return true;
}

int32_t AVBDDirectSpaceState3D::_intersect_point(const Vector3 &p_position, uint32_t p_collision_mask,
        bool p_collide_with_bodies, bool p_collide_with_areas,
        PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) {
    if (server == nullptr || r_results == nullptr || p_max_results <= 0) {
        return 0;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return 0;
    }
    const Vec3 point = to_modern(p_position);
    const WorldAabb probe{point, point};

    int32_t count = 0;
    for (const QueryTarget &target : collect_targets(probe.lo, probe.hi, p_collision_mask,
            p_collide_with_bodies, p_collide_with_areas)) {
        if (!point_in_shape(target.shape, target.transform, point)) {
            continue;
        }
        r_results[count].rid = target.rid;
        r_results[count].collider_id = target.collider_id;
        r_results[count].collider = nullptr;
        r_results[count].shape = target.shape_index;
        ++count;
        if (count >= p_max_results)
            break;
    }
    return count;
}

int32_t AVBDDirectSpaceState3D::_intersect_shape(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) {
    if (server == nullptr || r_results == nullptr || p_max_results <= 0) {
        return 0;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return 0;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return 0;
    }
    collision::Shape query;
    if (!AVBDPhysicsServer3D::query_shape_from_data(*data, p_transform, query)) {
        return 0;
    }
    Transform3D xform = p_transform;
    xform.origin += p_motion;
    WorldAabb qaabb = shape_aabb(query, xform);
    const Scalar_t margin = std::max(0.0f, p_margin);
    qaabb.lo -= Vec3{margin};
    qaabb.hi += Vec3{margin};

    int32_t count = 0;
    for (const QueryTarget &target : collect_targets(qaabb.lo, qaabb.hi, p_collision_mask,
            p_collide_with_bodies, p_collide_with_areas)) {
        WorldAabb taabb = shape_aabb(target.shape, target.transform);
        taabb.lo -= Vec3{margin};
        taabb.hi += Vec3{margin};
        if (!overlaps(qaabb, taabb)) {
            continue;
        }
        r_results[count].rid = target.rid;
        r_results[count].collider_id = target.collider_id;
        r_results[count].collider = nullptr;
        r_results[count].shape = target.shape_index;
        ++count;
        if (count >= p_max_results)
            break;
    }
    return count;
}

bool AVBDDirectSpaceState3D::_cast_motion(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, float *r_closest_safe, float *r_closest_unsafe,
        PhysicsServer3DExtensionShapeRestInfo *r_info) {
    if (r_closest_safe != nullptr)
        *r_closest_safe = 1.0f;
    if (r_closest_unsafe != nullptr)
        *r_closest_unsafe = 1.0f;
    if (server == nullptr || (!p_collide_with_bodies && !p_collide_with_areas)) {
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return false;
    }
    collision::Shape base;
    if (!AVBDPhysicsServer3D::query_shape_from_data(*data, p_transform, base)) {
        return false;
    }

    const Vec3 motion = to_modern(p_motion);
    const Scalar_t margin = std::max(0.0f, p_margin);
    WorldAabb swept = shape_aabb(base, p_transform);
    swept.lo = glm::min(swept.lo, swept.lo + motion) - Vec3{margin};
    swept.hi = glm::max(swept.hi, swept.hi + motion) + Vec3{margin};

    bool found = false;
    Scalar_t best_safe = 1;
    Scalar_t best_depth = 0;
    Vec3 best_normal{0};
    QueryTarget best_target{};

    for (const QueryTarget &target : collect_targets(swept.lo, swept.hi, p_collision_mask,
            p_collide_with_bodies, p_collide_with_areas)) {
        WorldAabb obstacle = shape_aabb(target.shape, target.transform);
        obstacle.lo -= Vec3{margin};
        obstacle.hi += Vec3{margin};
        Scalar_t safe = 0, depth = 0;
        Vec3 normal{0};
        if (!sweep_aabb(shape_aabb(base, p_transform), motion, obstacle, safe, depth, normal)) {
            continue;
        }
        if (!found || safe < best_safe) {
            found = true;
            best_safe = safe;
            best_depth = depth;
            best_normal = normal;
            best_target = target;
        }
    }
    if (!found) {
        return false;
    }

    if (r_closest_safe != nullptr)
        *r_closest_safe = best_safe;
    if (r_closest_unsafe != nullptr)
        *r_closest_unsafe = best_safe;
    if (r_info != nullptr) {
        Vec3 point = to_modern(p_transform.origin);
        point += motion * best_safe;
        r_info->point = to_godot(point);
        r_info->normal = to_godot(best_normal);
        r_info->rid = best_target.rid;
        r_info->collider_id = best_target.collider_id;
        r_info->shape = best_target.shape_index;
        r_info->linear_velocity = best_target.linear_velocity;
    }
    (void)best_depth;
    return true;
}

bool AVBDDirectSpaceState3D::_collide_shape(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, void *r_results, int32_t p_max_results, int32_t *r_result_count) {
    if (r_result_count != nullptr)
        *r_result_count = 0;
    if (server == nullptr || r_results == nullptr || r_result_count == nullptr) {
        return false;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return false;
    }
    collision::Shape query;
    if (!AVBDPhysicsServer3D::query_shape_from_data(*data, p_transform, query)) {
        return false;
    }
    Transform3D xform = p_transform;
    xform.origin += p_motion;

    WorldAabb qaabb = shape_aabb(query, xform);
    const Scalar_t margin = std::max(0.0f, p_margin);
    qaabb.lo -= Vec3{margin};
    qaabb.hi += Vec3{margin};

    Vector3 *points = static_cast<Vector3 *>(r_results);
    const int32_t max_pairs = p_max_results / 2;
    int32_t pairs = 0;
    for (const QueryTarget &target : collect_targets(qaabb.lo, qaabb.hi, p_collision_mask,
            p_collide_with_bodies, p_collide_with_areas)) {
        WorldAabb taabb = shape_aabb(target.shape, target.transform);
        taabb.lo -= Vec3{margin};
        taabb.hi += Vec3{margin};
        if (!overlaps(qaabb, taabb)) {
            continue;
        }
        if (pairs >= max_pairs)
            break;
        points[pairs * 2 + 0] = to_godot(closest_point(query, xform, to_modern(target.transform.origin)));
        points[pairs * 2 + 1] = to_godot(
                closest_point(target.shape, target.transform, to_modern(xform.origin)));
        ++pairs;
    }
    *r_result_count = pairs * 2;
    return pairs > 0;
}

bool AVBDDirectSpaceState3D::_rest_info(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, PhysicsServer3DExtensionShapeRestInfo *r_rest_info) {
    if (server == nullptr || r_rest_info == nullptr) {
        return false;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return false;
    }
    collision::Shape query;
    if (!AVBDPhysicsServer3D::query_shape_from_data(*data, p_transform, query)) {
        return false;
    }
    Transform3D xform = p_transform;
    xform.origin += p_motion;

    WorldAabb qaabb = shape_aabb(query, xform);
    const Scalar_t margin = std::max(0.0f, p_margin);
    qaabb.lo -= Vec3{margin};
    qaabb.hi += Vec3{margin};

    bool found = false;
    Scalar_t best_depth = -1;
    Vec3 best_point{0};
    Vec3 best_normal{0};
    QueryTarget best_target{};

    for (const QueryTarget &target : collect_targets(qaabb.lo, qaabb.hi, p_collision_mask,
            p_collide_with_bodies, p_collide_with_areas)) {
        const WorldAabb taabb = shape_aabb(target.shape, target.transform);
        if (!overlaps(qaabb, taabb))
            continue;
        // Minimum overlap axis = separating direction.
        const Vec3 overlap{std::min(qaabb.hi.x, taabb.hi.x) - std::max(qaabb.lo.x, taabb.lo.x),
                std::min(qaabb.hi.y, taabb.hi.y) - std::max(qaabb.lo.y, taabb.lo.y),
                std::min(qaabb.hi.z, taabb.hi.z) - std::max(qaabb.lo.z, taabb.lo.z)};
        int axis = 0;
        if (overlap.y < overlap[axis])
            axis = 1;
        if (overlap.z < overlap[axis])
            axis = 2;
        if (overlap[axis] <= best_depth) {
            continue;
        }
        const Vec3 center_a = to_modern(xform.origin);
        const Vec3 center_b = to_modern(target.transform.origin);
        Vec3 normal{0};
        normal[axis] = center_a[axis] >= center_b[axis] ? Scalar_t(1) : Scalar_t(-1);
        best_depth = overlap[axis];
        best_normal = normal;
        best_point = closest_point(query, xform, center_b);
        best_target = target;
        found = true;
    }
    if (!found) {
        return false;
    }

    r_rest_info->point = to_godot(best_point);
    r_rest_info->normal = to_godot(best_normal);
    r_rest_info->rid = best_target.rid;
    r_rest_info->collider_id = best_target.collider_id;
    r_rest_info->shape = best_target.shape_index;
    r_rest_info->linear_velocity = best_target.linear_velocity;
    return true;
}

Vector3 AVBDDirectSpaceState3D::_get_closest_point_to_object_volume(const RID &p_object, const Vector3 &p_point) const {
    if (server == nullptr) {
        return p_point;
    }
    const Vec3 point = to_modern(p_point);
    const auto consider = [&](const QueryTarget &target, Vec3 &best, Scalar_t &best_distance, bool &found) {
        const Vec3 candidate = closest_point(target.shape, target.transform, point);
        const Scalar_t distance = glm::dot(candidate - point, candidate - point);
        if (!found || distance < best_distance) {
            best_distance = distance;
            best = candidate;
            found = true;
        }
    };

    bool found = false;
    Scalar_t best_distance = std::numeric_limits<Scalar_t>::max();
    Vec3 best = point;

    if (const AVBDPhysicsServer3D::BodyData *body = server->find_body(p_object)) {
        for (std::size_t slot = 0; slot < body->shapes.size(); ++slot) {
            if (slot < body->shape_disabled.size() && body->shape_disabled[slot])
                continue;
            const Transform3D posed = server->body_shape_transform(*body, static_cast<int32_t>(slot));
            collision::Shape shape;
            if (!AVBDPhysicsServer3D::query_shape_from_data(body->shapes[slot], posed, shape))
                continue;
            consider(QueryTarget{shape, posed}, best, best_distance, found);
        }
    }
    if (const AVBDPhysicsServer3D::AreaData *area = server->find_area(p_object)) {
        for (std::size_t i = 0; i < area->shapes.size(); ++i) {
            if (i < area->shape_disabled.size() && area->shape_disabled[i])
                continue;
            const Transform3D posed = area->transform * area->shape_transforms[i];
            collision::Shape shape;
            if (!AVBDPhysicsServer3D::query_shape_from_data(area->shapes[i], posed, shape))
                continue;
            consider(QueryTarget{shape, posed}, best, best_distance, found);
        }
    }
    return found ? to_godot(best) : p_point;
}
