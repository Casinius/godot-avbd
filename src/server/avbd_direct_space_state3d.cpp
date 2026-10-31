/*
 * AVBDDirectSpaceState3D implementation. See the header.
 *
 * All coordinates arriving here are Godot world space; the queries convert to the
 * solver's Z-up frame, walk the space's bodies and answer from collideShapes() - the
 * same pair test the solver's broad phase uses, so query results and contact behaviour
 * cannot disagree.
 *
 * The engine reaches this through World3D.direct_space_state (RayCast3D, intersect_*,
 * ShapeCast3D, move_and_collide's rest info, ...). Queries are answered from the solver's
 * current state, i.e. the poses as of the last step, exactly like the built-in server.
 */

#include "server/avbd_direct_space_state3d.hpp"

#include <algorithm>
#include <limits>
#include <vector>

#include "avbd/query/query_shape.hpp"
#include "nodes/godot_convert.hpp"
#include "server/avbd_physics_server3d.hpp"

using namespace godot;
using namespace avbd_godot;

void AVBDDirectSpaceState3D::bind(AVBDPhysicsServer3D *p_server, const RID &p_space) {
    server = p_server;
    space = p_space;
}

namespace {

// The solver-side Shape for one body (the same flatten Manifold::collide performs).
avbd::Shape shape_of(const avbd::Rigid &b) {
    return AVBDPhysicsServer3D::query_shape_from_rigid(b);
}

// The solver-side Shape for one server ShapeData record, posed in sim space. The
// conversion lives on the server so bodies, areas and queries cannot drift apart.
bool query_shape_of(const AVBDPhysicsServer3D::ShapeData &data, const Transform3D &p_transform,
        avbd::Shape &r_shape) {
    return AVBDPhysicsServer3D::query_shape_from_data(data, p_transform, r_shape);
}

// Tight world-space (sim frame) AABB of a query shape - the same bound the broad phase
// builds its leaves with, so the candidate list prunes the query.
void shape_aabb(const avbd::Shape &s, avbd::float3 &r_min, avbd::float3 &r_max) {
    avbd::query::shapeAabb(s, r_min, r_max);
}

// // World-space (solver frame) position of body A's contact anchor.
// avbd::float3 contact_world_a(const avbd::Rigid &a, const avbd::Manifold::Contact &c) {
//     return a.positionLin + qrot(a.positionAng, c.rA);
// }

// The normal from A to B as a Godot-space vector.
Vector3 contact_normal_godot(const avbd::float3x3 &basis) {
    // Row 0 of the basis is the contact normal pointing from B to A; Godot's convention
    // (PhysicsServer3DExtensionRayResult.normal) is the surface normal facing the query,
    // which is from B towards A here as well - it is what the built-in server reports.
    return to_godot(avbd::float3{basis[0][0], basis[1][0], basis[2][0]});
}

} // namespace

std::vector<AVBDDirectSpaceState3D::QueryTarget> AVBDDirectSpaceState3D::collect_targets(
        const avbd::float3 &p_min, const avbd::float3 &p_max, uint32_t p_mask, bool p_bodies,
        bool p_areas, bool p_pick_ray) const {
    std::vector<QueryTarget> out;
    if (server == nullptr) {
        return out;
    }

    if (p_bodies) {
        for (const AVBDPhysicsServer3D::QueryCandidate &candidate :
                server->query_candidates_aabb(space, p_min, p_max, p_mask)) {
            const RID rid = server->solver_pick_rid(space, candidate.rigid);
            if (is_body_excluded_from_query(rid)) {
                continue;
            }
            if (p_pick_ray && !candidate.ray_pickable) {
                continue;
            }
            out.push_back(QueryTarget{shape_of(*candidate.rigid), rid,
                    server->body_instance_id(space, candidate.rigid), 0,
                    server->body_velocity_of(space, candidate.rigid), false});
        }
    }

    if (p_areas) {
        for (const AVBDPhysicsServer3D::QueryAreaCandidate &candidate :
                server->query_area_candidates_aabb(space, p_min, p_max, p_mask)) {
            if (is_body_excluded_from_query(candidate.rid)) {
                continue;
            }
            if (p_pick_ray && !candidate.ray_pickable) {
                continue;
            }
            out.push_back(QueryTarget{candidate.shape, candidate.rid, candidate.instance_id,
                    candidate.shape_index, Vector3(), true});
        }
    }
    return out;
}

bool AVBDDirectSpaceState3D::_intersect_ray(const Vector3 &p_from, const Vector3 &p_to, uint32_t p_collision_mask,
        bool p_collide_with_bodies, bool p_collide_with_areas, bool p_hit_from_inside,
        bool p_hit_back_faces, bool p_pick_ray, PhysicsServer3DExtensionRayResult *r_result) {
    (void)p_hit_back_faces; // every query shape is a solid, so faces are hit from both sides
    if (server == nullptr || r_result == nullptr) {
        return false;
    }
    if (!p_collide_with_bodies && !p_collide_with_areas) {
        return false;
    }

    const avbd::float3 origin = to_sim(p_from);
    const avbd::float3 end = to_sim(p_to);
    avbd::float3 dir = end - origin;
    const float max_distance = avbd::length(dir);
    if (max_distance <= 1.0e-6f) {
        return false;
    }
    dir = dir * (1.0f / max_distance);

    // Candidate broad phase: the segment's AABB. The exact ray test decides.
    const float pad = 1.0e-4f;
    const avbd::float3 seg_min{std::min(origin.x, end.x) - pad, std::min(origin.y, end.y) - pad,
            std::min(origin.z, end.z) - pad};
    const avbd::float3 seg_max{std::max(origin.x, end.x) + pad, std::max(origin.y, end.y) + pad,
            std::max(origin.z, end.z) + pad};

    struct Hit {
        float t = 0.0f;
        avbd::float3 point{0, 0, 0};
        avbd::float3 normal{0, 0, 1};
        RID rid;
        uint64_t collider_id = 0;
        int shape = 0;
        bool inside = false;
        bool valid = false;
    } best;

    const avbd::query::Ray ray{origin, dir, max_distance};
    const auto consider = [&](const avbd::Shape &shape, const RID &p_rid, uint64_t p_collider_id,
            int p_shape_index, bool p_ray_pickable) {
        if (p_pick_ray && !p_ray_pickable) {
            return;
        }
        avbd::query::RayHit hit{};
        if (!avbd::query::intersectRay(shape, ray, hit)) {
            return;
        }
        if (hit.inside && !p_hit_from_inside) {
            return;
        }
        if (hit.inside) {
            // The engine wants the origin itself when the ray starts inside a solid.
            best.t = 0.0f;
            best.point = origin;
            best.normal = dir * -1.0f;
            best.rid = p_rid;
            best.collider_id = p_collider_id;
            best.shape = p_shape_index;
            best.inside = true;
            best.valid = true;
            return;
        }
        if (hit.t > max_distance || (best.valid && hit.t >= best.t)) {
            return;
        }
        best.t = hit.t;
        best.point = hit.point;
        best.normal = hit.normal;
        best.rid = p_rid;
        best.collider_id = p_collider_id;
        best.shape = p_shape_index;
        best.inside = false;
        best.valid = true;
    };

    if (p_collide_with_bodies) {
        for (const AVBDPhysicsServer3D::QueryCandidate &candidate :
                server->query_candidates_aabb(space, seg_min, seg_max, p_collision_mask)) {
            const RID rid = server->make_rid(candidate.body_id);
            if (is_body_excluded_from_query(rid)) {
                continue;
            }
            consider(shape_of(*candidate.rigid), rid, server->body_instance_id(space, candidate.rigid), 0,
                    candidate.ray_pickable);
            if (best.valid && best.inside) {
                break;
            }
        }
    }

    if (p_collide_with_areas && !(best.valid && best.inside)) {
        for (const AVBDPhysicsServer3D::QueryAreaCandidate &candidate :
                server->query_area_candidates_aabb(space, seg_min, seg_max, p_collision_mask)) {
            if (is_body_excluded_from_query(candidate.rid)) {
                continue;
            }
            consider(candidate.shape, candidate.rid, candidate.instance_id, candidate.shape_index,
                    candidate.ray_pickable);
            if (best.valid && best.inside) {
                break;
            }
        }
    }

    if (!best.valid) {
        return false;
    }

    r_result->position = to_godot(best.point);
    Vector3 normal = to_godot(best.normal);
    if (!best.inside && normal.dot(p_to - p_from) > 0.0) {
        normal = -normal; // face the query origin
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

    const avbd::float3 point = to_sim(p_position);
    const float pad = 1.0e-4f;
    const avbd::float3 probe_min = point - avbd::float3{pad, pad, pad};
    const avbd::float3 probe_max = point + avbd::float3{pad, pad, pad};

    int32_t count = 0;
    for (const QueryTarget &target :
            collect_targets(probe_min, probe_max, p_collision_mask, p_collide_with_bodies,
                    p_collide_with_areas)) {
        if (!avbd::query::containsPoint(target.shape, point)) {
            continue;
        }
        r_results[count].rid = target.rid;
        r_results[count].collider_id = target.collider_id;
        r_results[count].collider = nullptr;
        r_results[count].shape = target.shape_index;
        count++;
        if (count >= p_max_results) {
            break;
        }
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
    avbd::Shape query;
    if (!query_shape_of(*data, p_transform, query)) {
        return 0;
    }
    query.center = query.center + to_sim(p_motion); // the query pose is transform moved by motion

    avbd::float3 qmin, qmax;
    shape_aabb(query, qmin, qmax);
    const float margin = std::max(0.0f, p_margin);
    if (margin > 0.0f) {
        const avbd::float3 m{margin, margin, margin};
        qmin -= m;
        qmax += m;
    }

    int32_t count = 0;
    for (const QueryTarget &target :
            collect_targets(qmin, qmax, p_collision_mask, p_collide_with_bodies,
                    p_collide_with_areas)) {
        avbd::Manifold::Contact contacts[8];
        avbd::float3x3 basis{};
        if (avbd::query::shapeContacts(query, target.shape,
                    std::span<avbd::Manifold::Contact>(contacts, 8), basis, margin) <= 0) {
            continue;
        }
        r_results[count].rid = target.rid;
        r_results[count].collider_id = target.collider_id;
        r_results[count].collider = nullptr;
        r_results[count].shape = target.shape_index;
        count++;
        if (count >= p_max_results) {
            break;
        }
    }
    return count;
}

bool AVBDDirectSpaceState3D::_cast_motion(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, float *r_closest_safe, float *r_closest_unsafe,
        PhysicsServer3DExtensionShapeRestInfo *r_info) {
    if (r_closest_safe != nullptr) {
        *r_closest_safe = 1.0f;
    }
    if (r_closest_unsafe != nullptr) {
        *r_closest_unsafe = 1.0f;
    }
    if (server == nullptr || (!p_collide_with_bodies && !p_collide_with_areas)) {
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return false;
    }
    avbd::Shape base;
    if (!query_shape_of(*data, p_transform, base)) {
        return false;
    }

    const avbd::float3 motion = to_sim(p_motion);
    const float margin = std::max(0.0f, p_margin);

    avbd::float3 swept_min, swept_max;
    shape_aabb(base, swept_min, swept_max);
    for (int i = 0; i < 3; i++) {
        const float m = motion[i];
        swept_min[i] = std::min(swept_min[i], swept_min[i] + m);
        swept_max[i] = std::max(swept_max[i], swept_max[i] + m);
    }
    if (margin > 0.0f) {
        const avbd::float3 pad{margin, margin, margin};
        swept_min -= pad;
        swept_max += pad;
    }

    bool found = false;
    avbd::query::SweepHit best;
    QueryTarget best_target{};
    for (const QueryTarget &target :
            collect_targets(swept_min, swept_max, p_collision_mask, p_collide_with_bodies,
                    p_collide_with_areas)) {
        const avbd::query::SweepHit hit = avbd::query::sweepShape(base, motion, target.shape, margin);
        if (!hit.hit) {
            continue;
        }
        if (!found || hit.safe < best.safe) {
            found = true;
            best = hit;
            best_target = target;
        }
    }

    if (!found) {
        return false;
    }

    if (r_closest_safe != nullptr) {
        *r_closest_safe = best.safe;
    }
    if (r_closest_unsafe != nullptr) {
        *r_closest_unsafe = best.unsafe;
    }
    if (r_info != nullptr) {
        r_info->point = to_godot(best.point);
        r_info->normal = to_godot(best.normal);
        r_info->rid = best_target.rid;
        r_info->collider_id = best_target.collider_id;
        r_info->shape = best_target.shape_index;
        r_info->linear_velocity = best_target.linear_velocity;
    }
    return true;
}

bool AVBDDirectSpaceState3D::_collide_shape(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, void *r_results, int32_t p_max_results, int32_t *r_result_count) {
    if (r_result_count != nullptr) {
        *r_result_count = 0;
    }
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
    avbd::Shape query;
    if (!query_shape_of(*data, p_transform, query)) {
        return false;
    }
    query.center = query.center + to_sim(p_motion);

    avbd::float3 qmin, qmax;
    shape_aabb(query, qmin, qmax);
    const float margin = std::max(0.0f, p_margin);
    if (margin > 0.0f) {
        const avbd::float3 m{margin, margin, margin};
        qmin -= m;
        qmax += m;
    }

    // Godot's convention: pairs of Vector3, first on the query shape, then on the collider.
    Vector3 *points = static_cast<Vector3 *>(r_results);
    const int32_t max_pairs = p_max_results / 2;
    int32_t pairs = 0;
    bool any = false;
    for (const QueryTarget &target :
            collect_targets(qmin, qmax, p_collision_mask, p_collide_with_bodies,
                    p_collide_with_areas)) {
        avbd::Manifold::Contact contacts[8];
        avbd::float3x3 basis{};
        const int count = avbd::query::shapeContacts(query, target.shape,
                std::span<avbd::Manifold::Contact>(contacts, 8), basis, margin);
        if (count <= 0) {
            continue;
        }
        any = true;
        for (int i = 0; i < count && pairs < max_pairs; i++) {
            const avbd::Manifold::Contact &contact = contacts[i];
            const avbd::float3 world_a = query.center + avbd::qrot(query.rotation, contact.rA);
            const avbd::float3 world_b = target.shape.center + avbd::qrot(target.shape.rotation, contact.rB);
            points[pairs * 2 + 0] = to_godot(world_a);
            points[pairs * 2 + 1] = to_godot(world_b);
            pairs++;
        }
        if (pairs >= max_pairs) {
            break;
        }
    }
    *r_result_count = pairs * 2;
    return any;
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
    avbd::Shape query;
    if (!query_shape_of(*data, p_transform, query)) {
        return false;
    }
    query.center = query.center + to_sim(p_motion);

    avbd::float3 qmin, qmax;
    shape_aabb(query, qmin, qmax);
    const float margin = std::max(0.0f, p_margin);
    if (margin > 0.0f) {
        const avbd::float3 m{margin, margin, margin};
        qmin -= m;
        qmax += m;
    }

    bool found = false;
    float best_depth = -1.0f;
    avbd::float3 best_point{0, 0, 0};
    avbd::float3x3 best_basis{};
    QueryTarget best_target{};
    for (const QueryTarget &target :
            collect_targets(qmin, qmax, p_collision_mask, p_collide_with_bodies,
                    p_collide_with_areas)) {
        avbd::Manifold::Contact contacts[8];
        avbd::float3x3 basis{};
        const int count = avbd::query::shapeContacts(query, target.shape,
                std::span<avbd::Manifold::Contact>(contacts, 8), basis, margin);
        for (int i = 0; i < count; i++) {
            const avbd::float3 normal = avbd::row(basis, 0);
            const float depth = avbd::query::penetrationDepth(query, target.shape, normal);
            if (depth <= best_depth) {
                continue;
            }
            best_depth = depth;
            best_basis = basis;
            best_point = query.center + avbd::qrot(query.rotation, contacts[i].rA);
            best_target = target;
            found = true;
        }
    }

    if (!found) {
        return false;
    }
    r_rest_info->point = to_godot(best_point);
    r_rest_info->normal = contact_normal_godot(best_basis);
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
    const avbd::float3 point = to_sim(p_point);
    bool found = false;
    float best_distance = std::numeric_limits<float>::infinity();
    avbd::float3 best_point = point;

    const auto consider = [&](const avbd::Shape &shape) {
        const avbd::float3 candidate = avbd::query::closestPoint(shape, point);
        const float distance = avbd::length2(candidate - point);
        if (distance < best_distance) {
            best_distance = distance;
            best_point = candidate;
            found = true;
        }
    };

    if (const AVBDPhysicsServer3D::BodyData *body = server->find_body(p_object)) {
        if (body->rigid != nullptr) {
            consider(AVBDPhysicsServer3D::query_shape_from_rigid(*body->rigid));
        } else if (!body->shapes.empty()) {
            avbd::Shape shape;
            if (AVBDPhysicsServer3D::query_shape_from_data(body->shapes.front(), body->transform, shape)) {
                consider(shape);
            }
        }
    }
    if (const AVBDPhysicsServer3D::AreaData *area = server->find_area(p_object)) {
        for (size_t i = 0; i < area->shapes.size(); i++) {
            if (area->shape_disabled[i]) {
                continue;
            }
            const Transform3D posed = area->transform * area->shape_transforms[i];
            avbd::Shape shape;
            if (AVBDPhysicsServer3D::query_shape_from_data(area->shapes[i], posed, shape)) {
                consider(shape);
            }
        }
    }
    return found ? to_godot(best_point) : p_point;
}
