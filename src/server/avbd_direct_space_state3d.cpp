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
#include <ranges>
#include <vector>

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
    avbd::Shape s;
    s.type = b.shape;
    s.center = b.positionLin;
    s.rotation = b.positionAng;
    s.half = b.size * 0.5f;
    s.radius = b.size.x();
    s.halfHeight = b.size.z() * 0.5f;
    s.axis = b.positionAng * avbd::float3{0, 0, 1};
    return s;
}

// The solver-side Shape for one server ShapeData record, posed in sim space.
avbd::Shape query_shape_of(const AVBDPhysicsServer3D::ShapeData &data, const Transform3D &p_transform) {
    avbd::Shape s;
    switch (data.type) {
        case PhysicsServer3D::SHAPE_SPHERE:
            s.type = avbd::ShapeType::Sphere;
            s.radius = static_cast<float>(data.radius);
            break;
        case PhysicsServer3D::SHAPE_CYLINDER: {
            s.type = avbd::ShapeType::Cylinder;
            s.radius = static_cast<float>(data.radius);
            s.halfHeight = static_cast<float>(data.height) * 0.5f;
            break;
        }
        default:
            s.type = avbd::ShapeType::Box;
            s.half = to_sim_extents(data.extents) * 0.5f;
            break;
    }
    s.center = to_sim(p_transform.origin);
    s.rotation = to_sim(p_transform.basis.get_rotation_quaternion());
    s.axis = s.rotation * avbd::float3{0, 0, 1};
    return s;
}

// Deepest contact of a pair, or false when separated. Returns the deepest contact by
// penetration depth - the one rest info and collide shape want.
bool deepest_contact(const avbd::Shape &a, const avbd::Shape &b, avbd::Manifold::Contact &r_contact,
        avbd::float3x3 &r_basis) {
    avbd::Manifold::Contact contacts[8];
    const int count = avbd::collideShapes(a, b, contacts, r_basis);
    if (count <= 0) {
        return false;
    }
    // Local anchors span the penetration; the longest span is the deepest point.
    // max_element keeps the first of ties, matching the strict > scan it replaced.
    const std::span span(contacts, static_cast<size_t>(count));
    const auto best = std::ranges::max_element(span, {},
            [](const avbd::Manifold::Contact &c) { return static_cast<avbd::float3>(c.rA - c.rB).norm(); });
    r_contact = *best;
    return true;
}

// World-space (solver frame) position of body A's contact anchor.
avbd::float3 contact_world_a(const avbd::Rigid &a, const avbd::Manifold::Contact &c) {
    return a.positionLin + a.positionAng * c.rA;
}

// The normal from A to B as a Godot-space vector.
Vector3 contact_normal_godot(const avbd::float3x3 &basis) {
    // Row 0 of the basis is the contact normal pointing from B to A; Godot's convention
    // (PhysicsServer3DExtensionRayResult.normal) is the surface normal facing the query,
    // which is from B towards A here as well - it is what the built-in server reports.
    return to_godot(avbd::float3{basis(0, 0), basis(0, 1), basis(0, 2)});
}

} // namespace

bool AVBDDirectSpaceState3D::_intersect_ray(const Vector3 &p_from, const Vector3 &p_to, uint32_t p_collision_mask,
        bool p_collide_with_bodies, bool p_collide_with_areas, bool p_hit_from_inside,
        bool p_hit_back_faces, bool p_pick_ray, PhysicsServer3DExtensionRayResult *r_result) {
    if (server == nullptr || r_result == nullptr || !p_collide_with_bodies) {
        return false;
    }
    (void)p_collide_with_areas; // areas are not in the solver; documented
    (void)p_hit_back_faces;     // every solver shape is a solid, so faces are hit from both sides
    (void)p_pick_ray;

    avbd::float3 origin = to_sim(p_from);
    avbd::float3 end = to_sim(p_to);
    avbd::float3 dir = end - origin;
    const float max_distance = dir.norm();
    if (max_distance <= 1.0e-6f) {
        return false;
    }
    dir = dir * (1.0f / max_distance);

    avbd::float3 local{};
    // Solver::pick already walks every body and tests layer & mask; it returns the closest
    // hit along the whole ray, so clamp by the segment length afterwards.
    avbd::Rigid *hit = server->solver_pick(space, origin, dir, p_collision_mask, local);
    if (hit == nullptr) {
        return false;
    }

    const avbd::float3 world_hit = hit->positionLin + hit->positionAng * local;
    const float distance = (world_hit - origin).norm();
    if (distance > max_distance) {
        return false;
    }

    // hit_from_inside: Godot reports the origin body itself when the ray starts inside a
    // shape. Probe with a 1 mm sphere at the origin; the deepest match wins.
    avbd::Manifold::Contact contact{};
    avbd::float3x3 basis{};
    Vector3 normal(0, 1, 0);
    if (p_hit_from_inside) {
        avbd::Shape probe;
        probe.type = avbd::ShapeType::Sphere;
        probe.radius = 1.0e-3f;
        probe.center = origin;
        probe.rotation = avbd::quat::Identity();
        probe.axis = avbd::float3{0, 0, 1};
        for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
            if (!deepest_contact(probe, shape_of(*candidate.rigid), contact, basis)) {
                continue;
            }
            // The ray starts inside this body: report the origin as the hit.
            r_result->position = p_from;
            r_result->normal = -to_godot(dir);
            r_result->rid = server->solver_pick_rid(space, candidate.rigid);
            r_result->collider_id = server->body_instance_id(space, candidate.rigid);
            r_result->collider = nullptr;
            r_result->shape = 0;
            return true;
        }
    }

    // Normal: probe a 1 mm sphere pushed into the surface at the hit point, which recovers
    // the exact surface frame collideShapes would produce there.
    avbd::Shape probe;
    probe.type = avbd::ShapeType::Sphere;
    probe.radius = 1.0e-3f;
    probe.center = world_hit - dir * 1.0e-3f;
    probe.rotation = hit->positionAng;
    probe.axis = hit->positionAng * avbd::float3{0, 0, 1};
    if (deepest_contact(probe, shape_of(*hit), contact, basis)) {
        normal = to_godot(avbd::float3{-basis(0, 0), -basis(0, 1), -basis(0, 2)});
        if (normal.dot(p_to - p_from) > 0.0) {
            normal = -normal; // face the query origin
        }
    }

    r_result->position = to_godot(world_hit);
    r_result->normal = normal;
    r_result->rid = server->solver_pick_rid(space, hit);
    r_result->collider_id = server->body_instance_id(space, hit);
    r_result->collider = nullptr;
    r_result->shape = 0;
    return true;
}

int32_t AVBDDirectSpaceState3D::_intersect_point(const Vector3 &p_position, uint32_t p_collision_mask,
        bool p_collide_with_bodies, bool p_collide_with_areas,
        PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) {
    (void)p_collide_with_areas;
    if (server == nullptr || r_results == nullptr || !p_collide_with_bodies || p_max_results <= 0) {
        return 0;
    }

    // A point query is a tiny sphere sweep with zero radius margin: use a 1 mm probe.
    avbd::Shape probe;
    probe.type = avbd::ShapeType::Sphere;
    probe.radius = 1.0e-3f;
    probe.center = to_sim(p_position);
    probe.rotation = avbd::quat::Identity();
    probe.axis = avbd::float3{0, 0, 1};

    int32_t count = 0;
    for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
        avbd::Manifold::Contact contact{};
        avbd::float3x3 basis{};
        if (!deepest_contact(probe, shape_of(*candidate.rigid), contact, basis)) {
            continue;
        }
        r_results[count].rid = server->solver_pick_rid(space, candidate.rigid);
        r_results[count].collider_id = server->body_instance_id(space, candidate.rigid);
        r_results[count].collider = nullptr;
        r_results[count].shape = 0;
        count++;
        if (count >= p_max_results) {
            break;
        }
    }
    return count;
}

int32_t AVBDDirectSpaceState3D::_intersect_shape(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, double p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) {
    (void)p_collide_with_areas;
    if (server == nullptr || r_results == nullptr || !p_collide_with_bodies || p_max_results <= 0) {
        return 0;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return 0;
    }

    avbd::Shape query = query_shape_of(*data, p_transform);
    query.center = query.center + to_sim(p_motion); // the query pose is transform moved by motion
    (void)p_margin; // the solver's COLLISION_MARGIN already keeps contacts stable

    int32_t count = 0;
    for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
        avbd::Manifold::Contact contact{};
        avbd::float3x3 basis{};
        if (!deepest_contact(query, shape_of(*candidate.rigid), contact, basis)) {
            continue;
        }
        r_results[count].rid = server->solver_pick_rid(space, candidate.rigid);
        r_results[count].collider_id = server->body_instance_id(space, candidate.rigid);
        r_results[count].collider = nullptr;
        r_results[count].shape = 0;
        count++;
        if (count >= p_max_results) {
            break;
        }
    }
    return count;
}

bool AVBDDirectSpaceState3D::_cast_motion(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, double p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, float *r_closest_safe, float *r_closest_unsafe,
        PhysicsServer3DExtensionShapeRestInfo *r_info) {
    (void)p_margin;
    (void)p_collide_with_areas;
    (void)r_info;
    if (server == nullptr || !p_collide_with_bodies) {
        if (r_closest_safe != nullptr) {
            *r_closest_safe = 1.0f;
        }
        if (r_closest_unsafe != nullptr) {
            *r_closest_unsafe = 1.0f;
        }
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        if (r_closest_safe != nullptr) {
            *r_closest_safe = 1.0f;
        }
        if (r_closest_unsafe != nullptr) {
            *r_closest_unsafe = 1.0f;
        }
        return false;
    }

    const avbd::float3 motion = to_sim(p_motion);
    const avbd::Shape base = query_shape_of(*data, p_transform);

    // Binary search the motion fraction for the first colliding pose (8 iterations is
    // 1/256 precision - finer than any gameplay use).
    const auto collides = [&](float t) {
        avbd::Shape posed = base;
        posed.center = posed.center + motion * t;
        for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
            avbd::Manifold::Contact contact{};
            avbd::float3x3 basis{};
            if (deepest_contact(posed, shape_of(*candidate.rigid), contact, basis)) {
                return true;
            }
        }
        return false;
    };

    if (!collides(1.0f)) {
        // The whole motion is safe.
        if (r_closest_safe != nullptr) {
            *r_closest_safe = 1.0f;
        }
        if (r_closest_unsafe != nullptr) {
            *r_closest_unsafe = 1.0f;
        }
        return false; // false = no collision along the motion
    }

    float low = 0.0f;
    float high = 1.0f;
    for (int i = 0; i < 8; i++) {
        const float mid = (low + high) * 0.5f;
        if (collides(mid)) {
            high = mid;
        } else {
            low = mid;
        }
    }
    if (r_closest_safe != nullptr) {
        *r_closest_safe = low;
    }
    if (r_closest_unsafe != nullptr) {
        *r_closest_unsafe = high;
    }
    return true; // true = the motion would collide
}

bool AVBDDirectSpaceState3D::_collide_shape(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, double p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, void *r_results, int32_t p_max_results, int32_t *r_result_count) {
    (void)p_collide_with_areas;
    
    if (server == nullptr || r_result_count == nullptr || !p_collide_with_bodies) {
        if (r_result_count != nullptr) {
            *r_result_count = 0;
        }
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        *r_result_count = 0;
        return false;
    }

    avbd::Shape query = query_shape_of(*data, p_transform);
    query.center = query.center + to_sim(p_motion);
    (void)p_margin;

    // The results buffer is pairs of Vector3 (Godot's convention: contact point on the
    // query shape, then on the collider, alternating).
    Vector3 *points = static_cast<Vector3 *>(r_results);
    const int32_t max_pairs = p_max_results / 2;

    int32_t pairs = 0;
    bool any = false;
    for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
        avbd::Manifold::Contact contact{};
        avbd::float3x3 basis{};
        if (!deepest_contact(query, shape_of(*candidate.rigid), contact, basis)) {
            continue;
        }
        any = true;
        if (pairs < max_pairs) {
            const avbd::float3 world_a = query.center + query.rotation * contact.rA;
            const avbd::Rigid &b = *candidate.rigid;
            const avbd::float3 world_b = b.positionLin + b.positionAng * contact.rB;
            points[pairs * 2 + 0] = to_godot(world_a);
            points[pairs * 2 + 1] = to_godot(world_b);
            pairs++;
        }
        break; // the engine takes the first colliding pair set
    }
    *r_result_count = pairs * 2;
    return any;
}

bool AVBDDirectSpaceState3D::_rest_info(const RID &p_shape_rid, const Transform3D &p_transform,
        const Vector3 &p_motion, double p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies,
        bool p_collide_with_areas, PhysicsServer3DExtensionShapeRestInfo *r_rest_info) {
    (void)p_collide_with_areas;
    if (server == nullptr || r_rest_info == nullptr || !p_collide_with_bodies) {
        return false;
    }
    const AVBDPhysicsServer3D::ShapeData *data = server->shape_data(p_shape_rid);
    if (data == nullptr) {
        return false;
    }

    avbd::Shape query = query_shape_of(*data, p_transform);
    query.center = query.center + to_sim(p_motion);
    (void)p_margin;

    // The deepest contact over all bodies wins.
    bool found = false;
    float best_depth = -1.0f;
    for (const AVBDPhysicsServer3D::QueryCandidate &candidate : server->query_candidates(space, p_collision_mask)) {
        avbd::Manifold::Contact contact{};
        avbd::float3x3 basis{};
        if (!deepest_contact(query, shape_of(*candidate.rigid), contact, basis)) {
            continue;
        }
        const float depth = (contact.rA - contact.rB).norm();
        if (depth <= best_depth) {
            continue;
        }
        best_depth = depth;
        found = true;
        const avbd::float3 world_a = query.center + query.rotation * contact.rA;
        r_rest_info->point = to_godot(world_a);
        r_rest_info->normal = contact_normal_godot(basis);
        r_rest_info->rid = server->solver_pick_rid(space, candidate.rigid);
        r_rest_info->collider_id = server->body_instance_id(space, candidate.rigid);
        r_rest_info->shape = 0;
        r_rest_info->linear_velocity = server->body_velocity_of(space, candidate.rigid);
    }
    return found;
}

Vector3 AVBDDirectSpaceState3D::_get_closest_point_to_object_volume(const RID &p_object, const Vector3 &p_point) const {
    // Approximation used by the editor's selection; not needed for gameplay queries.
    (void)p_object;
    return p_point;
}
