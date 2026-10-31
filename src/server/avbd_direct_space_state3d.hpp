/*
 * AVBDDirectSpaceState3D - the PhysicsDirectSpaceState3D behind space queries.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * The engine reaches this through World3D.direct_space_state (RayCast3D, intersect_*,
 * ShapeCast3D, ...). One instance per space, created lazily by the server.
 *
 * Queries walk the space's solver bodies and areas, coarse-filter by the query AABB and
 * answer with the same collideShapes() the solver uses, so query results and contact
 * behaviour cannot disagree. Exclusion lists (`params.exclude`) are honoured through
 * PhysicsDirectSpaceState3DExtension::is_body_excluded_from_query(); `ray_pickable` is
 * read for body and area hits. `cast_motion` sweeps a shape instead of testing only its
 * final pose; `_body_test_motion` uses the same sweep for CharacterBody3D/test_move.
 */

#ifndef AVBD_DIRECT_SPACE_STATE3D_HPP
#define AVBD_DIRECT_SPACE_STATE3D_HPP

#include <vector>

#include <godot_cpp/classes/physics_direct_space_state3d_extension.hpp>
#include <godot_cpp/variant/rid.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDPhysicsServer3D;

class AVBDDirectSpaceState3D : public PhysicsDirectSpaceState3DExtension {
    GDCLASS(AVBDDirectSpaceState3D, PhysicsDirectSpaceState3DExtension)

    AVBDPhysicsServer3D *server = nullptr;
    RID space;

    // One body or area shape a query may hit, resolved from RID to a pure avbd::Shape.
    struct QueryTarget {
        avbd::Shape shape;
        RID rid;
        uint64_t collider_id = 0;
        int shape_index = 0;
        Vector3 linear_velocity;
        bool is_area = false;
    };
    std::vector<QueryTarget> collect_targets(const avbd::float3 &p_min, const avbd::float3 &p_max,
            uint32_t p_mask, bool p_bodies, bool p_areas, bool p_pick_ray = false) const;

protected:
    static void _bind_methods() {}

public:
    AVBDDirectSpaceState3D() = default;
    ~AVBDDirectSpaceState3D() override = default;

    void bind(AVBDPhysicsServer3D *p_server, const RID &p_space);

    bool _intersect_ray(const Vector3 &p_from, const Vector3 &p_to, uint32_t p_collision_mask,
            bool p_collide_with_bodies, bool p_collide_with_areas, bool p_hit_from_inside,
            bool p_hit_back_faces, bool p_pick_ray, PhysicsServer3DExtensionRayResult *r_result) override;
    int32_t _intersect_point(const Vector3 &p_position, uint32_t p_collision_mask,
            bool p_collide_with_bodies, bool p_collide_with_areas,
            PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) override;
    int32_t _intersect_shape(const RID &p_shape_rid, const Transform3D &p_transform, const Vector3 &p_motion,
            float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies, bool p_collide_with_areas,
            PhysicsServer3DExtensionShapeResult *r_results, int32_t p_max_results) override;
    bool _cast_motion(const RID &p_shape_rid, const Transform3D &p_transform, const Vector3 &p_motion,
            float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies, bool p_collide_with_areas,
            float *r_closest_safe, float *r_closest_unsafe,
            PhysicsServer3DExtensionShapeRestInfo *r_info) override;
    bool _collide_shape(const RID &p_shape_rid, const Transform3D &p_transform, const Vector3 &p_motion,
            float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies, bool p_collide_with_areas,
            void *r_results, int32_t p_max_results, int32_t *r_result_count) override;
    bool _rest_info(const RID &p_shape_rid, const Transform3D &p_transform, const Vector3 &p_motion,
            float p_margin, uint32_t p_collision_mask, bool p_collide_with_bodies, bool p_collide_with_areas,
            PhysicsServer3DExtensionShapeRestInfo *r_rest_info) override;
    Vector3 _get_closest_point_to_object_volume(const RID &p_object, const Vector3 &p_point) const override;
};

} // namespace godot

#endif // AVBD_DIRECT_SPACE_STATE3D_HPP
