/*
 * AVBDDirectSpaceState3D - the PhysicsDirectSpaceState3D behind space queries.
 *
 * Part of godot-avbd.
 *
 * The engine reaches this through World3D.direct_space_state (RayCast3D, intersect_*,
 * ShapeCast3D, ...). One instance per space, created lazily by the server.
 *
 * Queries walk the space's modern body/area shapes and answer with analytic
 * ray/point/AABB tests. This is the MVP query path: box/sphere/cylinder are exact for
 * ray and point tests, while shape-shape queries use conservative AABB overlap (the
 * same contract the swept motion query uses). Exclusion lists (`params.exclude`) are
 * honoured through PhysicsDirectSpaceState3DExtension::is_body_excluded_from_query();
 * `ray_pickable` is read for body and area hits.
 */

#ifndef AVBD_DIRECT_SPACE_STATE3D_HPP
#define AVBD_DIRECT_SPACE_STATE3D_HPP

#include <vector>

#include <godot_cpp/classes/physics_direct_space_state3d_extension.hpp>
#include <godot_cpp/variant/rid.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "collision/shape.hpp"
#include "phys_core/state.hpp"

namespace godot {

class AVBDPhysicsServer3D;

class AVBDDirectSpaceState3D : public PhysicsDirectSpaceState3DExtension {
    GDCLASS(AVBDDirectSpaceState3D, PhysicsDirectSpaceState3DExtension)

    AVBDPhysicsServer3D *server = nullptr;
    RID space;

    // One body or area shape a query may hit, resolved to a local shape plus its world pose.
    struct QueryTarget {
        collision::Shape shape;
        Transform3D transform;
        RID rid;
        uint64_t collider_id = 0;
        int32_t shape_index = 0;
        Vector3 linear_velocity;
        bool is_area = false;
    };
    std::vector<QueryTarget> collect_targets(const Vec3 &p_min, const Vec3 &p_max, uint32_t p_mask,
            bool p_bodies, bool p_areas, bool p_pick_ray = false) const;

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
