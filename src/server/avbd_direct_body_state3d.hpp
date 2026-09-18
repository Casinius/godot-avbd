/*
 * AVBDDirectBodyState3D - the PhysicsDirectBodyState3D the server hands the engine.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * Godot's body state sync callback contract (verified against Godot 4.7,
 * scene/3d/physics/rigid_body_3d.cpp `_body_state_changed` and
 * modules/godot_physics_3d/godot_body_3d.cpp `call_queries`) is a *single argument*:
 * a PhysicsDirectBodyState3D. The same instance is also the first argument of the
 * force integration callback, followed by the stored userdata.
 *
 * One instance per space is cached by the server and re-pointed at each body before
 * the callbacks run, so the engine never holds a stale view. Every getter reads the
 * avbd::Rigid live; setters write straight through to it, mirroring what the engine's
 * own `_integrate_forces` users expect.
 *
 * The converter lives in src/nodes/godot_convert.hpp and maps Godot's Y-up scene
 * space onto the solver's Z-up space.
 */

#ifndef AVBD_DIRECT_BODY_STATE3D_HPP
#define AVBD_DIRECT_BODY_STATE3D_HPP

#include <godot_cpp/classes/physics_direct_body_state3d_extension.hpp>
#include <godot_cpp/variant/basis.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

namespace godot {

class AVBDPhysicsServer3D;

class AVBDDirectBodyState3D : public PhysicsDirectBodyState3DExtension {
    GDCLASS(AVBDDirectBodyState3D, PhysicsDirectBodyState3DExtension)

    AVBDPhysicsServer3D *server = nullptr;
    avbd::Rigid *rigid = nullptr; // the body this instance currently views
    RID body_rid;

protected:
    static void _bind_methods() {}

public:
    AVBDDirectBodyState3D() = default;
    ~AVBDDirectBodyState3D() override = default;

    // Re-aim the cached instance at one body before its callbacks run.
    void bind(AVBDPhysicsServer3D *p_server, avbd::Rigid *p_rigid, const RID &p_body_rid);

    // --- reads/writes onto the bound body ------------------------------------
    Vector3 _get_total_gravity() const override;
    double _get_total_linear_damp() const override;
    double _get_total_angular_damp() const override;
    Vector3 _get_center_of_mass() const override;
    Vector3 _get_center_of_mass_local() const override;
    Basis _get_principal_inertia_axes() const override;
    double _get_inverse_mass() const override;
    Vector3 _get_inverse_inertia() const override;
    Basis _get_inverse_inertia_tensor() const override;
    void _set_linear_velocity(const Vector3 &p_velocity) override;
    Vector3 _get_linear_velocity() const override;
    void _set_angular_velocity(const Vector3 &p_velocity) override;
    Vector3 _get_angular_velocity() const override;
    void _set_transform(const Transform3D &p_transform) override;
    Transform3D _get_transform() const override;
    Vector3 _get_velocity_at_local_position(const Vector3 &p_local_position) const override;
    void _apply_central_impulse(const Vector3 &p_impulse) override;
    void _apply_impulse(const Vector3 &p_impulse, const Vector3 &p_position) override;
    void _apply_torque_impulse(const Vector3 &p_impulse) override;
    void _apply_central_force(const Vector3 &p_force) override;
    void _apply_force(const Vector3 &p_force, const Vector3 &p_position) override;
    void _apply_torque(const Vector3 &p_torque) override;
    void _add_constant_central_force(const Vector3 &p_force) override;
    void _add_constant_force(const Vector3 &p_force, const Vector3 &p_position) override;
    void _add_constant_torque(const Vector3 &p_torque) override;
    void _set_constant_force(const Vector3 &p_force) override;
    Vector3 _get_constant_force() const override;
    void _set_constant_torque(const Vector3 &p_torque) override;
    Vector3 _get_constant_torque() const override;
    void _set_sleep_state(bool p_enabled) override;
    bool _is_sleeping() const override;
    int32_t _get_contact_count() const override;
    Vector3 _get_contact_local_position(int32_t p_contact_idx) const override;
    Vector3 _get_contact_local_normal(int32_t p_contact_idx) const override;
    Vector3 _get_contact_impulse(int32_t p_contact_idx) const override;
    int32_t _get_contact_local_shape(int32_t p_contact_idx) const override;
    Vector3 _get_contact_local_velocity_at_position(int32_t p_contact_idx) const override;
    RID _get_contact_collider(int32_t p_contact_idx) const override;
    Vector3 _get_contact_collider_position(int32_t p_contact_idx) const override;
    uint64_t _get_contact_collider_id(int32_t p_contact_idx) const override;
    Object *_get_contact_collider_object(int32_t p_contact_idx) const override;
    int32_t _get_contact_collider_shape(int32_t p_contact_idx) const override;
    Vector3 _get_contact_collider_velocity_at_position(int32_t p_contact_idx) const override;
    double _get_step() const override;
    void _integrate_forces() override;
    PhysicsDirectSpaceState3D *_get_space_state() override;
};

} // namespace godot

#endif // AVBD_DIRECT_BODY_STATE3D_HPP
