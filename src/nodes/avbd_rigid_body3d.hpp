/*
 * AVBDRigidBody3D - a box-shaped rigid body simulated by AVBD.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * The body's transform at the moment it joins the simulation defines its initial
 * pose; from then on the solver owns it and writes the result back every physics
 * tick. Use teleport() to move a dynamic body afterwards. Bodies flagged
 * `static_body` are immovable but still collide; their node transform stays
 * authoritative, so moving one in the editor or by script moves the collider.
 */

#ifndef AVBD_RIGID_BODY3D_HPP
#define AVBD_RIGID_BODY3D_HPP

#include <godot_cpp/classes/mesh_instance3d.hpp>
#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/variant/color.hpp>
#include <godot_cpp/variant/quaternion.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDWorld3D;

class AVBDRigidBody3D : public Node3D {
    GDCLASS(AVBDRigidBody3D, Node3D)

    Vector3 size = Vector3(1, 1, 1);
    double density = 1.0;
    double friction = 0.5;
    bool static_body = false;
    Vector3 initial_velocity;
    Vector3 initial_angular_velocity;
    bool visualize_shape = false;
    Color shape_color = Color(0.55, 0.62, 0.75);

    AVBDWorld3D *world = nullptr;
    avbd::Rigid *rigid = nullptr;
    MeshInstance3D *shape_visual = nullptr;

protected:
    static void _bind_methods();
    void _update_shape_visual();

public:
    AVBDRigidBody3D() = default;
    ~AVBDRigidBody3D() override = default;

    // Box extents in metres (full widths, Godot axes).
    void set_size(const Vector3 &p_size);
    Vector3 get_size() const;

    // Mass density in kg/m^3. Mass and inertia are derived from size and density;
    // a static body has infinite mass regardless of this value.
    void set_density(double p_density);
    double get_density() const;

    // Contact friction coefficient in [0, 1]. Two bodies in contact use the
    // geometric mean of their coefficients.
    void set_friction(double p_friction);
    double get_friction() const;

    // Immovable collider (ground, walls, ramps). Never integrated, never pushed.
    void set_static_body(bool p_static);
    bool is_static_body() const;

    // Velocity applied when the body enters the simulation.
    void set_initial_velocity(const Vector3 &p_velocity);
    Vector3 get_initial_velocity() const;
    void set_initial_angular_velocity(const Vector3 &p_velocity);
    Vector3 get_initial_angular_velocity() const;

    double get_mass() const;

    Vector3 get_linear_velocity() const;
    void set_linear_velocity(const Vector3 &p_velocity);

    // Angular velocity is reported and applied in world (Godot) axes.
    Vector3 get_angular_velocity() const;
    void set_angular_velocity(const Vector3 &p_velocity);

    // Impulse in world axes. `p_position_offset` is the point of application
    // relative to the body's centre, also in world axes; the linear and angular
    // contributions follow from mass and inertia.
    void apply_impulse(const Vector3 &p_impulse, const Vector3 &p_position_offset = Vector3());

    // Move the body to an absolute pose in Godot world space and stop it. Only
    // valid for dynamic bodies; static bodies follow their node transform.
    void teleport(const Vector3 &p_position, const Quaternion &p_rotation = Quaternion());

    // Debug/editor convenience: draw a box matching `size` as a child mesh. The
    // mesh is created at runtime only, so it never ends up in the scene file.
    void set_visualize_shape(bool p_visualize);
    bool is_visualize_shape() const;
    void set_shape_color(const Color &p_color);
    Color get_shape_color() const;

    void _ready() override;

    // --- simulation interface, used by AVBDWorld3D ---------------------------
    void _bind_simulation(AVBDWorld3D *p_world, avbd::Rigid *p_rigid) {
        world = p_world;
        rigid = p_rigid;
    }
    const avbd::Rigid *_get_rigid() const { return rigid; }
    bool _is_ready() const { return rigid != nullptr; }
};

} // namespace godot

#endif // AVBD_RIGID_BODY3D_HPP
