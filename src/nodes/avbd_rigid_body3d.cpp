/*
 * AVBDRigidBody3D - a box-shaped rigid body simulated by AVBD.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_rigid_body3d.hpp"

#include <godot_cpp/classes/box_mesh.hpp>
#include <godot_cpp/classes/standard_material3d.hpp>
#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_world3d.hpp"
#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

void AVBDRigidBody3D::set_size(const Vector3 &p_size) {
    size = p_size;
    _update_shape_visual();
}

Vector3 AVBDRigidBody3D::get_size() const {
    return size;
}

void AVBDRigidBody3D::set_density(double p_density) {
    density = p_density;
}

double AVBDRigidBody3D::get_density() const {
    return density;
}

void AVBDRigidBody3D::set_friction(double p_friction) {
    friction = p_friction;
}

double AVBDRigidBody3D::get_friction() const {
    return friction;
}

void AVBDRigidBody3D::set_static_body(bool p_static) {
    static_body = p_static;
}

bool AVBDRigidBody3D::is_static_body() const {
    return static_body;
}

void AVBDRigidBody3D::set_initial_velocity(const Vector3 &p_velocity) {
    initial_velocity = p_velocity;
}

Vector3 AVBDRigidBody3D::get_initial_velocity() const {
    return initial_velocity;
}

void AVBDRigidBody3D::set_initial_angular_velocity(const Vector3 &p_velocity) {
    initial_angular_velocity = p_velocity;
}

Vector3 AVBDRigidBody3D::get_initial_angular_velocity() const {
    return initial_angular_velocity;
}

double AVBDRigidBody3D::get_mass() const {
    if (rigid != nullptr) {
        return rigid->mass;
    }
    return static_body ? 0.0 : size.x * size.y * size.z * density;
}

Vector3 AVBDRigidBody3D::get_linear_velocity() const {
    if (rigid == nullptr) {
        return initial_velocity;
    }
    return direction_to_godot(world->get_global_transform(), rigid->velocityLin);
}

void AVBDRigidBody3D::set_linear_velocity(const Vector3 &p_velocity) {
    if (rigid == nullptr) {
        initial_velocity = p_velocity;
        return;
    }
    rigid->velocityLin = direction_to_sim(world->get_global_transform(), p_velocity);
}

// The solver integrates angular velocity in the world frame (AVBD composes
// orientation updates on the left of the quaternion), so both directions of this
// accessor are plain axis changes.
Vector3 AVBDRigidBody3D::get_angular_velocity() const {
    if (rigid == nullptr) {
        return initial_angular_velocity;
    }
    return direction_to_godot(world->get_global_transform(), rigid->velocityAng);
}

void AVBDRigidBody3D::set_angular_velocity(const Vector3 &p_velocity) {
    if (rigid == nullptr) {
        initial_angular_velocity = p_velocity;
        return;
    }
    rigid->velocityAng = direction_to_sim(world->get_global_transform(), p_velocity);
}

void AVBDRigidBody3D::apply_impulse(const Vector3 &p_impulse, const Vector3 &p_position_offset) {
    if (rigid == nullptr || rigid->mass <= 0.0f) {
        return;
    }
    const Transform3D world_global = world->get_global_transform();
    const avbd::float3 impulse = direction_to_sim(world_global, p_impulse);
    const avbd::float3 offset = direction_to_sim(world_global, p_position_offset);

    rigid->velocityLin += impulse / rigid->mass;

    // World-frame inverse inertia: I^-1 = R * Ilocal^-1 * R^T.
    // (Qualified: Node3D::rotate would otherwise hide the free function.)
    const avbd::float3 torque = avbd::cross(offset, impulse);
    const avbd::float3 local_torque = avbd::rotate(avbd::conjugate(rigid->positionAng), torque);
    const avbd::float3 local_delta{local_torque.x / rigid->moment.x, local_torque.y / rigid->moment.y,
            local_torque.z / rigid->moment.z};
    rigid->velocityAng += avbd::rotate(rigid->positionAng, local_delta);
}

void AVBDRigidBody3D::teleport(const Vector3 &p_position, const Quaternion &p_rotation) {
    const Transform3D pose(Basis(p_rotation), p_position);
    if (rigid == nullptr || static_body) {
        // Static bodies are driven by their node transform.
        set_global_transform(pose);
        return;
    }

    const Transform3D world_global = world->get_global_transform();
    const Quaternion world_rotation = world_global.basis.get_rotation_quaternion();
    rigid->positionLin = point_to_sim(world_global, p_position);
    rigid->positionAng = to_sim(world_rotation.inverse() * p_rotation);
    rigid->velocityLin = avbd::float3{0, 0, 0};
    rigid->velocityAng = avbd::float3{0, 0, 0};
    set_global_transform(pose);
}

void AVBDRigidBody3D::set_visualize_shape(bool p_visualize) {
    visualize_shape = p_visualize;
    _update_shape_visual();
}

bool AVBDRigidBody3D::is_visualize_shape() const {
    return visualize_shape;
}

void AVBDRigidBody3D::set_shape_color(const Color &p_color) {
    shape_color = p_color;
    _update_shape_visual();
}

Color AVBDRigidBody3D::get_shape_color() const {
    return shape_color;
}

void AVBDRigidBody3D::_ready() {
    _update_shape_visual();
}

void AVBDRigidBody3D::_update_shape_visual() {
    if (!is_inside_tree()) {
        return;
    }
    if (!visualize_shape) {
        if (shape_visual != nullptr) {
            shape_visual->queue_free();
            shape_visual = nullptr;
        }
        return;
    }
    if (shape_visual == nullptr) {
        shape_visual = memnew(MeshInstance3D);
        shape_visual->set_name("AVBDShape");
        add_child(shape_visual);
    }
    Ref<BoxMesh> mesh = shape_visual->get_mesh();
    if (mesh.is_null()) {
        mesh.instantiate();
        shape_visual->set_mesh(mesh);
    }
    mesh->set_size(size);
    Ref<StandardMaterial3D> material = mesh->get_material();
    if (material.is_null()) {
        material.instantiate();
        mesh->set_material(material);
    }
    material->set_albedo(shape_color);
}

void AVBDRigidBody3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_size", "size"), &AVBDRigidBody3D::set_size);
    ClassDB::bind_method(D_METHOD("get_size"), &AVBDRigidBody3D::get_size);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "size", PROPERTY_HINT_NONE, "suffix:m"), "set_size", "get_size");

    ClassDB::bind_method(D_METHOD("set_density", "density"), &AVBDRigidBody3D::set_density);
    ClassDB::bind_method(D_METHOD("get_density"), &AVBDRigidBody3D::get_density);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "density", PROPERTY_HINT_RANGE, "0,10000,0.01,or_greater"),
            "set_density", "get_density");

    ClassDB::bind_method(D_METHOD("set_friction", "friction"), &AVBDRigidBody3D::set_friction);
    ClassDB::bind_method(D_METHOD("get_friction"), &AVBDRigidBody3D::get_friction);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "friction", PROPERTY_HINT_RANGE, "0,1,0.01,or_greater"), "set_friction",
            "get_friction");

    ClassDB::bind_method(D_METHOD("set_static_body", "static_body"), &AVBDRigidBody3D::set_static_body);
    ClassDB::bind_method(D_METHOD("is_static_body"), &AVBDRigidBody3D::is_static_body);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "static_body"), "set_static_body", "is_static_body");

    ClassDB::bind_method(D_METHOD("set_initial_velocity", "velocity"), &AVBDRigidBody3D::set_initial_velocity);
    ClassDB::bind_method(D_METHOD("get_initial_velocity"), &AVBDRigidBody3D::get_initial_velocity);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "initial_velocity", PROPERTY_HINT_NONE, "suffix:m/s"),
            "set_initial_velocity", "get_initial_velocity");

    ClassDB::bind_method(D_METHOD("set_initial_angular_velocity", "velocity"),
            &AVBDRigidBody3D::set_initial_angular_velocity);
    ClassDB::bind_method(D_METHOD("get_initial_angular_velocity"), &AVBDRigidBody3D::get_initial_angular_velocity);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "initial_angular_velocity", PROPERTY_HINT_NONE, "suffix:rad/s"),
            "set_initial_angular_velocity", "get_initial_angular_velocity");

    ClassDB::bind_method(D_METHOD("get_mass"), &AVBDRigidBody3D::get_mass);

    ClassDB::bind_method(D_METHOD("set_linear_velocity", "velocity"), &AVBDRigidBody3D::set_linear_velocity);
    ClassDB::bind_method(D_METHOD("get_linear_velocity"), &AVBDRigidBody3D::get_linear_velocity);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "linear_velocity", PROPERTY_HINT_NONE, "suffix:m/s"),
            "set_linear_velocity", "get_linear_velocity");

    ClassDB::bind_method(D_METHOD("set_angular_velocity", "velocity"), &AVBDRigidBody3D::set_angular_velocity);
    ClassDB::bind_method(D_METHOD("get_angular_velocity"), &AVBDRigidBody3D::get_angular_velocity);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "angular_velocity", PROPERTY_HINT_NONE, "suffix:rad/s"),
            "set_angular_velocity", "get_angular_velocity");

    ClassDB::bind_method(D_METHOD("apply_impulse", "impulse", "position_offset"),
            &AVBDRigidBody3D::apply_impulse, DEFVAL(Vector3()));
    ClassDB::bind_method(D_METHOD("teleport", "position", "rotation"), &AVBDRigidBody3D::teleport,
            DEFVAL(Quaternion()));

    ClassDB::bind_method(D_METHOD("set_visualize_shape", "visualize"), &AVBDRigidBody3D::set_visualize_shape);
    ClassDB::bind_method(D_METHOD("is_visualize_shape"), &AVBDRigidBody3D::is_visualize_shape);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "visualize_shape"), "set_visualize_shape", "is_visualize_shape");

    ClassDB::bind_method(D_METHOD("set_shape_color", "color"), &AVBDRigidBody3D::set_shape_color);
    ClassDB::bind_method(D_METHOD("get_shape_color"), &AVBDRigidBody3D::get_shape_color);
    ADD_PROPERTY(PropertyInfo(Variant::COLOR, "shape_color"), "set_shape_color", "get_shape_color");
}
