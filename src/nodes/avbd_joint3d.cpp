/*
 * AVBDJoint3D - ball-and-socket constraint between two bodies or a body and space.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_joint3d.hpp"

#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_rigid_body3d.hpp"
#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

void AVBDJoint3D::set_anchor_a(const Vector3 &p_anchor) {
    anchor_a = p_anchor;
}

Vector3 AVBDJoint3D::get_anchor_a() const {
    return anchor_a;
}

void AVBDJoint3D::set_anchor_b(const Vector3 &p_anchor) {
    anchor_b = p_anchor;
}

Vector3 AVBDJoint3D::get_anchor_b() const {
    return anchor_b;
}

void AVBDJoint3D::set_linear_stiffness(double p_stiffness) {
    linear_stiffness = p_stiffness;
}

double AVBDJoint3D::get_linear_stiffness() const {
    return linear_stiffness;
}

void AVBDJoint3D::set_angular_stiffness(double p_stiffness) {
    angular_stiffness = p_stiffness;
}

double AVBDJoint3D::get_angular_stiffness() const {
    return angular_stiffness;
}

void AVBDJoint3D::set_fracture_force(double p_force) {
    fracture_force = p_force;
}

double AVBDJoint3D::get_fracture_force() const {
    return fracture_force;
}

bool AVBDJoint3D::is_broken() const {
    if (broken) {
        return true;
    }
    // The solver may have broken it during the last step, before the world noticed.
    return force != nullptr && force->isBroken();
}

void AVBDJoint3D::set_broken(bool p_broken) {
    broken = p_broken;
    if (p_broken && force != nullptr) {
        static_cast<avbd::Joint *>(force)->breakNow();
    }
}

void AVBDJoint3D::break_joint() {
    set_broken(true);
}

// A joint the solver deleted has broken: remember it, because the world will build a
// fresh, unbroken joint the next time it recreates the simulation.
void AVBDJoint3D::_on_force_removed() {
    broken = true;
}

avbd::Force *AVBDJoint3D::create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies,
        const Transform3D &p_world_global) const {
    AVBDRigidBody3D *body_a = get_body_a();
    AVBDRigidBody3D *body_b = get_body_b();
    if (body_b == nullptr) {
        return nullptr;
    }

    avbd::Rigid *rigid_a = nullptr;
    if (body_a != nullptr) {
        const auto found = p_bodies.find(body_a);
        if (found == p_bodies.end()) {
            return nullptr;
        }
        rigid_a = found->second;
    }
    const auto found_b = p_bodies.find(body_b);
    if (found_b == p_bodies.end()) {
        return nullptr;
    }

    const avbd::float3 r_a = (rigid_a != nullptr) ? to_sim(anchor_a) : point_to_sim(p_world_global, anchor_a);
    const avbd::float3 r_b = to_sim(anchor_b);

    return new avbd::Joint(&p_solver, rigid_a, found_b->second, r_a, r_b, stiffness_from_property(linear_stiffness),
            stiffness_from_property(angular_stiffness), stiffness_from_property(fracture_force));
}

void AVBDJoint3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_anchor_a", "anchor"), &AVBDJoint3D::set_anchor_a);
    ClassDB::bind_method(D_METHOD("get_anchor_a"), &AVBDJoint3D::get_anchor_a);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "anchor_a", PROPERTY_HINT_NONE, "suffix:m"), "set_anchor_a",
            "get_anchor_a");

    ClassDB::bind_method(D_METHOD("set_anchor_b", "anchor"), &AVBDJoint3D::set_anchor_b);
    ClassDB::bind_method(D_METHOD("get_anchor_b"), &AVBDJoint3D::get_anchor_b);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "anchor_b", PROPERTY_HINT_NONE, "suffix:m"), "set_anchor_b",
            "get_anchor_b");

    ClassDB::bind_method(D_METHOD("set_linear_stiffness", "stiffness"), &AVBDJoint3D::set_linear_stiffness);
    ClassDB::bind_method(D_METHOD("get_linear_stiffness"), &AVBDJoint3D::get_linear_stiffness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "linear_stiffness", PROPERTY_HINT_RANGE, "-1,1000000,1,or_greater"),
            "set_linear_stiffness", "get_linear_stiffness");

    ClassDB::bind_method(D_METHOD("set_angular_stiffness", "stiffness"), &AVBDJoint3D::set_angular_stiffness);
    ClassDB::bind_method(D_METHOD("get_angular_stiffness"), &AVBDJoint3D::get_angular_stiffness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "angular_stiffness", PROPERTY_HINT_RANGE, "-1,100000,1,or_greater"),
            "set_angular_stiffness", "get_angular_stiffness");

    ClassDB::bind_method(D_METHOD("set_fracture_force", "force"), &AVBDJoint3D::set_fracture_force);
    ClassDB::bind_method(D_METHOD("get_fracture_force"), &AVBDJoint3D::get_fracture_force);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "fracture_force", PROPERTY_HINT_RANGE, "-1,100000,1,or_greater"),
            "set_fracture_force", "get_fracture_force");

    ClassDB::bind_method(D_METHOD("is_broken"), &AVBDJoint3D::is_broken);
    ClassDB::bind_method(D_METHOD("set_broken", "broken"), &AVBDJoint3D::set_broken);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "broken"), "set_broken", "is_broken");

    ClassDB::bind_method(D_METHOD("break_joint"), &AVBDJoint3D::break_joint);
}
