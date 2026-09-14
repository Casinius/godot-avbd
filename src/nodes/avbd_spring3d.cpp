/*
 * AVBDSpring3D - linear spring between two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_spring3d.hpp"

#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_rigid_body3d.hpp"
#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

void AVBDSpring3D::set_anchor_a(const Vector3 &p_anchor) {
    anchor_a = p_anchor;
}

Vector3 AVBDSpring3D::get_anchor_a() const {
    return anchor_a;
}

void AVBDSpring3D::set_anchor_b(const Vector3 &p_anchor) {
    anchor_b = p_anchor;
}

Vector3 AVBDSpring3D::get_anchor_b() const {
    return anchor_b;
}

void AVBDSpring3D::set_stiffness(double p_stiffness) {
    stiffness = p_stiffness;
}

double AVBDSpring3D::get_stiffness() const {
    return stiffness;
}

void AVBDSpring3D::set_rest_length(double p_length) {
    rest_length = p_length;
}

double AVBDSpring3D::get_rest_length() const {
    return rest_length;
}

// The world transform only matters for anchors given in world space; both of these use
// body-local anchors, so it is unused.
avbd::Force *AVBDSpring3D::create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies, const Transform3D &) const {
    AVBDRigidBody3D *body_a = get_body_a();
    AVBDRigidBody3D *body_b = get_body_b();
    if (body_a == nullptr || body_b == nullptr) {
        return nullptr;
    }
    const auto found_a = p_bodies.find(body_a);
    const auto found_b = p_bodies.find(body_b);
    if (found_a == p_bodies.end() || found_b == p_bodies.end()) {
        return nullptr;
    }

    return new avbd::Spring(&p_solver, found_a->second, found_b->second, to_sim(anchor_a), to_sim(anchor_b),
            (float)stiffness, (float)rest_length);
}

void AVBDSpring3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_anchor_a", "anchor"), &AVBDSpring3D::set_anchor_a);
    ClassDB::bind_method(D_METHOD("get_anchor_a"), &AVBDSpring3D::get_anchor_a);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "anchor_a", PROPERTY_HINT_NONE, "suffix:m"), "set_anchor_a",
            "get_anchor_a");

    ClassDB::bind_method(D_METHOD("set_anchor_b", "anchor"), &AVBDSpring3D::set_anchor_b);
    ClassDB::bind_method(D_METHOD("get_anchor_b"), &AVBDSpring3D::get_anchor_b);
    ADD_PROPERTY(PropertyInfo(Variant::VECTOR3, "anchor_b", PROPERTY_HINT_NONE, "suffix:m"), "set_anchor_b",
            "get_anchor_b");

    ClassDB::bind_method(D_METHOD("set_stiffness", "stiffness"), &AVBDSpring3D::set_stiffness);
    ClassDB::bind_method(D_METHOD("get_stiffness"), &AVBDSpring3D::get_stiffness);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "stiffness", PROPERTY_HINT_RANGE, "0,1000000,1,or_greater"), "set_stiffness",
            "get_stiffness");

    ClassDB::bind_method(D_METHOD("set_rest_length", "length"), &AVBDSpring3D::set_rest_length);
    ClassDB::bind_method(D_METHOD("get_rest_length"), &AVBDSpring3D::get_rest_length);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "rest_length", PROPERTY_HINT_RANGE, "-1,1000,0.01,or_greater"),
            "set_rest_length", "get_rest_length");
}
