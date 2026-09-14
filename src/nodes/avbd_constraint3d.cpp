/*
 * AVBDConstraint3D - shared base for every AVBD force that connects two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_constraint3d.hpp"

#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_rigid_body3d.hpp"

using namespace godot;

void AVBDConstraint3D::set_node_a(const NodePath &p_path) {
    node_a = p_path;
}

NodePath AVBDConstraint3D::get_node_a() const {
    return node_a;
}

void AVBDConstraint3D::set_node_b(const NodePath &p_path) {
    node_b = p_path;
}

NodePath AVBDConstraint3D::get_node_b() const {
    return node_b;
}

AVBDRigidBody3D *AVBDConstraint3D::get_body_a() const {
    if (node_a.is_empty()) {
        return nullptr;
    }
    return Object::cast_to<AVBDRigidBody3D>(get_node_or_null(node_a));
}

AVBDRigidBody3D *AVBDConstraint3D::get_body_b() const {
    if (node_b.is_empty()) {
        return nullptr;
    }
    return Object::cast_to<AVBDRigidBody3D>(get_node_or_null(node_b));
}

bool AVBDConstraint3D::is_simulated() const {
    return force != nullptr;
}

void AVBDConstraint3D::_bind_force(avbd::Force *p_force) {
    force = p_force;
    force_alive = p_force != nullptr;
    if (p_force != nullptr) {
        p_force->owner = this;
    }
}

// The solver deletes a force whose initialize() reported it inactive. The world walks
// the surviving forces once per tick, so a force that is no longer there was removed:
// for a joint that means it broke. Dropping the handle here keeps the node from ever
// dereferencing freed memory.
void AVBDConstraint3D::_resolve_force_removal() {
    if (force != nullptr && !force_alive) {
        force = nullptr;
        _on_force_removed();
    }
}

void AVBDConstraint3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_node_a", "path"), &AVBDConstraint3D::set_node_a);
    ClassDB::bind_method(D_METHOD("get_node_a"), &AVBDConstraint3D::get_node_a);
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "node_a", PROPERTY_HINT_NODE_PATH_VALID_TYPES, "AVBDRigidBody3D"),
            "set_node_a", "get_node_a");

    ClassDB::bind_method(D_METHOD("set_node_b", "path"), &AVBDConstraint3D::set_node_b);
    ClassDB::bind_method(D_METHOD("get_node_b"), &AVBDConstraint3D::get_node_b);
    ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "node_b", PROPERTY_HINT_NODE_PATH_VALID_TYPES, "AVBDRigidBody3D"),
            "set_node_b", "get_node_b");

    ClassDB::bind_method(D_METHOD("get_body_a"), &AVBDConstraint3D::get_body_a);
    ClassDB::bind_method(D_METHOD("get_body_b"), &AVBDConstraint3D::get_body_b);
    ClassDB::bind_method(D_METHOD("is_simulated"), &AVBDConstraint3D::is_simulated);
}
