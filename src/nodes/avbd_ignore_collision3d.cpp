/*
 * AVBDIgnoreCollision3D - suppresses contacts between two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_ignore_collision3d.hpp"

#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_rigid_body3d.hpp"

using namespace godot;

// The world transform only matters for anchors given in world space; both of these use
// body-local anchors, so it is unused.
avbd::Force *AVBDIgnoreCollision3D::create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies, const Transform3D &) const {
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
    return new avbd::IgnoreCollision(&p_solver, found_a->second, found_b->second);
}

void AVBDIgnoreCollision3D::_bind_methods() {
}
