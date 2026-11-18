/*
 * AVBDDirectBodyState3D implementation. See the header for the contract.
 *
 * Everything reads the server's BodyData record and the body's ModernSpace. Both are
 * Y-up, so vectors and quaternions pass through unchanged. Contact data is filled by
 * the server after each step; before a body is bound (or when the RID is unknown) the
 * getters answer neutrally instead of dereferencing anything.
 */

#include "server/avbd_direct_body_state3d.hpp"

#include "server/avbd_physics_server3d.hpp"
#include "server/modern_convert.hpp"

#include <godot_cpp/classes/physics_direct_space_state3d.hpp>
#include <godot_cpp/core/object.hpp>

using namespace godot;

void AVBDDirectBodyState3D::bind(AVBDPhysicsServer3D *p_server, const RID &p_body_rid) {
    server = p_server;
    body_rid = p_body_rid;
}

uint32_t AVBDDirectBodyState3D::_get_collision_layer() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? static_cast<uint32_t>(body->collision_layer) : 1;
}

void AVBDDirectBodyState3D::_set_collision_layer(uint32_t p_layer) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->collision_layer = p_layer;
    // The filter tables live in ModernSpace; refresh them before the next step.
    if (body->space.is_valid()) {
        server->space_rebuild_pending[server->id_of(body->space)] = true;
    }
}

uint32_t AVBDDirectBodyState3D::_get_collision_mask() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? static_cast<uint32_t>(body->collision_mask) : 1;
}

void AVBDDirectBodyState3D::_set_collision_mask(uint32_t p_mask) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->collision_mask = p_mask;
    if (body->space.is_valid()) {
        server->space_rebuild_pending[server->id_of(body->space)] = true;
    }
}

Vector3 AVBDDirectBodyState3D::_get_total_gravity() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || !body->space.is_valid()) {
        return Vector3();
    }
    const AVBDPhysicsServer3D::SpaceData *space = server->find_space(body->space);
    if (space == nullptr) {
        return Vector3();
    }
    return to_godot(space->gravity * static_cast<Scalar_t>(body->gravity_scale));
}

float AVBDDirectBodyState3D::_get_total_linear_damp() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? static_cast<float>(body->linear_damp) : 0.0f;
}

float AVBDDirectBodyState3D::_get_total_angular_damp() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? static_cast<float>(body->angular_damp) : 0.0f;
}

Vector3 AVBDDirectBodyState3D::_get_center_of_mass() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Vector3();
    }
    const ModernSpace *modern = server->modern_for(*body);
    if (modern == nullptr) {
        return body->transform.origin;
    }
    return to_godot(modern->state().pose.position(static_cast<std::size_t>(body->modern_index)));
}

Vector3 AVBDDirectBodyState3D::_get_center_of_mass_local() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Vector3();
    }
    const ModernSpace *modern = server->modern_for(*body);
    if (modern == nullptr) {
        return Vector3();
    }
    return to_godot(modern->com_local(static_cast<std::size_t>(body->modern_index)));
}

Basis AVBDDirectBodyState3D::_get_principal_inertia_axes() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? body->transform.basis.orthonormalized() : Basis();
}

float AVBDDirectBodyState3D::_get_inverse_mass() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return 0.0f;
    }
    const ModernSpace *modern = server->modern_for(*body);
    if (modern == nullptr) {
        return 0.0f;
    }
    return modern->state().inertial.inv_mass(static_cast<std::size_t>(body->modern_index));
}

Vector3 AVBDDirectBodyState3D::_get_inverse_inertia() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Vector3();
    }
    const ModernSpace *modern = server->modern_for(*body);
    if (modern == nullptr) {
        return Vector3();
    }
    const std::size_t i = static_cast<std::size_t>(body->modern_index);
    const Mat3 inv = modern->state().inertial.local_inv_inertia(i);
    return Vector3(inv[0][0], inv[1][1], inv[2][2]);
}

Basis AVBDDirectBodyState3D::_get_inverse_inertia_tensor() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Basis();
    }
    const ModernSpace *modern = server->modern_for(*body);
    if (modern == nullptr) {
        return Basis();
    }
    const std::size_t i = static_cast<std::size_t>(body->modern_index);
    // World inverse inertia, refreshed against the current pose.
    const Mat3 inv_world = modern->state().inertial.world_inv_inertia(i);
    Basis out;
    for (int col = 0; col < 3; ++col) {
        out.set_column(col, Vector3(inv_world[col][0], inv_world[col][1], inv_world[col][2]));
    }
    return out;
}

void AVBDDirectBodyState3D::_set_linear_velocity(const Vector3 &p_velocity) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->linear_velocity = p_velocity;
    if (ModernSpace *modern = server->modern_for(*body)) {
        modern->set_linear_velocity(static_cast<std::size_t>(body->modern_index), to_modern(p_velocity));
    }
}

Vector3 AVBDDirectBodyState3D::_get_linear_velocity() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? body->linear_velocity : Vector3();
}

void AVBDDirectBodyState3D::_set_angular_velocity(const Vector3 &p_velocity) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->angular_velocity = p_velocity;
    if (ModernSpace *modern = server->modern_for(*body)) {
        modern->set_angular_velocity(static_cast<std::size_t>(body->modern_index), to_modern(p_velocity));
    }
}

Vector3 AVBDDirectBodyState3D::_get_angular_velocity() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? body->angular_velocity : Vector3();
}

void AVBDDirectBodyState3D::_set_transform(const Transform3D &p_transform) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->transform = p_transform;
    if (ModernSpace *modern = server->modern_for(*body)) {
        modern->set_body_transform(static_cast<std::size_t>(body->modern_index),
                to_modern(p_transform.origin),
                to_modern(p_transform.basis.get_rotation_quaternion()));
    }
}

Transform3D AVBDDirectBodyState3D::_get_transform() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? body->transform : Transform3D();
}

Vector3 AVBDDirectBodyState3D::_get_velocity_at_local_position(const Vector3 &p_local_position) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Vector3();
    }
    return body->linear_velocity + body->angular_velocity.cross(p_local_position);
}

void AVBDDirectBodyState3D::_apply_central_impulse(const Vector3 &p_impulse) {
    if (server != nullptr) {
        server->_body_apply_central_impulse(body_rid, p_impulse);
    }
}

void AVBDDirectBodyState3D::_apply_impulse(const Vector3 &p_impulse, const Vector3 &p_position) {
    if (server != nullptr) {
        server->_body_apply_impulse(body_rid, p_impulse, p_position);
    }
}

void AVBDDirectBodyState3D::_apply_torque_impulse(const Vector3 &p_impulse) {
    if (server != nullptr) {
        server->_body_apply_torque_impulse(body_rid, p_impulse);
    }
}

void AVBDDirectBodyState3D::_apply_central_force(const Vector3 &p_force) {
    if (server != nullptr) {
        server->body_state_add_constant_force(body_rid, p_force, Vector3());
    }
}

void AVBDDirectBodyState3D::_apply_force(const Vector3 &p_force, const Vector3 &p_position) {
    if (server != nullptr) {
        server->body_state_add_constant_force(body_rid, p_force, p_position);
    }
}

void AVBDDirectBodyState3D::_apply_torque(const Vector3 &p_torque) {
    if (server != nullptr) {
        server->body_state_add_constant_torque(body_rid, p_torque);
    }
}

void AVBDDirectBodyState3D::_add_constant_central_force(const Vector3 &p_force) {
    if (server != nullptr) {
        server->body_state_add_constant_force(body_rid, p_force, Vector3());
    }
}

void AVBDDirectBodyState3D::_add_constant_force(const Vector3 &p_force, const Vector3 &p_position) {
    if (server != nullptr) {
        server->body_state_add_constant_force(body_rid, p_force, p_position);
    }
}

void AVBDDirectBodyState3D::_add_constant_torque(const Vector3 &p_torque) {
    if (server != nullptr) {
        server->body_state_add_constant_torque(body_rid, p_torque);
    }
}

void AVBDDirectBodyState3D::_set_constant_force(const Vector3 &p_force) {
    if (server != nullptr) {
        server->body_state_set_constant_force(body_rid, p_force);
    }
}

Vector3 AVBDDirectBodyState3D::_get_constant_force() const {
    return server != nullptr ? server->body_state_get_constant_force(body_rid) : Vector3();
}

void AVBDDirectBodyState3D::_set_constant_torque(const Vector3 &p_torque) {
    if (server != nullptr) {
        server->body_state_set_constant_torque(body_rid, p_torque);
    }
}

Vector3 AVBDDirectBodyState3D::_get_constant_torque() const {
    return server != nullptr ? server->body_state_get_constant_torque(body_rid) : Vector3();
}

void AVBDDirectBodyState3D::_set_sleep_state(bool p_enabled) {
    AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return;
    }
    body->sleeping = p_enabled;
    if (p_enabled) {
        body->linear_velocity = Vector3();
        body->angular_velocity = Vector3();
        if (ModernSpace *modern = server->modern_for(*body)) {
            modern->set_linear_velocity(static_cast<std::size_t>(body->modern_index), Vec3{0});
            modern->set_angular_velocity(static_cast<std::size_t>(body->modern_index), Vec3{0});
        }
    }
}

bool AVBDDirectBodyState3D::_is_sleeping() const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    return body != nullptr ? body->sleeping : false;
}

int32_t AVBDDirectBodyState3D::_get_contact_count() const {
    return server != nullptr ? server->body_state_contact_count(body_rid) : 0;
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_position(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_position(body_rid, p_contact_idx) : Vector3();
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_normal(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || p_contact_idx < 0 ||
            static_cast<std::size_t>(p_contact_idx) >= body->frame_contacts.size()) {
        return Vector3();
    }
    return body->frame_contacts[p_contact_idx].normal;
}

Vector3 AVBDDirectBodyState3D::_get_contact_impulse(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || p_contact_idx < 0 ||
            static_cast<std::size_t>(p_contact_idx) >= body->frame_contacts.size()) {
        return Vector3();
    }
    return body->frame_contacts[p_contact_idx].impulse;
}

int32_t AVBDDirectBodyState3D::_get_contact_local_shape(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || p_contact_idx < 0 ||
            static_cast<std::size_t>(p_contact_idx) >= body->frame_contacts.size()) {
        return 0;
    }
    return body->frame_contacts[p_contact_idx].local_shape;
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_velocity_at_position(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr) {
        return Vector3();
    }
    const Vector3 world = server->body_state_contact_position(body_rid, p_contact_idx);
    const Vector3 r = world - body->transform.origin;
    return body->linear_velocity + body->angular_velocity.cross(r);
}

RID AVBDDirectBodyState3D::_get_contact_collider(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_collider(body_rid, p_contact_idx) : RID();
}

Vector3 AVBDDirectBodyState3D::_get_contact_collider_position(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_position(body_rid, p_contact_idx) : Vector3();
}

uint64_t AVBDDirectBodyState3D::_get_contact_collider_id(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_collider_id(body_rid, p_contact_idx) : 0;
}

Object *AVBDDirectBodyState3D::_get_contact_collider_object(int32_t p_contact_idx) const {
    if (server == nullptr) {
        return nullptr;
    }
    const uint64_t id = server->body_state_contact_collider_id(body_rid, p_contact_idx);
    return id != 0 ? ObjectDB::get_instance(id) : nullptr;
}

int32_t AVBDDirectBodyState3D::_get_contact_collider_shape(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || p_contact_idx < 0 ||
            static_cast<std::size_t>(p_contact_idx) >= body->frame_contacts.size()) {
        return 0;
    }
    return body->frame_contacts[p_contact_idx].collider_shape;
}

Vector3 AVBDDirectBodyState3D::_get_contact_collider_velocity_at_position(int32_t p_contact_idx) const {
    const AVBDPhysicsServer3D::BodyData *body = server != nullptr ? server->find_body(body_rid) : nullptr;
    if (body == nullptr || p_contact_idx < 0 ||
            static_cast<std::size_t>(p_contact_idx) >= body->frame_contacts.size()) {
        return Vector3();
    }
    const auto other = server->find_body(server->body_state_contact_collider(body_rid, p_contact_idx));
    if (other == nullptr) {
        return Vector3();
    }
    const Vector3 world = server->body_state_contact_position(body_rid, p_contact_idx);
    const Vector3 r = world - other->transform.origin;
    return other->linear_velocity + other->angular_velocity.cross(r);
}

float AVBDDirectBodyState3D::_get_step() const {
    return server != nullptr ? static_cast<float>(server->body_state_step()) : 0.0f;
}

void AVBDDirectBodyState3D::_integrate_forces() {
    // The solver already integrated this step; nothing to do.
}

PhysicsDirectSpaceState3D *AVBDDirectBodyState3D::_get_space_state() {
    return server != nullptr ? server->body_state_space_state(body_rid) : nullptr;
}
