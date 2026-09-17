/*
 * AVBDDirectBodyState3D implementation. See the header for the contract.
 *
 * Frame conversion: the solver is Z-up, Godot is Y-up; to_sim/to_godot in
 * src/nodes/godot_convert.hpp do the mapping for vectors, quaternions and
 * transforms. Everything the engine reads comes straight off the bound
 * avbd::Rigid; writes go straight back onto it.
 *
 * Contact data is filled by the server after each step (frame_contacts on the
 * BodyData); until the server binds a body, this object reports neutral values
 * rather than dereferencing a null Rigid.
 */

#include "server/avbd_direct_body_state3d.hpp"

#include "nodes/godot_convert.hpp"
#include "server/avbd_physics_server3d.hpp"

#include <godot_cpp/classes/physics_direct_space_state3d.hpp>
#include <godot_cpp/core/object.hpp>

using namespace godot;
using namespace avbd_godot;

void AVBDDirectBodyState3D::bind(AVBDPhysicsServer3D *p_server, avbd::Rigid *p_rigid, const RID &p_body_rid) {
    server = p_server;
    rigid = p_rigid;
    body_rid = p_body_rid;
}

Vector3 AVBDDirectBodyState3D::_get_total_gravity() const {
    // The solver applies gravity along its -Z with magnitude |gravity|; express that
    // back in Godot's frame (straight down).
    if (rigid == nullptr) {
        return Vector3();
    }
    return Vector3(0, -std::fabs(rigid->gravity), 0);
}

double AVBDDirectBodyState3D::_get_total_linear_damp() const {
    return 0.0; // damping is not simulated (documented in the README)
}

double AVBDDirectBodyState3D::_get_total_angular_damp() const {
    return 0.0;
}

Vector3 AVBDDirectBodyState3D::_get_center_of_mass() const {
    // Every solver body's origin is its centre of mass.
    return rigid != nullptr ? to_godot(rigid->positionLin) : Vector3();
}

Vector3 AVBDDirectBodyState3D::_get_center_of_mass_local() const {
    return Vector3();
}

Basis AVBDDirectBodyState3D::_get_principal_inertia_axes() const {
    // The moment vector is already diagonal in the body's local frame.
    return rigid != nullptr ? to_godot_basis(rigid->positionAng) : Basis();
}

double AVBDDirectBodyState3D::_get_inverse_mass() const {
    return rigid != nullptr && rigid->mass > 0.0f ? 1.0 / rigid->mass : 0.0;
}

Vector3 AVBDDirectBodyState3D::_get_inverse_inertia() const {
    if (rigid == nullptr || rigid->mass <= 0.0f) {
        return Vector3();
    }
    const avbd::float3 moment = rigid->moment;
    return to_godot(avbd::float3{1.0f / moment.x, 1.0f / moment.y, 1.0f / moment.z});
}

Basis AVBDDirectBodyState3D::_get_inverse_inertia_tensor() const {
    if (rigid == nullptr) {
        return Basis();
    }
    const Basis rot = to_godot_basis(rigid->positionAng);
    const Vector3 inv = _get_inverse_inertia();
    return rot.scaled(Vector3(inv.x, inv.y, inv.z)) * rot.transposed();
}

void AVBDDirectBodyState3D::_set_linear_velocity(const Vector3 &p_velocity) {
    if (rigid != nullptr) {
        rigid->velocityLin = to_sim(p_velocity);
    }
}

Vector3 AVBDDirectBodyState3D::_get_linear_velocity() const {
    return rigid != nullptr ? to_godot(rigid->velocityLin) : Vector3();
}

void AVBDDirectBodyState3D::_set_angular_velocity(const Vector3 &p_velocity) {
    if (rigid != nullptr) {
        rigid->velocityAng = to_sim(p_velocity);
    }
}

Vector3 AVBDDirectBodyState3D::_get_angular_velocity() const {
    return rigid != nullptr ? to_godot(rigid->velocityAng) : Vector3();
}

void AVBDDirectBodyState3D::_set_transform(const Transform3D &p_transform) {
    if (rigid != nullptr) {
        rigid->positionLin = to_sim(p_transform.origin);
        rigid->positionAng = to_sim(p_transform.basis.get_rotation_quaternion());
    }
}

Transform3D AVBDDirectBodyState3D::_get_transform() const {
    if (rigid == nullptr) {
        return Transform3D();
    }
    return sim_transform(rigid->positionLin, rigid->positionAng);
}

Vector3 AVBDDirectBodyState3D::_get_velocity_at_local_position(const Vector3 &p_local_position) const {
    if (rigid == nullptr) {
        return Vector3();
    }
    const avbd::float3 r = to_sim(p_local_position);
    const avbd::float3 velocity = rigid->velocityLin + avbd::cross(rigid->velocityAng, r);
    return to_godot(velocity);
}

void AVBDDirectBodyState3D::_apply_central_impulse(const Vector3 &p_impulse) {
    if (rigid != nullptr && rigid->mass > 0.0f) {
        rigid->velocityLin += to_sim(p_impulse) / rigid->mass;
    }
}

void AVBDDirectBodyState3D::_apply_impulse(const Vector3 &p_impulse, const Vector3 &p_position) {
    if (rigid == nullptr || rigid->mass <= 0.0f) {
        return;
    }
    const avbd::float3 impulse = to_sim(p_impulse);
    const avbd::float3 r = to_sim(p_position);
    rigid->velocityLin += impulse / rigid->mass;
    const avbd::float3 torque = avbd::cross(r, impulse);
    const avbd::float3 local = avbd::rotate(avbd::conjugate(rigid->positionAng), torque);
    const avbd::float3 delta{local.x / rigid->moment.x, local.y / rigid->moment.y, local.z / rigid->moment.z};
    rigid->velocityAng += avbd::rotate(rigid->positionAng, delta);
}

void AVBDDirectBodyState3D::_apply_torque_impulse(const Vector3 &p_impulse) {
    if (rigid == nullptr || rigid->mass <= 0.0f) {
        return;
    }
    const avbd::float3 local = avbd::rotate(avbd::conjugate(rigid->positionAng), to_sim(p_impulse));
    const avbd::float3 delta{local.x / rigid->moment.x, local.y / rigid->moment.y, local.z / rigid->moment.z};
    rigid->velocityAng += avbd::rotate(rigid->positionAng, delta);
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
    if (rigid != nullptr) {
        rigid->sleeping = p_enabled;
        if (p_enabled) {
            rigid->velocityLin = avbd::float3{0, 0, 0};
            rigid->velocityAng = avbd::float3{0, 0, 0};
        }
    }
}

bool AVBDDirectBodyState3D::_is_sleeping() const {
    return rigid != nullptr ? rigid->sleeping : false;
}

int32_t AVBDDirectBodyState3D::_get_contact_count() const {
    return server != nullptr ? server->body_state_contact_count(body_rid) : 0;
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_position(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_position(body_rid, p_contact_idx) : Vector3();
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_normal(int32_t) const {
    return Vector3(); // per-contact normals are not tracked this round
}

Vector3 AVBDDirectBodyState3D::_get_contact_impulse(int32_t) const {
    return Vector3(); // impulses are not accumulated per contact
}

int32_t AVBDDirectBodyState3D::_get_contact_local_shape(int32_t) const {
    return 0; // one shape per solver body
}

Vector3 AVBDDirectBodyState3D::_get_contact_local_velocity_at_position(int32_t p_contact_idx) const {
    if (rigid == nullptr) {
        return Vector3();
    }
    const Vector3 world = server != nullptr ? server->body_state_contact_position(body_rid, p_contact_idx) : Vector3();
    const avbd::float3 r = to_sim(world - to_godot(rigid->positionLin));
    return to_godot(rigid->velocityLin + avbd::cross(rigid->velocityAng, r));
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

int32_t AVBDDirectBodyState3D::_get_contact_collider_shape(int32_t) const {
    return 0; // one shape per solver body
}

Vector3 AVBDDirectBodyState3D::_get_contact_collider_velocity_at_position(int32_t p_contact_idx) const {
    return server != nullptr ? server->body_state_contact_collider_velocity(body_rid, p_contact_idx) : Vector3();
}

double AVBDDirectBodyState3D::_get_step() const {
    return server != nullptr ? server->body_state_step() : 0.0;
}

void AVBDDirectBodyState3D::_integrate_forces() {
    // The solver already integrated this step; nothing to do.
}

PhysicsDirectSpaceState3D *AVBDDirectBodyState3D::_get_space_state() {
    return server != nullptr ? server->body_state_space_state(body_rid) : nullptr;
}
