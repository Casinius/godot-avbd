/*
 * AVBDWorld3D - root of an AVBD simulation.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#include "nodes/avbd_world3d.hpp"

#include <thread>

#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/classes/time.hpp>
#include <godot_cpp/core/class_db.hpp>

#include "nodes/avbd_constraint3d.hpp"
#include "nodes/avbd_rigid_body3d.hpp"
#include "nodes/avbd_soft_body3d.hpp"
#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

// The solver's parameters are default-initialised on the members; the node's properties are
// pushed into it once per step (see _apply_params), so no setup is needed here.
AVBDWorld3D::AVBDWorld3D() = default;

AVBDWorld3D::~AVBDWorld3D() {
    solver.clear();
}

// ---------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------

void AVBDWorld3D::set_gravity(double p_gravity) {
    gravity = p_gravity;
}

double AVBDWorld3D::get_gravity() const {
    return gravity;
}

void AVBDWorld3D::set_iterations(int p_iterations) {
    iterations = p_iterations;
}

int AVBDWorld3D::get_iterations() const {
    return iterations;
}

void AVBDWorld3D::set_alpha(double p_alpha) {
    alpha = p_alpha;
}

double AVBDWorld3D::get_alpha() const {
    return alpha;
}

void AVBDWorld3D::set_beta_linear(double p_beta) {
    beta_linear = p_beta;
}

double AVBDWorld3D::get_beta_linear() const {
    return beta_linear;
}

void AVBDWorld3D::set_beta_angular(double p_beta) {
    beta_angular = p_beta;
}

double AVBDWorld3D::get_beta_angular() const {
    return beta_angular;
}

void AVBDWorld3D::set_gamma(double p_gamma) {
    gamma = p_gamma;
}

double AVBDWorld3D::get_gamma() const {
    return gamma;
}

void AVBDWorld3D::set_substeps(int p_substeps) {
    substeps = p_substeps;
}

int AVBDWorld3D::get_substeps() const {
    return substeps;
}

void AVBDWorld3D::set_paused(bool p_paused) {
    paused = p_paused;
}

bool AVBDWorld3D::is_paused() const {
    return paused;
}

void AVBDWorld3D::set_threads(int p_threads) {
    threads = p_threads;
}

int AVBDWorld3D::get_threads() const {
    return threads;
}

int AVBDWorld3D::get_thread_count() const {
    if (solver.threadCount() > 1) {
        return solver.threadCount();
    }
    // No worker pool yet (or a single-worker one): report what the setting asks for.
    return threads <= 0 ? static_cast<int>(std::thread::hardware_concurrency()) : threads;
}

// ---------------------------------------------------------------------------
// Introspection
// ---------------------------------------------------------------------------

int AVBDWorld3D::get_body_count() const {
    int count = 0;
    for (const avbd::Rigid *body = solver.bodies; body != nullptr; body = body->next) {
        count++;
    }
    return count;
}

int AVBDWorld3D::get_force_count() const {
    int count = 0;
    for (const avbd::Force *force = solver.forces; force != nullptr; force = force->next) {
        count++;
    }
    return count;
}

int AVBDWorld3D::get_contact_count() const {
    int count = 0;
    for (const avbd::Force *force = solver.forces; force != nullptr; force = force->next) {
        if (force->contactPointCount() > 0) {
            count++;
        }
    }
    return count;
}

int AVBDWorld3D::get_contact_point_count() const {
    int count = 0;
    for (const avbd::Force *force = solver.forces; force != nullptr; force = force->next) {
        count += force->contactPointCount();
    }
    return count;
}

uint64_t AVBDWorld3D::get_step_time_usec() const {
    return last_step_usec;
}

uint64_t AVBDWorld3D::get_step_count() const {
    return step_count;
}

avbd::Rigid *AVBDWorld3D::rigid_for(const AVBDRigidBody3D *p_body) const {
    for (const BodyEntry &entry : bodies) {
        if (entry.node == p_body) {
            return entry.rigid;
        }
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// Scene graph -> solver
// ---------------------------------------------------------------------------

void AVBDWorld3D::_collect(Node *p_node) {
    const int child_count = p_node->get_child_count();
    for (int i = 0; i < child_count; i++) {
        Node *child = p_node->get_child(i);
        if (child == nullptr || child->is_queued_for_deletion()) {
            continue;
        }
        if (Object::cast_to<AVBDWorld3D>(child) != nullptr) {
            continue; // a nested world owns its own subtree
        }
        if (AVBDRigidBody3D *body = Object::cast_to<AVBDRigidBody3D>(child)) {
            scan_bodies.push_back(body);
        } else if (AVBDSoftBody3D *soft = Object::cast_to<AVBDSoftBody3D>(child)) {
            scan_soft_bodies.push_back(soft);
        } else if (AVBDConstraint3D *constraint = Object::cast_to<AVBDConstraint3D>(child)) {
            scan_constraints.push_back(constraint);
        }
        _collect(child);
    }
}

void AVBDWorld3D::_scan() {
    scan_bodies.clear();
    scan_constraints.clear();
    scan_soft_bodies.clear();
    _collect(this);
}

bool AVBDWorld3D::_scan_changed() {
    _scan();
    if (scan_bodies.size() != bodies.size() || scan_constraints.size() != constraints.size() ||
            scan_soft_bodies.size() != soft_bodies.size()) {
        return true;
    }
    for (size_t i = 0; i < scan_bodies.size(); i++) {
        if (scan_bodies[i] != bodies[i].node) {
            return true;
        }
    }
    for (size_t i = 0; i < scan_constraints.size(); i++) {
        if (scan_constraints[i] != constraints[i]) {
            return true;
        }
    }
    for (size_t i = 0; i < scan_soft_bodies.size(); i++) {
        if (scan_soft_bodies[i] != soft_bodies[i].node) {
            return true;
        }
    }
    return false;
}

void AVBDWorld3D::rebuild() {
    _scan();
    _rebuild();
}

void AVBDWorld3D::_release() {
    for (BodyEntry &entry : bodies) {
        if (entry.node != nullptr) {
            entry.node->_bind_simulation(nullptr, nullptr);
        }
    }
    for (AVBDConstraint3D *node : constraints) {
        node->_bind_force(nullptr); // before the forces are freed
    }
    for (SoftBodyEntry &entry : soft_bodies) {
        entry.node->_clear();
    }
    solver.clear();
    bodies.clear();
    constraints.clear();
    soft_bodies.clear();
}

void AVBDWorld3D::_rebuild() {
    // Snapshot the current pose of every simulated body before the solver drops it,
    // so adding or removing one node does not reset the rest of the scene.
    std::vector<BodyEntry> previous = std::move(bodies);
    for (BodyEntry &entry : previous) {
        if (entry.rigid == nullptr) {
            continue;
        }
        entry.positionLin = entry.rigid->positionLin;
        entry.positionAng = entry.rigid->positionAng;
        entry.velocityLin = entry.rigid->velocityLin;
        entry.velocityAng = entry.rigid->velocityAng;
    }
    for (BodyEntry &entry : previous) {
        if (entry.node != nullptr) {
            entry.node->_bind_simulation(nullptr, nullptr);
        }
    }
    for (SoftBodyEntry &entry : soft_bodies) {
        entry.node->_store_state();
        entry.node->_clear();
    }

    solver.clear();
    bodies.clear();
    constraints.clear();
    soft_bodies.clear();

    const Transform3D world_global = get_global_transform();

    auto previous_entry = [&previous](const AVBDRigidBody3D *p_node) -> const BodyEntry * {
        for (const BodyEntry &entry : previous) {
            if (entry.node == p_node) {
                return &entry;
            }
        }
        return nullptr;
    };

    AVBDBodyMap map;
    bodies.reserve(scan_bodies.size());
    for (AVBDRigidBody3D *node : scan_bodies) {
        avbd::float3 position;
        avbd::quat rotation;
        avbd::float3 velocity;
        avbd::float3 angular_velocity;

        const BodyEntry *old = previous_entry(node);
        if (old != nullptr) {
            position = old->positionLin;
            rotation = old->positionAng;
            velocity = old->velocityLin;
            angular_velocity = old->velocityAng;
        } else {
            const Transform3D in_sim = node_in_sim(this, node);
            position = to_sim(in_sim.get_origin());
            rotation = to_sim(in_sim.get_basis().get_rotation_quaternion());
            velocity = direction_to_sim(world_global, node->get_initial_velocity());
            angular_velocity = direction_to_sim(world_global, node->get_initial_angular_velocity());
        }

        const float density = node->is_static_body() ? 0.0f : static_cast<float>(node->get_density());
        avbd::Rigid *rigid = new avbd::Rigid(&solver, to_sim_extents(node->get_size()), density,
                static_cast<float>(node->get_friction()), position, velocity);
        rigid->velocityAng = angular_velocity;
        node->_bind_simulation(this, rigid);

        BodyEntry entry;
        entry.node = node;
        entry.rigid = rigid;
        entry.positionLin = position;
        entry.positionAng = rotation;
        entry.velocityLin = velocity;
        entry.velocityAng = angular_velocity;
        bodies.push_back(entry);
        map[node] = rigid;
    }

    constraints = scan_constraints;
    for (AVBDConstraint3D *node : constraints) {
        if (node->_is_inactive()) {
            // E.g. a joint that already broke: no force, so its bodies collide again.
            node->_bind_force(nullptr);
            continue;
        }
        avbd::Force *created = node->create_force(solver, map, world_global);
        node->_bind_force(created);
        if (created == nullptr) {
            WARN_PRINT(vformat("AVBD: constraint '%s' could not be resolved and was skipped.",
                    String(node->get_path())));
        }
    }

    soft_bodies.reserve(scan_soft_bodies.size());
    for (AVBDSoftBody3D *node : scan_soft_bodies) {
        const Transform3D in_sim = node_in_sim(this, node);
        node->_build(solver, in_sim);
        SoftBodyEntry entry;
        entry.node = node;
        soft_bodies.push_back(entry);
    }
}

// ---------------------------------------------------------------------------
// Stepping
// ---------------------------------------------------------------------------

void AVBDWorld3D::_apply_params() {
    // The solver is Z-up with gravity along -Z; the node property is a magnitude along
    // Godot's -Y, which maps onto it exactly.
    solver.gravity = -static_cast<float>(gravity);
    solver.iterations = std::max(1, iterations);
    solver.alpha = static_cast<float>(alpha);
    solver.betaLin = static_cast<float>(beta_linear);
    solver.betaAng = static_cast<float>(beta_angular);
    solver.gamma = static_cast<float>(gamma);
    solver.threads = threads;
}

void AVBDWorld3D::_step(double p_delta) {
    _apply_params();

    const int sub = std::max(1, substeps);
    solver.dt = static_cast<float>(p_delta / sub);

    const uint64_t start = Time::get_singleton()->get_ticks_usec();
    for (int i = 0; i < sub; i++) {
        solver.step();
    }
    last_step_usec = Time::get_singleton()->get_ticks_usec() - start;
    step_count++;
}

// Static bodies are driven by their node transform: if a user (or an animation)
// moved one, adopt the new pose and stop it. Dynamic bodies are owned by the
// solver instead, and are moved with AVBDRigidBody3D::teleport().
void AVBDWorld3D::_sync_in_static_bodies() {
    for (BodyEntry &entry : bodies) {
        if (entry.rigid == nullptr || !entry.node->is_static_body()) {
            continue;
        }
        const Transform3D current = entry.node->get_global_transform();
        if (entry.pushed && current.is_equal_approx(entry.pushed_transform)) {
            continue;
        }
        const Transform3D in_sim = node_in_sim(this, entry.node);
        entry.rigid->positionLin = to_sim(in_sim.get_origin());
        entry.rigid->positionAng = to_sim(in_sim.get_basis().get_rotation_quaternion());
        entry.rigid->velocityLin = avbd::float3{0, 0, 0};
        entry.rigid->velocityAng = avbd::float3{0, 0, 0};
        entry.pushed_transform = current;
        entry.pushed = true;
    }
}

void AVBDWorld3D::_sync_out() {
    const Transform3D world_global = get_global_transform();

    for (BodyEntry &entry : bodies) {
        avbd::Rigid *rigid = entry.rigid;
        if (rigid == nullptr) {
            continue;
        }
        entry.positionLin = rigid->positionLin;
        entry.positionAng = rigid->positionAng;
        entry.velocityLin = rigid->velocityLin;
        entry.velocityAng = rigid->velocityAng;
        if (!entry.node->is_static_body()) {
            const Transform3D pose = world_global * sim_transform(rigid->positionLin, rigid->positionAng);
            entry.node->set_global_transform(pose);
            entry.pushed_transform = pose;
            entry.pushed = true;
        }
    }

    for (SoftBodyEntry &entry : soft_bodies) {
        entry.node->_store_state();
        const Transform3D sim_in_node = node_in_sim(this, entry.node).affine_inverse();
        entry.node->_sync_visuals(sim_in_node);
    }
}

void AVBDWorld3D::_physics_process(double p_delta) {
    if (_scan_changed()) {
        _rebuild();
    }
    if (paused) {
        return;
    }
    if (bodies.empty() && soft_bodies.empty()) {
        return;
    }
    _sync_in_static_bodies();
    _step(p_delta);
    _resolve_removed_forces();
    _sync_out();
}

// The solver deletes a force whose initialize() reported it inactive - for a joint,
// that happens when it breaks. Walk the surviving forces once, then let each
// constraint node notice that its force is gone, so it never holds a freed pointer.
void AVBDWorld3D::_resolve_removed_forces() {
    for (AVBDConstraint3D *node : constraints) {
        node->_mark_force_alive(false);
    }
    for (avbd::Force *force = solver.forces; force != nullptr; force = force->next) {
        if (force->owner != nullptr) {
            static_cast<AVBDConstraint3D *>(force->owner)->_mark_force_alive(true);
        }
    }
    for (AVBDConstraint3D *node : constraints) {
        node->_resolve_force_removal();
    }
}

void AVBDWorld3D::_ready() {
    set_physics_process(true);
}

void AVBDWorld3D::_notification(int p_what) {
    switch (p_what) {
        case NOTIFICATION_EXIT_TREE:
        case NOTIFICATION_PREDELETE:
            _release();
            break;
        default:
            break;
    }
}

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------

// Which node (and lattice cell) a solver body belongs to, or empty if none does.
AVBDWorld3D::Hit AVBDWorld3D::identify(const avbd::Rigid *p_rigid) const {
    for (const BodyEntry &entry : bodies) {
        if (entry.rigid == p_rigid) {
            return Hit{entry.node, Vector3(), 0.0, -1};
        }
    }
    for (const SoftBodyEntry &entry : soft_bodies) {
        const std::size_t cells = static_cast<std::size_t>(entry.node->get_cell_count());
        for (std::size_t i = 0; i < cells; i++) {
            if (entry.node->_get_cell_rigid(i) == p_rigid) {
                return Hit{entry.node, Vector3(), 0.0, static_cast<int>(i)};
            }
        }
    }
    return {};
}

// Trace a ray and report the first body it meets. The solver works in its own Z-up space, so
// the ray goes in and the hit comes back converted; everything in between is plain C++.
AVBDWorld3D::Hit AVBDWorld3D::trace(const Vector3 &p_origin, const Vector3 &p_direction, double p_max_distance) {
    Hit hit;
    const Vector3 direction = p_direction.normalized();
    if (direction.length_squared() <= 0.0) {
        return hit;
    }

    const Transform3D world_global = get_global_transform();
    avbd::float3 local{0, 0, 0};
    const avbd::Rigid *rigid = solver.pick(point_to_sim(world_global, p_origin),
            direction_to_sim(world_global, direction), local);
    if (rigid == nullptr) {
        return hit;
    }

    const avbd::float3 position = rigid->positionLin + avbd::rotate(rigid->positionAng, local);
    hit.position = world_global.xform(to_godot(position));
    hit.distance = p_origin.distance_to(hit.position);
    if (hit.distance > p_max_distance) {
        return Hit{};
    }

    const Hit owner = identify(rigid);
    hit.node = owner.node;
    hit.cell = owner.cell;
    return hit;
}

Dictionary AVBDWorld3D::raycast(const Vector3 &p_origin, const Vector3 &p_direction, double p_max_distance) {
    const Hit hit = trace(p_origin, p_direction, p_max_distance);
    Dictionary result;
    if (hit.node == nullptr) {
        return result;
    }
    result["body"] = hit.node;
    if (hit.cell >= 0) {
        result["cell"] = hit.cell;
    }
    result["position"] = hit.position;
    result["distance"] = hit.distance;
    return result;
}

// ---------------------------------------------------------------------------
// Bindings
// ---------------------------------------------------------------------------

void AVBDWorld3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_gravity", "gravity"), &AVBDWorld3D::set_gravity);
    ClassDB::bind_method(D_METHOD("get_gravity"), &AVBDWorld3D::get_gravity);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "gravity", PROPERTY_HINT_RANGE, "0,100,0.01,or_greater,suffix:m/s^2"),
            "set_gravity", "get_gravity");

    ClassDB::bind_method(D_METHOD("set_iterations", "iterations"), &AVBDWorld3D::set_iterations);
    ClassDB::bind_method(D_METHOD("get_iterations"), &AVBDWorld3D::get_iterations);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "iterations", PROPERTY_HINT_RANGE, "1,200,1"), "set_iterations",
            "get_iterations");

    ClassDB::bind_method(D_METHOD("set_alpha", "alpha"), &AVBDWorld3D::set_alpha);
    ClassDB::bind_method(D_METHOD("get_alpha"), &AVBDWorld3D::get_alpha);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "alpha", PROPERTY_HINT_RANGE, "0,1,0.001"), "set_alpha", "get_alpha");

    ClassDB::bind_method(D_METHOD("set_beta_linear", "beta"), &AVBDWorld3D::set_beta_linear);
    ClassDB::bind_method(D_METHOD("get_beta_linear"), &AVBDWorld3D::get_beta_linear);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "beta_linear", PROPERTY_HINT_RANGE, "0,100000000,1,or_greater"),
            "set_beta_linear", "get_beta_linear");

    ClassDB::bind_method(D_METHOD("set_beta_angular", "beta"), &AVBDWorld3D::set_beta_angular);
    ClassDB::bind_method(D_METHOD("get_beta_angular"), &AVBDWorld3D::get_beta_angular);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "beta_angular", PROPERTY_HINT_RANGE, "0,1000000,1,or_greater"),
            "set_beta_angular", "get_beta_angular");

    ClassDB::bind_method(D_METHOD("set_gamma", "gamma"), &AVBDWorld3D::set_gamma);
    ClassDB::bind_method(D_METHOD("get_gamma"), &AVBDWorld3D::get_gamma);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "gamma", PROPERTY_HINT_RANGE, "0,1,0.001"), "set_gamma", "get_gamma");

    ClassDB::bind_method(D_METHOD("set_substeps", "substeps"), &AVBDWorld3D::set_substeps);
    ClassDB::bind_method(D_METHOD("get_substeps"), &AVBDWorld3D::get_substeps);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "substeps", PROPERTY_HINT_RANGE, "1,8,1"), "set_substeps", "get_substeps");

    ClassDB::bind_method(D_METHOD("set_paused", "paused"), &AVBDWorld3D::set_paused);
    ClassDB::bind_method(D_METHOD("is_paused"), &AVBDWorld3D::is_paused);
    ADD_PROPERTY(PropertyInfo(Variant::BOOL, "paused"), "set_paused", "is_paused");

    ClassDB::bind_method(D_METHOD("set_threads", "threads"), &AVBDWorld3D::set_threads);
    ClassDB::bind_method(D_METHOD("get_threads"), &AVBDWorld3D::get_threads);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "threads", PROPERTY_HINT_RANGE, "0,256,1"), "set_threads", "get_threads");

    ClassDB::bind_method(D_METHOD("get_thread_count"), &AVBDWorld3D::get_thread_count);

    ClassDB::bind_method(D_METHOD("get_body_count"), &AVBDWorld3D::get_body_count);
    ClassDB::bind_method(D_METHOD("get_force_count"), &AVBDWorld3D::get_force_count);
    ClassDB::bind_method(D_METHOD("get_contact_count"), &AVBDWorld3D::get_contact_count);
    ClassDB::bind_method(D_METHOD("get_contact_point_count"), &AVBDWorld3D::get_contact_point_count);
    ClassDB::bind_method(D_METHOD("get_step_time_usec"), &AVBDWorld3D::get_step_time_usec);
    ClassDB::bind_method(D_METHOD("get_step_count"), &AVBDWorld3D::get_step_count);

    ClassDB::bind_method(D_METHOD("rebuild"), &AVBDWorld3D::rebuild);
    ClassDB::bind_method(D_METHOD("raycast", "origin", "direction", "max_distance"), &AVBDWorld3D::raycast,
            DEFVAL(10000.0));
}
