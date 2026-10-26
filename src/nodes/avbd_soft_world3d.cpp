/*
 * AVBDSoftWorld3D implementation. See the header.
 *
 * The rebuild / step / sync pattern follows the soft-body branches of the retired
 * AVBDWorld3D: snapshot each lattice's cell poses before the solver drops them so a
 * rebuild does not reset the scene, rebuild only when the node set changes, and push
 * cell poses back through _sync_visuals every simulated frame.
 */

#include "nodes/avbd_soft_world3d.hpp"

#include <algorithm>
#include <ranges>

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/object.hpp>

#include "avbd/list_range.hpp"
#include "nodes/avbd_soft_body3d.hpp"
#include "nodes/godot_convert.hpp"

using namespace godot;
using namespace avbd_godot;

AVBDSoftWorld3D::~AVBDSoftWorld3D() {
    solver.clear();
}

void AVBDSoftWorld3D::set_gravity(double p_gravity) {
    gravity = p_gravity;
}

double AVBDSoftWorld3D::get_gravity() const {
    return gravity;
}

void AVBDSoftWorld3D::set_threads(int p_threads) {
    threads = p_threads;
}

int AVBDSoftWorld3D::get_threads() const {
    return threads;
}

int AVBDSoftWorld3D::get_body_count() const {
    return static_cast<int>(std::ranges::distance(avbd::next_range(solver.bodies)));
}

void AVBDSoftWorld3D::set_tick_rate(double p_tick_rate) {
    // 0 means "no clock"; anything else is a rate the solver can actually divide by.
    tick_rate = p_tick_rate > 0.0 ? std::clamp(p_tick_rate, 1.0, 1000.0) : 0.0;
}

double AVBDSoftWorld3D::get_tick_rate() const {
    return tick_rate;
}

void AVBDSoftWorld3D::set_max_catchup_ticks(int p_ticks) {
    max_catchup_ticks = std::clamp(p_ticks, 1, 64);
}

int AVBDSoftWorld3D::get_max_catchup_ticks() const {
    return max_catchup_ticks;
}

void AVBDSoftWorld3D::rebuild() {
    _collect(this);
    _rebuild();
}

void AVBDSoftWorld3D::_collect(Node *p_node) {
    const int child_count = p_node->get_child_count();
    for (int i = 0; i < child_count; i++) {
        Node *child = p_node->get_child(i);
        if (child == nullptr || child->is_queued_for_deletion()) {
            continue;
        }
        if (Object::cast_to<AVBDSoftWorld3D>(child) != nullptr) {
            continue; // a nested soft world owns its own subtree
        }
        if (AVBDSoftBody3D *soft = Object::cast_to<AVBDSoftBody3D>(child)) {
            scan_soft_bodies.push_back(soft);
        } else if (StaticBody3D *ground = Object::cast_to<StaticBody3D>(child)) {
            // A ground body must carry exactly one box shape; that is the demo's needs and
            // keeps this node minimal. Anything else is skipped with a warning.
            GroundEntry entry;
            entry.node = ground;
            if (const BoxShape3D *box = Object::cast_to<const BoxShape3D>(_first_shape(ground))) {
                entry.extents = box->get_size();
            }
            if (entry.extents == Vector3()) {
                WARN_PRINT("AVBDSoftWorld3D: a StaticBody3D has no BoxShape3D; skipped as ground.");
            } else {
                scan_grounds.push_back(entry);
            }
        }
        _collect(child);
    }
}

Shape3D *AVBDSoftWorld3D::_first_shape(StaticBody3D *p_body) {
	// Scan children directly: shape-owner registration is timing-sensitive during
	// scene instantiation, and the ground contract is "a CollisionShape3D child
	// carrying one BoxShape3D".
	const int child_count = p_body->get_child_count();
	for (int i = 0; i < child_count; i++) {
		CollisionShape3D *collision = Object::cast_to<CollisionShape3D>(p_body->get_child(i));
		if (collision != nullptr && collision->get_shape().is_valid()) {
			return collision->get_shape().ptr();
		}
	}
	return nullptr;
}

bool AVBDSoftWorld3D::_scan_changed() {
    scan_soft_bodies.clear();
    scan_grounds.clear();
    _collect(this);
    if (scan_soft_bodies.size() != soft_bodies.size() || scan_grounds.size() != grounds.size()) {
        return true;
    }
    if (!std::ranges::equal(scan_soft_bodies, soft_bodies)) {
        return true;
    }
    if (!std::ranges::equal(scan_grounds, grounds, {}, &GroundEntry::node, &GroundEntry::node)) {
        return true;
    }
    return false;
}

void AVBDSoftWorld3D::_release() {
    for (AVBDSoftBody3D *node : soft_bodies) {
        node->_clear();
    }
    solver.clear();
    soft_bodies.clear();
    grounds.clear();
}

void AVBDSoftWorld3D::_rebuild() {
    // Snapshot the lattices' current poses before the solver drops them, so adding or
    // removing one soft body does not reset the others.
    for (AVBDSoftBody3D *node : soft_bodies) {
        node->_store_state();
        node->_clear();
    }
    solver.clear();
    soft_bodies.clear();
    grounds.clear();

    grounds.reserve(scan_grounds.size());
    for (GroundEntry &entry : scan_grounds) {
        const Transform3D in_sim = node_in_sim(this, entry.node);
        avbd::Rigid *rigid = new avbd::Rigid(&solver, to_sim_extents(entry.extents), 0.0f, 0.6f,
                to_sim(in_sim.origin));
        rigid->positionAng = to_sim(in_sim.basis.get_rotation_quaternion());
        entry.rigid = rigid;
        grounds.push_back(entry);
    }

    soft_bodies.reserve(scan_soft_bodies.size());
    for (AVBDSoftBody3D *node : scan_soft_bodies) {
        const Transform3D in_sim = node_in_sim(this, node);
        node->_build(solver, in_sim);
        soft_bodies.push_back(node);
    }
}

void AVBDSoftWorld3D::_physics_process(double p_delta) {
    // Solver gravity is the signed magnitude along the simulation's -Z; the exposed
    // property is a positive magnitude. Set BEFORE the rebuild: new rigids snapshot
    // the solver's gravity at construction.
    solver.gravity = -static_cast<float>(gravity);

    if (_scan_changed()) {
        _rebuild();
    }
    if (soft_bodies.empty()) {
        return;
    }

    solver.threads = threads;

    // The same two clocks the physics server offers. With a rate the solver steps on its own fixed
    // dt and the scene tree's delta only says how much of it has been paid for; without one the
    // delta is the step, which is what this node always did - and what made a lattice soft or stiff
    // depend on the game's tick rate.
    int ticks = 1;
    double alpha = 1.0;
    if (tick_rate > 0.0) {
        clock.tick_dt = 1.0 / tick_rate;
        clock.max_ticks = max_catchup_ticks;
        ticks = clock.consume(p_delta, alpha);
        solver.dt = static_cast<float>(clock.tick_dt);
    } else {
        solver.dt = static_cast<float>(p_delta);
    }

    // Static grounds follow their node: if a user or animation moved one, adopt the pose.
    for (GroundEntry &entry : grounds) {
        if (entry.rigid == nullptr) {
            continue;
        }
        const Transform3D in_sim = node_in_sim(this, entry.node);
        entry.rigid->positionLin = to_sim(in_sim.origin);
        entry.rigid->positionAng = to_sim(in_sim.basis.get_rotation_quaternion());
        entry.rigid->velocityLin = avbd::float3{0, 0, 0};
        entry.rigid->velocityAng = avbd::float3{0, 0, 0};
    }

    for (int tick = 0; tick < ticks; tick++) {
        // The pose entering the last tick is the half the visuals blend from, so snapshot it there
        // rather than after the step. A notification that ran no tick leaves the pair alone, which
        // is right: nothing new was simulated, and alpha has moved closer to the live pose.
        if (tick + 1 == ticks) {
            for (AVBDSoftBody3D *node : soft_bodies) {
                node->_store_state();
            }
        }
        solver.step();
    }

    const Transform3D sim_in_node = node_in_sim(this, this).affine_inverse();
    for (AVBDSoftBody3D *node : soft_bodies) {
        node->_sync_visuals(sim_in_node, alpha);
    }
}

void AVBDSoftWorld3D::_ready() {
    set_physics_process(true);
}

void AVBDSoftWorld3D::_notification(int p_what) {
    switch (p_what) {
        case NOTIFICATION_EXIT_TREE:
        case NOTIFICATION_PREDELETE:
            _release();
            break;
        default:
            break;
    }
}

void AVBDSoftWorld3D::_bind_methods() {
    ClassDB::bind_method(D_METHOD("set_gravity", "gravity"), &AVBDSoftWorld3D::set_gravity);
    ClassDB::bind_method(D_METHOD("get_gravity"), &AVBDSoftWorld3D::get_gravity);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "gravity", PROPERTY_HINT_RANGE, "0,100,0.01,or_greater,suffix:m/s^2"),
            "set_gravity", "get_gravity");

    ClassDB::bind_method(D_METHOD("set_threads", "threads"), &AVBDSoftWorld3D::set_threads);
    ClassDB::bind_method(D_METHOD("get_threads"), &AVBDSoftWorld3D::get_threads);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "threads", PROPERTY_HINT_RANGE, "0,256,1"), "set_threads", "get_threads");

    ClassDB::bind_method(D_METHOD("get_body_count"), &AVBDSoftWorld3D::get_body_count);
    ClassDB::bind_method(D_METHOD("rebuild"), &AVBDSoftWorld3D::rebuild);

    ClassDB::bind_method(D_METHOD("set_tick_rate", "tick_rate"), &AVBDSoftWorld3D::set_tick_rate);
    ClassDB::bind_method(D_METHOD("get_tick_rate"), &AVBDSoftWorld3D::get_tick_rate);
    ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "tick_rate", PROPERTY_HINT_RANGE, "0,1000,0.1,or_greater,suffix:Hz"),
            "set_tick_rate", "get_tick_rate");

    ClassDB::bind_method(D_METHOD("set_max_catchup_ticks", "ticks"), &AVBDSoftWorld3D::set_max_catchup_ticks);
    ClassDB::bind_method(D_METHOD("get_max_catchup_ticks"), &AVBDSoftWorld3D::get_max_catchup_ticks);
    ADD_PROPERTY(PropertyInfo(Variant::INT, "max_catchup_ticks", PROPERTY_HINT_RANGE, "1,64,1"),
            "set_max_catchup_ticks", "get_max_catchup_ticks");
}
