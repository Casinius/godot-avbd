/*
 * AVBD as a Godot physics server. See the header for the design and the unsupported list.
 */

#include "server/avbd_physics_server3d.hpp"

#include "server/avbd_direct_body_state3d.hpp"
#include "server/avbd_direct_space_state3d.hpp"
#include "server/modern_convert.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <array>
#include <ranges>
#include <span>
#include <tuple>

#include <godot_cpp/classes/box_shape3d.hpp>
#include <godot_cpp/classes/cylinder_shape3d.hpp>
#include <godot_cpp/classes/sphere_shape3d.hpp>
#include <godot_cpp/classes/project_settings.hpp>
#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>
#include <godot_cpp/variant/array.hpp>

using namespace godot;
// Singleton pointer declaration.
AVBDPhysicsServer3D *AVBDPhysicsServer3D::SpaceData::singleton = nullptr;

// --- step profiling (AVBD_PROFILE=1) -----------------------------------------
//
// Env-gated aggregate timers for headless runs that need to attribute per-tick
// cost to a phase: rebuild, project-settings reads, force integration, the
// solver itself, contact collection, sync dispatch. Off by default; the on
// path is a handful of steady_clock reads per step.

namespace {
using ProfClock = std::chrono::steady_clock;

double seconds_since(const ProfClock::time_point &p_since) {
    return std::chrono::duration<double, std::nano>(ProfClock::now() - p_since).count();
}

struct StepProfile {
    bool enabled = false;
    bool decided = false;
    uint64_t ticks = 0;
    double total_ns = 0.0;
    double solver_ns = 0.0;
    double sync_ns = 0.0;
    double contacts_ns = 0.0;
    double rebuild_ns = 0.0;
    double params_ns = 0.0;
    double forces_ns = 0.0;
    // Spike attribution: per-phase maxima within the reporting window. The average
    // hides exactly what a player feels - the one tick in a hundred that stalls the
    // frame (a deep-impact tick, a broad-phase pair storm, a rebuild).
    double max_total_ns = 0.0;
    double max_solver_ns = 0.0;
    double max_sync_ns = 0.0;
    double max_contacts_ns = 0.0;
    double max_rebuild_ns = 0.0;
    double max_params_ns = 0.0;
    double max_forces_ns = 0.0;
    uint32_t over_2ms_ticks = 0; // ticks above a 2 ms budget in the window
    uint32_t rebuilds = 0;
    int threads_used = 1;
    int substeps = 0;
    int iterations = 0;
};

StepProfile &step_profile() {
    static StepProfile profile;
    if (!profile.decided) {
        profile.decided = true;
        const char *env = std::getenv("AVBD_PROFILE");
        profile.enabled = env != nullptr && env[0] != '\0' && env[0] != '0';
    }
    return profile;
}
} // namespace

// Godot 4.3 added per-axis springs to the generic 6-DOF joint and grew the server enums; the
// godot-cpp headers this extension builds against predate that, so the new slots are named
// here with the engine's values (servers/physics_server_3d.h, G6DOFJointAxisParam/Flag).
namespace g6dof_spring {
constexpr int LINEAR_STIFFNESS = 7;
constexpr int LINEAR_DAMPING = 8;
constexpr int LINEAR_EQUILIBRIUM = 9;
constexpr int ANGULAR_STIFFNESS = 19;
constexpr int ANGULAR_DAMPING = 20;
constexpr int ANGULAR_EQUILIBRIUM = 21;
constexpr int FLAG_ENABLE_ANGULAR = 2;
constexpr int FLAG_ENABLE_LINEAR = 3;
} // namespace g6dof_spring

// -----------------------------------------------------------------------------
// Handles
//
// The extension interface offers no way to allocate a RID, so the server mints its own: an
// integer id written into the RID's opaque bytes through the pointer godot-cpp does expose. The
// engine only ever hands these values back to us, so the encoding is ours to choose - but it must
// round-trip, which is what id_of() checks by reading them back.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::make_rid(uint64_t p_id) const {
    RID rid;
    std::memcpy(rid._native_ptr(), &p_id, sizeof(p_id));
    return rid;
}

uint64_t AVBDPhysicsServer3D::id_of(const RID &p_rid) const {
    uint64_t id = 0;
    std::memcpy(&id, p_rid._native_ptr(), sizeof(id));
    return id;
}

AVBDPhysicsServer3D::BodyData *AVBDPhysicsServer3D::find_body(const RID &p_rid) {
    const auto found = bodies.find(id_of(p_rid));
    return found == bodies.end() ? nullptr : &found->second;
}

const AVBDPhysicsServer3D::BodyData *AVBDPhysicsServer3D::find_body(const RID &p_rid) const {
    const auto found = bodies.find(id_of(p_rid));
    return found == bodies.end() ? nullptr : &found->second;
}

AVBDPhysicsServer3D::ShapeData *AVBDPhysicsServer3D::find_shape(const RID &p_rid) {
    const auto found = shapes.find(id_of(p_rid));
    return found == shapes.end() ? nullptr : &found->second;
}

AVBDPhysicsServer3D::SpaceData *AVBDPhysicsServer3D::find_space(const RID &p_rid) {
    const auto found = spaces.find(id_of(p_rid));
    return found == spaces.end() ? nullptr : found->second.get();
}

const AVBDPhysicsServer3D::SpaceData *AVBDPhysicsServer3D::find_space(const RID &p_rid) const {
    const auto found = spaces.find(id_of(p_rid));
    return found == spaces.end() ? nullptr : found->second.get();
}

ModernSpace *AVBDPhysicsServer3D::modern_for(const BodyData &p_body) {
    if (p_body.modern_index < 0) {
        return nullptr;
    }
    const auto found = spaces.find(id_of(p_body.space));
    return found == spaces.end() ? nullptr : found->second->modern.get();
}

const ModernSpace *AVBDPhysicsServer3D::modern_for(const BodyData &p_body) const {
    if (p_body.modern_index < 0) {
        return nullptr;
    }
    const auto found = spaces.find(id_of(p_body.space));
    return found == spaces.end() ? nullptr : found->second->modern.get();
}

// -----------------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------------

// Spaces may outlive _finish in edge paths (a space freed without the engine stopping
// first): the destructor is the last chance to release the cached state objects.
AVBDPhysicsServer3D::~AVBDPhysicsServer3D() {
    for (auto &[id, space_ptr] : spaces) {
        (void)id;
        SpaceData &space = *space_ptr;
        if (space.direct_state != nullptr) {
            memdelete(space.direct_state);
            space.direct_state = nullptr;
        }
    }
    if (body_state != nullptr) {
        memdelete(body_state);
        body_state = nullptr;
    }
}

void AVBDPhysicsServer3D::_finish() {
    // The engine is shutting the server down; drop everything it owns. The cached state
    // objects are raw pointers inside the records, so they go first.
    for (auto &[id, space_ptr] : spaces) {
        (void)id;
        SpaceData &space = *space_ptr;
        if (space.direct_state != nullptr) {
            memdelete(space.direct_state);
            space.direct_state = nullptr;
        }
    }
    spaces.clear();
    bodies.clear();
    shapes.clear();
    joints.clear();
    areas.clear();
    if (body_state != nullptr) {
        memdelete(body_state);
        body_state = nullptr;
    }
}

void AVBDPhysicsServer3D::_init() {
    // Spaces created after this moment read gravity themselves; this covers the ones the engine
    // creates before any of our code runs.
    for (auto &[id, space_ptr] : spaces) {
        read_project_gravity(*space_ptr);
        read_project_params(*space_ptr);
    }
    if (ProjectSettings *settings = ProjectSettings::get_singleton()) {
        // Solver parameters are re-read only when a setting actually changes (the signal
        // is deferred), instead of hashing the setting names on every tick.
        settings->connect("settings_changed", Callable(this, "_on_project_settings_changed"));
    }
}

void AVBDPhysicsServer3D::_on_project_settings_changed() {
    for (auto &[id, space_ptr] : spaces) {
        (void)id;
        space_ptr->params_dirty = true;
    }
}

void AVBDPhysicsServer3D::_flush_queries() {
    flushing_queries = true;
    // The engine flushes before stepping, so this is where the last step's poses reach
    // the scene: the per-body state-sync callable, then scripted integration. Both read
    // and write the solver's live state, so forces a script applies in _integrate_forces
    // land in the step about to run (the GodotPhysics contract, and godot-jolt's shape).
    // The first flush of a scene has no step behind it yet: nothing to deliver.
    if (has_stepped) {
        StepProfile &prof = step_profile();
        const ProfClock::time_point prof_t_sync = ProfClock::now();

        static thread_local Array one_arg = []() {
            Array array;
            array.resize(1);
            return array;
        }();
        static thread_local Array two_args = []() {
            Array array;
            array.resize(2);
            return array;
        }();

        for (auto &[body_id, body] : bodies) {
            (void)body_id;
            if (!body.sync.is_valid() && !body.force_integration_callback.is_valid()) {
                continue;
            }
            body_state->bind(this, make_rid(body_id));
            if (body.sync.is_valid()) {
                one_arg[0] = body_state;
                body.sync.callv(one_arg);
            }
            if (body.force_integration_callback.is_valid()) {
                if (body.force_integration_userdata.get_type() != Variant::NIL) {
                    two_args[0] = body_state;
                    two_args[1] = body.force_integration_userdata;
                    body.force_integration_callback.callv(two_args);
                } else {
                    one_arg[0] = body_state;
                    body.force_integration_callback.callv(one_arg);
                }
            }
        }
        if (prof.enabled) {
            const double sync_elapsed = seconds_since(prof_t_sync);
            prof.sync_ns += sync_elapsed;
            prof.max_sync_ns = std::max(prof.max_sync_ns, sync_elapsed);
        }
    }
    flushing_queries = false;
}

bool AVBDPhysicsServer3D::_is_flushing_queries() const {
    return flushing_queries;
}

int32_t AVBDPhysicsServer3D::_get_process_info(PhysicsServer3D::ProcessInfo p_info) {
    switch (p_info) {
        case PhysicsServer3D::INFO_ACTIVE_OBJECTS:
            return static_cast<int32_t>(std::ranges::count_if(bodies | std::views::values,
                    [](const BodyData &body) { return body.modern_index >= 0; }));
        case PhysicsServer3D::INFO_ISLAND_COUNT:
            return static_cast<int32_t>(spaces.size());
        case PhysicsServer3D::INFO_COLLISION_PAIRS: {
            // A manifold with at least one contact is a colliding pair.
            int32_t count = 0;
            for (const auto &[id, space] : spaces) {
                (void)id;
                if (space->modern == nullptr)
                    continue;
                for (const collision::Manifold &m : space->modern->manifolds())
                    if (m.count > 0)
                        count++;
            }
            return count;
        }
        default:
            return 0;
    }
}

void AVBDPhysicsServer3D::_set_active(bool p_active) {
    // SceneTree::set_pause calls PhysicsServer3D.set_active(!paused): this is the point
    // where simulation must follow the game. _step refuses to run while false.
    server_active = p_active;
}

void AVBDPhysicsServer3D::_sync() {
    // Called from the physics thread before stepping. The bodies were already written by
    // _body_set_state, so there is nothing to pull in here.
}

void AVBDPhysicsServer3D::_end_sync() {
    // ...and after stepping, which is when a real implementation would push results back. Ours
    // pushes them at the end of _step, where the solver's output is in hand.
}

void AVBDPhysicsServer3D::_free_rid(const RID &p_rid) {
    const uint64_t id = id_of(p_rid);

    const auto space_found = spaces.find(id);
    if (space_found != spaces.end()) {
        if (space_found->second->direct_state != nullptr) {
            memdelete(space_found->second->direct_state);
        }
        spaces.erase(space_found);
        return;
    }
    if (bodies.erase(id) > 0) {
        return;
    }
    if (shapes.erase(id) > 0) {
        return;
    }
    if (joints.erase(id) > 0) {
        return;
    }
    areas.erase(id);
}

// -----------------------------------------------------------------------------
// Spaces
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_space_create() {
    const uint64_t id = next_id++;
    // Solver owns intrusive lists and is deliberately non-movable, so SpaceData lives behind a
    // unique_ptr: the flat map relocates the handle, never the solver it owns.
    auto [space_it, space_inserted] = spaces.emplace(id, std::make_unique<SpaceData>());
    (void)space_inserted;
    SpaceData &space = *space_it->second;
    read_project_gravity(space);
    read_project_params(space);
    return make_rid(id);
}

// Solver parameters come from project settings the same way gravity does, so a project can
// tune the solver without a node layer. Keys absent from a project fall back to the solver's
// defaults.
void AVBDPhysicsServer3D::read_project_params(SpaceData &p_space) {
    const ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    // Policy selection. B-PXPBD is the default now; "avbd" remains an opt-in legacy
    // policy. The defaults (20 iterations, zeta 0.4, 2e-4 contact compliance) are the
    // parity configuration the solver/constraint suites pass with, so an absent
    // setting gets those instead of the old AVBD values.
    const String solver_mode = settings->get_setting("physics/avbd/solver", "bpxpbd");
    const bool bpxpbd = solver_mode == "bpxpbd";
    p_space.solver = bpxpbd ? 1 : 0;
    p_space.threads = static_cast<int>(settings->get_setting("physics/avbd/threads", 0));
    p_space.iterations = static_cast<int>(
            settings->get_setting("physics/avbd/iterations", bpxpbd ? 20 : 8));
    p_space.bpxpbd_zeta = settings->get_setting("physics/avbd/bpxpbd_zeta", 0.4);
    p_space.bpxpbd_contact_compliance =
            settings->get_setting("physics/avbd/bpxpbd_contact_compliance", 2e-4);
    p_space.alpha = settings->get_setting("physics/avbd/alpha", 0.99);
    p_space.beta_linear = settings->get_setting("physics/avbd/beta_linear", 1e4);
    p_space.beta_angular = settings->get_setting("physics/avbd/beta_angular", 1e2);
    p_space.gamma = settings->get_setting("physics/avbd/gamma", 0.999);
    p_space.adaptive_beta = settings->get_setting("physics/avbd/adaptive_beta", false);
    p_space.newton_ratio = settings->get_setting("physics/avbd/newton_ratio", 1.0);
    p_space.stiffness_decay = settings->get_setting("physics/avbd/stiffness_decay", 1.0);
    p_space.substeps = static_cast<int>(settings->get_setting("physics/avbd/substeps", 2));
    p_space.jacobian_rebuild_distance = settings->get_setting("physics/avbd/jacobian_rebuild_distance", 0.0);
    p_space.interleaved_substeps = settings->get_setting("physics/avbd/interleaved_substeps", false);
    p_space.rounds_per_substep = static_cast<int>(settings->get_setting("physics/avbd/rounds_per_substep", 0.0));
    const double internal_hz = std::max(static_cast<double>(settings->get_setting(
            "physics/avbd/internal_ticks_per_second", 0.0)), 0.0);
    if (internal_hz != p_space.internal_ticks_per_second) {
        // A rate change is a schedule change: do not carry a stale backlog into it.
        p_space.internal_accumulator = 0.0;
    }
    p_space.internal_ticks_per_second = internal_hz;
    p_space.max_internal_steps = std::max(static_cast<int>(settings->get_setting(
            "physics/avbd/max_internal_steps", 4.0)), 1);
    p_space.position_projection = settings->get_setting("physics/avbd/position_projection", true);
    p_space.projection_iterations = static_cast<int>(settings->get_setting(
            "physics/avbd/projection_iterations", 0.0));
    p_space.projection_beta = settings->get_setting("physics/avbd/projection_beta", 1.0);
    p_space.projection_slop = settings->get_setting("physics/avbd/projection_slop", 0.008);
    p_space.projection_max_correction = settings->get_setting(
            "physics/avbd/projection_max_correction", 0.8);
    p_space.sleep_linear_threshold = settings->get_setting(
            "physics/avbd/sleep_linear_threshold", 0.05);
    p_space.sleep_angular_threshold = settings->get_setting(
            "physics/avbd/sleep_angular_threshold", 0.05);
    p_space.sleep_time = settings->get_setting("physics/avbd/sleep_time", 1.0);
    p_space.contact_grace_ticks = static_cast<int>(settings->get_setting(
            "physics/avbd/contact_grace_ticks", 0.0));
    p_space.contact_gs = settings->get_setting("physics/avbd/contact_gs", false);
    p_space.speculative_margin = settings->get_setting("physics/avbd/speculative_margin", 0.01);
}

void AVBDPhysicsServer3D::_space_set_active(const RID &p_space, bool p_active) {
    if (SpaceData *space = find_space(p_space)) {
        space->active = p_active;
    }
}

bool AVBDPhysicsServer3D::_space_is_active(const RID &p_space) const {
    const auto found = spaces.find(id_of(p_space));
    return found == spaces.end() ? false : found->second->active;
}

void AVBDPhysicsServer3D::_space_set_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param, float p_value) {
    // The space parameters are contact tolerances and sleep thresholds, none of which map onto
    // AVBD's model. They are accepted and ignored rather than rejected, so a project's settings do
    // not error out on load. Gravity does not come through here at all: Godot keeps it in the
    // project settings, and a server is expected to read it (see read_project_gravity).
    (void)p_space;
    (void)p_param;
    (void)p_value;
}

// Godot stores gravity in the project settings rather than passing it to the server, so the server
// reads it - the same thing the built-in one does.
void AVBDPhysicsServer3D::read_project_gravity(SpaceData &p_space) {
    const ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    const double magnitude = settings->get_setting("physics/3d/default_gravity", 9.8);
    const Vector3 direction = settings->get_setting("physics/3d/default_gravity_vector", Vector3(0, -1, 0));

    // Both the engine and the modern core are Y-up, so the vector passes through as-is.
    // `default_gravity_vector` is expected to be normalised (Godot normalises it when set);
    // normalise defensively so a hand-edited project setting cannot scale gravity twice.
    const Vector3 unit = direction.is_zero_approx() ? Vector3(0, -1, 0) : direction.normalized();
    p_space.gravity = to_modern(unit * magnitude);
}

float AVBDPhysicsServer3D::_space_get_param(const RID &p_space, PhysicsServer3D::SpaceParameter p_param) const {
    (void)p_space;
    (void)p_param;
    return 0.0;
}

// One state object per space, created on first ask. PhysicsDirectSpaceState3D is
// reference-counted by the engine, so the server holds a raw pointer to the object it
// owns for the space's lifetime.
PhysicsDirectSpaceState3D *AVBDPhysicsServer3D::_space_get_direct_state(const RID &p_space) {
    SpaceData *space = find_space(p_space);
    if (space == nullptr) {
        return nullptr;
    }
    if (space->direct_state == nullptr) {
        space->direct_state = memnew(AVBDDirectSpaceState3D);
        space->direct_state->bind(this, p_space);
    }
    return space->direct_state;
}

void AVBDPhysicsServer3D::_space_set_debug_contacts(const RID &p_space, int32_t p_max_contacts) {
    if (SpaceData *space = find_space(p_space)) {
        space->debug_contacts_requested = p_max_contacts > 0;
        space->debug_contacts_max = p_max_contacts;
    }
}

PackedVector3Array AVBDPhysicsServer3D::_space_get_contacts(const RID &p_space) const {
    PackedVector3Array out;
    const auto found = spaces.find(id_of(p_space));
    if (found == spaces.end()) {
        return out;
    }
    const SpaceData &space = *found->second;
    if (space.modern == nullptr) {
        return out;
    }
    // The engine's debug contacts are world-space points of recent contacts; serve every
    // manifold's A-side contact points seen this step.
    const State &state = space.modern->state();
    for (const collision::Manifold &m : space.modern->manifolds()) {
        for (std::size_t i = 0; i < m.count; ++i) {
            const Vec3 world = state.pose.position(m.body_a) +
                    rotate(state.pose.rotation(m.body_a), m.contacts[i].rA);
            out.push_back(to_godot(world));
            if (out.size() >= static_cast<int64_t>(space.debug_contacts_max)) {
                return out;
            }
        }
    }
    return out;
}

int32_t AVBDPhysicsServer3D::_space_get_contact_count(const RID &p_space) const {
    const auto found = spaces.find(id_of(p_space));
    if (found == spaces.end()) {
        return 0;
    }
    if (found->second->modern == nullptr) {
        return 0;
    }
    int32_t points = 0;
    for (const collision::Manifold &m : found->second->modern->manifolds())
        points += static_cast<int32_t>(m.count);
    return points;
}

// -----------------------------------------------------------------------------
// Shapes
//
// The engine sends a shape's parameters as a Variant from the matching *Shape3D resource, which
// is also what it can read back, so the getters rebuild the same object.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_box_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_BOX;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_sphere_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_SPHERE;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_cylinder_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CYLINDER;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_capsule_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CAPSULE;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_convex_polygon_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CONVEX_POLYGON;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_concave_polygon_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CONCAVE_POLYGON;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_heightmap_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_HEIGHTMAP;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_world_boundary_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_WORLD_BOUNDARY;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_separation_ray_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_SEPARATION_RAY;
    return make_rid(id);
}

RID AVBDPhysicsServer3D::_custom_shape_create() {
    const uint64_t id = next_id++;
    ShapeData &shape = shapes[id];
    shape.type = PhysicsServer3D::SHAPE_CUSTOM;
    return make_rid(id);
}

void AVBDPhysicsServer3D::_shape_set_custom_solver_bias(const RID &p_shape, float p_bias) {
    if (ShapeData *shape = find_shape(p_shape)) {
        shape->custom_bias = p_bias;
    }
}

float AVBDPhysicsServer3D::_shape_get_custom_solver_bias(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? 0.0f : found->second.custom_bias;
}

void AVBDPhysicsServer3D::_shape_set_margin(const RID &p_shape, float p_margin) {
    if (ShapeData *shape = find_shape(p_shape)) {
        shape->margin = p_margin;
    }
}

float AVBDPhysicsServer3D::_shape_get_margin(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? 0.0f : found->second.margin;
}

// -----------------------------------------------------------------------------
// Soft-body server API: documented unsupported. Godot only drives these through
// SoftBody3D nodes; every method answers neutrally and nothing is simulated.
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_soft_body_create() {
    return make_rid(next_id++);
}
void AVBDPhysicsServer3D::_soft_body_set_space(const RID &p_body, const RID &p_space) { (void)p_body; (void)p_space; }
RID AVBDPhysicsServer3D::_soft_body_get_space(const RID &p_body) const { (void)p_body; return RID(); }
void AVBDPhysicsServer3D::_soft_body_set_mesh(const RID &p_body, const RID &p_mesh) { (void)p_body; (void)p_mesh; }
void AVBDPhysicsServer3D::_soft_body_set_transform(const RID &p_body, const Transform3D &p_transform) { (void)p_body; (void)p_transform; }
void AVBDPhysicsServer3D::_soft_body_set_collision_layer(const RID &p_body, uint32_t p_layer) { (void)p_body; (void)p_layer; }
uint32_t AVBDPhysicsServer3D::_soft_body_get_collision_layer(const RID &p_body) const { (void)p_body; return 1; }
void AVBDPhysicsServer3D::_soft_body_set_collision_mask(const RID &p_body, uint32_t p_mask) { (void)p_body; (void)p_mask; }
uint32_t AVBDPhysicsServer3D::_soft_body_get_collision_mask(const RID &p_body) const { (void)p_body; return 1; }
void AVBDPhysicsServer3D::_soft_body_add_collision_exception(const RID &p_body, const RID &p_body_b) { (void)p_body; (void)p_body_b; }
void AVBDPhysicsServer3D::_soft_body_remove_collision_exception(const RID &p_body, const RID &p_body_b) { (void)p_body; (void)p_body_b; }
TypedArray<RID> AVBDPhysicsServer3D::_soft_body_get_collision_exceptions(const RID &p_body) const { (void)p_body; return {}; }
void AVBDPhysicsServer3D::_soft_body_set_total_mass(const RID &p_body, float p_total_mass) { (void)p_body; (void)p_total_mass; }
float AVBDPhysicsServer3D::_soft_body_get_total_mass(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_linear_stiffness(const RID &p_body, float p_linear_stiffness) { (void)p_body; (void)p_linear_stiffness; }
float AVBDPhysicsServer3D::_soft_body_get_linear_stiffness(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_pressure_coefficient(const RID &p_body, float p_pressure_coefficient) { (void)p_body; (void)p_pressure_coefficient; }
float AVBDPhysicsServer3D::_soft_body_get_pressure_coefficient(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_damping_coefficient(const RID &p_body, float p_damping_coefficient) { (void)p_body; (void)p_damping_coefficient; }
float AVBDPhysicsServer3D::_soft_body_get_damping_coefficient(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_drag_coefficient(const RID &p_body, float p_drag_coefficient) { (void)p_body; (void)p_drag_coefficient; }
float AVBDPhysicsServer3D::_soft_body_get_drag_coefficient(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_simulation_precision(const RID &p_body, int32_t p_simulation_precision) { (void)p_body; (void)p_simulation_precision; }
int32_t AVBDPhysicsServer3D::_soft_body_get_simulation_precision(const RID &p_body) const { (void)p_body; return 0; }
void AVBDPhysicsServer3D::_soft_body_set_shrinking_factor(const RID &p_body, float p_shrinking_factor) { (void)p_body; (void)p_shrinking_factor; }
float AVBDPhysicsServer3D::_soft_body_get_shrinking_factor(const RID &p_body) const { (void)p_body; return 0.0f; }
void AVBDPhysicsServer3D::_soft_body_set_state(const RID &p_body, PhysicsServer3D::BodyState p_state, const Variant &p_variant) { (void)p_body; (void)p_state; (void)p_variant; }
Variant AVBDPhysicsServer3D::_soft_body_get_state(const RID &p_body, PhysicsServer3D::BodyState p_state) const { (void)p_body; (void)p_state; return {}; }
void AVBDPhysicsServer3D::_soft_body_apply_central_force(const RID &p_body, const Vector3 &p_force) { (void)p_body; (void)p_force; }
void AVBDPhysicsServer3D::_soft_body_apply_central_impulse(const RID &p_body, const Vector3 &p_impulse) { (void)p_body; (void)p_impulse; }
void AVBDPhysicsServer3D::_soft_body_apply_point_force(const RID &p_body, int32_t p_point_index, const Vector3 &p_force) { (void)p_body; (void)p_point_index; (void)p_force; }
void AVBDPhysicsServer3D::_soft_body_apply_point_impulse(const RID &p_body, int32_t p_point_index, const Vector3 &p_impulse) { (void)p_body; (void)p_point_index; (void)p_impulse; }
void AVBDPhysicsServer3D::_soft_body_move_point(const RID &p_body, int32_t p_point_index, const Vector3 &p_global_position) { (void)p_body; (void)p_point_index; (void)p_global_position; }
void AVBDPhysicsServer3D::_soft_body_pin_point(const RID &p_body, int32_t p_point_index, bool p_pin) { (void)p_body; (void)p_point_index; (void)p_pin; }
bool AVBDPhysicsServer3D::_soft_body_is_point_pinned(const RID &p_body, int32_t p_point_index) const { (void)p_body; (void)p_point_index; return false; }
Vector3 AVBDPhysicsServer3D::_soft_body_get_point_global_position(const RID &p_body, int32_t p_point_index) const { (void)p_body; (void)p_point_index; return Vector3(); }
void AVBDPhysicsServer3D::_soft_body_remove_all_pinned_points(const RID &p_body) { (void)p_body; }
void AVBDPhysicsServer3D::_soft_body_set_ray_pickable(const RID &p_body, bool p_enable) { (void)p_body; (void)p_enable; }
AABB AVBDPhysicsServer3D::_soft_body_get_bounds(const RID &p_body) const { (void)p_body; return AABB(); }
void AVBDPhysicsServer3D::_soft_body_update_rendering_server(const RID &p_body,
        PhysicsServer3DRenderingServerHandler *p_rendering_server_handler) { (void)p_body; (void)p_rendering_server_handler; }

// Motion queries are implemented against the modern shapes' world AABBs rather than a
// full swept narrowphase: this is the MVP contract. It gets CharacterBody3D/move_and_slide
// a conservative, monotone answer (no tunnelling, correct stop distance for boxes and
// spheres); the exact shape-vs-shape sweep is Phase 3 work.
static void shape_world_aabb(const collision::Shape &p_shape, const Transform3D &p_xform,
        Vec3 &r_lo, Vec3 &r_hi) {
    const Vec3 half = p_shape.local_half_extents();
    const Mat3 R = rot_matrix(to_modern(p_xform.basis.get_rotation_quaternion()));
    const Vec3 world_half{std::abs(R[0][0]) * half.x + std::abs(R[1][0]) * half.y + std::abs(R[2][0]) * half.z,
            std::abs(R[0][1]) * half.x + std::abs(R[1][1]) * half.y + std::abs(R[2][1]) * half.z,
            std::abs(R[0][2]) * half.x + std::abs(R[1][2]) * half.y + std::abs(R[2][2]) * half.z};
    const Vec3 center = to_modern(p_xform.origin);
    r_lo = center - world_half;
    r_hi = center + world_half;
}

bool AVBDPhysicsServer3D::_body_test_motion(const RID &p_body, const Transform3D &p_from, const Vector3 &p_motion,
        float p_margin, int32_t p_max_collisions, bool p_collide_separation_ray, bool p_recovery_as_collision,
        PhysicsServer3DExtensionMotionResult *p_result) const {
    (void)p_collide_separation_ray;
    (void)p_recovery_as_collision;

    if (p_result != nullptr) {
        p_result->travel = p_motion;
        p_result->remainder = Vector3();
        p_result->collision_depth = 0.0f;
        p_result->collision_safe_fraction = 1.0f;
        p_result->collision_unsafe_fraction = 1.0f;
        p_result->collision_count = 0;
    }

    const auto moving_it = bodies.find(id_of(p_body));
    if (moving_it == bodies.end() || !moving_it->second.space.is_valid()) {
        return false;
    }
    const BodyData &moving = moving_it->second;
    const uint64_t moving_id = id_of(p_body);
    const SpaceData *space = find_space(moving.space);
    if (space == nullptr || space->modern == nullptr) {
        return false;
    }

    if (moving.shapes.empty()) {
        return false;
    }
    collision::Shape base;
    if (!query_shape_from_data(moving.shapes.front(), p_from, base)) {
        return false;
    }

    Vec3 base_lo, base_hi;
    shape_world_aabb(base, p_from, base_lo, base_hi);
    const Vec3 motion = to_modern(p_motion);
    const Scalar_t margin = std::max(0.0f, p_margin);

    struct Hit {
        Scalar_t safe = 1;
        Scalar_t depth = 0;
        Vector3 point;
        Vector3 normal;
        uint64_t body_id = 0;
        uint64_t instance_id = 0;
        Vector3 linear_velocity;
        Vector3 angular_velocity;
    };
    std::vector<Hit> hits;

    const uint64_t space_id = id_of(moving.space);
    for (auto &[other_id, other] : bodies) {
        if (other_id == moving_id || id_of(other.space) != space_id) {
            continue;
        }
        if ((other.collision_layer & moving.collision_mask) == 0) {
            continue;
        }
        if (body_is_excepted(moving_id, other_id)) {
            continue;
        }
        for (std::size_t slot = 0; slot < other.shapes.size(); ++slot) {
            if (slot < other.shape_disabled.size() && other.shape_disabled[slot]) {
                continue;
            }
            collision::Shape shape;
            if (!query_shape_from_data(other.shapes[slot], body_shape_transform(other, static_cast<int32_t>(slot)), shape)) {
                continue;
            }
            Vec3 obs_lo, obs_hi;
            shape_world_aabb(shape, body_shape_transform(other, static_cast<int32_t>(slot)), obs_lo, obs_hi);
            obs_lo -= Vec3{margin};
            obs_hi += Vec3{margin};

            // Swept AABB (Minkowski): the interval of t where the moving box overlaps.
            Scalar_t t_enter = 0;
            Scalar_t t_exit = 1;
            int enter_axis = -1;
            Scalar_t enter_sign = 0;
            bool separated = false;
            for (int axis = 0; axis < 3; ++axis) {
                const Scalar_t m = motion[axis];
                if (std::abs(m) < Scalar_t(1e-9)) {
                    if (base_hi[axis] < obs_lo[axis] || base_lo[axis] > obs_hi[axis]) {
                        separated = true;
                        break;
                    }
                    continue;
                }
                Scalar_t t1, t2;
                if (m > 0) {
                    t1 = (obs_lo[axis] - base_hi[axis]) / m;
                    t2 = (obs_hi[axis] - base_lo[axis]) / m;
                } else {
                    t1 = (obs_hi[axis] - base_lo[axis]) / m;
                    t2 = (obs_lo[axis] - base_hi[axis]) / m;
                }
                if (t1 > t_enter) {
                    t_enter = t1;
                    enter_axis = axis;
                    enter_sign = m > 0 ? Scalar_t(-1) : Scalar_t(1);
                }
                t_exit = std::min(t_exit, t2);
            }
            if (separated || enter_axis < 0 || t_enter > t_exit || t_enter > 1 || t_exit < 0) {
                continue;
            }
            Hit hit;
            hit.safe = std::max<Scalar_t>(t_enter, 0);
            hit.depth = std::max<Scalar_t>(0, std::min(base_hi[enter_axis], obs_hi[enter_axis]) -
                                                  std::max(base_lo[enter_axis], obs_lo[enter_axis]));
            Vec3 point = to_modern(p_from.origin);
            point[enter_axis] += motion[enter_axis] * hit.safe;
            hit.point = to_godot(point);
            Vec3 normal{0};
            normal[enter_axis] = enter_sign;
            hit.normal = to_godot(normal);
            hit.body_id = other_id;
            hit.instance_id = other.instance_id;
            hit.linear_velocity = other.linear_velocity;
            hit.angular_velocity = other.angular_velocity;
            hits.push_back(hit);
        }
    }

    if (hits.empty()) {
        return false;
    }
    std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) { return a.safe < b.safe; });

    const Hit &nearest = hits.front();
    if (p_result != nullptr) {
        p_result->travel = p_motion * nearest.safe;
        p_result->remainder = p_motion * (1.0f - nearest.safe);
        p_result->collision_depth = nearest.depth;
        p_result->collision_safe_fraction = nearest.safe;
        p_result->collision_unsafe_fraction = nearest.safe;
        const int32_t max_collisions = std::clamp(p_max_collisions, 0, 32);
        const int32_t count = std::min<int32_t>(static_cast<int32_t>(hits.size()), max_collisions);
        p_result->collision_count = count;
        for (int32_t i = 0; i < count; i++) {
            const Hit &hit = hits[static_cast<std::size_t>(i)];
            PhysicsServer3DExtensionMotionCollision &collision = p_result->collisions[i];
            collision.position = hit.point;
            collision.normal = hit.normal;
            collision.collider_velocity = hit.linear_velocity;
            collision.collider_angular_velocity = hit.angular_velocity;
            collision.depth = hit.depth;
            collision.local_shape = 0;
            collision.collider_id = hit.instance_id;
            collision.collider = make_rid(hit.body_id);
            collision.collider_shape = 0;
        }
    }
    return nearest.safe < 1.0f || nearest.depth > 0.0f;
}

void AVBDPhysicsServer3D::_shape_set_data(const RID &p_shape, const Variant &p_data) {
    ShapeData *shape = find_shape(p_shape);
    if (shape == nullptr) {
        return;
    }
    // Godot 4.x payload per shape type (scene/resources/3d/*_shape_3d.cpp::_update_shape):
    //   box      Vector3 - HALF extents (the resource sends `size / 2`)
    //   sphere   float radius
    //   cylinder Dictionary { radius: float, height: float } - full height
    switch (shape->type) {
        case PhysicsServer3D::SHAPE_SPHERE: {
            if (p_data.get_type() == Variant::FLOAT) {
                shape->radius = p_data;
            } else if (p_data.get_type() == Variant::DICTIONARY && ((const Dictionary)p_data).has("radius")) {
                shape->radius = ((const Dictionary)p_data)["radius"];
            }
            break;
        }
        case PhysicsServer3D::SHAPE_CYLINDER: {
            if (p_data.get_type() == Variant::DICTIONARY) {
                const Dictionary data = p_data;
                if (data.has("radius")) {
                    shape->radius = data["radius"];
                }
                if (data.has("height")) {
                    shape->height = data["height"];
                }
            }
            break;
        }
        case PhysicsServer3D::SHAPE_CAPSULE: {
            if (p_data.get_type() == Variant::DICTIONARY) {
                const Dictionary data = p_data;
                if (data.has("radius")) {
                    shape->radius = data["radius"];
                }
                if (data.has("height")) {
                    shape->height = data["height"];
                }
            }
            break;
        }
        case PhysicsServer3D::SHAPE_CONVEX_POLYGON: {
            // Proxy contract: the convex hull is simulated as the local AABB of its
            // points, stored as full extents (the solver rotates it as an OBB).
            if (p_data.get_type() == Variant::PACKED_VECTOR3_ARRAY) {
                const PackedVector3Array points = p_data;
                Vector3 lo = points.is_empty() ? Vector3() : points[0];
                Vector3 hi = lo;
                for (const Vector3 &pt : points) {
                    lo = lo.min(pt);
                    hi = hi.max(pt);
                }
                shape->extents = hi - lo;
            }
            break;
        }
        default: {
            if (p_data.get_type() == Variant::VECTOR3) {
                // Half extents in, full extents stored (the resource's `size / 2`).
                const Vector3 half = p_data;
                shape->extents = half * 2.0;
            } else if (p_data.get_type() == Variant::DICTIONARY && ((const Dictionary)p_data).has("size")) {
                shape->extents = ((const Dictionary)p_data)["size"];
            }
            break;
        }
    }
    // Round-trip for _shape_get_data (the engine reads back what it sent).
    shape->data = p_data;
}

PhysicsServer3D::ShapeType AVBDPhysicsServer3D::_shape_get_type(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? PhysicsServer3D::SHAPE_BOX : found->second.type;
}

Variant AVBDPhysicsServer3D::_shape_get_data(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    if (found == shapes.end()) {
        return Variant();
    }
    const ShapeData &shape = found->second;

    // The exact form the engine's resources send: box = half extents Vector3,
    // sphere = radius float, cylinder = Dictionary.
    switch (shape.type) {
        case PhysicsServer3D::SHAPE_SPHERE:
            return Variant(shape.radius);
        case PhysicsServer3D::SHAPE_CYLINDER: {
            Dictionary data;
            data["radius"] = shape.radius;
            data["height"] = shape.height;
            return data;
        }
        case PhysicsServer3D::SHAPE_CAPSULE: {
            Dictionary data;
            data["radius"] = shape.radius;
            data["height"] = shape.height;
            return data;
        }
        default:
            // Proxied and excluded shapes round-trip the engine's own payload.
            return shape.data;
    }
}

static bool modern_shape_from_data(const AVBDPhysicsServer3D::ShapeData &p_data,
        collision::Shape &r_shape);

bool AVBDPhysicsServer3D::query_shape_from_data(const ShapeData &p_data, const Transform3D &p_transform,
        collision::Shape &r_shape) {
    // The transform is carried by the query alongside the shape; the modern Shape itself
    // is purely local geometry.
    (void)p_transform;
    return modern_shape_from_data(p_data, r_shape);
}

// -----------------------------------------------------------------------------
// Shape conversion + space rebuild
// -----------------------------------------------------------------------------

// Godot shape record -> modern collision shape. The map is deliberately narrow:
// box/sphere/cylinder are the MVP set, capsule is simulated as the same cylinder the
// legacy server used, and convex points become their local AABB box (a proxy, not the
// hull). Everything else is refused by the caller.
static bool modern_shape_from_data(const AVBDPhysicsServer3D::ShapeData &p_data,
        collision::Shape &r_shape) {
    switch (p_data.type) {
        case PhysicsServer3D::SHAPE_BOX: {
            // ShapeData.extents holds Godot's full extents for boxes.
            const Vec3 half{static_cast<Scalar_t>(p_data.extents.x * 0.5),
                    static_cast<Scalar_t>(p_data.extents.y * 0.5),
                    static_cast<Scalar_t>(p_data.extents.z * 0.5)};
            r_shape = collision::Shape::box(half);
            return true;
        }
        case PhysicsServer3D::SHAPE_SPHERE:
            r_shape = collision::Shape::sphere(static_cast<Scalar_t>(p_data.radius));
            return p_data.radius > 0.0;
        case PhysicsServer3D::SHAPE_CYLINDER:
            // Godot's cylinder height is full height; the modern shape takes it the same way
            // and both axes are local +Y, so no permutation.
            r_shape = collision::Shape::cylinder(static_cast<Scalar_t>(p_data.radius),
                    static_cast<Scalar_t>(p_data.height));
            return p_data.radius > 0.0 && p_data.height > 0.0;
        case PhysicsServer3D::SHAPE_CAPSULE:
            // Documented proxy: a cylinder over the capsule's radius and full height.
            r_shape = collision::Shape::cylinder(static_cast<Scalar_t>(p_data.radius),
                    static_cast<Scalar_t>(p_data.height));
            return p_data.radius > 0.0 && p_data.height > 0.0;
        case PhysicsServer3D::SHAPE_CONVEX_POLYGON: {
            // Documented proxy: the points' local AABB as a box.
            const Vec3 half{static_cast<Scalar_t>(p_data.extents.x * 0.5),
                    static_cast<Scalar_t>(p_data.extents.y * 0.5),
                    static_cast<Scalar_t>(p_data.extents.z * 0.5)};
            if (half.x <= 0 || half.y <= 0 || half.z <= 0)
                return false;
            r_shape = collision::Shape::box(half);
            return true;
        }
        default:
            return false;
    }
}

void AVBDPhysicsServer3D::rebuild_space(const RID &p_space, std::source_location where) {
    SpaceData *space = find_space(p_space);
    if (space == nullptr) {
        return;
    }

    const uint64_t space_id = id_of(p_space);
    ++rebuild_total;
    const std::string reason = std::string(where.function_name()) + ":" + std::to_string(where.line());
    ++rebuild_reasons[reason];
    rebuild_last_step = step_total;
    rebuild_last_reason = String::utf8(reason.c_str());

    for (auto &[id, body] : bodies) {
        (void)id;
        if (id_of(body.space) == space_id) {
            body.modern_index = -1;
        }
    }

    if (space->modern == nullptr) {
        space->modern = std::make_unique<ModernSpace>();
    } else {
        space->modern->clear();
    }

    // Deterministic order: ids are minted in creation order, so sorting by id makes the
    // state index - and with it the Gauss-Seidel order and every contact digest - a
    // function of the scene rather than of unordered_map bucket layout.
    std::vector<uint64_t> space_body_ids;
    for (const auto &[id, body] : bodies) {
        if (id_of(body.space) == space_id) {
            space_body_ids.push_back(id);
        }
    }
    std::ranges::sort(space_body_ids);

    std::vector<ModernBodyDesc> descs;
    descs.reserve(space_body_ids.size());

    for (const uint64_t id : space_body_ids) {
        BodyData &body = bodies.find(id)->second;

        // Refresh the shape snapshots from the engine's records: resource data may have
        // arrived after the shape was attached to the body.
        for (std::size_t i = 0; i < body.shapes.size() && i < body.shape_rids.size(); i++) {
            if (const ShapeData *fresh = shape_data(body.shape_rids[i])) {
                body.shapes[i] = *fresh;
            }
        }
        if (body.shape_disabled.size() < body.shapes.size()) {
            body.shape_disabled.resize(body.shapes.size(), false);
        }

        ModernBodyDesc desc;
        desc.id = id;
        const bool immovable_mode = body.mode == PhysicsServer3D::BODY_MODE_STATIC ||
                body.mode == PhysicsServer3D::BODY_MODE_KINEMATIC;
        desc.is_static = immovable_mode || !(body.mass > 0.0);
        desc.mass = body.mass;
        desc.position = to_modern(body.transform.origin);
        desc.rotation = to_modern(body.transform.basis.get_rotation_quaternion());
        desc.linear_velocity = to_modern(body.linear_velocity);
        desc.angular_velocity = to_modern(body.angular_velocity);
        desc.layer = body.collision_layer;
        desc.mask = body.collision_mask;
        desc.axis_locks = body.axis_locks;
        desc.exceptions = body.exceptions;

        bool unsupported = false;
        for (std::size_t i = 0; i < body.shapes.size(); ++i) {
            if (body.shape_disabled[i]) {
                continue;
            }
            collision::Shape shape;
            if (!modern_shape_from_data(body.shapes[i], shape)) {
                unsupported = true;
                break;
            }
            // Friction is carried by the shape material (the core combines the two
            // bodies' values as sqrt(fA*fB), Godot's own rule). The engine stores the
            // per-body value in BodyData, so pass it through here.
            shape.material.friction = static_cast<Scalar_t>(body.friction);
            collision::ShapeInstance inst;
            inst.shape = shape;
            const Transform3D &local = i < body.shape_transforms.size()
                    ? body.shape_transforms[i]
                    : Transform3D();
            inst.local_position = to_modern(local.origin);
            inst.local_rotation = to_modern(local.basis.get_rotation_quaternion());
            desc.shapes.push_back(inst);
        }

        if (unsupported) {
            if (!body.unsupported_warned) {
                body.unsupported_warned = true;
                WARN_PRINT("AVBD physics server: a body uses a shape AVBD cannot simulate "
                           "(concave/heightmap/world boundary/separation ray/custom); the body "
                           "is excluded from simulation. Use box/sphere/cylinder shapes.");
            }
            continue;
        }
        if (desc.shapes.empty()) {
            // No shape yet: the engine attaches one before it matters, and until then the
            // body has neither geometry nor inertia.
            continue;
        }

        descs.push_back(std::move(desc));
    }

    space->modern->rebuild(std::move(descs));

    // Map indices back onto the records and stamp per-body gravity/damping, which the
    // rebuild does not take from the desc (they are engine parameters, not geometry).
    for (std::size_t i = 0; i < space->modern->size(); ++i) {
        const uint64_t id = space->modern->body_id(i);
        const auto found = bodies.find(id);
        if (found == bodies.end()) {
            continue;
        }
        BodyData &body = found->second;
        body.modern_index = static_cast<int32_t>(i);
        space->modern->set_gravity_scale(i, body.gravity_scale);
        space->modern->set_damping(i, body.linear_damp, body.angular_damp);
    }

    // Joints last: they need the body indices the loop above just wrote.
    rebuild_joints_for_space(*space, space_id);
}

// -----------------------------------------------------------------------------
// Joints
//
// Every engine joint is translated into modern JointAxisConstraint rows:
//   Pin      - 3 linear locked, angular free.
//   Hinge    - 3 linear locked; angular free about the hinge axis (joint frame Z),
//              locked about the other two. HINGE_JOINT_FLAG_USE_LIMIT turns the hinge
//              axis into a Limited row (Godot's hinge angle is B-relative-to-A, which is
//              exactly what the row measures).
//   6DOF     - each axis Free / Locked / Limited from its enable flag and limits. Godot
//              spells "locked" as an enabled limit with lower == upper.
//   Slider   - linear/angular about the joint frame X; the orthogonal axes locked.
//   ConeTwist- not mapped yet (documented gap); it has no row.
// Spring and motor parameters are still record-only: one row carries one (K, lambda, C)
// pair, so a hard constraint and a spring cannot share an axis in this representation.
//
// A joint edit calls rebuild_space, which rebuilds the body set and then reinstalls the
// rows from the records; that is what makes parameter changes take effect.
// -----------------------------------------------------------------------------

namespace {

double joint_param(const AVBDPhysicsServer3D::JointData &p_joint, int p_key) {
    const auto found = p_joint.joint_params.find(p_key);
    return found == p_joint.joint_params.end() ? 0.0 : found->second;
}

bool joint_flag(const AVBDPhysicsServer3D::JointData &p_joint, int p_key) {
    const auto found = p_joint.joint_flags.find(p_key);
    return found != p_joint.joint_flags.end() && found->second;
}

// Godot's generic-6DOF convention: an enabled limit with lower == upper means "locked"
// (the legacy mapping used the same rule).
ModernJointAxis axis_from_limit(bool p_enabled, double p_lower, double p_upper) {
    ModernJointAxis axis;
    if (!p_enabled) {
        axis.mode = JointAxisMode::Free;
        return axis;
    }
    if (std::fabs(p_upper - p_lower) < 1.0e-9) {
        axis.mode = JointAxisMode::Locked;
        return axis;
    }
    axis.mode = JointAxisMode::Limited;
    axis.lower = static_cast<Scalar_t>(p_lower);
    axis.upper = static_cast<Scalar_t>(p_upper);
    return axis;
}

Quat basis_rotation(const Basis &p_basis) {
    return to_modern(p_basis.get_rotation_quaternion());
}

} // namespace

void AVBDPhysicsServer3D::rebuild_joints_for_space(SpaceData &p_space, uint64_t p_space_id) {
    if (p_space.modern == nullptr) {
        return;
    }

    std::vector<uint64_t> space_joint_ids;
    space_joint_ids.reserve(joints.size());
    for (const auto &[id, joint] : joints) {
        (void)joint;
        space_joint_ids.push_back(id);
    }
    std::ranges::sort(space_joint_ids);

    std::vector<ModernJointDesc> descs;
    for (const uint64_t joint_id : space_joint_ids) {
        const JointData &joint = joints.find(joint_id)->second;

        const BodyData *body_a = joint.body_a.is_valid() ? find_body(joint.body_a) : nullptr;
        const BodyData *body_b = find_body(joint.body_b);
        if (body_b == nullptr || body_b->modern_index < 0 || id_of(body_b->space) != p_space_id) {
            continue;
        }
        if (joint.body_a.is_valid() &&
                (body_a == nullptr || body_a->modern_index < 0 ||
                        id_of(body_a->space) != p_space_id)) {
            continue;
        }

        ModernJointDesc desc;
        desc.id = joint_id;
        desc.body_a = joint.body_a.is_valid() ? body_a->modern_index : -1;
        desc.body_b = body_b->modern_index;
        desc.anchor_a = to_modern(joint.ref_a.origin);
        desc.anchor_b = to_modern(joint.ref_b.origin);
        desc.frame_a = basis_rotation(joint.ref_a.basis);
        desc.frame_b = basis_rotation(joint.ref_b.basis);

        for (int axis = 0; axis < 3; ++axis) {
            desc.linear[axis] = ModernJointAxis{}; // Free
            desc.angular[axis] = ModernJointAxis{};
        }

        switch (joint.type) {
            case PhysicsServer3D::JOINT_TYPE_PIN: {
                // Godot's pin is a ball joint: the anchor is fixed, rotation is free.
                for (int axis = 0; axis < 3; ++axis)
                    desc.linear[axis].mode = JointAxisMode::Locked;
                break;
            }
            case PhysicsServer3D::JOINT_TYPE_HINGE: {
                for (int axis = 0; axis < 3; ++axis)
                    desc.linear[axis].mode = JointAxisMode::Locked;
                desc.angular[0].mode = JointAxisMode::Locked;
                desc.angular[1].mode = JointAxisMode::Locked;
                const bool use_limit = joint_flag(joint, PhysicsServer3D::HINGE_JOINT_FLAG_USE_LIMIT);
                const double upper = joint_param(joint, PhysicsServer3D::HINGE_JOINT_LIMIT_UPPER);
                const double lower = joint_param(joint, PhysicsServer3D::HINGE_JOINT_LIMIT_LOWER);
                desc.angular[2] = axis_from_limit(use_limit, lower, upper);
                break;
            }
            case PhysicsServer3D::JOINT_TYPE_SLIDER: {
                const double lin_upper =
                        joint_param(joint, PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_UPPER);
                const double lin_lower =
                        joint_param(joint, PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_LOWER);
                desc.linear[0] = axis_from_limit(true, lin_lower, lin_upper);
                desc.linear[1].mode = JointAxisMode::Locked;
                desc.linear[2].mode = JointAxisMode::Locked;

                const double ang_upper =
                        joint_param(joint, PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_UPPER);
                const double ang_lower =
                        joint_param(joint, PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_LOWER);
                desc.angular[0] = axis_from_limit(true, ang_lower, ang_upper);
                desc.angular[1].mode = JointAxisMode::Locked;
                desc.angular[2].mode = JointAxisMode::Locked;
                break;
            }
            case PhysicsServer3D::JOINT_TYPE_6DOF: {
                for (int axis = 0; axis < 3; ++axis) {
                    const double lin_lower =
                            joint.params[axis][PhysicsServer3D::G6DOF_JOINT_LINEAR_LOWER_LIMIT];
                    const double lin_upper =
                            joint.params[axis][PhysicsServer3D::G6DOF_JOINT_LINEAR_UPPER_LIMIT];
                    const bool lin_enabled = joint.flags[axis]
                            [PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_LINEAR_LIMIT];
                    desc.linear[axis] = axis_from_limit(lin_enabled, lin_lower, lin_upper);

                    // Angular params are stored with the angular enum family rebased to 0
                    // (see _generic_6dof_joint_set_param): the angular lower/upper slots
                    // are 0/1 on row 3+axis, not the enum's 10/11.
                    constexpr int kAngularParamBase =
                            PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT;
                    const double ang_lower =
                            joint.params[axis + 3][PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT -
                                                   kAngularParamBase];
                    const double ang_upper =
                            joint.params[axis + 3][PhysicsServer3D::G6DOF_JOINT_ANGULAR_UPPER_LIMIT -
                                                   kAngularParamBase];
                    const bool ang_enabled = joint.flags[axis + 3]
                            [PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT];
                    desc.angular[axis] = axis_from_limit(ang_enabled, ang_lower, ang_upper);
                }
                break;
            }
            case PhysicsServer3D::JOINT_TYPE_CONE_TWIST:
            default:
                // Cone-twist has no row yet: its swing cone is not a per-axis interval,
                // and shipping a wrong approximation is worse than no constraint.
                continue;
        }

        descs.push_back(std::move(desc));
    }

    p_space.modern->rebuild_joints(std::move(descs));
}

void AVBDPhysicsServer3D::_rebuild_for_joint(const RID &p_joint) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    const JointData &joint = found->second;
    if (BodyData *body = find_body(joint.body_b)) {
        rebuild_space(body->space);
    } else if (BodyData *body = find_body(joint.body_a)) {
        rebuild_space(body->space);
    }
}

// Run every active space for one step, then tell each body what happened.
//
// Stepping follows the game: while the engine has the server inactive (SceneTree paused)
// nothing steps. The modern core owns the actual integration; this function only turns
// engine records into per-tick force injection, drives one ModernSpace step, and copies
// the result back into the host records.
void AVBDPhysicsServer3D::_step(float p_step) {
    if (!server_active) {
        return;
    }
    ++step_total;
    diag_max_penetration = 0.0;
    diag_deepest_a = 0;
    diag_deepest_b = 0;

    // One state object serves every callback this step; re-aimed per body before each call.
    if (body_state == nullptr) {
        body_state = memnew(AVBDDirectBodyState3D);
    }

    StepProfile &prof = step_profile();
    const bool profiling = prof.enabled;
    const ProfClock::time_point prof_t0 = ProfClock::now();
    double prof_rebuild = 0.0;
    double prof_params = 0.0;
    double prof_forces = 0.0;
    double prof_solver = 0.0;
    double prof_contacts = 0.0;
    int prof_threads = 1;
    int prof_iterations = 0;
    int prof_rebuild_count = 0;
    bool stepped_any = false;

    for (auto &[id, space_ptr] : spaces) {
        SpaceData &space = *space_ptr;
        if (!space.active) {
            continue;
        }

        // A body joined/left this space since the last step; rebuild now that the engine
        // has finished sending its setup (shapes, mass, material).
        if (space_rebuild_pending.erase(id) > 0) {
            const ProfClock::time_point prof_t_rebuild = ProfClock::now();
            rebuild_space(make_rid(id));
            prof_rebuild += seconds_since(prof_t_rebuild);
            prof_rebuild_count++;
        }

        if (space.modern == nullptr || space.modern->empty()) {
            continue;
        }

        const ProfClock::time_point prof_t_params = ProfClock::now();
        if (space.params_dirty) {
            read_project_params(space);
            read_project_gravity(space);
            space.params_dirty = false;
        }
        prof_params += seconds_since(prof_t_params);

        const int iterations = space.iterations > 0 ? space.iterations : 8;
        // Policy selection first: set_policy_stiffness is type-aware, and switching the
        // solver resets dual state (no force<->multiplier cross-injection). set_solver is
        // a cheap no-op when the kind has not changed.
        space.modern->set_solver(
                space.solver == 1 ? SolverKind::BPXPBD : SolverKind::AVBD);
        // Contact stiffness. The legacy beta_linear setting was the per-contact penalty
        // scale; the modern AVBD policy exposes exactly that as its per-step k0.
        // B-PXPBD uses zeta/contact-compliance instead (see ModernSpace).
        space.modern->set_policy_stiffness(
                static_cast<Scalar_t>(std::clamp(space.beta_linear, 1.0e4, 1.0e6)), Scalar_t(1e-2),
                static_cast<Scalar_t>(std::clamp(space.bpxpbd_zeta, 0.0, 1.0)),
                static_cast<Scalar_t>(std::max(space.bpxpbd_contact_compliance, 0.0)));
        // Velocity-neutral position projection. The legacy projection_* settings were
        // dead after the modern-core migration; they now drive B-PXPBD's step-end
        // projection. AVBD ignores them (its implicit contact already depenetrates),
        // so the AVBD path stays bit-identical.
        space.modern->set_contact_projection(space.position_projection,
                space.projection_iterations,
                static_cast<Scalar_t>(std::clamp(space.projection_beta, 0.0, 1.0)),
                static_cast<Scalar_t>(std::max(space.projection_slop, 0.0)),
                static_cast<Scalar_t>(std::max(space.projection_max_correction, 0.0)));
        space.modern->set_contact_margin(Scalar_t(0.002));

        // Godot's physics tick is the sampling rate. If an internal rate is configured,
        // accumulate engine time and advance the solver on its own fixed tick; the engine
        // tick then only samples the state. 0 keeps one solver tick per engine tick.
        // Substepping (`physics/avbd/substeps`) is folded into this: the modern core's
        // per-iteration solve replaces the legacy round split, and running the whole
        // collision+force pass `substeps` times per tick would only multiply cost.
        double dt = p_step;
        int tickCount = 1;
        if (space.internal_ticks_per_second > 0.0) {
            dt = 1.0 / space.internal_ticks_per_second;
            space.internal_accumulator += p_step;
            tickCount = 0;
            while (space.internal_accumulator + 1.0e-9 >= dt && tickCount < space.max_internal_steps) {
                space.internal_accumulator -= dt;
                ++tickCount;
            }
            if (space.internal_accumulator >= dt) {
                space.internal_accumulator = std::fmod(space.internal_accumulator, dt);
            }
        }

        ModernSpace &modern = *space.modern;
        // Contact generation distance from project settings. Only push it when it
        // changes: set_speculative_distance is cheap, but the value participates in
        // the collision broadphase/narrowphase options.
        const Scalar_t want_speculative = static_cast<Scalar_t>(
                std::max(space.speculative_margin, 0.0));
        if (modern.speculative_distance() != want_speculative) {
            modern.set_speculative_distance(want_speculative);
        }
        for (int tick = 0; tick < tickCount; ++tick) {
            const Scalar_t tickDt = static_cast<Scalar_t>(dt);

            const ProfClock::time_point prof_t_forces = ProfClock::now();
            for (auto &[body_id, body] : bodies) {
                (void)body_id;
                if (id_of(body.space) == id && body.modern_index >= 0) {
                    modern.set_sleeping(static_cast<std::size_t>(body.modern_index), body.sleeping);
                }
            }
            modern.begin_tick(space.gravity, tickDt);

            // Forces accumulated since the last solver tick (persistent constant_* plus
            // one-shot frame ones). The modern core converts a force into acceleration
            // through the body's inverse mass, exactly like the legacy path.
            for (auto &[body_id, body] : bodies) {
                (void)body_id;
                if (id_of(body.space) != id || body.modern_index < 0) {
                    continue;
                }
                if (body.omit_force_integration) {
                    // A body with scripted integration applies its own forces in its callback.
                    body.frame_force = Vector3();
                    body.frame_torque = Vector3();
                    continue;
                }
                const Vector3 force = body.constant_force + body.frame_force;
                const Vector3 torque = body.constant_torque + body.frame_torque;
                body.frame_force = Vector3();
                body.frame_torque = Vector3();
                if (force.length_squared() > 0.0) {
                    modern.add_force(static_cast<std::size_t>(body.modern_index), to_modern(force));
                }
                if (torque.length_squared() > 0.0) {
                    modern.add_torque(static_cast<std::size_t>(body.modern_index), to_modern(torque));
                }
            }
            prof_forces += seconds_since(prof_t_forces);

            const ProfClock::time_point prof_t_solver = ProfClock::now();
            modern.step(iterations, space.threads);
            prof_solver += seconds_since(prof_t_solver);
        }
        if (tickCount > 0) {
            stepped_any = true;
        }

        // --- diagnostics: deepest contact of the step -----------------------------
        // ManifoldContact.depth is the geometric overlap measured by narrowphase at the
        // pose the last solve worked against. Mapping back through the space's stable
        // index gives the Godot instance ids the headless stress suites resolve.
        for (const collision::Manifold &m : modern.manifolds()) {
            for (std::size_t c = 0; c < m.count; ++c) {
                const double depth = static_cast<double>(m.contacts[c].depth);
                if (depth <= diag_max_penetration) {
                    continue;
                }
                diag_max_penetration = depth;
                const auto rec_a = bodies.find(modern.body_id(m.body_a));
                const auto rec_b = bodies.find(modern.body_id(m.body_b));
                diag_deepest_a = rec_a != bodies.end() ? rec_a->second.instance_id : 0;
                diag_deepest_b = rec_b != bodies.end() ? rec_b->second.instance_id : 0;
            }
        }

        // Area monitoring runs once per engine step, after the solver's output is in hand
        // and before the scene reads anything back.
        _monitor_areas();
        const ProfClock::time_point prof_t_contacts = ProfClock::now();

        // --- per-body contact reporting (RigidBody3D contact monitoring) -----------
        bool any_reports = false;
        for (auto &[collect_id, collect_body] : bodies) {
            (void)collect_id;
            if (id_of(collect_body.space) != id) {
                continue;
            }
            if (collect_body.max_contacts_reported > 0) {
                collect_body.frame_contacts.clear();
                any_reports = true;
            }
        }
        if (any_reports) {
            const State &state = modern.state();
            for (const collision::Manifold &m : modern.manifolds()) {
                if (m.count == 0) {
                    continue;
                }
                const auto found_a = bodies.find(modern.body_id(m.body_a));
                const auto found_b = bodies.find(modern.body_id(m.body_b));
                if (found_a == bodies.end() || found_b == bodies.end()) {
                    continue;
                }
                BodyData *data[2] = {&found_a->second, &found_b->second};
                const std::size_t index[2] = {m.body_a, m.body_b};
                const Vec3 normal = m.basis.normal;

                for (int side = 0; side < 2; side++) {
                    BodyData *fill = data[side];
                    if (fill->max_contacts_reported <= 0) {
                        continue;
                    }
                    const std::size_t self_index = index[side];
                    for (std::size_t c = 0;
                            c < m.count &&
                            static_cast<int32_t>(fill->frame_contacts.size()) < fill->max_contacts_reported;
                            ++c) {
                        const collision::ManifoldContact &contact = m.contacts[c];
                        if (static_cast<double>(contact.depth) < fill->contacts_depth_threshold) {
                            continue;
                        }
                        const Vec3 r = side == 0 ? contact.rA : contact.rB;
                        const Vec3 world = state.pose.position(self_index) +
                                rotate(state.pose.rotation(self_index), r);
                        BodyData::FrameContact fc;
                        fc.position = to_godot(world);
                        // basis.normal points from B to A; each side reports the normal as
                        // pointing out of the other body towards itself.
                        fc.normal = to_godot(normal * (side == 0 ? Scalar_t(1) : Scalar_t(-1)));
                        fc.impulse = to_godot(normal * (contact.lambda.x *
                                (side == 0 ? Scalar_t(1) : Scalar_t(-1))));
                        fc.collider_body_id = modern.body_id(index[side == 0 ? 1 : 0]);
                        // The manifold stores proxy ids; the slot needs the shape table,
                        // which the MVP contact path does not expose. Single-shape bodies
                        // (the common case) report slot 0.
                        fc.local_shape = 0;
                        fc.collider_shape = 0;
                        fill->frame_contacts.push_back(fc);
                    }
                }
            }
        }
        prof_contacts += seconds_since(prof_t_contacts);

        // --- write back ------------------------------------------------------------
        for (auto &[body_id, body] : bodies) {
            (void)body_id;
            if (id_of(body.space) != id || body.modern_index < 0) {
                continue;
            }
            const std::size_t i = static_cast<std::size_t>(body.modern_index);
            Transform3D transform;
            transform.origin = to_godot(modern.body_origin(i));
            transform.basis = Basis(to_godot(modern.body_rotation(i)).normalized());
            body.transform = transform;
            body.linear_velocity = to_godot(modern.linear_velocity(i));
            body.angular_velocity = to_godot(modern.angular_velocity(i));
            // body.sleeping is owned by the engine (explicit sleep only; the modern core
            // has no automatic idle detection yet).
        }

        prof_threads = std::max(prof_threads, modern.threads_used());
        prof_iterations = std::max(prof_iterations, iterations);
    }

    if (stepped_any) {
        has_stepped = true;
    }

    if (profiling) {
        const double prof_total = seconds_since(prof_t0);
        prof.total_ns += prof_total;
        prof.solver_ns += prof_solver;
        prof.contacts_ns += prof_contacts;
        prof.rebuild_ns += prof_rebuild;
        prof.params_ns += prof_params;
        prof.forces_ns += prof_forces;
        prof.max_total_ns = std::max(prof.max_total_ns, prof_total);
        prof.max_solver_ns = std::max(prof.max_solver_ns, prof_solver);
        prof.max_contacts_ns = std::max(prof.max_contacts_ns, prof_contacts);
        prof.max_rebuild_ns = std::max(prof.max_rebuild_ns, prof_rebuild);
        prof.max_params_ns = std::max(prof.max_params_ns, prof_params);
        prof.max_forces_ns = std::max(prof.max_forces_ns, prof_forces);
        prof.rebuilds += prof_rebuild_count;
        if (prof_total > 2.0e6)
            prof.over_2ms_ticks++;
        prof.threads_used = std::max(prof.threads_used, prof_threads);
        prof.iterations = prof_iterations;
        prof.ticks++;
        if (prof.ticks % 240 == 0) {
            constexpr double to_us = 0.001;
            std::fprintf(stdout,
                    "[avbd-profile] ticks=%llu total=%.1fus max=%.1fus solver=%.1fus(max %.1fus) sync=%.1fus(max %.1fus) "
                    "contacts=%.1fus rebuild=%.1fus(max %.1fus x%u) params=%.1fus forces=%.1fus over2ms=%u "
                    "threads=%d substeps=%d iterations=%d\n",
                    static_cast<unsigned long long>(prof.ticks), prof.total_ns * to_us,
                    prof.max_total_ns * to_us,
                    prof.solver_ns * to_us, prof.max_solver_ns * to_us,
                    prof.sync_ns * to_us, prof.max_sync_ns * to_us,
                    prof.contacts_ns * to_us,
                    prof.rebuild_ns * to_us, prof.max_rebuild_ns * to_us, prof.rebuilds,
                    prof.params_ns * to_us, prof.forces_ns * to_us, prof.over_2ms_ticks,
                    prof.threads_used, prof.substeps, prof.iterations);
            std::fflush(stdout);
            prof.ticks = 0;
            prof.total_ns = 0.0;
            prof.solver_ns = 0.0;
            prof.sync_ns = 0.0;
            prof.contacts_ns = 0.0;
            prof.rebuild_ns = 0.0;
            prof.params_ns = 0.0;
            prof.forces_ns = 0.0;
            prof.max_total_ns = 0.0;
            prof.max_solver_ns = 0.0;
            prof.max_sync_ns = 0.0;
            prof.max_contacts_ns = 0.0;
            prof.max_rebuild_ns = 0.0;
            prof.max_params_ns = 0.0;
            prof.max_forces_ns = 0.0;
            prof.over_2ms_ticks = 0;
            prof.rebuilds = 0;
            prof.threads_used = 1;
        }
    }
}

RID AVBDPhysicsServer3D::_joint_create() {
    const uint64_t id = next_id++;
    joints.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

void AVBDPhysicsServer3D::_joint_clear(const RID &p_joint) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    // Clearing resets the record to an untyped joint, as the engine's contract states; bodies
    // the joint constrained get a rebuild so the solver constraint goes away.
    const RID body_a = found->second.body_a;
    const RID body_b = found->second.body_b;
    found->second = JointData();
    rebuild_space(find_body(body_a) ? find_body(body_a)->space : RID());
    if (BodyData *b = find_body(body_b)) {
        rebuild_space(b->space);
    }
}

PhysicsServer3D::JointType AVBDPhysicsServer3D::_joint_get_type(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? PhysicsServer3D::JOINT_TYPE_6DOF : found->second.type;
}

void AVBDPhysicsServer3D::_joint_set_solver_priority(const RID &p_joint, int32_t p_priority) {
    const auto found = joints.find(id_of(p_joint));
    if (found != joints.end()) {
        found->second.priority = p_priority; // AVBD has one solver; ordering is not modelled
    }
}

int32_t AVBDPhysicsServer3D::_joint_get_solver_priority(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? 1 : found->second.priority;
}

void AVBDPhysicsServer3D::_joint_disable_collisions_between_bodies(const RID &p_joint, bool p_disable) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    // A solver joint keeps its bodies from colliding no matter what (they share a force), so
    // Godot's default exclusion comes free; asking for contact between the pair cannot be
    // honoured and the record only reflects the request.
    found->second.exclude_collision = p_disable;
}

bool AVBDPhysicsServer3D::_joint_is_disabled_collisions_between_bodies(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? true : found->second.exclude_collision;
}

void AVBDPhysicsServer3D::_joint_make_generic_6dof(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_6DOF;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;

    if (BodyData *body_b = find_body(p_body_b)) {
        rebuild_space(body_b->space);
    }
}

void AVBDPhysicsServer3D::_joint_make_pin(const RID &p_joint, const RID &p_body_a, const Vector3 &p_local_a,
        const RID &p_body_b, const Vector3 &p_local_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_PIN;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.pin_local_a = p_local_a;
    joint.pin_local_b = p_local_b;
    joint.ref_a = Transform3D(Basis(), p_local_a);
    joint.ref_b = Transform3D(Basis(), p_local_b);
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_pin_joint_set_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param,
        float p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    // BIAS/DAMPING/IMPULSE_CLAMP are recorded only; the solver's stabilisation is its own.
}

float AVBDPhysicsServer3D::_pin_joint_get_param(const RID &p_joint, PhysicsServer3D::PinJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_pin_joint_set_local_a(const RID &p_joint, const Vector3 &p_local_a) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.pin_local_a = p_local_a;
    found->second.ref_a = Transform3D(Basis(), p_local_a);
    _rebuild_for_joint(p_joint);
}

Vector3 AVBDPhysicsServer3D::_pin_joint_get_local_a(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? Vector3() : found->second.pin_local_a;
}

void AVBDPhysicsServer3D::_pin_joint_set_local_b(const RID &p_joint, const Vector3 &p_local_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.pin_local_b = p_local_b;
    found->second.ref_b = Transform3D(Basis(), p_local_b);
    _rebuild_for_joint(p_joint);
}

Vector3 AVBDPhysicsServer3D::_pin_joint_get_local_b(const RID &p_joint) const {
    const auto found = joints.find(id_of(p_joint));
    return found == joints.end() ? Vector3() : found->second.pin_local_b;
}

void AVBDPhysicsServer3D::_joint_make_hinge(const RID &p_joint, const RID &p_body_a, const Transform3D &p_hinge_a,
        const RID &p_body_b, const Transform3D &p_hinge_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_HINGE;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_hinge_a;
    joint.ref_b = p_hinge_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_joint_make_hinge_simple(const RID &p_joint, const RID &p_body_a,
        const Vector3 &p_pivot_a, const Vector3 &p_axis_a, const RID &p_body_b, const Vector3 &p_pivot_b,
        const Vector3 &p_axis_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_HINGE;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    // Store the hinge axis in the Z column of the stored frames, matching the full form's
    // convention (the hinge axis is the frame's Z, per godot_hinge_joint_3d.cpp:160).
    const Vector3 axis_a = p_axis_a.normalized();
    const Vector3 axis_b = p_axis_b.normalized();
    Vector3 ortho_a = std::fabs(axis_a.z) < 0.9 ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
    Vector3 ortho_b = std::fabs(axis_b.z) < 0.9 ? Vector3(0, 0, 1) : Vector3(1, 0, 0);
    joint.ref_a = Transform3D(Basis(ortho_a.cross(axis_a).normalized(), axis_a.cross(ortho_a.cross(axis_a)).normalized(),
            axis_a), p_pivot_a);
    joint.ref_b = Transform3D(Basis(ortho_b.cross(axis_b).normalized(), axis_b.cross(ortho_b.cross(axis_b)).normalized(),
            axis_b), p_pivot_b);
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_hinge_joint_set_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param,
        float p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::HINGE_JOINT_LIMIT_LOWER || p_param == PhysicsServer3D::HINGE_JOINT_LIMIT_UPPER) {
        _rebuild_for_joint(p_joint);
    }
}

float AVBDPhysicsServer3D::_hinge_joint_get_param(const RID &p_joint, PhysicsServer3D::HingeJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_hinge_joint_set_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag,
        bool p_enabled) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_flags[static_cast<int>(p_flag)] = p_enabled;
    if (p_flag == PhysicsServer3D::HINGE_JOINT_FLAG_USE_LIMIT) {
        _rebuild_for_joint(p_joint);
    }
}

bool AVBDPhysicsServer3D::_hinge_joint_get_flag(const RID &p_joint, PhysicsServer3D::HingeJointFlag p_flag) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return false;
    }
    const auto flag = found->second.joint_flags.find(static_cast<int>(p_flag));
    return flag != found->second.joint_flags.end() && flag->second;
}

void AVBDPhysicsServer3D::_joint_make_slider(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_SLIDER;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_slider_joint_set_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param,
        float p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_UPPER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_LINEAR_LIMIT_LOWER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_UPPER ||
            p_param == PhysicsServer3D::SLIDER_JOINT_ANGULAR_LIMIT_LOWER) {
        _rebuild_for_joint(p_joint);
    }
}

float AVBDPhysicsServer3D::_slider_joint_get_param(const RID &p_joint, PhysicsServer3D::SliderJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_joint_make_cone_twist(const RID &p_joint, const RID &p_body_a,
        const Transform3D &p_local_ref_a, const RID &p_body_b, const Transform3D &p_local_ref_b) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    JointData &joint = found->second;
    joint.type = PhysicsServer3D::JOINT_TYPE_CONE_TWIST;
    joint.body_a = p_body_a;
    joint.body_b = p_body_b;
    joint.ref_a = p_local_ref_a;
    joint.ref_b = p_local_ref_b;
    _rebuild_for_joint(p_joint);
}

void AVBDPhysicsServer3D::_cone_twist_joint_set_param(const RID &p_joint, PhysicsServer3D::ConeTwistJointParam p_param,
        float p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return;
    }
    found->second.joint_params[static_cast<int>(p_param)] = p_value;
    if (p_param == PhysicsServer3D::CONE_TWIST_JOINT_SWING_SPAN ||
            p_param == PhysicsServer3D::CONE_TWIST_JOINT_TWIST_SPAN) {
        _rebuild_for_joint(p_joint);
    }
}

float AVBDPhysicsServer3D::_cone_twist_joint_get_param(const RID &p_joint,
        PhysicsServer3D::ConeTwistJointParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end()) {
        return 0.0;
    }
    const auto param = found->second.joint_params.find(static_cast<int>(p_param));
    return param == found->second.joint_params.end() ? 0.0 : param->second;
}

void AVBDPhysicsServer3D::_generic_6dof_joint_set_param(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisParam p_param, float p_value) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_param) >= JointData::kMaxParams) {
        return;
    }
    // Rows 0..2 are linear axes, 3..5 angular. The param enum numbers linear and angular
    // families separately; the row offset puts them side by side.
    const bool angular = static_cast<int>(p_param) >= PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT;
    const int param_index = angular
            ? static_cast<int>(p_param) - PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT
            : static_cast<int>(p_param);
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    found->second.params[row][param_index] = p_value;
    if (BodyData *body_b = find_body(found->second.body_b)) {
        rebuild_space(body_b->space);
    } else if (BodyData *body_a = find_body(found->second.body_a)) {
        rebuild_space(body_a->space);
    }
}

float AVBDPhysicsServer3D::_generic_6dof_joint_get_param(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisParam p_param) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_param) >= JointData::kMaxParams) {
        return 0.0;
    }
    const bool angular = static_cast<int>(p_param) >= PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT;
    const int param_index = angular
            ? static_cast<int>(p_param) - PhysicsServer3D::G6DOF_JOINT_ANGULAR_LOWER_LIMIT
            : static_cast<int>(p_param);
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    return found->second.params[row][param_index];
}

void AVBDPhysicsServer3D::_generic_6dof_joint_set_flag(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisFlag p_flag, bool p_enable) {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_flag) >= JointData::kMaxFlags) {
        return;
    }
    const bool angular = static_cast<int>(p_flag) == PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT ||
            static_cast<int>(p_flag) == g6dof_spring::FLAG_ENABLE_ANGULAR;
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    found->second.flags[row][static_cast<int>(p_flag)] = p_enable;
    if (BodyData *body_b = find_body(found->second.body_b)) {
        rebuild_space(body_b->space);
    } else if (BodyData *body_a = find_body(found->second.body_a)) {
        rebuild_space(body_a->space);
    }
}

bool AVBDPhysicsServer3D::_generic_6dof_joint_get_flag(const RID &p_joint, Vector3::Axis p_axis,
        PhysicsServer3D::G6DOFJointAxisFlag p_flag) const {
    const auto found = joints.find(id_of(p_joint));
    if (found == joints.end() || p_axis < 0 || p_axis > Vector3::AXIS_Z ||
            static_cast<int>(p_flag) >= JointData::kMaxFlags) {
        return false;
    }
    const bool angular = static_cast<int>(p_flag) == PhysicsServer3D::G6DOF_JOINT_FLAG_ENABLE_ANGULAR_LIMIT ||
            static_cast<int>(p_flag) == g6dof_spring::FLAG_ENABLE_ANGULAR;
    const int row = (angular ? 3 : 0) + static_cast<int>(p_axis);
    return found->second.flags[row][static_cast<int>(p_flag)];
}

// -----------------------------------------------------------------------------
// Bodies
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_body_create() {
    const uint64_t id = next_id++;
    bodies.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

void AVBDPhysicsServer3D::_body_set_space(const RID &p_body, const RID &p_space) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    if (body->space == p_space) {
        return; // pointless, as the built-in server says
    }
    (void)0;
    const RID previous = body->space;
    body->space = p_space;
    // The engine calls this the moment a node enters the tree - before its shapes and
    // physics material arrive. Rebuilding now would freeze an empty body into the
    // solver; deferring to the next step lets the whole setup land first.
    space_rebuild_pending[id_of(p_space)] = true;
    if (previous.is_valid()) {
        space_rebuild_pending[id_of(previous)] = true;
    }
}

RID AVBDPhysicsServer3D::_body_get_space(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? RID() : found->second.space;
}

void AVBDPhysicsServer3D::_body_set_mode(const RID &p_body, PhysicsServer3D::BodyMode p_mode) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->mode = p_mode;
    rebuild_space(body->space);
}

PhysicsServer3D::BodyMode AVBDPhysicsServer3D::_body_get_mode(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? PhysicsServer3D::BODY_MODE_RIGID : found->second.mode;
}

void AVBDPhysicsServer3D::_body_set_state(const RID &p_body, PhysicsServer3D::BodyState p_state, const Variant &p_value) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }

    switch (p_state) {
        case PhysicsServer3D::BODY_STATE_TRANSFORM:
            body->transform = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                // A transform set from the scene is the body's new pose right now, so the
                // solver is told to match rather than left to interpolate towards it.
                const Transform3D &t = body->transform;
                modern->set_body_transform(static_cast<std::size_t>(body->modern_index),
                        to_modern(t.origin), to_modern(t.basis.get_rotation_quaternion()));
            }
            body->sleeping = false;
            break;

        case PhysicsServer3D::BODY_STATE_LINEAR_VELOCITY:
            body->linear_velocity = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                modern->set_linear_velocity(static_cast<std::size_t>(body->modern_index),
                        to_modern(static_cast<Vector3>(body->linear_velocity)));
            }
            body->sleeping = false;
            break;

        case PhysicsServer3D::BODY_STATE_ANGULAR_VELOCITY:
            body->angular_velocity = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                modern->set_angular_velocity(static_cast<std::size_t>(body->modern_index),
                        to_modern(static_cast<Vector3>(body->angular_velocity)));
            }
            body->sleeping = false;
            break;

        case PhysicsServer3D::BODY_STATE_SLEEPING: {
            body->sleeping = p_value;
            if (p_value) {
                // Godot zeroes velocity when a body is put to sleep explicitly. There is
                // no sleep in the modern core yet, so this only records the state and
                // parks the velocity the engine can see.
                body->linear_velocity = Vector3();
                body->angular_velocity = Vector3();
                if (ModernSpace *modern = modern_for(*body)) {
                    modern->set_linear_velocity(static_cast<std::size_t>(body->modern_index), Vec3{0});
                    modern->set_angular_velocity(static_cast<std::size_t>(body->modern_index), Vec3{0});
                }
            }
            break;
        }

        case PhysicsServer3D::BODY_STATE_CAN_SLEEP:
            // Rebuild applies sleep_mode = can_sleep ? SLEEP : NEVER; the setter stores the
            // flag and lets the next rebuild stamp it. Sleeping bodies wake to re-evaluate.
            body->can_sleep = p_value;
            rebuild_space(body->space);
            break;

        default:
            break;
    }
}

Variant AVBDPhysicsServer3D::_body_get_state(const RID &p_body, PhysicsServer3D::BodyState p_state) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return Variant();
    }
    const BodyData &body = found->second;

    switch (p_state) {
        case PhysicsServer3D::BODY_STATE_TRANSFORM:
            return body.transform;
        case PhysicsServer3D::BODY_STATE_LINEAR_VELOCITY:
            return body.linear_velocity;
        case PhysicsServer3D::BODY_STATE_ANGULAR_VELOCITY:
            return body.angular_velocity;
        case PhysicsServer3D::BODY_STATE_SLEEPING:
            return body.sleeping;
        case PhysicsServer3D::BODY_STATE_CAN_SLEEP:
            return body.can_sleep;
        default:
            return false;
    }
}

void AVBDPhysicsServer3D::_body_set_state_sync_callback(const RID &p_body, const Callable &p_callable) {
    if (BodyData *body = find_body(p_body)) {
        body->sync = p_callable;
    }
}

void AVBDPhysicsServer3D::_body_attach_object_instance_id(const RID &p_body, uint64_t p_id) {
    if (BodyData *body = find_body(p_body)) {
        body->instance_id = p_id;
    }
}

uint64_t AVBDPhysicsServer3D::_body_get_object_instance_id(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.instance_id;
}

void AVBDPhysicsServer3D::_body_set_ray_pickable(const RID &p_body, bool p_enable) {
    // Space queries do not read this yet; storing keeps the getter honest.
    if (BodyData *body = find_body(p_body)) {
        body->ray_pickable = p_enable;
    }
}

// --- DirectBodyState accessors -------------------------------------------------
//
// The state object is the engine's view of one body during (or after) a step; these
// helpers do the BodyData lookups it cannot do itself. All of them tolerate an unknown
// RID: the engine may ask about a body this server never saw.

void AVBDPhysicsServer3D::body_state_add_constant_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force += p_force;
        if (p_position.length_squared() > 0.0) {
            body->constant_torque += p_position.cross(p_force);
        }
    }
}

void AVBDPhysicsServer3D::body_state_add_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque += p_torque;
    }
}

void AVBDPhysicsServer3D::body_state_set_constant_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

Vector3 AVBDPhysicsServer3D::body_state_get_constant_force(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_force;
}

void AVBDPhysicsServer3D::body_state_set_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque = p_torque;
    }
}

Vector3 AVBDPhysicsServer3D::body_state_get_constant_torque(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_torque;
}

int32_t AVBDPhysicsServer3D::body_state_contact_count(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : static_cast<int32_t>(found->second.frame_contacts.size());
}

Vector3 AVBDPhysicsServer3D::body_state_contact_position(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return Vector3();
    }
    return found->second.frame_contacts[p_index].position;
}

RID AVBDPhysicsServer3D::body_state_contact_collider(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return RID();
    }
    return make_rid(found->second.frame_contacts[p_index].collider_body_id);
}

uint64_t AVBDPhysicsServer3D::body_state_contact_collider_id(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return 0;
    }
    const BodyData *other = find_body(make_rid(found->second.frame_contacts[p_index].collider_body_id));
    return other != nullptr ? other->instance_id : 0;
}

Vector3 AVBDPhysicsServer3D::body_state_contact_collider_velocity(const RID &p_body, int32_t p_index) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_index < 0 || static_cast<size_t>(p_index) >= found->second.frame_contacts.size()) {
        return Vector3();
    }
    const BodyData *other = find_body(make_rid(found->second.frame_contacts[p_index].collider_body_id));
    return other != nullptr ? other->linear_velocity : Vector3();
}

double AVBDPhysicsServer3D::body_state_step() const {
    // Every space steps with the engine's fixed timestep; report the common value.
    for (const auto &[id, space] : spaces) {
        (void)id;
        if (space->modern != nullptr)
            return static_cast<double>(space->modern->state().dt);
    }
    return 0.0;
}

// The modern core reports one total per step; the legacy phase split (broad phase /
// colour graph / solve / finish) does not exist in it, so those getters stay 0 rather
// than returning numbers that would be made up.
double AVBDPhysicsServer3D::get_solver_broad_phase_time() const {
    return 0.0;
}

double AVBDPhysicsServer3D::get_solver_colour_graph_time() const {
    return 0.0;
}

double AVBDPhysicsServer3D::get_solver_solve_time() const {
    double total = 0.0;
    for (const auto &[id, space] : spaces) {
        (void)id;
        if (space->modern != nullptr)
            total = std::max(total, space->modern->last_step_ms());
    }
    return total;
}

double AVBDPhysicsServer3D::get_solver_finish_time() const {
    return 0.0;
}

// The deepest normal overlap across every live manifold, as snapshotted at the end of the
// last step. A contact's C0 is stamped at the start of the step (Manifold::initialize), so
// this is the geometric penetration the last solve was working against, not a
// solver-internal residual.
double AVBDPhysicsServer3D::get_solver_max_penetration() const {
    return diag_max_penetration;
}

Array AVBDPhysicsServer3D::get_solver_deepest_contact() const {
    Array out;
    if (diag_max_penetration > 0.0) {
        out.push_back(static_cast<int64_t>(diag_deepest_a));
        out.push_back(static_cast<int64_t>(diag_deepest_b));
        out.push_back(diag_max_penetration);
    }
    return out;
}

// Sleep/wake diagnostics. The modern core has no sleeping yet, so every movable body is
// reported awake and no manifold is frozen; the counters stay meaningful for the stress
// suites (`sleeping == 0` is itself the statement that sleep is not implemented).
Dictionary AVBDPhysicsServer3D::get_solver_sleep_stats() const {
    int bodies_count = 0;
    int movable = 0;
    int manifolds = 0;

    for (const auto &[id, space] : spaces) {
        (void)id;
        if (!space->active || space->modern == nullptr) {
            continue;
        }
        const ModernSpace &modern = *space->modern;
        bodies_count += static_cast<int>(modern.size());
        for (std::size_t i = 0; i < modern.size(); ++i) {
            if (!modern.is_static(i)) {
                ++movable;
            }
        }
        manifolds += static_cast<int>(modern.manifolds().size());
    }

    Dictionary out;
    out["bodies"] = bodies_count;
    out["movable"] = movable;
    out["sleeping"] = 0;
    out["sleep_mode_never"] = 0;
    out["max_still_frames"] = 0;
    out["update_set"] = movable;
    out["manifolds"] = manifolds;
    out["frozen_manifolds"] = 0;
    out["sleep_frames"] = 0;
    out["rebuilds"] = static_cast<int64_t>(rebuild_total);
    out["steps"] = static_cast<int64_t>(step_total);
    out["rebuild_last_step"] = static_cast<int64_t>(rebuild_last_step);
    out["rebuild_last_reason"] = rebuild_last_reason;
    Dictionary reasons;
    for (const auto &[reason, count] : rebuild_reasons) {
        reasons[String::utf8(reason.c_str())] = static_cast<int64_t>(count);
    }
    out["rebuild_reasons"] = reasons;
    return out;
}

AVBDPhysicsServer3D *AVBDPhysicsServer3D::get_singleton() {
    return SpaceData::singleton;
}PhysicsDirectSpaceState3D *AVBDPhysicsServer3D::body_state_space_state(const RID &p_body) {
    const BodyData *body = find_body(p_body);
    if (body == nullptr || body->space.is_valid() == false) {
        return nullptr;
    }
    return _space_get_direct_state(body->space);
}

PhysicsDirectBodyState3D *AVBDPhysicsServer3D::_body_get_direct_state(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return nullptr;
    }
    if (body_state == nullptr) {
        body_state = memnew(AVBDDirectBodyState3D);
    }
    body_state->bind(this, p_body);
    return body_state;
}

// --- query helpers -------------------------------------------------------------

Transform3D AVBDPhysicsServer3D::body_shape_transform(const BodyData &p_body, int32_t p_shape_idx) const {
    if (p_shape_idx < 0 || static_cast<std::size_t>(p_shape_idx) >= p_body.shapes.size()) {
        return p_body.transform;
    }
    const Transform3D local = static_cast<std::size_t>(p_shape_idx) < p_body.shape_transforms.size()
            ? p_body.shape_transforms[p_shape_idx]
            : Transform3D();
    return p_body.transform * local;
}

std::vector<AVBDPhysicsServer3D::QueryShapeHit> AVBDPhysicsServer3D::query_body_shapes(
        const RID &p_space, uint32_t p_mask) const {
    std::vector<QueryShapeHit> out;
    const uint64_t space_id = id_of(p_space);
    if (!spaces.contains(space_id)) {
        return out;
    }
    for (const auto &[body_id, body] : bodies) {
        if (id_of(body.space) != space_id) {
            continue;
        }
        if ((body.collision_layer & p_mask) == 0) {
            continue;
        }
        for (std::size_t slot = 0; slot < body.shapes.size(); ++slot) {
            if (slot < body.shape_disabled.size() && body.shape_disabled[slot]) {
                continue;
            }
            const Transform3D posed = body_shape_transform(body, static_cast<int32_t>(slot));
            collision::Shape shape;
            if (!query_shape_from_data(body.shapes[slot], posed, shape)) {
                continue;
            }
            QueryShapeHit hit;
            hit.shape = shape;
            hit.transform = posed;
            hit.rid = make_rid(body_id);
            hit.collider_id = body.instance_id;
            hit.body_id = body_id;
            hit.shape_index = static_cast<int32_t>(slot);
            hit.linear_velocity = body.linear_velocity;
            hit.is_area = false;
            hit.ray_pickable = body.ray_pickable;
            out.push_back(hit);
        }
    }
    return out;
}

std::vector<AVBDPhysicsServer3D::QueryShapeHit> AVBDPhysicsServer3D::query_area_shapes(
        const RID &p_space, uint32_t p_mask) const {
    std::vector<QueryShapeHit> out;
    const uint64_t space_id = id_of(p_space);
    for (const auto &[area_id, area] : areas) {
        if (!area.space.is_valid() || id_of(area.space) != space_id) {
            continue;
        }
        if ((area.collision_layer & p_mask) == 0) {
            continue;
        }
        for (std::size_t i = 0; i < area.shapes.size(); ++i) {
            if (i < area.shape_disabled.size() && area.shape_disabled[i]) {
                continue;
            }
            const Transform3D posed = area.transform * area.shape_transforms[i];
            collision::Shape shape;
            if (!query_shape_from_data(area.shapes[i], posed, shape)) {
                continue;
            }
            QueryShapeHit hit;
            hit.shape = shape;
            hit.transform = posed;
            hit.rid = make_rid(area_id);
            hit.collider_id = area.instance_id;
            hit.body_id = 0;
            hit.shape_index = static_cast<int32_t>(i);
            hit.is_area = true;
            hit.ray_pickable = area.ray_pickable;
            out.push_back(hit);
        }
    }
    return out;
}

bool AVBDPhysicsServer3D::body_is_excepted(uint64_t p_body_a, uint64_t p_body_b) const {
    const auto a = bodies.find(p_body_a);
    const auto b = bodies.find(p_body_b);
    if (a == bodies.end() || b == bodies.end()) {
        return false;
    }
    const auto lists = [](const BodyData &from, uint64_t other) {
        return std::find(from.exceptions.begin(), from.exceptions.end(), other) != from.exceptions.end();
    };
    return lists(a->second, p_body_b) || lists(b->second, p_body_a);
}

const AVBDPhysicsServer3D::ShapeData *AVBDPhysicsServer3D::shape_data(const RID &p_shape) const {
    const auto found = shapes.find(id_of(p_shape));
    return found == shapes.end() ? nullptr : &found->second;
}

// -----------------------------------------------------------------------------
// Body parameters and forces
// -----------------------------------------------------------------------------

void AVBDPhysicsServer3D::_body_set_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param,
        const Variant &p_value) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }

    switch (p_param) {
        case PhysicsServer3D::BODY_PARAM_MASS:
            body->mass = p_value;
            rebuild_space(body->space);
            break;
        case PhysicsServer3D::BODY_PARAM_FRICTION:
            body->friction = p_value;
            // Friction is carried by the shape material, so a change needs a rebuild.
            space_rebuild_pending[id_of(body->space)] = true;
            break;
        case PhysicsServer3D::BODY_PARAM_BOUNCE:
            // The AVBD contact model is friction-only: a bounce value is remembered for
            // _body_get_param but no restitution is simulated.
            body->bounce = p_value;
            break;
        case PhysicsServer3D::BODY_PARAM_GRAVITY_SCALE: {
            body->gravity_scale = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                modern->set_gravity_scale(static_cast<std::size_t>(body->modern_index),
                        static_cast<double>(body->gravity_scale));
            }
            break;
        }
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP:
            body->linear_damp = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                modern->set_damping(static_cast<std::size_t>(body->modern_index),
                        body->linear_damp, body->angular_damp);
            }
            break;
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP:
            body->angular_damp = p_value;
            if (ModernSpace *modern = modern_for(*body)) {
                modern->set_damping(static_cast<std::size_t>(body->modern_index),
                        body->linear_damp, body->angular_damp);
            }
            break;
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP_MODE:
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP_MODE:
            break; // damp modes are meaningless without damping
        case PhysicsServer3D::BODY_PARAM_INERTIA:
        case PhysicsServer3D::BODY_PARAM_CENTER_OF_MASS:
            break; // inertia follows the shape; the centre of mass is the body origin
        default:
            break;
    }
}

Variant AVBDPhysicsServer3D::_body_get_param(const RID &p_body, PhysicsServer3D::BodyParameter p_param) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return Variant();
    }
    const BodyData &body = found->second;

    switch (p_param) {
        case PhysicsServer3D::BODY_PARAM_MASS:
            return body.mass;
        case PhysicsServer3D::BODY_PARAM_FRICTION:
            return body.friction;
        case PhysicsServer3D::BODY_PARAM_BOUNCE:
            return body.bounce;
        case PhysicsServer3D::BODY_PARAM_GRAVITY_SCALE:
            return body.gravity_scale;
        case PhysicsServer3D::BODY_PARAM_LINEAR_DAMP:
            return body.linear_damp;
        case PhysicsServer3D::BODY_PARAM_ANGULAR_DAMP:
            return body.angular_damp;
        default:
            return 0.0;
    }
}

// Impulses and forces arrive in Godot's world frame, which is the modern core's frame
// too; the modern space converts an impulse into a velocity change with the body's
// inverse mass and world inverse inertia.

void AVBDPhysicsServer3D::_body_apply_central_impulse(const RID &p_body, const Vector3 &p_impulse) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    if (ModernSpace *modern = modern_for(*body)) {
        modern->apply_central_impulse(static_cast<std::size_t>(body->modern_index),
                to_modern(p_impulse));
    }
}

void AVBDPhysicsServer3D::_body_apply_impulse(const RID &p_body, const Vector3 &p_impulse,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    if (ModernSpace *modern = modern_for(*body)) {
        modern->apply_impulse(static_cast<std::size_t>(body->modern_index),
                to_modern(p_impulse), to_modern(p_position));
    }
}

void AVBDPhysicsServer3D::_body_apply_torque_impulse(const RID &p_body, const Vector3 &p_impulse) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    if (ModernSpace *modern = modern_for(*body)) {
        modern->apply_torque_impulse(static_cast<std::size_t>(body->modern_index),
                to_modern(p_impulse));
    }
}

void AVBDPhysicsServer3D::_body_apply_central_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->frame_force += p_force;
    }
}

void AVBDPhysicsServer3D::_body_apply_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->frame_force += p_force;
    if (p_position.length_squared() > 0.0) {
        body->frame_torque += p_position.cross(p_force);
    }
}

void AVBDPhysicsServer3D::_body_apply_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->frame_torque += p_torque;
    }
}

void AVBDPhysicsServer3D::_body_add_constant_central_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

void AVBDPhysicsServer3D::_body_add_constant_force(const RID &p_body, const Vector3 &p_force,
        const Vector3 &p_position) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->constant_force += p_force;
    if (p_position.length_squared() > 0.0) {
        body->constant_torque += p_position.cross(p_force);
    }
}

void AVBDPhysicsServer3D::_body_add_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque += p_torque;
    }
}

// Godot BodyAxis bits (linear X/Y/Z = 1/2/4, angular X/Y/Z = 8/16/32) to the solver's
// sim-axis bits (0..2 = sim x/y/z). Godot Y is sim Z, Godot Z is sim -Y; the lock is a
// directionless bit so the sign does not matter, only the permutation.
void AVBDPhysicsServer3D::_body_set_axis_lock(const RID &p_body, PhysicsServer3D::BodyAxis p_axis, bool p_lock) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    const uint32_t bit = static_cast<uint32_t>(p_axis);
    if (p_lock) {
        body->axis_locks |= bit;
    } else {
        body->axis_locks &= ~bit;
    }
    // Axis locks change the inverse inertia (and, for a full linear lock, immovability),
    // which is baked at rebuild time.
    if (body->space.is_valid()) {
        space_rebuild_pending[id_of(body->space)] = true;
    }
}

bool AVBDPhysicsServer3D::_body_is_axis_locked(const RID &p_body, PhysicsServer3D::BodyAxis p_axis) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return false;
    }
    return (found->second.axis_locks & static_cast<uint32_t>(p_axis)) != 0;
}

void AVBDPhysicsServer3D::_body_add_collision_exception(const RID &p_body, const RID &p_excepted_body) {
    BodyData *body = find_body(p_body);
    const uint64_t other = id_of(p_excepted_body);
    if (body == nullptr || other == id_of(p_body)) {
        return;
    }
    if (std::find(body->exceptions.begin(), body->exceptions.end(), other) == body->exceptions.end()) {
        body->exceptions.push_back(other);
    }
    if (body->space.is_valid()) {
        rebuild_space(body->space);
    }
}

void AVBDPhysicsServer3D::_body_remove_collision_exception(const RID &p_body, const RID &p_excepted_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    std::vector<uint64_t> &list = body->exceptions;
    list.erase(std::remove(list.begin(), list.end(), id_of(p_excepted_body)), list.end());
    if (body->space.is_valid()) {
        rebuild_space(body->space);
    }
}

TypedArray<RID> AVBDPhysicsServer3D::_body_get_collision_exceptions(const RID &p_body) const {
    TypedArray<RID> out;
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end()) {
        return out;
    }
    for (const uint64_t other_id : found->second.exceptions) {
        out.push_back(make_rid(other_id));
    }
    return out;
}

void AVBDPhysicsServer3D::_body_set_max_contacts_reported(const RID &p_body, int32_t p_amount) {
    if (BodyData *body = find_body(p_body)) {
        body->max_contacts_reported = p_amount;
    }
}

int32_t AVBDPhysicsServer3D::_body_get_max_contacts_reported(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.max_contacts_reported;
}

void AVBDPhysicsServer3D::_body_set_contacts_reported_depth_threshold(const RID &p_body, float p_threshold)  {
    if (BodyData *body = find_body(p_body)) {
        body->contacts_depth_threshold = p_threshold; // recorded; not applied (contact depths come from the solver)
    }
}

float AVBDPhysicsServer3D::_body_get_contacts_reported_depth_threshold(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0.0 : found->second.contacts_depth_threshold;
}

void AVBDPhysicsServer3D::_body_set_omit_force_integration(const RID &p_body, bool p_enable) {
    if (BodyData *body = find_body(p_body)) {
        body->omit_force_integration = p_enable;
    }
}

bool AVBDPhysicsServer3D::_body_is_omitting_force_integration(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? false : found->second.omit_force_integration;
}

void AVBDPhysicsServer3D::_body_set_force_integration_callback(const RID &p_body, const Callable &p_callable,
        const Variant &p_userdata) {
    if (BodyData *body = find_body(p_body)) {
        body->force_integration_callback = p_callable;
        body->force_integration_userdata = p_userdata;
    }
}

void AVBDPhysicsServer3D::_body_set_constant_force(const RID &p_body, const Vector3 &p_force) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_force = p_force;
    }
}

Vector3 AVBDPhysicsServer3D::_body_get_constant_force(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_force;
}

void AVBDPhysicsServer3D::_body_set_constant_torque(const RID &p_body, const Vector3 &p_torque) {
    if (BodyData *body = find_body(p_body)) {
        body->constant_torque = p_torque;
    }
}

Vector3 AVBDPhysicsServer3D::_body_get_constant_torque(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? Vector3() : found->second.constant_torque;
}

void AVBDPhysicsServer3D::_body_set_axis_velocity(const RID &p_body, const Vector3 &p_axis_velocity) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    // Godot semantics: each nonzero component replaces the velocity along that axis.
    Vector3 velocity = body->linear_velocity;
    if (p_axis_velocity.x != 0.0) {
        velocity.x = p_axis_velocity.x;
    }
    if (p_axis_velocity.y != 0.0) {
        velocity.y = p_axis_velocity.y;
    }
    if (p_axis_velocity.z != 0.0) {
        velocity.z = p_axis_velocity.z;
    }
    body->linear_velocity = velocity;
    if (ModernSpace *modern = modern_for(*body)) {
        modern->set_linear_velocity(static_cast<std::size_t>(body->modern_index), to_modern(velocity));
    }
}

void AVBDPhysicsServer3D::_body_set_enable_continuous_collision_detection(const RID &p_body, bool p_enable) {
    if (BodyData *body = find_body(p_body)) {
        body->ccd = p_enable; // recorded; CCD is not simulated (documented)
    }
}

bool AVBDPhysicsServer3D::_body_is_continuous_collision_detection_enabled(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? false : found->second.ccd;
}

void AVBDPhysicsServer3D::_body_set_collision_layer(const RID &p_body, uint32_t p_layer) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_layer = p_layer;
        // The layer/mask tables live in ModernSpace; refresh them before the next step.
        space_rebuild_pending[id_of(body->space)] = true;
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_collision_layer(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1 : found->second.collision_layer;
}

void AVBDPhysicsServer3D::_body_set_collision_mask(const RID &p_body, uint32_t p_mask) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_mask = p_mask;
        space_rebuild_pending[id_of(body->space)] = true;
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_collision_mask(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1 : found->second.collision_mask;
}

void AVBDPhysicsServer3D::_body_set_collision_priority(const RID &p_body, float p_priority) {
    if (BodyData *body = find_body(p_body)) {
        body->collision_priority = p_priority; // recorded; solver has no priority weighting
    }
}

float AVBDPhysicsServer3D::_body_get_collision_priority(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 1.0 : found->second.collision_priority;
}

void AVBDPhysicsServer3D::_body_set_user_flags(const RID &p_body, uint32_t p_flags) {
    if (BodyData *body = find_body(p_body)) {
        body->user_flags = p_flags;
    }
}

uint32_t AVBDPhysicsServer3D::_body_get_user_flags(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : found->second.user_flags;
}

void AVBDPhysicsServer3D::_body_reset_mass_properties(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    // Godot's contract for this call is "restore the default inertia and centre of mass to
    // cancel custom values". Both are always derived from the shape here, so the rebuild is
    // the whole of the reset - it must not touch the engine's mass.
    rebuild_space(body->space);
}

void AVBDPhysicsServer3D::_body_add_shape(const RID &p_body, const RID &p_shape, const Transform3D &p_transform,
        bool p_disabled) {
    BodyData *body = find_body(p_body);
    ShapeData *shape = find_shape(p_shape);
    if (body == nullptr || shape == nullptr) {
        return;
    }
    // Keep the RID: the engine may set the shape resource's data after attaching it, and
    // the deferred rebuild refreshes the snapshot from this RID.
    body->shape_rids.push_back(p_shape);
    body->shapes.push_back(*shape);
    body->shape_transforms.push_back(p_transform);
    body->shape_disabled.push_back(p_disabled);
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_set_shape(const RID &p_body, int32_t p_shape_idx, const RID &p_shape) {
    BodyData *body = find_body(p_body);
    ShapeData *shape = find_shape(p_shape);
    if (body == nullptr || shape == nullptr || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= body->shapes.size()) {
        return;
    }
    body->shape_rids[p_shape_idx] = p_shape;
    body->shapes[p_shape_idx] = *shape;
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_set_shape_transform(const RID &p_body, int32_t p_shape_idx,
        const Transform3D &p_transform) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= body->shape_transforms.size()) {
        return;
    }
    body->shape_transforms[p_shape_idx] = p_transform;
}

void AVBDPhysicsServer3D::_body_set_shape_disabled(const RID &p_body, int32_t p_shape_idx, bool p_disabled) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0 || static_cast<std::size_t>(p_shape_idx) >= body->shapes.size()) {
        return;
    }
    if (body->shape_disabled.size() < body->shapes.size()) {
        body->shape_disabled.resize(body->shapes.size(), false);
    }
    body->shape_disabled[p_shape_idx] = p_disabled;
    space_rebuild_pending[id_of(body->space)] = true;
}

int32_t AVBDPhysicsServer3D::_body_get_shape_count(const RID &p_body) const {
    const auto found = bodies.find(id_of(p_body));
    return found == bodies.end() ? 0 : static_cast<int32_t>(found->second.shapes.size());
}

RID AVBDPhysicsServer3D::_body_get_shape(const RID &p_body, int32_t p_shape_idx) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_shape_idx < 0 ||
            static_cast<std::size_t>(p_shape_idx) >= found->second.shape_rids.size()) {
        return RID();
    }
    return found->second.shape_rids[p_shape_idx];
}

Transform3D AVBDPhysicsServer3D::_body_get_shape_transform(const RID &p_body, int32_t p_shape_idx) const {
    const auto found = bodies.find(id_of(p_body));
    if (found == bodies.end() || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= found->second.shape_transforms.size()) {
        return Transform3D();
    }
    return found->second.shape_transforms[p_shape_idx];
}

void AVBDPhysicsServer3D::_body_remove_shape(const RID &p_body, int32_t p_shape_idx) {
    BodyData *body = find_body(p_body);
    if (body == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= body->shapes.size()) {
        return;
    }
    body->shapes.erase(body->shapes.begin() + p_shape_idx);
    body->shape_transforms.erase(body->shape_transforms.begin() + p_shape_idx);
    body->shape_rids.erase(body->shape_rids.begin() + p_shape_idx);
    if (static_cast<std::size_t>(p_shape_idx) < body->shape_disabled.size()) {
        body->shape_disabled.erase(body->shape_disabled.begin() + p_shape_idx);
    }
    space_rebuild_pending[id_of(body->space)] = true;
}

void AVBDPhysicsServer3D::_body_clear_shapes(const RID &p_body) {
    BodyData *body = find_body(p_body);
    if (body == nullptr) {
        return;
    }
    body->shapes.clear();
    body->shape_transforms.clear();
    body->shape_rids.clear();
    body->shape_disabled.clear();
    space_rebuild_pending[id_of(body->space)] = true;
}

// -----------------------------------------------------------------------------
// Areas
//
// Every storage virtual keeps the engine's description in an AreaData record, and
// _monitor_areas diffs each area against the bodies of its space every step, calling
// the monitor callbacks with Godot's four-argument contract (verified against 4.7
// scene/3d/physics/area_3d.cpp): (status, body_rid, instance_id, body_shape, area_shape)
// - status 0 = entered (AREA_BODY_ADDED), 1 = exited (AREA_BODY_REMOVED).
// -----------------------------------------------------------------------------

RID AVBDPhysicsServer3D::_area_create() {
    const uint64_t id = next_id++;
    areas.emplace(std::piecewise_construct, std::forward_as_tuple(id), std::forward_as_tuple());
    return make_rid(id);
}

AVBDPhysicsServer3D::AreaData *AVBDPhysicsServer3D::find_area(const RID &p_rid) {
    const auto found = areas.find(id_of(p_rid));
    return found == areas.end() ? nullptr : &found->second;
}

void AVBDPhysicsServer3D::_area_set_space(const RID &p_area, const RID &p_space) {
    if (AreaData *area = find_area(p_area)) {
        area->space = p_space;
        area->prev_inside_bodies.clear();
        area->prev_inside_areas.clear();
    }
}

RID AVBDPhysicsServer3D::_area_get_space(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? RID() : found->second.space;
}

void AVBDPhysicsServer3D::_area_set_transform(const RID &p_area, const Transform3D &p_transform) {
    if (AreaData *area = find_area(p_area)) {
        area->transform = p_transform;
    }
}

Transform3D AVBDPhysicsServer3D::_area_get_transform(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? Transform3D() : found->second.transform;
}

void AVBDPhysicsServer3D::_area_add_shape(const RID &p_area, const RID &p_shape, const Transform3D &p_transform,
        bool p_disabled) {
    AreaData *area = find_area(p_area);
    ShapeData *shape = find_shape(p_shape);
    if (area == nullptr || shape == nullptr) {
        return;
    }
    area->shapes.push_back(*shape);
    area->shape_transforms.push_back(p_transform);
    area->shape_disabled.push_back(p_disabled);
}

void AVBDPhysicsServer3D::_area_set_shape(const RID &p_area, int32_t p_shape_idx, const RID &p_shape) {
    AreaData *area = find_area(p_area);
    ShapeData *shape = find_shape(p_shape);
    if (area == nullptr || shape == nullptr || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= area->shapes.size()) {
        return;
    }
    area->shapes[p_shape_idx] = *shape;
}

void AVBDPhysicsServer3D::_area_set_shape_transform(const RID &p_area, int32_t p_shape_idx,
        const Transform3D &p_transform) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shape_transforms.size()) {
        return;
    }
    area->shape_transforms[p_shape_idx] = p_transform;
}

void AVBDPhysicsServer3D::_area_set_shape_disabled(const RID &p_area, int32_t p_shape_idx, bool p_disabled) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shape_disabled.size()) {
        return;
    }
    area->shape_disabled[p_shape_idx] = p_disabled;
}

int32_t AVBDPhysicsServer3D::_area_get_shape_count(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 0 : static_cast<int32_t>(found->second.shapes.size());
}

RID AVBDPhysicsServer3D::_area_get_shape(const RID &, int32_t) const {
    // As with bodies: the engine owns the shape RIDs; this server only keeps parameters.
    return RID();
}

Transform3D AVBDPhysicsServer3D::_area_get_shape_transform(const RID &p_area, int32_t p_shape_idx) const {
    const auto found = areas.find(id_of(p_area));
    if (found == areas.end() || p_shape_idx < 0 ||
            static_cast<size_t>(p_shape_idx) >= found->second.shape_transforms.size()) {
        return Transform3D();
    }
    return found->second.shape_transforms[p_shape_idx];
}

void AVBDPhysicsServer3D::_area_remove_shape(const RID &p_area, int32_t p_shape_idx) {
    AreaData *area = find_area(p_area);
    if (area == nullptr || p_shape_idx < 0 || static_cast<size_t>(p_shape_idx) >= area->shapes.size()) {
        return;
    }
    area->shapes.erase(area->shapes.begin() + p_shape_idx);
    area->shape_transforms.erase(area->shape_transforms.begin() + p_shape_idx);
    area->shape_disabled.erase(area->shape_disabled.begin() + p_shape_idx);
}

void AVBDPhysicsServer3D::_area_clear_shapes(const RID &p_area) {
    if (AreaData *area = find_area(p_area)) {
        area->shapes.clear();
        area->shape_transforms.clear();
        area->shape_disabled.clear();
    }
}

void AVBDPhysicsServer3D::_area_set_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param,
        const Variant &p_value) {
    if (AreaData *area = find_area(p_area)) {
        // Recorded so the getter replays it; gravity overrides and the other area physics
        // parameters do not act on the solver this round (documented).
        area->params[p_param] = p_value;
    }
}

Variant AVBDPhysicsServer3D::_area_get_param(const RID &p_area, PhysicsServer3D::AreaParameter p_param) const {
    const auto found = areas.find(id_of(p_area));
    if (found == areas.end()) {
        return Variant();
    }
    const auto param = found->second.params.find(p_param);
    return param == found->second.params.end() ? Variant() : param->second;
}

void AVBDPhysicsServer3D::_area_set_collision_layer(const RID &p_area, uint32_t p_layer) {
    if (AreaData *area = find_area(p_area)) {
        area->collision_layer = p_layer;
    }
}

uint32_t AVBDPhysicsServer3D::_area_get_collision_layer(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 1 : found->second.collision_layer;
}

void AVBDPhysicsServer3D::_area_set_collision_mask(const RID &p_area, uint32_t p_mask) {
    if (AreaData *area = find_area(p_area)) {
        area->collision_mask = p_mask;
    }
}

uint32_t AVBDPhysicsServer3D::_area_get_collision_mask(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 1 : found->second.collision_mask;
}

void AVBDPhysicsServer3D::_area_set_monitorable(const RID &p_area, bool p_monitorable) {
    if (AreaData *area = find_area(p_area)) {
        area->monitorable = p_monitorable;
    }
}

void AVBDPhysicsServer3D::_area_set_ray_pickable(const RID &p_area, bool p_enable) {
    if (AreaData *area = find_area(p_area)) {
        area->ray_pickable = p_enable; // space queries do not hit areas this round
    }
}

void AVBDPhysicsServer3D::_area_attach_object_instance_id(const RID &p_area, uint64_t p_id) {
    if (AreaData *area = find_area(p_area)) {
        area->instance_id = p_id;
    }
}

uint64_t AVBDPhysicsServer3D::_area_get_object_instance_id(const RID &p_area) const {
    const auto found = areas.find(id_of(p_area));
    return found == areas.end() ? 0 : found->second.instance_id;
}

void AVBDPhysicsServer3D::_area_set_monitor_callback(const RID &p_area, const Callable &p_callback) {
    if (AreaData *area = find_area(p_area)) {
        area->monitor_callback = p_callback;
        area->monitor_callback_set = p_callback.is_valid();
    }
}

void AVBDPhysicsServer3D::_area_set_area_monitor_callback(const RID &p_area, const Callable &p_callback) {
    if (AreaData *area = find_area(p_area)) {
        area->area_monitor_callback = p_callback;
        area->area_monitor_callback_set = p_callback.is_valid();
    }
}

// One step of area monitoring: diff every area's inside-set against the previous one and
// fire the callbacks on the transitions. Godot's body-shape/area-shape arguments are the
// first overlapping pair (this server keeps one solver shape per body, so 0 on the body
// side and the area's first enabled shape on its side).
// One step of area monitoring: diff every area's inside-set against the previous one and
// fire the callbacks on the transitions. Overlap is AABB-based on the modern shapes: the
// area path has no contact solver of its own, and a conservative box test is the right
// contract for enter/exit notifications.
void AVBDPhysicsServer3D::_monitor_areas() {
    struct Placed {
        collision::Shape shape;
        Transform3D transform;
    };
    const auto overlap = [](const Placed &a, const Placed &b) {
        Vec3 a_lo, a_hi, b_lo, b_hi;
        shape_world_aabb(a.shape, a.transform, a_lo, a_hi);
        shape_world_aabb(b.shape, b.transform, b_lo, b_hi);
        return a_lo.x <= b_hi.x && a_hi.x >= b_lo.x && a_lo.y <= b_hi.y && a_hi.y >= b_lo.y &&
               a_lo.z <= b_hi.z && a_hi.z >= b_lo.z;
    };

    for (auto &[area_id, area] : areas) {
        (void)area_id;
        if (!area.space.is_valid()) {
            continue;
        }
        if (area.shape_disabled.size() < area.shapes.size()) {
            area.shape_disabled.resize(area.shapes.size(), false);
        }

        std::vector<Placed> area_shapes;
        for (std::size_t i = 0; i < area.shapes.size(); ++i) {
            if (area.shape_disabled[i]) {
                continue;
            }
            const Transform3D posed = area.transform * area.shape_transforms[i];
            collision::Shape shape;
            if (query_shape_from_data(area.shapes[i], posed, shape)) {
                area_shapes.push_back(Placed{shape, posed});
            }
        }
        if (area_shapes.empty()) {
            continue;
        }

        // --- bodies: the area reports what enters its mask ---
        std::vector<uint64_t> inside;
        if (area.monitor_callback_set || area.monitor_callback.is_valid()) {
            for (const auto &[body_id, body] : bodies) {
                if (id_of(body.space) != id_of(area.space) || body.modern_index < 0) {
                    continue;
                }
                // Godot's area test: the body's layer must meet the area's mask.
                if ((body.collision_layer & area.collision_mask) == 0) {
                    continue;
                }
                bool hits = false;
                for (std::size_t slot = 0; slot < body.shapes.size() && !hits; ++slot) {
                    if (slot < body.shape_disabled.size() && body.shape_disabled[slot]) {
                        continue;
                    }
                    const Transform3D posed = body_shape_transform(body, static_cast<int32_t>(slot));
                    collision::Shape shape;
                    if (!query_shape_from_data(body.shapes[slot], posed, shape)) {
                        continue;
                    }
                    const Placed body_shape{shape, posed};
                    for (const Placed &shape_area : area_shapes) {
                        if (overlap(shape_area, body_shape)) {
                            hits = true;
                            break;
                        }
                    }
                }
                if (hits) {
                    inside.push_back(body_id);
                }
            }

            auto was_inside = [&area](uint64_t id) -> bool {
                return std::find(area.prev_inside_bodies.begin(), area.prev_inside_bodies.end(), id)
                        != area.prev_inside_bodies.end();
            };
            for (const uint64_t body_id : inside) {
                if (was_inside(body_id) || !area.monitor_callback.is_valid()) {
                    continue;
                }
                const BodyData *body = find_body(make_rid(body_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_ADDED));
                args.push_back(make_rid(body_id));
                args.push_back(body != nullptr ? body->instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.monitor_callback.callv(args);
            }
            for (const uint64_t body_id : area.prev_inside_bodies) {
                if (std::find(inside.begin(), inside.end(), body_id) != inside.end()) {
                    continue;
                }
                if (!area.monitor_callback.is_valid()) {
                    break;
                }
                const BodyData *body = find_body(make_rid(body_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_REMOVED));
                args.push_back(make_rid(body_id));
                args.push_back(body != nullptr ? body->instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.monitor_callback.callv(args);
            }
            area.prev_inside_bodies = inside;
        } else {
            area.prev_inside_bodies.clear();
        }

        // --- areas monitoring areas (only when the other area is monitorable) ---
        if (area.area_monitor_callback.is_valid()) {
            std::vector<uint64_t> areas_inside;
            for (auto &[other_id, other] : areas) {
                if (other_id == area_id || !other.space.is_valid() ||
                        id_of(other.space) != id_of(area.space)) {
                    continue;
                }
                if (!other.monitorable) {
                    continue;
                }
                if ((other.collision_layer & area.collision_mask) == 0) {
                    continue;
                }
                if (other.shape_disabled.size() < other.shapes.size()) {
                    other.shape_disabled.resize(other.shapes.size(), false);
                }
                bool hits = false;
                for (std::size_t i = 0; i < other.shapes.size() && !hits; ++i) {
                    if (other.shape_disabled[i]) {
                        continue;
                    }
                    const Transform3D posed = other.transform * other.shape_transforms[i];
                    collision::Shape shape;
                    if (!query_shape_from_data(other.shapes[i], posed, shape)) {
                        continue;
                    }
                    const Placed other_shape{shape, posed};
                    for (const Placed &shape_area : area_shapes) {
                        if (overlap(shape_area, other_shape)) {
                            hits = true;
                            break;
                        }
                    }
                }
                if (hits) {
                    areas_inside.push_back(other_id);
                }
            }

            auto area_was_inside = [&areas_inside, &area](uint64_t id) -> bool {
                return std::find(area.prev_inside_areas.begin(), area.prev_inside_areas.end(), id)
                        != area.prev_inside_areas.end();
            };
            for (const uint64_t other_id : areas_inside) {
                if (area_was_inside(other_id) || !area.area_monitor_callback.is_valid()) {
                    continue;
                }
                const AreaData *other = find_area(make_rid(other_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_ADDED));
                args.push_back(make_rid(other_id));
                args.push_back(other != nullptr ? other->instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.area_monitor_callback.callv(args);
            }
            for (const uint64_t other_id : area.prev_inside_areas) {
                if (std::find(areas_inside.begin(), areas_inside.end(), other_id) != areas_inside.end()) {
                    continue;
                }
                if (!area.area_monitor_callback.is_valid()) {
                    break;
                }
                const AreaData *other = find_area(make_rid(other_id));
                Array args;
                args.push_back(static_cast<int>(PhysicsServer3D::AREA_BODY_REMOVED));
                args.push_back(make_rid(other_id));
                args.push_back(other != nullptr ? other->instance_id : 0);
                args.push_back(0);
                args.push_back(0);
                area.area_monitor_callback.callv(args);
            }
            area.prev_inside_areas = areas_inside;
        } else {
            area.prev_inside_areas.clear();
        }
    }
}
