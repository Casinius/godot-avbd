/*
 * AVBD (Augmented Vertex Block Descent) for Godot 4 - GDExtension entry point.
 *
 * Part of godot-avbd. The solver core in src/avbd is ported from the authors'
 * reference implementation (github.com/savant117/avbd-demo3d, MIT).
 *
 * Rigid bodies, joints, areas and space queries run through the AVBD physics server
 * (standard Godot nodes, physics/3d/physics_engine = "AVBD"). The only node classes
 * left here are the soft-body pair: AVBDSoftBody3D (no server representation yet)
 * and AVBDSoftWorld3D, which owns their solver.
 */

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/godot.hpp>

#include <godot_cpp/classes/project_settings.hpp>
#include <godot_cpp/classes/physics_server3d_manager.hpp>
#include <godot_cpp/variant/variant.hpp>

#include "nodes/avbd_soft_body3d.hpp"
#include "nodes/avbd_soft_world3d.hpp"
#include "server/avbd_direct_body_state3d.hpp"
#include "server/avbd_direct_space_state3d.hpp"
#include "server/avbd_physics_server3d.hpp"
#include "server/physics_server_factory.hpp"

using namespace godot;

// Offer AVBD to the engine as a selectable 3D physics engine. This is what lets a scene built
// from Godot's own nodes - RigidBody3D, CollisionShape3D and so on - be simulated by AVBD, without
// the project adopting any AVBD node classes. Selecting it is a project setting:
//
//   [physics]
//   3d/physics_engine="AVBD"
// The physics server factory the manager's creation Callable points at; freed at
// deinitialization because nothing else owns it.
static PhysicsServerFactory *avbd_factory = nullptr;

static void register_physics_server() {
    PhysicsServer3DManager *manager = PhysicsServer3DManager::get_singleton();
    if (manager == nullptr) {
        return;
    }

    // The manager wants a callable that produces the server, so a small factory object holds one.
    // The manager only keeps the Callable, which does not own the object - the factory is
    // freed at deinitialization (see uninitialize_avbd_module).
    PhysicsServerFactory *factory = memnew(PhysicsServerFactory);
    avbd_factory = factory;
    manager->register_server("AVBD", Callable(factory, "create_server"));
}

// Register every physics/avbd/* solver setting with ProjectSettings, so the editor's
// Project Settings dialog shows an "AVBD" group instead of the settings existing only
// after a script writes them. Defaults must match read_project_params in the server.
static void register_avbd_settings() {
    ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    const std::pair<const char *, Variant> entries[] = {
            {"physics/avbd/threads", Variant(0)},
            {"physics/avbd/iterations", Variant(4)},
            {"physics/avbd/alpha", Variant(0.99)},
            {"physics/avbd/beta_linear", Variant(10000.0)},
            {"physics/avbd/beta_angular", Variant(100.0)},
            {"physics/avbd/gamma", Variant(0.999)},
            {"physics/avbd/substeps", Variant(2)},
            // 0 = the engine step is the clock, split across `substeps`. A project that wants the
            // solver's own fixed tick sets a rate here (the suite measurements below were taken
            // at 60 and 120) and physics_ticks_per_second becomes a sampling rate for it.
            {"physics/avbd/tick_rate", Variant(0.0)},
            // Ceiling on the ticks one engine step may run when the clock has fallen behind.
            // Time past this is dropped, which is what keeps a stall from becoming a longer
            // internal step.
            {"physics/avbd/max_catchup_ticks", Variant(4)},
            // Delivery smoothing, not a solver knob: the scene is handed the blend of the last two
            // deliveries instead of the newest pose. Off keeps every existing project's numbers.
            {"physics/avbd/interpolate_poses", Variant(false)},
            {"physics/avbd/newton_ratio", Variant(1.0)},
            {"physics/avbd/stiffness_decay", Variant(0.8)},
            {"physics/avbd/jacobian_rebuild_distance", Variant(0.0)},
            {"physics/avbd/chebyshev_rho", Variant(0.0)},
            {"physics/avbd/adaptive_max_iterations", Variant(16)},
            {"physics/avbd/adaptive_beta", Variant(false)},
            // The position projection layer: a positional pass over the constraints that declares
            // itself hard (contacts, locked and limited joint axes), so a delivered pose is legal
            // even where a penalty method can only approach one. Rounds is the switch and 0 keeps
            // a project exactly as it was. What it is worth turning it on for is split by scope,
            // hence projection_joints and projection_contacts: the joint half holds up chains that
            // the plain solver loses (Bridge, Rope), the contact half buys a lower residual at the
            // cost of squirrelling deep piles. See SolverParams::Projection for the numbers.
            {"physics/avbd/projection_rounds", Variant(0)},
            {"physics/avbd/projection_scale", Variant(0.9)},
            // Ceiling on how far one body may be moved by one pass, metres. It is a per-body
            // ceiling rather than a per-contact one, so a body pinched between many contacts is
            // moved once per pass and not once per pinch. Keep it small: a pass that hits the
            // ceiling is a teleport, and a teleport is what breaks a jointed chain.
            {"physics/avbd/projection_max_correction", Variant(0.003)},
            // 0 = the projected move is given back out of the step's velocity, so the layer costs
            // no momentum; 1 = standard PBD, where the move *is* the velocity. Measured: 0 leaves
            // the pose and the velocity disagreeing and the jointed scenes diverge.
            {"physics/avbd/projection_velocity_transfer", Variant(1.0)},
            // Charge the correction to the constraint's force multiplier so the next round's
            // penalty ramp does not buy the same lift twice.
            {"physics/avbd/projection_dual_feedback", Variant(true)},
            {"physics/avbd/projection_joints", Variant(true)},
            {"physics/avbd/projection_contacts", Variant(true)},
    };
    for (const auto &entry : entries) {
        if (settings->has_setting(entry.first)) {
            continue;
        }
        settings->set_setting(entry.first, entry.second);
        settings->set_initial_value(entry.first, entry.second);
        // The knobs people reach for first are basic; the rest live under "Advanced".
        const StringName name(entry.first);
        settings->set_as_basic(name,
                name == StringName("physics/avbd/adaptive_beta") ||
                        name == StringName("physics/avbd/projection_rounds") ||
                        name == StringName("physics/avbd/substeps") ||
                        name == StringName("physics/avbd/tick_rate") ||
                        name == StringName("physics/avbd/iterations"));
    }
}

// Everything registers at the SERVERS level (as godot-jolt does): the physics server
// object must keep its ClassDB entry until PhysicsServer3D::finish() has called
// _finish on it at exit. Registered at SCENE level the classes deregister with the
// main loop first, and teardown fails with "Required virtual method ::_finish must be
// overridden" plus a leaked instance.
static void uninitialize_avbd_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SERVERS) {
        return;
    }
    if (avbd_factory != nullptr) {
        memdelete(avbd_factory);
        avbd_factory = nullptr;
    }
}

static void initialize_avbd_module(ModuleInitializationLevel p_level) {
    // SERVERS: the engine-facing classes (the physics server must exist before the
    // engine picks its physics engine, and its ClassDB entry must outlive
    // PhysicsServer3D::finish()).
    if (p_level == MODULE_INITIALIZATION_LEVEL_SERVERS) {
        GDREGISTER_CLASS(AVBDPhysicsServer3D);
        GDREGISTER_CLASS(PhysicsServerFactory);
        // The engine instantiates these two itself (body state sync, direct space state).
        GDREGISTER_CLASS(AVBDDirectBodyState3D);
        GDREGISTER_CLASS(AVBDDirectSpaceState3D);
        register_avbd_settings();
        register_physics_server();
    }
    // SCENE: the soft-body pair derive from Node3D/MultiMeshInstance3D, which do not
    // exist yet at the SERVERS level.
    if (p_level == MODULE_INITIALIZATION_LEVEL_SCENE) {
        GDREGISTER_CLASS(AVBDSoftBody3D);
        GDREGISTER_CLASS(AVBDSoftWorld3D);
    }
}

extern "C" {
GDExtensionBool GDE_EXPORT avbd_library_init(
        GDExtensionInterfaceGetProcAddress p_get_proc_address,
        const GDExtensionClassLibraryPtr p_library,
        GDExtensionInitialization *r_initialization) {
    GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);
    init_obj.register_initializer(initialize_avbd_module);
    init_obj.register_terminator(uninitialize_avbd_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SERVERS);
    return init_obj.init();
}
}
