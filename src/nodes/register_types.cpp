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
static void register_physics_server() {
    PhysicsServer3DManager *manager = PhysicsServer3DManager::get_singleton();
    if (manager == nullptr) {
        return;
    }

    // The manager wants a callable that produces the server, so a small factory object holds one.
    PhysicsServerFactory *factory = memnew(PhysicsServerFactory);
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
            {"physics/avbd/newton_ratio", Variant(1.0)},
            {"physics/avbd/stiffness_decay", Variant(0.8)},
            {"physics/avbd/jacobian_rebuild_distance", Variant(0.0)},
            {"physics/avbd/chebyshev_rho", Variant(0.0)},
            {"physics/avbd/adaptive_max_iterations", Variant(16)},
            {"physics/avbd/adaptive_beta", Variant(false)},
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
                        name == StringName("physics/avbd/substeps") ||
                        name == StringName("physics/avbd/iterations"));
    }
}

static void uninitialize_avbd_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
}

static void initialize_avbd_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
    GDREGISTER_CLASS(AVBDSoftBody3D);
    GDREGISTER_CLASS(AVBDSoftWorld3D);
    GDREGISTER_CLASS(AVBDPhysicsServer3D);
    GDREGISTER_CLASS(PhysicsServerFactory);
    // The engine instantiates these two itself (body state sync, direct space state).
    GDREGISTER_CLASS(AVBDDirectBodyState3D);
    GDREGISTER_CLASS(AVBDDirectSpaceState3D);
    register_avbd_settings();
    register_physics_server();
}

extern "C" {
GDExtensionBool GDE_EXPORT avbd_library_init(
        GDExtensionInterfaceGetProcAddress p_get_proc_address,
        const GDExtensionClassLibraryPtr p_library,
        GDExtensionInitialization *r_initialization) {
    GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);
    init_obj.register_initializer(initialize_avbd_module);
    init_obj.register_terminator(uninitialize_avbd_module);
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);
    return init_obj.init();
}
}
