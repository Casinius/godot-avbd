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

// Offer AVBD as a selectable 3D engine ([physics] 3d/physics_engine="AVBD") so scenes built
// from Godot's own nodes simulate without AVBD classes. avbd_factory is freed at deinit.
static PhysicsServerFactory *avbd_factory = nullptr;

static void register_physics_server() {
    PhysicsServer3DManager *manager = PhysicsServer3DManager::get_singleton();
    if (manager == nullptr) {
        return;
    }

    // The manager keeps only the Callable, not the factory object, so we free the latter at
    // deinitialization (see uninitialize_avbd_module).
    PhysicsServerFactory *factory = memnew(PhysicsServerFactory);
    avbd_factory = factory;
    manager->register_server("AVBD", Callable(factory, "create_server"));
}

// Register every physics/avbd/* setting so the editor shows an "AVBD" group; defaults
// must match read_project_params in the server.
static void register_avbd_settings() {
    ProjectSettings *settings = ProjectSettings::get_singleton();
    if (settings == nullptr) {
        return;
    }

    const std::pair<const char *, Variant> entries[] = {
            {"physics/avbd/threads", Variant(0)},
            {"physics/avbd/iterations", Variant(8)},
            {"physics/avbd/alpha", Variant(0.99)},
            {"physics/avbd/beta_linear", Variant(10000.0)},
            {"physics/avbd/beta_angular", Variant(100.0)},
            {"physics/avbd/gamma", Variant(0.999)},
            {"physics/avbd/adaptive_beta", Variant(false)},
            {"physics/avbd/newton_ratio", Variant(1.0)},
            {"physics/avbd/stiffness_decay", Variant(1.0)},
            {"physics/avbd/substeps", Variant(2)},
            {"physics/avbd/jacobian_rebuild_distance", Variant(0.0)},
            {"physics/avbd/interleaved_substeps", Variant(false)},
            {"physics/avbd/rounds_per_substep", Variant(0)},
            {"physics/avbd/internal_ticks_per_second", Variant(0.0)},
            {"physics/avbd/max_internal_steps", Variant(4)},
            {"physics/avbd/position_projection", Variant(true)},
            {"physics/avbd/projection_iterations", Variant(4)},
            {"physics/avbd/projection_beta", Variant(1.0)},
            {"physics/avbd/projection_slop", Variant(0.008)},
            {"physics/avbd/projection_max_correction", Variant(0.0)},
            {"physics/avbd/sleep_linear_threshold", Variant(0.05)},
            {"physics/avbd/sleep_angular_threshold", Variant(0.05)},
            {"physics/avbd/sleep_time", Variant(1.0)},
            {"physics/avbd/contact_grace_ticks", Variant(0)},
            {"physics/avbd/speculative_margin", Variant(0.01)},
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
                name == StringName("physics/avbd/interleaved_substeps") ||
                        name == StringName("physics/avbd/internal_ticks_per_second") ||
                        name == StringName("physics/avbd/position_projection") ||
                        name == StringName("physics/avbd/substeps") ||
                        name == StringName("physics/avbd/iterations"));
    }
}

// Register at SERVERS level (as godot-jolt does): the ClassDB entry must outlive
// PhysicsServer3D::finish(); SCENE level deregisters early and breaks teardown/leaks.
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
    // SERVERS: engine-facing classes; the server must exist before the engine picks it,
    // and its ClassDB entry must outlive PhysicsServer3D::finish().
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
