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

#include "nodes/avbd_soft_body3d.hpp"
#include "nodes/avbd_soft_world3d.hpp"
#include "server/avbd_direct_body_state3d.hpp"
#include "server/avbd_direct_space_state3d.hpp"
#include "server/avbd_physics_server3d.hpp"
#include "server/physics_server_factory.hpp"

#include <godot_cpp/classes/physics_server3d_manager.hpp>

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
