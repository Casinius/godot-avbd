/*
 * AVBD (Augmented Vertex Block Descent) for Godot 4 - GDExtension entry point.
 *
 * Part of godot-avbd. The solver core in src/avbd is ported from the authors'
 * reference implementation (github.com/savant117/avbd-demo3d, MIT).
 */

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/godot.hpp>

#include "nodes/avbd_constraint3d.hpp"
#include "nodes/avbd_ignore_collision3d.hpp"
#include "nodes/avbd_joint3d.hpp"
#include "nodes/avbd_rigid_body3d.hpp"
#include "nodes/avbd_soft_body3d.hpp"
#include "nodes/avbd_spring3d.hpp"
#include "nodes/avbd_world3d.hpp"

using namespace godot;

static void initialize_avbd_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
    }
    GDREGISTER_CLASS(AVBDWorld3D);
    GDREGISTER_CLASS(AVBDRigidBody3D);
    GDREGISTER_CLASS(AVBDSoftBody3D);
    GDREGISTER_ABSTRACT_CLASS(AVBDConstraint3D);
    GDREGISTER_CLASS(AVBDJoint3D);
    GDREGISTER_CLASS(AVBDSpring3D);
    GDREGISTER_CLASS(AVBDIgnoreCollision3D);
}

static void uninitialize_avbd_module(ModuleInitializationLevel p_level) {
    if (p_level != MODULE_INITIALIZATION_LEVEL_SCENE) {
        return;
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
    init_obj.set_minimum_library_initialization_level(MODULE_INITIALIZATION_LEVEL_SCENE);
    return init_obj.init();
}
}
