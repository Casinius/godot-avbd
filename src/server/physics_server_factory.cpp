/*
 * Factory for the AVBD physics server. See the header.
 */

#include "server/physics_server_factory.hpp"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/memory.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include "server/avbd_physics_server3d.hpp"

using namespace godot;

PhysicsServer3D* PhysicsServerFactory::create_server() {
    UtilityFunctions::print("AVBD physics server: factory invoked - the engine selected AVBD");
    return memnew(AVBDPhysicsServer3D);
}

void PhysicsServerFactory::_bind_methods() {
    // The engine invokes this through a Callable, so the method has to be bound to
    // ClassDB - without this the callable is invalid, PhysicsServer3DManager silently
    // falls back to the default engine and "AVBD" runs Godot Physics.
    ClassDB::bind_method(D_METHOD("create_server"), &PhysicsServerFactory::create_server);
}
