/*
 * Factory for the AVBD physics server. See the header.
 */

#include "server/physics_server_factory.hpp"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/memory.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include "server/avbd_physics_server3d.hpp"

using namespace godot;

PhysicsServer3D *PhysicsServerFactory::create_server() {
    UtilityFunctions::print("AVBD physics server: factory invoked - the engine selected AVBD");
    return memnew(AVBDPhysicsServer3D);
}

void PhysicsServerFactory::_bind_methods() {
    // The manager calls the registered callable to instantiate the server; without this
    // binding the call fails and Godot silently falls back to the built-in engine.
    ClassDB::bind_method(D_METHOD("create_server"), &PhysicsServerFactory::create_server);
}
