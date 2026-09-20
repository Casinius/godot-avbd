/*
 * Factory for the AVBD physics server. See the header.
 */

#include "server/physics_server_factory.hpp"

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/memory.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include "server/avbd_physics_server3d.hpp"

using namespace godot;

Object *PhysicsServerFactory::create_server() {
    UtilityFunctions::print("AVBD physics server: factory invoked - the engine selected AVBD");
    return memnew(AVBDPhysicsServer3D);
}

void PhysicsServerFactory::_bind_methods() {
    // ClassDB::bind_method(D_METHOD("create_server"), &PhysicsServerFactory::create_server);
}
