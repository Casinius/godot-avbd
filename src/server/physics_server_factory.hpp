/*
 * The physics server manager expects a Callable that produces a server instance, so this is the
 * object that holds that method. It exists only for registration.
 */

#ifndef AVBD_PHYSICS_SERVER_FACTORY_HPP
#define AVBD_PHYSICS_SERVER_FACTORY_HPP

#include <godot_cpp/classes/object.hpp>
#include <godot_cpp/classes/physics_server3d.hpp>

namespace godot {

class AVBDPhysicsServer3D;

class PhysicsServerFactory : public Object {
    GDCLASS(PhysicsServerFactory, Object)

public:
    // Returns PhysicsServer3D* to match Godot's PhysicsServer3DManager interface. The engine
    // uses the returned value as a PhysicsServer3D*; AVBDPhysicsServer3D inherits from
    // PhysicsServer3DExtension which inherits from PhysicsServer3D.
    PhysicsServer3D* create_server();

protected:
    static void _bind_methods();
};

} // namespace godot

#endif // AVBD_PHYSICS_SERVER_FACTORY_HPP
