/*
 * The physics server manager expects a Callable that produces a server instance, so this is the
 * object that holds that method. It exists only for registration.
 */

#ifndef AVBD_PHYSICS_SERVER_FACTORY_HPP
#define AVBD_PHYSICS_SERVER_FACTORY_HPP

#include <godot_cpp/classes/object.hpp>
#include <godot_cpp/core/class_db.hpp>

namespace godot {

class AVBDPhysicsServer3D;

class PhysicsServerFactory : public Object {
    GDCLASS(PhysicsServerFactory, Object)

protected:
    static void _bind_methods();

public:
    // Returns Object* so the method binder has an unambiguous type; the engine casts it to the
    // server interface it asked for.
    Object *create_server();
};

} // namespace godot

#endif // AVBD_PHYSICS_SERVER_FACTORY_HPP
