extends Node

# Test script to verify AVBD physics server initialization
# Expected: No errors; physics simulation starts successfully

var start_time: float = 0.0
var steps: int = 0
var max_steps: int = 600  # 10 seconds at 60 FPS

func _ready():
    print("AVBD physics server initialization test")
    print("=" * 60)
    
    # Get the physics server
    var physics_server = PhysicsServer3D.get_singleton()
    if physics_server == null:
        print("ERROR: PhysicsServer3D not available")
        return
        
    print("Physics server obtained: %s" % physics_server.get_class())
    
    # Create a space
    var space_rid = physics_server.space_create()
    print("Space created (RID: %s)" % space_rid)
    
    # Create a rigid body
    var body_rid = physics_server.body_create()
    print("Body created (RID: %s)" % body_rid)
    
    # Set body type to rigid
    physics_server.body_set_mode(body_rid, PhysicsServer3D.BODY_MODE_RIGID)
    print("Body mode set to RIGID")
    
    # Set position and linear velocity
    physics_server.body_set_state(body_rid, PhysicsServer3D.BODY_STATE_TRANSFORM, Transform3D(Vector3.ZERO, Vector3.ONE))
    physics_server.body_set_state(body_rid, PhysicsServer3D.BODY_STATE_LINEAR_VELOCITY, Vector3.ZERO)
    print("Body initial state set")
    
    # Create a simple shape (box)
    var shape_rid = physics_server.shape_create(PhysicsServer3D.SHAPE_BOX)
    var box_size = Vector3(1.0, 1.0, 1.0)
    physics_server.shape_set_data(shape_rid, box_size)
    print("Shape created (RID: %s, type: BOX)" % shape_rid)
    
    # Attach shape to body
    physics_server.body_add_shape(body_rid, shape_rid, Transform3D.IDENTITY, Vector3.ZERO)
    print("Shape attached to body")
    
    # Set gravity scale
    physics_server.body_set_param(body_rid, PhysicsServer3D.BODY_PARAM_GRAVITY_SCALE, 1.0)
    print("Gravity scale set to 1.0")
    
    # Activate the space
    physics_server.space_set_active(space_rid, true)
    print("Space activated")
    
    # Test stepping
    print("\nStarting physics steps...")
    start_time = Time.get_ticks_msec()
    
    while steps < max_steps:
        # Step the simulation by a fixed timestep
        physics_server.step(1.0 / 60.0)
        
        steps += 1
        
        # Progress indicator
        if steps % 60 == 0:
            var elapsed = Time.get_ticks_msec() - start_time
            var fps = steps / (elapsed / 1000.0)
            print("Step %d / %d (FPS: %.1f)" % [steps, max_steps, fps])
    
    var elapsed = Time.get_ticks_msec() - start_time
    var avg_fps = steps / (elapsed / 1000.0)
    
    print("\n" + "=" * 60)
    print("TEST COMPLETE")
    print("Steps executed: %d" % steps)
    print("Total time: %.2f seconds" % (elapsed / 1000.0))
    print("Average FPS: %.2f" % avg_fps)
    print("AVBD physics server initialized successfully!")
    print("=" * 60)
    
    # Cleanup
    physics_server.free_rid(shape_rid)
    physics_server.free_rid(body_rid)
    physics_server.free_rid(space_rid)
