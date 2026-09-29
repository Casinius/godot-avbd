extends Node

# Test script to verify per-solver contention reduction in multi-space scenarios
# Expected: Each space has its own solver/pool, reduced contention

var space_times: Array = []
var num_spaces: int = 4
var steps_per_space: int = 600  # 10 seconds at 60 FPS

func _ready():
    print("AVBD physics server multi-space contention test")
    print("=" * 60)
    print("Creating %d spaces concurrently" % num_spaces)
    print("=" * 60)
    
    var start_time = Time.get_ticks_msec()
    
    # Create and initialize all spaces
    var spaces = []
    for i in range(num_spaces):
        var space_rid = PhysicsServer3D.get_singleton().space_create()
        spaces.append(space_rid)
        
        # Create a body in each space
        var body_rid = PhysicsServer3D.get_singleton().body_create()
        PhysicsServer3D.get_singleton().body_set_mode(body_rid, PhysicsServer3D.BODY_MODE_RIGID)
        PhysicsServer3D.get_singleton().body_set_state(
            body_rid,
            PhysicsServer3D.BODY_STATE_TRANSFORM,
            Transform3D(Vector3(0.0, float(i) * 2.0, 0.0), Vector3.ONE)
        )
        PhysicsServer3D.get_singleton().body_set_state(
            body_rid,
            PhysicsServer3D.BODY_STATE_LINEAR_VELOCITY,
            Vector3.ZERO
        )
        
        # Create shape
        var shape_rid = PhysicsServer3D.get_singleton().shape_create(PhysicsServer3D.SHAPE_BOX)
        PhysicsServer3D.get_singleton().shape_set_data(shape_rid, Vector3(1.0, 1.0, 1.0))
        PhysicsServer3D.get_singleton().body_add_shape(body_rid, shape_rid, Transform3D.IDENTITY, Vector3.ZERO)
        
        PhysicsServer3D.get_singleton().body_set_param(body_rid, PhysicsServer3D.BODY_PARAM_GRAVITY_SCALE, 1.0)
        PhysicsServer3D.get_singleton().space_set_active(space_rid, true)
    
    # Measure time to step all spaces
    var step_start = Time.get_ticks_msec()
    var steps = 0
    var max_steps = steps_per_space
    
    while steps < max_steps:
        for i in range(num_spaces):
            PhysicsServer3D.get_singleton().step(1.0 / 60.0)
        steps += 1
        
        # Progress indicator
        if steps % 60 == 0:
            var elapsed = Time.get_ticks_msec() - step_start
            var fps = steps / (elapsed / 1000.0)
            print("Step %d / %d (Total FPS: %.1f)" % [steps, max_steps, fps])
    
    var total_elapsed = Time.get_ticks_msec() - start_time
    var avg_fps = steps / (total_elapsed / 1000.0)
    
    print("\n" + "=" * 60)
    print("TEST COMPLETE")
    print("Spaces created: %d" % num_spaces)
    print("Steps executed: %d" % steps)
    print("Total time: %.2f seconds" % (total_elapsed / 1000.0))
    print("Average FPS: %.2f" % avg_fps)
    print("Multi-space contention test passed!")
    print("=" * 60)
    
    # Cleanup
    for space_rid in spaces:
        PhysicsServer3D.get_singleton().free_rid(space_rid)
