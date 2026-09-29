extends Node

# Multithreading throughput benchmark
# Expected: Linear scalability with thread count, approaching hardware concurrency

var thread_configs = [1, 2, 4, 8, 16]
var test_duration: float = 10.0  # seconds
var step_interval: float = 1.0 / 60.0

func _ready():
    print("AVBD physics server multithreading throughput benchmark")
    print("=" * 60)
    
    for threads in thread_configs:
        print("\n" + "=" * 60)
        print("Testing with %d threads" % threads)
        print("=" * 60)
        
        # Set threads for this benchmark
        PhysicsServer3D.get_singleton().set_threads(threads)
        
        # Create space and bodies
        var space_rid = PhysicsServer3D.get_singleton().space_create()
        PhysicsServer3D.get_singleton().space_set_active(space_rid, true)
        
        # Create 100 bodies
        var body_rids = []
        for i in range(100):
            var body_rid = PhysicsServer3D.get_singleton().body_create()
            PhysicsServer3D.get_singleton().body_set_mode(body_rid, PhysicsServer3D.BODY_MODE_RIGID)
            
            var y_pos = float(i) * 0.5
            PhysicsServer3D.get_singleton().body_set_state(
                body_rid,
                PhysicsServer3D.BODY_STATE_TRANSFORM,
                Transform3D(Vector3.ZERO, Vector3(0.1, 0.1, 0.1))
            )
            PhysicsServer3D.get_singleton().body_set_state(
                body_rid,
                PhysicsServer3D.BODY_STATE_LINEAR_VELOCITY,
                Vector3.ZERO
            )
            
            var shape_rid = PhysicsServer3D.get_singleton().shape_create(PhysicsServer3D.SHAPE_BOX)
            PhysicsServer3D.get_singleton().shape_set_data(shape_rid, Vector3(0.5, 0.5, 0.5))
            PhysicsServer3D.get_singleton().body_add_shape(body_rid, shape_rid, Transform3D.IDENTITY, Vector3.ZERO)
            PhysicsServer3D.get_singleton().body_set_param(body_rid, PhysicsServer3D.BODY_PARAM_GRAVITY_SCALE, 1.0)
            
            body_rids.append(body_rid)
        
        # Benchmark
        var steps = 0
        var start_time = Time.get_ticks_msec()
        
        while (Time.get_ticks_msec() - start_time) / 1000.0 < test_duration:
            PhysicsServer3D.get_singleton().step(step_interval)
            steps += 1
        
        var elapsed = (Time.get_ticks_msec() - start_time) / 1000.0
        var fps = steps / elapsed
        
        print("Results:")
        print("  Threads: %d" % threads)
        print("  Steps: %d" % steps)
        print("  Time: %.2f seconds" % elapsed)
        print("  FPS: %.2f" % fps)
        print("  Bodies: 100")
        
        # Cleanup
        for body_rid in body_rids:
            PhysicsServer3D.get_singleton().free_rid(body_rid)
        PhysicsServer3D.get_singleton().free_rid(space_rid)
    
    print("\n" + "=" * 60)
    print("BENCHMARK COMPLETE")
    print("=" * 60)
