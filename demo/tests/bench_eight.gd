extends SceneTree

func _initialize():
	# Round-budget sweep from the shell, so the same scene can be timed at the demo's own settings
	# and at the cheaper ones without editing project.godot:
	#   AVBD_BENCH_ITERATIONS=8 godot --headless --path demo --script res://tests/bench_eight.gd
	var iterations := OS.get_environment("AVBD_BENCH_ITERATIONS")
	if iterations != "":
		ProjectSettings.set_setting("physics/avbd/iterations", int(iterations))
		print("bench: iterations overridden to ", iterations)
	# The projection layer and its scope, same idea: the round budget and the layer are a pair, and the
	# cheap configuration on this scene is a low budget plus the layer.
	# Rest test tolerance: how far a body may travel over its window and still count as idle. Loosening
	# it is what makes a jittering pile actually sleep - and a sleeping body is one whose state callback
	# the server can skip.
	var sleep_tol := OS.get_environment("AVBD_BENCH_SLEEP_TOL")
	if sleep_tol != "":
		ProjectSettings.set_setting("physics/avbd/sleep_position_tolerance", float(sleep_tol))
		ProjectSettings.set_setting("physics/avbd/sleep_rotation_tolerance", float(sleep_tol) * 4.0)
		ProjectSettings.set_setting("physics/avbd/sleep_hysteresis_ratio", 6.0)
		print("bench: sleep tolerance overridden to ", sleep_tol)
	# Sampling rate (engine ticks per second). The loop below waits for `physics_frame`, so a scene that
	# keeps up reports exactly this period and one that does not reports its own cost: raising the rate
	# until the reported period exceeds 1000/rate is how much headroom the step actually has.
	var rate := OS.get_environment("AVBD_BENCH_SAMPLE_RATE")
	if rate != "":
		ProjectSettings.set_setting("physics/common/physics_ticks_per_second", int(rate))
		print("bench: physics_ticks_per_second overridden to ", rate)
	# Physics on the calling thread instead of Godot's physics thread. The project default runs the
	# solver on a separate thread, and then `physics_frame` paces on whatever the main thread is doing -
	# which is not the solver, so solver settings barely move the number. Set this before reading a
	# timing as "the cost of a solver step".
	if OS.get_environment("AVBD_BENCH_MAIN_THREAD") != "":
		ProjectSettings.set_setting("physics/3d/run_on_separate_thread", false)
		print("bench: physics on the main thread")
	# The solver's own clock (`avbd/tick_rate`), which is independent of the engine's sampling rate
	# (`physics/common/physics_ticks_per_second`): 0 follows the engine step, N pins the solver to N Hz
	# and turns the engine tick into a sampling rate with interpolation between solver steps.
	var tick := OS.get_environment("AVBD_BENCH_TICK_RATE")
	if tick != "":
		ProjectSettings.set_setting("physics/avbd/tick_rate", float(tick))
		print("bench: tick_rate overridden to ", tick)
	var threads := OS.get_environment("AVBD_BENCH_THREADS")
	if threads != "":
		ProjectSettings.set_setting("physics/avbd/threads", int(threads))
		print("bench: threads overridden to ", threads)
	var projection := OS.get_environment("AVBD_BENCH_PROJECTION")
	if projection != "":
		ProjectSettings.set_setting("physics/avbd/projection_rounds", int(projection))
		print("bench: projection rounds overridden to ", projection)
	var scene = load("res://scenes/pyramid_eight.tscn")
	if not scene:
		print("NO_SCENE")
		quit()
		return
	var inst = scene.instantiate()
	# The demo scenes attach `scenes/demo.gd`, which prints a stats line per physics frame - a full
	# traversal of the pyramid plus three formatted writes. That is real work and it is not the solver,
	# so a timing that includes it measures the print loop once the solver is cheap. Mute it to read the
	# solver's cost.
	if OS.get_environment("AVBD_BENCH_MUTE") != "":
		inst.set_script(null)
		print("bench: scene script muted")
	root.add_child(inst)
	await process_frame
	await physics_frame
	
	var steps = int(OS.get_environment("AVBD_BENCH_STEPS")) if OS.get_environment("AVBD_BENCH_STEPS") != "" else 200
	var t0 = Time.get_ticks_usec()
	for i in steps:
		await physics_frame
	var t1 = Time.get_ticks_usec()
	var avg_us = (t1 - t0) / float(steps)
	print("PYRAMID_EIGHT_AVG_US_PER_PHYSICS_STEP:", avg_us)
	print("PYRAMID_EIGHT_AVG_MS_PER_PHYSICS_STEP:", avg_us / 1000.0)

	# A fast number is worth nothing if the pile came apart to get it: report the reach of the scene,
	# which for this pyramid is tens of metres when it holds and hundreds when it does not.
	var worst := 0.0
	var escaped := 0
	var stack := [root]
	while not stack.is_empty():
		var node = stack.pop_back()
		if node is RigidBody3D:
			var reach: float = node.global_position.length()
			worst = maxf(worst, reach)
			if reach > 200.0:
				escaped += 1
		stack.append_array(node.get_children())
	print("PYRAMID_EIGHT_REACH:", worst, " ESCAPED:", escaped)
	quit()
