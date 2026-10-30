extends SceneTree

func _initialize():
	# Round-budget sweep from the shell, so the same scene can be timed at the demo's own settings
	# and at the cheaper ones without editing project.godot:
	#   AVBD_BENCH_ITERATIONS=8 godot --headless --path demo --script res://tests/bench_eight.gd
	var iterations := OS.get_environment("AVBD_BENCH_ITERATIONS")
	if iterations != "":
		ProjectSettings.set_setting("physics/avbd/iterations", int(iterations))
		print("bench: iterations overridden to ", iterations)
	var scene = load("res://scenes/pyramid_eight.tscn")
	if not scene:
		print("NO_SCENE")
		quit()
		return
	var inst = scene.instantiate()
	root.add_child(inst)
	await process_frame
	await physics_frame
	
	var steps = 200
	var t0 = Time.get_ticks_usec()
	for i in steps:
		await physics_frame
	var t1 = Time.get_ticks_usec()
	var avg_us = (t1 - t0) / float(steps)
	print("PYRAMID_EIGHT_AVG_US_PER_PHYSICS_STEP:", avg_us)
	print("PYRAMID_EIGHT_AVG_MS_PER_PHYSICS_STEP:", avg_us / 1000.0)
	quit()
