extends SceneTree

func _initialize():
	var args = OS.get_cmdline_user_args()
	var iters = 4
	var proj = 1
	var substeps = 1
	var interleaved = false
	var specs = 0.005
	var grace = 2
	var conv = 0.00001
	var threads = 0
	for a in args:
		if a.begins_with("--iters="): iters = int(a.split("=")[1])
		elif a.begins_with("--proj="): proj = int(a.split("=")[1])
		elif a.begins_with("--substeps="): substeps = int(a.split("=")[1])
		elif a.begins_with("--interleaved="): interleaved = a.split("=")[1] == "true"
		elif a.begins_with("--specs="): specs = float(a.split("=")[1])
		elif a.begins_with("--grace="): grace = int(a.split("=")[1])
		elif a.begins_with("--conv="): conv = float(a.split("=")[1])
		elif a.begins_with("--threads="): threads = int(a.split("=")[1])
	
	ProjectSettings.set_setting("physics/avbd/iterations", iters)
	ProjectSettings.set_setting("physics/avbd/projection_iterations", proj)
	ProjectSettings.set_setting("physics/avbd/substeps", substeps)
	ProjectSettings.set_setting("physics/avbd/interleaved_substeps", interleaved)
	ProjectSettings.set_setting("physics/avbd/speculative_margin", specs)
	ProjectSettings.set_setting("physics/avbd/contact_grace_ticks", grace)
	ProjectSettings.set_setting("physics/avbd/convergence_threshold", conv)
	ProjectSettings.set_setting("physics/avbd/threads", threads)
	
	var scene = load("res://scenes/pyramid_eight.tscn")
	var inst = scene.instantiate()
	root.add_child(inst)
	await process_frame
	await physics_frame
	var steps = 150
	var t0 = Time.get_ticks_usec()
	for i in steps:
		await physics_frame
	var t1 = Time.get_ticks_usec()
	var avg_us = (t1 - t0) / float(steps)
	print("PYR8: iters=%d proj=%d sub=%d inter=%s specs=%f grace=%d conv=%g thr=%d avg_ms=%.2f" % [iters, proj, substeps, str(interleaved), specs, grace, conv, threads, avg_us/1000.0])
	quit()
