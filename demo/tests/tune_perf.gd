extends SceneTree

func _initialize():
	var args = OS.get_cmdline_user_args()
	var iters = 6
	var proj = 2
	var conv = 0.00001
	var substeps = 1
	var interleaved = false
	var specs = 0.01
	var grace = 0
	
	for a in args:
		if a.begins_with("--iters="):
			iters = int(a.split("=")[1])
		elif a.begins_with("--proj="):
			proj = int(a.split("=")[1])
		elif a.begins_with("--conv="):
			conv = float(a.split("=")[1])
		elif a.begins_with("--substeps="):
			substeps = int(a.split("=")[1])
		elif a.begins_with("--interleaved="):
			interleaved = a.split("=")[1] == "true"
		elif a.begins_with("--specs="):
			specs = float(a.split("=")[1])
		elif a.begins_with("--grace="):
			grace = int(a.split("=")[1])
	
	ProjectSettings.set_setting("physics/avbd/iterations", iters)
	ProjectSettings.set_setting("physics/avbd/projection_iterations", proj)
	ProjectSettings.set_setting("physics/avbd/convergence_threshold", conv)
	ProjectSettings.set_setting("physics/avbd/substeps", substeps)
	ProjectSettings.set_setting("physics/avbd/interleaved_substeps", interleaved)
	ProjectSettings.set_setting("physics/avbd/speculative_margin", specs)
	ProjectSettings.set_setting("physics/avbd/contact_grace_ticks", grace)
	
	print("SETTINGS: iters=%d proj=%d conv=%f substeps=%d interleaved=%s specs=%f grace=%d" % [iters, proj, conv, substeps, interleaved, specs, grace])
	
	# Run load test
	var result = []
	var output = []
	
	# Just print and quit - we'll run separately
	quit()
