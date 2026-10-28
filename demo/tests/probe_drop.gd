extends SceneTree

# Single-impact probe: a 1x0.5x1 box dropped onto the ground plane with a controlled
# downward speed. Prints the per-tick penetration so the contact response of two solver
# branches can be compared tick by tick.
#
#   godot --headless --path demo --script res://tests/probe_drop.gd -- --speed=8

var _speed := 8.0

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--speed="):
			_speed = float(arg.split("=")[1])
		elif arg.begins_with("--internal="):
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", float(arg.split("=")[1]))
		elif arg.begins_with("--substeps="):
			ProjectSettings.set_setting("physics/avbd/substeps", int(arg.split("=")[1]))
		elif arg.begins_with("--spec="):
			ProjectSettings.set_setting("physics/avbd/speculative_margin", float(arg.split("=")[1]))
		elif arg.begins_with("--alpha="):
			ProjectSettings.set_setting("physics/avbd/alpha", float(arg.split("=")[1]))
		elif arg.begins_with("--proj="):
			ProjectSettings.set_setting("physics/avbd/position_projection", arg.split("=")[1] == "1")
		elif arg.begins_with("--beta="):
			ProjectSettings.set_setting("physics/avbd/beta_linear", float(arg.split("=")[1]))
	_run()

func _run() -> void:
	var ground := StaticBody3D.new()
	var gshape := CollisionShape3D.new()
	var gbox := BoxShape3D.new()
	gbox.size = Vector3(10, 1, 10)
	gshape.shape = gbox
	ground.add_child(gshape)
	ground.position = Vector3(0, -0.5, 0)
	root.add_child(ground)

	var body := RigidBody3D.new()
	var shape := CollisionShape3D.new()
	var box := BoxShape3D.new()
	box.size = Vector3(1, 0.5, 1)
	shape.shape = box
	body.add_child(shape)
	body.mass = 1.0
	body.can_sleep = false
	body.gravity_scale = 1.0
	body.position = Vector3(0, 0.26, 0)
	body.linear_velocity = Vector3(0, -_speed, 0)
	root.add_child(body)

	var max_pen := 0.0
	var max_v := 0.0
	for tick in 240:
		await physics_frame
		var y := body.global_position.y
		var pen: float = maxf(0.0, 0.25 - y)
		max_pen = maxf(max_pen, pen)
		max_v = maxf(max_v, body.linear_velocity.length())
		if tick < 24:
			print("tick=%d y=%.5f v=%.3f pen=%.5f" % [tick, y, body.linear_velocity.y, pen])
	print("RESULT speed=%.1f max_pen=%.5f max_v=%.2f final_y=%.5f" % [
		_speed, max_pen, max_v, body.global_position.y])
	quit()
