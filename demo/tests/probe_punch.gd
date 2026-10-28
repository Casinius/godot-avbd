extends "res://tests/avbd_harness.gd"

const SCENE := "res://scenes/pyramid_dbl.tscn"

var _steps_total := 600
var _warmup := 120
var _config := ""

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--steps="):
			_steps_total = int(arg.split("=")[1])
		elif arg.begins_with("--warmup="):
			_warmup = int(arg.split("=")[1])
		elif arg.begins_with("--config="):
			_config = arg.split("=")[1]
	match _config:
		"", "current":
			pass
		"proj":
			ProjectSettings.set_setting("physics/avbd/position_projection", true)
		"proj32":
			ProjectSettings.set_setting("physics/avbd/position_projection", true)
			ProjectSettings.set_setting("physics/avbd/projection_iterations", 32)
		"internal-240":
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", 0.0)
		"internal-240-s1":
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", 0.0)
			ProjectSettings.set_setting("physics/avbd/substeps", 1)
		"no-interleave":
			ProjectSettings.set_setting("physics/avbd/interleaved_substeps", false)
		"beta-1e6":
			ProjectSettings.set_setting("physics/avbd/beta_linear", 1000000.0)
		"beta-1e5":
			ProjectSettings.set_setting("physics/avbd/beta_linear", 100000.0)
		"beta-1e4":
			ProjectSettings.set_setting("physics/avbd/beta_linear", 10000.0)
		"beta-1e3":
			ProjectSettings.set_setting("physics/avbd/beta_linear", 1000.0)
		"no-spec":
			ProjectSettings.set_setting("physics/avbd/speculative_margin", 0.0)
		"alpha-09":
			ProjectSettings.set_setting("physics/avbd/alpha", 0.9)
		"alpha-08":
			ProjectSettings.set_setting("physics/avbd/alpha", 0.8)
		"alpha-05":
			ProjectSettings.set_setting("physics/avbd/alpha", 0.5)
		"iters-32":
			ProjectSettings.set_setting("physics/avbd/iterations", 32)
		"s1-iters64":
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", 0.0)
			ProjectSettings.set_setting("physics/avbd/substeps", 1)
			ProjectSettings.set_setting("physics/avbd/iterations", 64)
		"s1-iters12":
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", 0.0)
			ProjectSettings.set_setting("physics/avbd/substeps", 1)
			ProjectSettings.set_setting("physics/avbd/iterations", 12)
		_:
			push_error("unknown config")
	super._initialize()

func _scenarios() -> Array:
	return [scenario("impact", run_impact)]

func _collect(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect(child, out)

func run_impact() -> Variant:
	var scene := (load(SCENE) as PackedScene).instantiate()
	scene.set_script(null)
	root.add_child(scene)
	var bodies: Array[RigidBody3D] = []
	_collect(scene, bodies)
	for body in bodies:
		body.mass = 1.0
		body.sleeping = false
		body.can_sleep = false

	var pens: Array[float] = []
	var min_y := 1.0e9
	var max_abs_x := 0.0
	var max_speed := 0.0
	var speed_sum := 0.0
	for step in _steps_total:
		await steps(1)
		var pen: float = PhysicsServer3D.call("get_solver_max_penetration")
		if step >= _warmup:
			pens.append(pen)
		for b in bodies:
			min_y = minf(min_y, b.global_position.y)
			max_abs_x = maxf(max_abs_x, absf(b.global_position.x))
			var v := b.linear_velocity.length()
			speed_sum += v
			max_speed = maxf(max_speed, v)
	pens.sort()
	var n := pens.size()
	var p50 := pens[n / 2] if n > 0 else 0.0
	var p95 := pens[int(n * 0.95)] if n > 0 else 0.0
	var pmax := pens[n - 1] if n > 0 else 0.0
	var frames_over := 0
	for p in pens:
		if p > 0.05:
			frames_over += 1
	print("cfg=%-15s p50=%.4f p95=%.4f max=%.4f over5cm=%d/%d min_y=%.3f max|x|=%.1f max_v=%.2f avg_v=%.2f" % [
		_config, p50, p95, pmax, frames_over, n, min_y, max_abs_x, max_speed,
		speed_sum / float(maxi(1, _steps_total * bodies.size()))])
	scene.queue_free()
	await frames(2)
	return {}
