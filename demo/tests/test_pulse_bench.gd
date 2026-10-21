# Headless pulse benchmark: one scene, one pulse drive, any physics engine.
#
# Everything here goes through Godot's own RigidBody3D nodes, so the exact same run
# measures AVBD, Jolt Physics or Godot Physics - only project.godot's
# `physics/3d/physics_engine` decides which server is behind the nodes. That is what
# makes the two runs comparable: same bodies, same impulses, same tick rate, same
# measurement code.
#
#   godot --headless --path demo --script res://tests/test_pulse_bench.gd
#   godot --headless --path demo --script res://tests/test_pulse_bench.gd -- --ticks=120,240,480
#   ./scripts/pulse_bench.sh                          # AVBD vs Jolt, side by side
#
# Output is one `pulse-bench:` line per tick rate, in key=value form, so the two
# engines' runs can be diffed line by line.
#
# What is measured, per tick rate:
#   wall_ms_per_step  clock time for one physics step while the engine keeps pace.
#                     Below budget_ms this is the frame-pacing floor (headless runs at
#                     60 fps), above it the engine is the bottleneck.
#   physics_ms_per_frame  the engine's own physics time for one frame, from the
#                     Performance monitor - the same quantity the profiler shows.
#   realtime_factor   achieved steps per second / requested ticks. 1.0 means the
#                     engine kept up with the tick rate; below 1.0 it is dropping steps
#                     (Godot caps catch-up at max_physics_steps_per_frame).
#   plus stability: max speed, max angular speed, worst position, sleeping bodies,
#   active objects, collision pairs and a pose digest.
#
# Bodies never sleep (can_sleep = false) so all three engines carry the same load:
# their sleep thresholds disagree, and a sleeping body would make one run artificially
# cheap. The pulses are applied from the main thread between physics steps, which both
# servers queue the same way, so the drive does not favour either engine.
extends SceneTree

const PATTERNS := ["rest", "spin", "push", "both"]

var _opts := {
	"ticks": "60,120,240,360,480",
	"pattern": "both",
	"grid": 23,              # square grid side; grid^2 bodies
	"gap": 1.02,             # spacing between box centres (boxes are 1 m^3)
	"steps": 240,            # measured physics steps per tick rate
	"warmup": 60,            # physics steps before the clock starts
	"pulse_every": 2,        # 0 = no pulses at all
	"impulse": 5.0,          # N s, along +/-X, alternating
	"torque": 5.0,           # N m s, around +/-Y, alternating
	"iterations": 0,         # 0 = leave each engine on its own defaults
	"max_steps_per_frame": 0,  # 0 = use the project's setting
	"probe": true,           # cheap engine sanity probe (see _engine_probe)
}

var _ground: StaticBody3D = null
var _bodies: Array[RigidBody3D] = []
var _initial: Array[Vector3] = []


# --- argument handling ------------------------------------------------------

func _parse_args() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--ticks="):
			_opts.ticks = arg.split("=")[1]
		elif arg.begins_with("--pattern="):
			_opts.pattern = arg.split("=")[1]
			if not PATTERNS.has(_opts.pattern):
				push_error("unknown pattern %s (want one of %s)" % [_opts.pattern, str(PATTERNS)])
				quit(2)
		elif arg.begins_with("--grid="):
			_opts.grid = int(arg.split("=")[1])
		elif arg.begins_with("--gap="):
			_opts.gap = float(arg.split("=")[1])
		elif arg.begins_with("--steps="):
			_opts.steps = int(arg.split("=")[1])
		elif arg.begins_with("--warmup="):
			_opts.warmup = int(arg.split("=")[1])
		elif arg.begins_with("--pulse-every="):
			_opts.pulse_every = int(arg.split("=")[1])
		elif arg.begins_with("--impulse="):
			_opts.impulse = float(arg.split("=")[1])
		elif arg.begins_with("--torque="):
			_opts.torque = float(arg.split("=")[1])
		elif arg.begins_with("--iterations="):
			_opts.iterations = int(arg.split("=")[1])
		elif arg.begins_with("--max-steps="):
			_opts.max_steps_per_frame = int(arg.split("=")[1])
		elif arg == "--":
			pass
		elif arg == "--no-probe":
			_opts.probe = false
		else:
			push_error("unknown argument %s" % arg)
			quit(2)


# --- engine pinning ---------------------------------------------------------

# `--iterations=N` cannot mean the same thing to every engine: AVBD runs N full
# primal+dual rounds, Jolt splits its work into velocity and position steps, Godot
# Physics exposes a single solver iteration count. Pin each one to N rounds of work so
# a pinned run is a controlled comparison; unpinned runs use each engine's defaults,
# which is what a user gets out of the box.
func _pin_iterations(iterations: int) -> String:
	if iterations <= 0:
		return "defaults"
	var engine := String(ProjectSettings.get_setting("physics/3d/physics_engine"))
	match engine:
		"AVBD":
			ProjectSettings.set_setting("physics/avbd/iterations", iterations)
			return "avbd/iterations=%d" % iterations
		"Jolt Physics":
			ProjectSettings.set_setting("physics/jolt_physics_3d/simulation/velocity_steps", iterations)
			ProjectSettings.set_setting("physics/jolt_physics_3d/simulation/position_steps", 1)
			return "jolt velocity_steps=%d position_steps=1" % iterations
		_:
			ProjectSettings.set_setting("physics/3d/solver/solver_iterations", iterations)
			return "solver_iterations=%d" % iterations


# --- scene ------------------------------------------------------------------

func _build_scene() -> void:
	_ground = StaticBody3D.new()
	_ground.name = "Ground"
	var ground_shape := CollisionShape3D.new()
	var ground_box := BoxShape3D.new()
	ground_box.size = Vector3(200, 1, 200)
	ground_shape.shape = ground_box
	_ground.add_child(ground_shape)
	var ground_material := PhysicsMaterial.new()
	ground_material.friction = 0.5
	_ground.physics_material_override = ground_material
	_ground.position = Vector3(0, -0.5, 0)
	root.add_child(_ground)

	# A single layer of grid^2 boxes, one body per cell, spacing slightly larger than
	# the box size so every box starts exactly touching its neighbours and the ground.
	var side := int(_opts.grid)
	var gap := float(_opts.gap)
	var offset := -0.5 * gap * float(side - 1)
	for x in side:
		for z in side:
			var body := RigidBody3D.new()
			body.name = "Box%d_%d" % [x, z]
			var shape := CollisionShape3D.new()
			var box := BoxShape3D.new()
			box.size = Vector3.ONE
			shape.shape = box
			body.add_child(shape)
			var material := PhysicsMaterial.new()
			material.friction = 0.5
			body.physics_material_override = material
			body.mass = 1.0
			body.can_sleep = false
			body.position = Vector3(offset + gap * x, 0.5, offset + gap * z)
			root.add_child(body)
			_bodies.append(body)
			_initial.append(body.position)


func _free_scene() -> void:
	for body in _bodies:
		body.queue_free()
	_bodies.clear()
	_initial.clear()
	if _ground != null:
		_ground.queue_free()
		_ground = null
	await process_frame
	await process_frame


func _reset_scene() -> void:
	for i in _bodies.size():
		_bodies[i].position = _initial[i]
		_bodies[i].rotation = Vector3.ZERO
		_bodies[i].linear_velocity = Vector3.ZERO
		_bodies[i].angular_velocity = Vector3.ZERO
		_bodies[i].sleeping = false


# --- pulse drive ------------------------------------------------------------

func _apply_pulse(sign: float) -> void:
	var spin: bool = _opts.pattern == "spin" or _opts.pattern == "both"
	var push: bool = _opts.pattern == "push" or _opts.pattern == "both"
	for body in _bodies:
		if push:
			body.apply_impulse(Vector3(sign * float(_opts.impulse), 0.0, 0.0))
		if spin:
			body.apply_torque_impulse(Vector3(0.0, sign * float(_opts.torque), 0.0))


# --- measurement ------------------------------------------------------------

func _steps(count: int) -> void:
	for i in count:
		await physics_frame


func _digest() -> int:
	# Pose digest, quantised to a millimetre so a one-tick measurement offset between
	# runs does not flip it while real divergence (centimetres) does.
	var hash := 5381
	for body in _bodies:
		# Local transform is the same as the global one here: every body hangs directly
		# off the root viewport at the origin.
		var t := body.transform
		t.origin = (t.origin * 1000.0).round() / 1000.0
		for byte in var_to_bytes(t):
			hash = (hash * 33 + byte) & 0xFFFFFFFF
	return hash


func _measure(ticks: int) -> Dictionary:
	Engine.physics_ticks_per_second = ticks
	if _opts.max_steps_per_frame > 0:
		Engine.max_physics_steps_per_frame = _opts.max_steps_per_frame
	_reset_scene()
	await _steps(_opts.warmup)

	var t0 := Time.get_ticks_usec()
	var steps := 0
	while steps < _opts.steps:
		var every := int(_opts.pulse_every)
		if every > 0 and steps % every == 0:
			var sign := 1.0 if (steps / every) % 2 == 0 else -1.0
			_apply_pulse(sign)
		await physics_frame
		steps += 1
	var elapsed_us := Time.get_ticks_usec() - t0
	var physics_ms := Performance.get_monitor(Performance.TIME_PHYSICS_PROCESS)

	var max_speed := 0.0
	var max_angular := 0.0
	var max_position := 0.0
	var sleeping := 0
	var finite := true
	for body in _bodies:
		max_speed = maxf(max_speed, body.linear_velocity.length())
		max_angular = maxf(max_angular, body.angular_velocity.length())
		max_position = maxf(max_position, body.transform.origin.length())
		if body.sleeping:
			sleeping += 1
		if not body.transform.origin.is_finite() or not body.linear_velocity.is_finite():
			finite = false

	var elapsed_s := float(elapsed_us) / 1000000.0
	var budget_ms := 1000.0 / float(ticks)
	var achieved := steps / elapsed_s if elapsed_s > 0.0 else 0.0
	return {
		"engine": ProjectSettings.get_setting("physics/3d/physics_engine"),
		"pattern": _opts.pattern,
		"bodies": _bodies.size(),
		"ticks": ticks,
		"budget_ms": budget_ms,
		"wall_ms_per_step": 1000.0 * elapsed_s / float(steps),
		"physics_ms_per_frame": physics_ms,
		"realtime_factor": achieved / float(ticks),
		"max_speed": max_speed,
		"max_angular": max_angular,
		"max_position": max_position,
		"sleeping": sleeping,
		"finite": finite,
		"active_objects": PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_ACTIVE_OBJECTS),
		"collision_pairs": PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_COLLISION_PAIRS),
		"digest": _digest(),
	}


# Engine names that are not registered silently fall back to the default engine, so a
# typo would produce two identical runs that look like a tie. A short stack of boxes is
# enough to tell them apart: the engines resolve a resting contact very differently.
func _engine_probe() -> float:
	# Ground first: without it the boxes simply fall, and every engine reports the same
	# gravitational distance instead of its own contact resolution.
	var ground := StaticBody3D.new()
	var ground_shape := CollisionShape3D.new()
	var ground_box := BoxShape3D.new()
	ground_box.size = Vector3(50, 1, 50)
	ground_shape.shape = ground_box
	ground.add_child(ground_shape)
	ground.position = Vector3(0, -0.5, 0)
	root.add_child(ground)

	var probe := []
	for i in 6:
		var body := RigidBody3D.new()
		var shape := CollisionShape3D.new()
		var box := BoxShape3D.new()
		box.size = Vector3.ONE
		shape.shape = box
		body.add_child(shape)
		body.mass = 1.0
		body.can_sleep = false
		body.position = Vector3(0.0, 0.5 + i, 0.0)
		root.add_child(body)
		probe.append(body)
	await _steps(240)
	var drift := 0.0
	for i in probe.size():
		var p: Vector3 = probe[i].position
		drift = maxf(drift, absf(p.y - (0.5 + i)))
		drift = maxf(drift, Vector2(p.x, p.z).length())
	for body in probe:
		body.queue_free()
	ground.queue_free()
	await process_frame
	return drift


# --- driver -----------------------------------------------------------------

func _initialize() -> void:
	_parse_args()
	var engine := String(ProjectSettings.get_setting("physics/3d/physics_engine"))
	var pinned := _pin_iterations(int(_opts.iterations))

	print("pulse-bench: engine=%s pattern=%s grid=%d gap=%.2f bodies=%d iterations=%s probe_drift=%.6f" % [
		engine, _opts.pattern, _opts.grid, _opts.gap, _opts.grid * _opts.grid, pinned,
		await _engine_probe() if _opts.probe else 0.0])

	for tick_text in String(_opts.ticks).split(","):
		var ticks := int(tick_text.strip_edges())
		if ticks <= 0:
			continue
		_build_scene()
		var evidence := await _measure(ticks)
		var line := "pulse-bench:"
		for key in evidence:
			line += " %s=%s" % [key, evidence[key]]
		print(line)
		_free_scene()

	quit(0)
