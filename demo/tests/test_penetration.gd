# Long-run penetration detection for the double-pyramid stress scene.
#
# pyramid_dbl stacks two 20-row pyramids: the lower one settles on the slab and the upper
# one starts ~1.4 m above it, so the pile goes through a *secondary* collision
# (upper-on-lower) after the first ground contacts. Deep shaking/tunnelling would show up
# there first, so the scenario watches penetration for a long wall of simulated time.
#
# Penetration is read from the solver's own contact error through the diagnostic added to
# the AVBD server (PhysicsServer3D.call("get_solver_max_penetration")): a manifold's C0 is
# the geometric overlap captured at the start of the step, so this is the real interpenetration
# and stays correct once the pile rotates (an AABB/transform heuristic would not).
#
# Tracked every physics frame:
#   - the deepest contact anywhere, and the instance-id pair when it is sampled
#   - "secondary collisions": a body that had come to rest is set moving again by a new
#     impact (speed crosses IMPACT_SPEED after REST_FRAMES frames under REST_SPEED)
#   - the deepest penetration seen while such a secondary impact is resolving
#   - sink-through: the lowest body centre against the slab top
#
#   Godot --headless --path demo --script res://tests/test_penetration.gd
#   ... -- --only=pyramid_dbl_penetration
#   ... -- --pen-seconds=180          (longer run)

extends "res://tests/avbd_harness.gd"

const SCENE := "res://scenes/pyramid_dbl.tscn"
const BODY_COUNT := 420
const RUN_SECONDS := 90.0  # logical seconds (60 Hz); override with --pen-seconds=
const REST_FRAMES := 3     # frames under REST_SPEED before a body counts as resting
const REST_SPEED := 0.25   # m/s
const IMPACT_SPEED := 0.5  # m/s that re-wakes a resting body = secondary collision
const SECONDARY_WINDOW := 12
# Bodies above this authored height belong to the dropped (upper) pyramid; the lower
# pyramid's apex is at y ~ 10.1 and the upper base is at y ~ 12.1.
const UPPER_PYRAMID_Y := 11.0

var _run_seconds := RUN_SECONDS


func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--pen-seconds="):
			_run_seconds = float(arg.split("=")[1])
	super._initialize()


func _scenarios() -> Array:
	return [scenario("pyramid_dbl_penetration", test_pyramid_dbl_penetration)]


func _collect_rigid_bodies(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect_rigid_bodies(child, out)


func _describe_instance(instance_id: int) -> String:
	if instance_id == 0:
		return "world"
	var obj := instance_from_id(instance_id)
	return str(obj.name) if obj != null else "body#%d" % instance_id


func test_pyramid_dbl_penetration() -> Variant:
	check_engine()
	check(PhysicsServer3D.has_method("get_solver_max_penetration"),
			"penetration diagnostic present", "PhysicsServer3D.get_solver_max_penetration")

	var packed: PackedScene = load(SCENE)
	var scene := packed.instantiate()
	# The demo driver adds random perturbation and prints stats; a regression run needs the
	# authored poses and no per-frame noise.
	scene.set_script(null)
	root.add_child(scene)

	var bodies: Array[RigidBody3D] = []
	_collect_rigid_bodies(scene, bodies)
	check(bodies.size() == BODY_COUNT, "body count",
			"found %d (expected %d)" % [bodies.size(), BODY_COUNT])

	var n := bodies.size()
	var initial_y := PackedFloat32Array()
	initial_y.resize(n)
	for i in n:
		var body := bodies[i]
		body.mass = 1.0
		body.sleeping = false
		body.can_sleep = false
		initial_y[i] = body.global_position.y if body.is_inside_tree() else body.position.y
	# Settle the flagged flag array sizing.
	var rest_streak := PackedInt32Array()
	rest_streak.resize(n)
	var secondary_armed := PackedInt32Array()
	secondary_armed.resize(n)
	var secondary_events := PackedInt32Array()
	secondary_events.resize(n)

	var worst_pen := 0.0
	var worst_pen_frame := 0
	var worst_pen_a := 0
	var worst_pen_b := 0
	var pen_at_secondary := 0.0
	var sec_upper := 0
	var sec_lower := 0
	var sec_upper_bodies := 0
	var sec_lower_bodies := 0
	var min_y := 1.0e30
	var max_y := -1.0e30
	var max_xz := 0.0
	var finite := true

	var ticks := int(_run_seconds * float(Engine.physics_ticks_per_second))
	var t0 := Time.get_ticks_msec()
	for frame in ticks:
		await physics_frame

		var frame_pen: float = PhysicsServer3D.call("get_solver_max_penetration")
		if frame_pen > worst_pen:
			worst_pen = frame_pen
			worst_pen_frame = frame
			var deepest: Array = PhysicsServer3D.call("get_solver_deepest_contact")
			if deepest.size() == 3:
				worst_pen_a = int(deepest[0])
				worst_pen_b = int(deepest[1])

		for i in n:
			var body := bodies[i]
			var speed := body.linear_velocity.length()
			if rest_streak[i] >= REST_FRAMES and speed > IMPACT_SPEED:
				secondary_events[i] += 1
				if secondary_events[i] == 1 and initial_y[i] > UPPER_PYRAMID_Y:
					sec_upper_bodies += 1
				elif secondary_events[i] == 1:
					sec_lower_bodies += 1
				if initial_y[i] > UPPER_PYRAMID_Y:
					sec_upper += 1
				else:
					sec_lower += 1
				secondary_armed[i] = SECONDARY_WINDOW
				rest_streak[i] = 0
			elif speed < REST_SPEED:
				rest_streak[i] = mini(rest_streak[i] + 1, REST_FRAMES)
			else:
				rest_streak[i] = 0

			if secondary_armed[i] > 0:
				secondary_armed[i] -= 1
				pen_at_secondary = maxf(pen_at_secondary, frame_pen)

			var pos := body.global_position
			min_y = minf(min_y, pos.y)
			max_y = maxf(max_y, pos.y)
			max_xz = maxf(max_xz, Vector2(pos.x, pos.z).length())
			finite = finite and is_finite_body(body)

	var elapsed := (Time.get_ticks_msec() - t0) / 1000.0
	var logical_time := float(worst_pen_frame) / float(Engine.physics_ticks_per_second)

	check(finite, "pile stays finite", "bodies=%d" % n)
	# A body whose centre drops below the slab surface by more than its own half-height has
	# gone through the ground (or fallen off the slab's edge, which the report makes visible).
	check(min_y >= -0.25, "no deep ground sink", "min centre y=%.4f m" % min_y)
	# Half a box height (0.25 m) is the hard ceiling for a surviving stack; a resting contact
	# should stay around the 1 mm collision margin. 5 cm already means the pile is being crushed.
	check(worst_pen <= 0.05, "penetration bounded",
			"worst=%.4f m at t=%.1f s (%s <-> %s)" % [
				worst_pen, logical_time, _describe_instance(worst_pen_a), _describe_instance(worst_pen_b)])

	var evidence := {
		"seconds": _run_seconds,
		"wall_seconds": elapsed,
		"worst_penetration": worst_pen,
		"worst_penetration_time": logical_time,
		"worst_pair": "%s <-> %s" % [_describe_instance(worst_pen_a), _describe_instance(worst_pen_b)],
		"penetration_at_secondary_collisions": pen_at_secondary,
		"secondary_events_upper": sec_upper,
		"secondary_events_lower": sec_lower,
		"secondary_bodies_upper": sec_upper_bodies,
		"secondary_bodies_lower": sec_lower_bodies,
		"min_y": min_y,
		"max_y": max_y,
		"max_xz": max_xz,
	}
	scene.queue_free()
	await frames(2)
	return evidence
