# Double-pyramid stress regression.
#
# pyramid_dbl.tscn stacks two 20-row, 210-box pyramids. The demo driver drops
# both with can_sleep=false, so this suite does the same: the pile must settle
# through the solver, not by freezing bodies. It runs long enough for the upper
# pyramid to fall and then measures a window for runaway drift / launched boxes.
#
#   Godot --headless --path demo --script res://tests/test_pyramid.gd
#   ... -- --list | --only=pyramid_dbl | --digest | --strict
extends "res://tests/avbd_harness.gd"

const SCENE := "res://scenes/pyramid_dbl.tscn"
const BODY_COUNT := 420
const SETTLE_STEPS := 1200 # 20 s at the logical 60 Hz tick
const WINDOW_STEPS := 600  # 10 s observation window

func _scenarios() -> Array:
	return [scenario("pyramid_dbl", test_pyramid_dbl)]

func _collect_rigid_bodies(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect_rigid_bodies(child, out)

func test_pyramid_dbl() -> Variant:
	var packed: PackedScene = load(SCENE)
	var scene := packed.instantiate()
	# The demo driver adds a random perturbation and prints stats; a regression
	# test needs the authored poses and no per-frame noise.
	scene.set_script(null)
	root.add_child(scene)

	var bodies: Array[RigidBody3D] = []
	_collect_rigid_bodies(scene, bodies)
	check(bodies.size() == BODY_COUNT, "body count", "found %d (expected %d)" % [bodies.size(), BODY_COUNT])
	# Match the demo: sleep is disabled, so a stable pile has to come from the solver.
	for body in bodies:
		body.mass = 1.0
		body.sleeping = false
		body.can_sleep = false

	await steps(SETTLE_STEPS)

	var p0: Array[Vector3] = []
	var com0 := Vector3.ZERO
	var max_y0 := -1.0e30
	for body in bodies:
		p0.append(body.global_position)
		com0 += body.global_position
		max_y0 = maxf(max_y0, body.global_position.y)
	com0 /= float(bodies.size())

	await steps(WINDOW_STEPS)

	var max_delta := 0.0
	var max_v := 0.0
	var max_ang := 0.0
	var max_xz := 0.0
	var max_y := -1.0e30
	var min_y := 1.0e30
	var com1 := Vector3.ZERO
	var finite := true
	for i in bodies.size():
		var body := bodies[i]
		var pos := body.global_position
		var d := pos - p0[i]
		max_delta = maxf(max_delta, d.length())
		max_v = maxf(max_v, body.linear_velocity.length())
		max_ang = maxf(max_ang, body.angular_velocity.length())
		max_xz = maxf(max_xz, Vector2(pos.x, pos.z).length())
		max_y = maxf(max_y, pos.y)
		min_y = minf(min_y, pos.y)
		com1 += pos
		finite = finite and is_finite_body(body)
	com1 /= float(bodies.size())
	var com_drift := (com1 - com0).length()

	check(finite, "pile stays finite", "bodies=%d" % bodies.size())
	# A body that is still climbing means the sleeping-mode/contact feedback bug
	# is back: the old runaway took boxes to y > 40 m with zero reported velocity.
	check(max_y <= 15.0, "no launched boxes", "max y=%.3f m (settled max was %.3f)" % [max_y, max_y0])
	check(min_y >= -2.0, "no boxes below the ground", "min y=%.3f m" % min_y)
	check(max_xz <= 40.0, "pile stays on the slab", "max |xz|=%.3f m" % max_xz)
	# A pile may still rearrange slowly, but 8 m of per-body motion in 10 s or
	# 0.5 m of centre-of-mass drift is a runaway, not settling.
	check(max_delta <= 8.0, "no long-range drift", "max body move=%.3f m in %.0f s" % [max_delta, WINDOW_STEPS / 60.0])
	check(com_drift <= 0.5, "centre of mass stays put", "COM drift=%.3f m in %.0f s" % [com_drift, WINDOW_STEPS / 60.0])
	check(max_v <= 2.5, "pile speed bounded", "max linear speed=%.3f m/s" % max_v)
	check(max_ang <= 8.0, "pile spin bounded", "max angular speed=%.3f rad/s" % max_ang)

	var evidence := {
		"max_delta": max_delta, "max_v": max_v, "max_ang": max_ang,
		"max_xz": max_xz, "max_y": max_y, "min_y": min_y, "com_drift": com_drift,
	}
	scene.queue_free()
	await frames(2)
	return evidence
