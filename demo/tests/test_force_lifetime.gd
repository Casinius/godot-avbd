# Regression tests for the lifetime of direct-state force APIs.
#
#   Godot --headless --path demo --script res://tests/test_force_lifetime.gd
extends "res://tests/avbd_harness.gd"


class ForceBody extends RigidBody3D:
	enum ForceMode { TRANSIENT_CENTRAL, TRANSIENT_OFFSET, ONE_SHOT, TORQUE, CONSTANT }

	var mode: ForceMode = ForceMode.TRANSIENT_CENTRAL
	var force := Vector3.ZERO
	var offset := Vector3.ZERO
	var torque := Vector3.ZERO
	var fired := false

	func _integrate_forces(state: PhysicsDirectBodyState3D) -> void:
		match mode:
			ForceMode.TRANSIENT_CENTRAL:
				state.apply_central_force(force)
			ForceMode.TRANSIENT_OFFSET:
				state.apply_force(force, offset)
			ForceMode.ONE_SHOT:
				if not fired:
					state.apply_central_force(force)
					fired = true
			ForceMode.TORQUE:
				state.apply_torque(torque)
			ForceMode.CONSTANT:
				pass


func _scenarios() -> Array:
	return [
		scenario("transient_force_lifetime", test_transient_force_lifetime),
		scenario("one_shot_force_lifetime", test_one_shot_force_lifetime),
		scenario("torque_lifetime", test_torque_lifetime),
		scenario("constant_force_lifetime", test_constant_force_lifetime),
	]


func _make_body(body_name: String, body_mode: ForceBody.ForceMode) -> ForceBody:
	var body := ForceBody.new()
	body.name = body_name
	body.mode = body_mode
	body.mass = 1.0
	body.gravity_scale = 0.0
	body.linear_damp = 0.0
	body.angular_damp = 0.0
	body.can_sleep = false
	body.collision_layer = 1
	body.collision_mask = 0
	body.position = Vector3.ZERO
	body.rotation = Vector3.ZERO

	var shape_node := CollisionShape3D.new()
	shape_node.name = "BoxShape"
	var shape := BoxShape3D.new()
	shape.size = Vector3.ONE
	shape_node.shape = shape
	body.add_child(shape_node)
	root.add_child(body)
	return body


func _velocity_samples(body: RigidBody3D, count: int) -> Array:
	var samples := []
	for i in count:
		await steps(1)
		samples.append(body.linear_velocity.x)
	return samples


func _angular_samples(body: RigidBody3D, count: int) -> Array:
	var samples := []
	for i in count:
		await steps(1)
		samples.append(body.angular_velocity.y)
	return samples


func _successive_deltas(samples: Array) -> Array:
	var deltas := []
	for i in range(1, samples.size()):
		deltas.append(samples[i] - samples[i - 1])
	return deltas


func test_transient_force_lifetime() -> Variant:
	check_engine()
	var central := _make_body("CentralForce", ForceBody.ForceMode.TRANSIENT_CENTRAL)
	central.force = Vector3(6.0, 0.0, 0.0)
	var offset := _make_body("OffsetForce", ForceBody.ForceMode.TRANSIENT_OFFSET)
	offset.force = Vector3(6.0, 0.0, 0.0)
	offset.offset = Vector3(0.0, 0.0, 0.25)

	var central_samples := await _velocity_samples(central, 6)
	var offset_samples := await _velocity_samples(offset, 6)
	var central_deltas := _successive_deltas(central_samples)
	var offset_deltas := _successive_deltas(offset_samples)
	var central_error := 0.0
	var offset_error := 0.0
	for delta in central_deltas:
		central_error = maxf(central_error, absf(delta - central_deltas[0]))
	for delta in offset_deltas:
		offset_error = maxf(offset_error, absf(delta - offset_deltas[0]))

	check(central_deltas[0] > 0.0, "central force accelerates", "dv=%.6f" % central_deltas[0])
	check(central_error <= 0.00001, "central force is transient", "max dv spread=%.8f" % central_error)
	check(offset_deltas[0] > 0.0, "offset force accelerates", "dv=%.6f" % offset_deltas[0])
	check(absf(offset.angular_velocity.y) > 0.0, "offset force produces torque", "wy=%.8f" % offset.angular_velocity.y)
	check(offset_error <= 0.00001, "offset force is transient", "max dv spread=%.8f" % offset_error)

	var evidence := {"central_dv": central_deltas[0], "offset_dv": offset_deltas[0]}
	central.queue_free()
	offset.queue_free()
	await frames(2)
	return evidence


func test_one_shot_force_lifetime() -> Variant:
	check_engine()
	var body := _make_body("OneShot", ForceBody.ForceMode.ONE_SHOT)
	body.force = Vector3(6.0, 0.0, 0.0)
	var samples := await _velocity_samples(body, 6)
	var deltas := _successive_deltas(samples)
	var late_change := 0.0
	for i in range(1, deltas.size()):
		late_change = maxf(late_change, absf(deltas[i]))

	check(deltas[0] > 0.0, "one-shot force applies once", "first dv=%.6f" % deltas[0])
	check(late_change <= 0.00001, "one-shot velocity stays constant", "late max dv=%.8f" % late_change)
	var evidence := {"first_dv": deltas[0], "late_change": late_change}
	body.queue_free()
	await frames(2)
	return evidence


func test_torque_lifetime() -> Variant:
	check_engine()
	var body := _make_body("Torque", ForceBody.ForceMode.TORQUE)
	body.torque = Vector3(0.0, 0.02, 0.0)
	var samples := await _angular_samples(body, 6)
	var deltas := _successive_deltas(samples)
	var spread := 0.0
	for delta in deltas:
		spread = maxf(spread, absf(delta - deltas[0]))

	check(deltas[0] > 0.0, "torque changes angular velocity", "dw=%.8f" % deltas[0])
	check(spread <= 0.00001, "torque is transient", "max dw spread=%.8f" % spread)
	check(body.rotation.is_equal_approx(Vector3.ZERO) == false or body.angular_velocity.y > 0.0,
			"torque starts without initial orientation", "rotation=%v" % body.rotation)
	var evidence := {"angular_dv": deltas[0], "angular_speed": body.angular_velocity.y}
	body.queue_free()
	await frames(2)
	return evidence


func test_constant_force_lifetime() -> Variant:
	check_engine()
	var body := _make_body("ConstantForce", ForceBody.ForceMode.CONSTANT)
	body.force = Vector3(6.0, 0.0, 0.0)
	body.constant_force = body.force
	var samples := await _velocity_samples(body, 6)
	var deltas := _successive_deltas(samples)
	var spread := 0.0
	for delta in deltas:
		spread = maxf(spread, absf(delta - deltas[0]))

	check(deltas[0] > 0.0, "constant force applies", "first dv=%.6f" % deltas[0])
	check(spread <= 0.00001, "constant force persists", "max dv spread=%.8f" % spread)
	var evidence := {"dv": deltas[0], "spread": spread}
	body.queue_free()
	await frames(2)
	return evidence
