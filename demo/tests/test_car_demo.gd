# Automated single-chassis raycast-car maneuver and wall-clock benchmark.
#   godot --headless --path demo --script res://tests/test_car_demo.gd -- --only=car_sensitivity_and_speed --strict
extends "res://tests/avbd_harness.gd"

const CAR_SCRIPT := preload("res://scenes/car_demo.gd")
const CHASSIS_SIZE := Vector3(1.8, 0.6, 3.6)
const CHASSIS_MASS := 900.0
const FRICTION := 0.9
const STEP_COUNT := 240


func _scenarios() -> Array:
	return [scenario("car_sensitivity_and_speed", car_sensitivity_and_speed)]


func _make_car() -> RigidBody3D:
	var car := RigidBody3D.new()
	car.name = "Chassis"
	car.set_script(CAR_SCRIPT)
	car.mass = CHASSIS_MASS
	car.physics_material_override = _material(FRICTION)
	car.position = Vector3(0, 1.05, 0)

	var shape_node := CollisionShape3D.new()
	shape_node.name = "Shape"
	var shape := BoxShape3D.new()
	shape.size = CHASSIS_SIZE
	shape_node.shape = shape
	car.add_child(shape_node)

	var offsets := [
		Vector3(-0.72, -0.15, -1.2), Vector3(0.72, -0.15, -1.2),
		Vector3(-0.72, -0.15, 1.2), Vector3(0.72, -0.15, 1.2)]
	for i in offsets.size():
		var ray := RayCast3D.new()
		ray.name = "WheelRay%d" % i
		ray.position = offsets[i]
		ray.target_position = Vector3(0, -0.85, 0)
		ray.collision_mask = 1
		ray.exclude_parent = true
		car.add_child(ray)
	root.add_child(car)
	return car


func _material(friction: float) -> PhysicsMaterial:
	var material := PhysicsMaterial.new()
	material.friction = friction
	return material


func _run_maneuver(steering: float, throttle: float, measure_time: bool) -> Dictionary:
	var ground := add_ground(root, 0.0, FRICTION, 50.0)
	var car := _make_car()
	# The demo controller's integration callback is the force application path used here.
	car.control_override = true
	car.override_throttle = 0.0
	car.override_brake = 0.0
	car.override_steering = 0.0
	var start_position := car.position
	var start_time := Time.get_ticks_usec()
	await steps(60)
	car.override_throttle = throttle
	await steps(120)
	car.override_throttle = 0.0
	car.override_steering = steering
	await steps(60)
	var end_time := Time.get_ticks_usec()

	var forward := car.global_transform.basis * Vector3(0, 0, -1)
	var yaw := atan2(forward.x, -forward.z)
	var displacement := car.position - start_position
	var elapsed_usec := maxi(1, end_time - start_time)
	var result := {
		"forward_displacement": -displacement.z,
		"lateral_displacement": displacement.x,
		"speed": car.linear_velocity.length(),
		"yaw_change": yaw,
		"finite": is_finite_body(car),
		"elapsed_usec": elapsed_usec,
	}
	print("car run: steer=%.2f drive=%.2f forward=%.4f m lateral=%.4f m speed=%.4f m/s yaw=%.6f rad finite=%s" % [
		steering, throttle, result.forward_displacement, result.lateral_displacement,
		result.speed, result.yaw_change, result.finite])
	if measure_time:
		var elapsed_seconds := elapsed_usec / 1000000.0
		print("car benchmark: elapsed=%.3f ms steps/s=%.3f ms/physics-step=%.6f" % [
			elapsed_usec / 1000.0, STEP_COUNT / elapsed_seconds, elapsed_usec / 1000.0 / STEP_COUNT])
	for node in [car, ground]:
		node.queue_free()
	await frames(2)
	return result


func car_sensitivity_and_speed() -> Variant:
	check_engine()
	var baseline := await _run_maneuver(1.0, 1.0, true)
	var control := await _run_maneuver(0.0, 0.0, false)
	check(baseline.finite and control.finite, "car states finite", "baseline=%s control=%s" % [baseline.finite, control.finite])
	check(baseline.forward_displacement > control.forward_displacement + 1.0,
			"driven forward response", "baseline=%.4f m control=%.4f m" % [
			baseline.forward_displacement, control.forward_displacement])
	check(absf(baseline.yaw_change) > absf(control.yaw_change) + 0.05,
			"driven steering response", "baseline=%.6f rad control=%.6f rad" % [
			baseline.yaw_change, control.yaw_change])

	var half_steer := await _run_maneuver(0.5, 1.0, false)
	var full_steer := await _run_maneuver(1.0, 1.0, false)
	var yaw_delta := absf(full_steer.yaw_change - half_steer.yaw_change)
	print("steering sensitivity: steer=0.50 yaw=%.6f, steer=1.00 yaw=%.6f, delta=%.6f rad" % [
		half_steer.yaw_change, full_steer.yaw_change, yaw_delta])
	check(yaw_delta >= 0.02, "steering sensitivity", "yaw delta=%.6f rad" % yaw_delta)
	return {
		"baseline": baseline,
		"control": control,
		"half_steer": half_steer,
		"full_steer": full_steer,
		"yaw_delta": yaw_delta,
	}
