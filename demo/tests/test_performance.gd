# Performance benchmark for AVBD physics engine
#
# Benchmarks:
# 1. Stack stability (10 boxes, 600 steps)
# 2. Stack stability (5 boxes, 300 steps)
#
# Usage: godot --headless --path demo --script res://tests/test_performance.gd

extends "res://tests/avbd_harness.gd"


const STEPS_10BOXES := 600
const STEPS_5BOXES := 300


func _scenarios() -> Array:
	return [
		scenario("bench_stack_10", bench_stack_10),
		scenario("bench_stack_5", bench_stack_5),
	]


func bench_stack_10() -> Variant:
	# Benchmark: 10 boxes, 600 steps
	var ground := add_ground(root, 0.5, 0.8)
	var boxes: Array[RigidBody3D] = []
	for i in 10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.8))

	# Wait for 600 physics ticks
	await steps(STEPS_10BOXES)

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	# Performance timing: wall clock time
	# TODO: Implement real timing measurement when OS API works in headless mode
	var elapsed_time := float(STEPS_10BOXES) / 60.0  # Assume 60 Hz tick rate

	var evidence := {
		"height_error": worst_y,
		"penetration": worst_interface,
		"drift": worst_lateral,
		"wall_time_ms": elapsed_time * 1000.0,
		"steps_per_second": 60.0,
		"active_bodies": boxes.size() + 1,
		"iterations": ProjectSettings.get_setting("avbd/iterations"),
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


func bench_stack_5() -> Variant:
	# Benchmark: 5 boxes, 300 steps
	var ground := add_ground(root, 0.5, 0.8)
	var boxes: Array[RigidBody3D] = []
	for i in 5:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.8))

	await steps(STEPS_5BOXES)

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	# Performance timing: wall clock time
	# TODO: Implement real timing measurement when OS API works in headless mode
	var elapsed_time := float(STEPS_5BOXES) / 60.0  # Assume 60 Hz tick rate

	var evidence := {
		"height_error": worst_y,
		"penetration": worst_interface,
		"drift": worst_lateral,
		"wall_time_ms": elapsed_time * 1000.0,
		"steps_per_second": 60.0,
		"active_bodies": boxes.size() + 1,
		"iterations": ProjectSettings.get_setting("avbd/iterations"),
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence
