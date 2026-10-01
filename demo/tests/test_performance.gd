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
	var start_time := OS.get_ticks_msec()
	await steps(STEPS_10BOXES)
	var end_time := OS.get_ticks_msec()
	var elapsed_time_ms := end_time - start_time
	var elapsed_time := elapsed_time_ms / 1000.0

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	var evidence := {
		"height_error": worst_y,
		"penetration": worst_interface,
		"drift": worst_lateral,
		"wall_time_ms": elapsed_time_ms,
		"wall_time_s": elapsed_time,
		"steps_per_second": STEPS_10BOXES / elapsed_time,
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

	var start_time := OS.get_ticks_msec()
	await steps(STEPS_5BOXES)
	var end_time := OS.get_ticks_msec()
	var elapsed_time_ms := end_time - start_time
	var elapsed_time := elapsed_time_ms / 1000.0

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	var evidence := {
		"height_error": worst_y,
		"penetration": worst_interface,
		"drift": worst_lateral,
		"wall_time_ms": elapsed_time_ms,
		"wall_time_s": elapsed_time,
		"steps_per_second": STEPS_5BOXES / elapsed_time,
		"active_bodies": boxes.size() + 1,
		"iterations": ProjectSettings.get_setting("avbd/iterations"),
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence
