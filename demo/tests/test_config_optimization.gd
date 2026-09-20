# Configuration optimization test for AVBD
# Tests different solver configurations to find best balance

extends "res://tests/avbd_harness.gd"


const TEST_STEPS := 300
const BOXES_10 := 10
const BOXES_5 := 5


func _scenarios() -> Array:
	return [
		scenario("config_default", test_config_default),
		scenario("config_high_friction", test_config_high_friction),
		scenario("config_high_iterations", test_config_high_iterations),
		scenario("config_both", test_config_both),
	]


func test_config_default() -> Variant:
	# Default: friction 0.5, iterations 10
	var ground := add_ground(root, 0.5, 0.5, 50.0)
	var boxes: Array[RigidBody3D] = []
	for i in BOXES_10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.5))

	await steps(TEST_STEPS)

	var worst_drift := 0.0
	for i in boxes.size():
		var drift := Vector2(boxes[i].global_position.x, boxes[i].global_position.z).length()
		worst_drift = maxf(worst_drift, drift)

	var ground_pos := ground.global_position
	var evidence := {
		"config": "default (friction=0.5, iterations=10)",
		"worst_drift": worst_drift,
		"max_y_error": absf(boxes[boxes.size() - 1].global_position.y - (1.0 + BOXES_10 - 1)),
		"iterations": ProjectSettings.get_setting("avbd/iterations"),
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


func test_config_high_friction() -> Variant:
	# High friction: friction 0.9, iterations 10
	var ground := add_ground(root, 0.5, 0.9, 50.0)
	var boxes: Array[RigidBody3D] = []
	for i in BOXES_10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.9))

	await steps(TEST_STEPS)

	var worst_drift := 0.0
	for i in boxes.size():
		var drift := Vector2(boxes[i].global_position.x, boxes[i].global_position.z).length()
		worst_drift = maxf(worst_drift, drift)

	var ground_pos := ground.global_position
	var evidence := {
		"config": "high friction (friction=0.9, iterations=10)",
		"worst_drift": worst_drift,
		"max_y_error": absf(boxes[boxes.size() - 1].global_position.y - (1.0 + BOXES_10 - 1)),
		"iterations": ProjectSettings.get_setting("avbd/iterations"),
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


func test_config_high_iterations() -> Variant:
	# High iterations: friction 0.5, iterations 20
	var ground := add_ground(root, 0.5, 0.5, 50.0)
	var boxes: Array[RigidBody3D] = []
	for i in BOXES_10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.5))

	# Temporarily increase iterations
	var old_iterations: int = ProjectSettings.get_setting("avbd/iterations", 10)
	ProjectSettings.set_setting("avbd/iterations", 20)

	await steps(TEST_STEPS)

	var worst_drift := 0.0
	for i in boxes.size():
		var drift := Vector2(boxes[i].global_position.x, boxes[i].global_position.z).length()
		worst_drift = maxf(worst_drift, drift)

	# Restore
	ProjectSettings.set_setting("avbd/iterations", old_iterations)

	var ground_pos := ground.global_position
	var evidence := {
		"config": "high iterations (friction=0.5, iterations=20)",
		"worst_drift": worst_drift,
		"max_y_error": absf(boxes[boxes.size() - 1].global_position.y - (1.0 + BOXES_10 - 1)),
		"iterations": 20,
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


func test_config_both() -> Variant:
	# High both: friction 0.9, iterations 20
	var ground := add_ground(root, 0.5, 0.9, 50.0)
	var boxes: Array[RigidBody3D] = []
	for i in BOXES_10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.9))

	# Temporarily increase iterations
	var old_iterations: int = ProjectSettings.get_setting("avbd/iterations", 10)
	ProjectSettings.set_setting("avbd/iterations", 20)

	await steps(TEST_STEPS)

	var worst_drift := 0.0
	for i in boxes.size():
		var drift := Vector2(boxes[i].global_position.x, boxes[i].global_position.z).length()
		worst_drift = maxf(worst_drift, drift)

	# Restore
	ProjectSettings.set_setting("avbd/iterations", old_iterations)

	var ground_pos := ground.global_position
	var evidence := {
		"config": "high both (friction=0.9, iterations=20)",
		"worst_drift": worst_drift,
		"max_y_error": absf(boxes[boxes.size() - 1].global_position.y - (1.0 + BOXES_10 - 1)),
		"iterations": 20,
		"beta_linear": ProjectSettings.get_setting("avbd/beta_linear"),
		"beta_angular": ProjectSettings.get_setting("avbd/beta_angular"),
	}

	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence
