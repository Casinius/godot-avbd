# Simple stack test for diagnosis
extends "res://tests/avbd_harness.gd"


func _scenarios() -> Array:
	return [
		scenario("stack_simple_5", test_stack_5),
		scenario("stack_simple_3", test_stack_3),
	]


func test_stack_5() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)

	var boxes: Array[RigidBody3D] = []
	for i in 5:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0)))

	print("=== Stack 5 layers, 600 steps ===")
	await steps(600)

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		print("Box%d: pos=(%.3f, %.3f, %.3f)" % [i, pos.x, pos.y, pos.z])
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	print("Results: height_error=%.5f, penetration=%.5f, drift=%.5f" % [worst_y, worst_interface, worst_lateral])
	print("Check: height<=0.2? %s, penetration<=0.05? %s, drift<=0.05? %s" % [
		worst_y <= 0.2, worst_interface <= 0.05, worst_lateral <= 0.05])

	check(worst_y <= 0.2, "stack height", "max |y - expected| = %.5f" % worst_y)
	check(worst_interface <= 0.05, "stack penetration", "worst interface = %.5f m" % worst_interface)
	check(worst_lateral <= 0.05, "stack drift", "max lateral = %.5f m" % worst_lateral)

	var evidence := {"height_error": worst_y, "penetration": worst_interface, "drift": worst_lateral}
	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


func test_stack_3() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)

	var boxes: Array[RigidBody3D] = []
	for i in 3:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0)))

	print("=== Stack 3 layers, 300 steps ===")
	await steps(300)

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		print("Box%d: pos=(%.3f, %.3f, %.3f)" % [i, pos.x, pos.y, pos.z])
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	print("Results: height_error=%.5f, penetration=%.5f, drift=%.5f" % [worst_y, worst_interface, worst_lateral])
	print("Check: height<=0.2? %s, penetration<=0.05? %s, drift<=0.05? %s" % [
		worst_y <= 0.2, worst_interface <= 0.05, worst_lateral <= 0.05])

	check(worst_y <= 0.2, "stack height", "max |y - expected| = %.5f" % worst_y)
	check(worst_interface <= 0.05, "stack penetration", "worst interface = %.5f m" % worst_interface)
	check(worst_lateral <= 0.05, "stack drift", "max lateral = %.5f m" % worst_lateral)

	var evidence := {"height_error": worst_y, "penetration": worst_interface, "drift": worst_lateral}
	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence
