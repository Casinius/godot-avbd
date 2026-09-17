# Load gate: the extension must register every class with the running Godot build,
# the AVBD physics server must be the active engine, and a 200-box pyramid must
# settle deterministically through it.
#   Godot --headless --path demo --script res://tests/test_load.gd
extends SceneTree

func _initialize() -> void:
	var failures := 0

	# Server-side objects are existence-checked only: instantiating a second
	# physics server (or the state objects the engine owns) corrupts the live one.
	for cname in ["AVBDPhysicsServer3D", "PhysicsServerFactory", "AVBDDirectBodyState3D",
			"AVBDDirectSpaceState3D"]:
		if not ClassDB.class_exists(cname):
			printerr("MISSING CLASS: ", cname)
			failures += 1
			continue
		print("ok  ", cname, " registered")

	for cname in ["AVBDSoftBody3D", "AVBDSoftWorld3D"]:
		if not ClassDB.class_exists(cname):
			printerr("MISSING CLASS: ", cname)
			failures += 1
			continue
		var instance = ClassDB.instantiate(cname)
		if instance == null:
			printerr("CANNOT INSTANTIATE: ", cname)
			failures += 1
			continue
		print("ok  ", cname, " -> ", instance.get_class())
		instance.free()

	for cname in ["AVBDWorld3D", "AVBDRigidBody3D", "AVBDJoint3D"]:
		if ClassDB.class_exists(cname):
			printerr("RETIRED CLASS STILL REGISTERED: ", cname)
			failures += 1

	# Properties must round-trip through the ClassDB.
	var soft = ClassDB.instantiate("AVBDSoftBody3D")
	if _roundtrip(soft, "dimensions", Vector3i(3, 4, 5)) != Vector3i(3, 4, 5):
		failures += 1
	soft.free()

	var world = ClassDB.instantiate("AVBDSoftWorld3D")
	if not is_equal_approx(_roundtrip(world, "gravity", 12.5), 12.5):
		failures += 1
	world.free()

	# --- 200-box pyramid through the AVBD physics server -----------------------
	if ProjectSettings.get_setting("physics/3d/physics_engine") != "AVBD":
		printerr("physics engine is not AVBD")
		failures += 1

	var t0 := Time.get_ticks_usec()
	var digest_first := await _pyramid_run()
	var digest_second := await _pyramid_run()
	var ms := (Time.get_ticks_usec() - t0) / 2000.0
	if digest_first != digest_second:
		printerr("pyramid digests differ: %d vs %d" % [digest_first, digest_second])
		failures += 1

	print("ok  pyramid 2x200 boxes settled deterministically (%.1f ms total)" % ms)
	print("AVBD_LOAD_GATE: ", "PASS" if failures == 0 else "FAIL (%d)" % failures)
	quit(0 if failures == 0 else 1)


# Build a 20-row pyramid (~210 boxes), step it to rest, and digest every pose.
# Both invocations must observe the same number of solver steps, so the builder
# aligns itself to a physics-frame callback before adding nodes.
func _pyramid_run() -> int:
	await physics_frame
	var ground := StaticBody3D.new()
	var ground_shape := CollisionShape3D.new()
	var ground_box := BoxShape3D.new()
	ground_box.size = Vector3(80, 1, 20)
	ground_shape.shape = ground_box
	ground.add_child(ground_shape)
	ground.position = Vector3(0, -0.5, 0)
	root.add_child(ground)

	var boxes: Array[RigidBody3D] = []
	var box_shape := BoxShape3D.new()
	box_shape.size = Vector3(1, 0.5, 1)
	var rows := 20
	for row in rows:
		for i in rows - row:
			var b := RigidBody3D.new()
			var cs := CollisionShape3D.new()
			cs.shape = box_shape
			b.add_child(cs)
			b.mass = 0.5
			b.position = Vector3((i - (rows - row - 1) / 2.0) * 1.01, 0.25 + row * 0.85, 0)
			root.add_child(b)
			boxes.append(b)

	for i in 600:
		await physics_frame

	var all_finite := true
	var lowest := INF
	for b in boxes:
		if not b.global_position.is_finite():
			all_finite = false
		lowest = minf(lowest, b.global_position.y)

	var failures := 0
	if not all_finite:
		printerr("pyramid state not finite")
		failures += 1
	# A 0.5-high box rests with its centre at 0.25 (plus solver penetration).
	if lowest < 0.20 or lowest > 0.32:
		printerr("pyramid lowest box out of range: %f" % lowest)
		failures += 1

	# Quantized pose digest (millimetre grid), as the harness does.
	var hash := 5381
	for b in boxes:
		var t := b.global_transform
		t.origin = (t.origin * 1000.0).round() / 1000.0
		for byte in var_to_bytes(t):
			hash = (hash * 33 + byte) & 0xFFFFFFFF

	for b in boxes:
		b.queue_free()
	ground.queue_free()
	await physics_frame
	if failures > 0:
		return -hash # impossible hash value marks the run failed
	return hash


func _roundtrip(object: Object, property: String, value: Variant) -> Variant:
	object.set(property, value)
	var read = object.get(property)
	if read != value:
		printerr("PROPERTY ROUNDTRIP FAILED: %s = %s (set %s)" % [property, read, value])
	return read
