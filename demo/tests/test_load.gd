# Load gate: the extension must register every class with the running Godot build.
#   Godot_v4.6.3 --headless --path demo --script res://tests/test_load.gd
extends SceneTree

const EXPECTED := [
	"AVBDWorld3D",
	"AVBDRigidBody3D",
	"AVBDSoftBody3D",
	"AVBDConstraint3D",
	"AVBDJoint3D",
	"AVBDSpring3D",
	"AVBDIgnoreCollision3D",
]


func _initialize() -> void:
	var failures := 0

	for cname in EXPECTED:
		if not ClassDB.class_exists(cname):
			printerr("MISSING CLASS: ", cname)
			failures += 1
			continue
		if not ClassDB.can_instantiate(cname):
			# AVBDConstraint3D is abstract: it must exist but not be constructible.
			print("ok  ", cname, " (abstract, not instantiable)")
			continue
		var instance = ClassDB.instantiate(cname)
		if instance == null:
			printerr("CANNOT INSTANTIATE: ", cname)
			failures += 1
			continue
		print("ok  ", cname, " -> ", instance.get_class())
		instance.free()

	# Properties must round-trip through the ClassDB, including some whose values
	# are stored in the solver's Z-up space.
	var world = ClassDB.instantiate("AVBDWorld3D")
	if not is_equal_approx(_roundtrip(world, "gravity", 12.5), 12.5):
		failures += 1
	if int(_roundtrip(world, "iterations", 24)) != 24:
		failures += 1
	if not is_equal_approx(_roundtrip(world, "beta_linear", 250000.0), 250000.0):
		failures += 1
	world.free()

	var body = ClassDB.instantiate("AVBDRigidBody3D")
	if not Vector3(1, 2, 3).is_equal_approx(_roundtrip(body, "size", Vector3(1, 2, 3))):
		failures += 1
	if not is_equal_approx(_roundtrip(body, "friction", 0.75), 0.75):
		failures += 1
	if bool(_roundtrip(body, "static_body", true)) != true:
		failures += 1
	body.free()

	var joint = ClassDB.instantiate("AVBDJoint3D")
	# -1 is the "infinitely stiff" sentinel; it must survive serialisation.
	if not is_equal_approx(_roundtrip(joint, "linear_stiffness", -1.0), -1.0):
		failures += 1
	if not is_equal_approx(_roundtrip(joint, "angular_stiffness", 250.0), 250.0):
		failures += 1
	# Fracture state and the constraint/force plumbing.
	for method in ["is_broken", "set_broken", "break_joint", "is_simulated"]:
		if not joint.has_method(method):
			printerr("MISSING METHOD: AVBDJoint3D.%s" % method)
			failures += 1
	if joint.is_broken():
		printerr("a fresh joint reports itself broken")
		failures += 1
	joint.break_joint()
	if not joint.is_broken():
		printerr("break_joint() did not break the joint")
		failures += 1
	joint.set_broken(false)
	if joint.is_broken():
		printerr("set_broken(false) did not mend the joint")
		failures += 1
	if joint.is_simulated():
		printerr("an unsimulated joint reports a solver force")
		failures += 1
	joint.free()

	for clazz in ["AVBDJoint3D", "AVBDSpring3D", "AVBDIgnoreCollision3D"]:
		if not ClassDB.class_has_method(clazz, "is_simulated"):
			printerr("MISSING METHOD: %s.is_simulated" % clazz)
			failures += 1

	var world_probe = ClassDB.instantiate("AVBDWorld3D")
	if not world_probe.has_method("get_step_count"):
		printerr("MISSING METHOD: AVBDWorld3D.get_step_count")
		failures += 1
	if not world_probe.has_method("get_thread_count"):
		printerr("MISSING METHOD: AVBDWorld3D.get_thread_count")
		failures += 1
	# The thread setting must round-trip, and 0 must mean "automatic".
	if int(_roundtrip(world_probe, "threads", 3)) != 3:
		failures += 1
	world_probe.set("threads", 0)
	var auto_workers := int(world_probe.call("get_thread_count"))
	if auto_workers <= 0:
		printerr("threads=0 should report a usable worker count, got %d" % auto_workers)
		failures += 1
	print("ok  AVBDWorld3D.threads=0 -> %d workers" % auto_workers)
	world_probe.free()

	var soft = ClassDB.instantiate("AVBDSoftBody3D")
	if _roundtrip(soft, "dimensions", Vector3i(3, 4, 5)) != Vector3i(3, 4, 5):
		failures += 1
	soft.free()

	print("AVBD_LOAD_GATE: ", "PASS" if failures == 0 else "FAIL (%d)" % failures)
	quit(0 if failures == 0 else 1)


func _roundtrip(object: Object, property: String, value: Variant) -> Variant:
	object.set(property, value)
	var read = object.get(property)
	if read != value:
		printerr("PROPERTY ROUNDTRIP FAILED: %s = %s (set %s)" % [property, read, value])
	return read
