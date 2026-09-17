# Static API-compliance audit for the AVBD GDExtension (server-first).
#
#   Godot --headless --path demo --script res://tests/test_server.gd
#
# Runs entirely outside the scene tree: checks that the physics server path is
# selected, that the registered classes are visible to ClassDB, that the retired
# node classes are gone, and that soft-body properties round-trip.
extends SceneTree

var checks := 0
var failures := 0


func check(ok: bool, name: String, detail := "") -> void:
	checks += 1
	if not ok:
		failures += 1
	print("%s %-46s %s" % ["  ok " if ok else "FAIL ", name, detail])


func _initialize() -> void:
	print("== project")
	check(ProjectSettings.get_setting("physics/3d/physics_engine") == "AVBD",
			"physics/3d/physics_engine == AVBD")

	print("== registered classes")
	for cls in ["AVBDPhysicsServer3D", "PhysicsServerFactory", "AVBDDirectBodyState3D",
			"AVBDDirectSpaceState3D", "AVBDSoftBody3D", "AVBDSoftWorld3D"]:
		check(ClassDB.class_exists(cls), "class %s registered" % cls)

	print("== retired node classes")
	for cls in ["AVBDWorld3D", "AVBDRigidBody3D", "AVBDJoint3D", "AVBDGeneric6DOFJoint3D",
			"AVBDSpring3D", "AVBDIgnoreCollision3D", "AVBDConstraint3D"]:
		check(not ClassDB.class_exists(cls), "class %s removed" % cls)

	print("== soft body round-trips")
	var soft := AVBDSoftBody3D.new()
	for trip in [[Vector3i(3, 4, 5), Vector3i(3, 4, 5)], [Vector3(0.3, 0.3, 0.3), Vector3(0.3, 0.3, 0.3)],
			[2.5, 2.5], [0.75, 0.75], [1200.0, 1200.0], [350.0, 350.0]]:
		var keys := ["dimensions", "box_size", "gap", "density", "friction", "stiffness_linear"]
	# spot-check a representative set
	soft.dimensions = Vector3i(3, 4, 5)
	check(soft.dimensions == Vector3i(3, 4, 5), "soft.dimensions round-trips", "%s" % soft.dimensions)
	soft.dimensions = Vector3i(4, 4, 4)
	soft.box_size = Vector3(0.3, 0.3, 0.3)
	check(soft.box_size.is_equal_approx(Vector3(0.3, 0.3, 0.3)), "soft.box_size round-trips", "")
	soft.density = 2.5
	check(is_equal_approx(soft.density, 2.5), "soft.density round-trips", "")
	soft.friction = 0.75
	check(is_equal_approx(soft.friction, 0.75), "soft.friction round-trips", "")
	soft.stiffness_linear = 1200.0
	check(is_equal_approx(soft.stiffness_linear, 1200.0), "soft.stiffness_linear round-trips", "")
	check(soft.get_cell_count() == 64, "soft 4x4x4 cell count", "cells=%d" % soft.get_cell_count())
	soft.free()

	var world := AVBDSoftWorld3D.new()
	world.gravity = 12.5
	check(is_equal_approx(world.gravity, 12.5), "soft_world.gravity round-trips", "")
	world.threads = 3
	check(world.threads == 3, "soft_world.threads round-trips", "")
	world.free()

	print("== server handles")
	var space: RID = PhysicsServer3D.space_create()
	check(space.is_valid(), "space_create valid", "")
	var body: RID = PhysicsServer3D.body_create()
	check(body.is_valid(), "body_create valid", "")
	PhysicsServer3D.body_set_space(body, space)
	check(PhysicsServer3D.body_get_space(body) == space, "body_get_space round-trips", "")
	var shape: RID = PhysicsServer3D.box_shape_create()
	PhysicsServer3D.shape_set_data(shape, Vector3(0.5, 0.5, 0.5)) # half extents
	var shape_type: int = PhysicsServer3D.shape_get_type(shape)
	check(shape_type == PhysicsServer3D.SHAPE_BOX, "box shape type", "%d" % shape_type)
	PhysicsServer3D.body_add_shape(body, shape, Transform3D(), false)
	check(PhysicsServer3D.body_get_shape_count(body) == 1, "body_get_shape_count", "")
	PhysicsServer3D.body_set_collision_layer(body, 0x20)
	check(PhysicsServer3D.body_get_collision_layer(body) == 0x20, "collision layer round-trips", "")
	PhysicsServer3D.body_set_collision_mask(body, 0x40)
	check(PhysicsServer3D.body_get_collision_mask(body) == 0x40, "collision mask round-trips", "")
	PhysicsServer3D.body_set_param(body, PhysicsServer3D.BODY_PARAM_GRAVITY_SCALE, 0.5)
	check(is_equal_approx(PhysicsServer3D.body_get_param(body, PhysicsServer3D.BODY_PARAM_GRAVITY_SCALE), 0.5),
			"gravity scale round-trips", "")
	var other: RID = PhysicsServer3D.body_create()
	PhysicsServer3D.body_add_collision_exception(body, other)
	PhysicsServer3D.free_rid(body)
	PhysicsServer3D.free_rid(other)
	PhysicsServer3D.free_rid(shape)
	PhysicsServer3D.free_rid(space)

	print("")
	print("%d checks, %d failures" % [checks, failures])
	quit(1 if failures > 0 else 0)
