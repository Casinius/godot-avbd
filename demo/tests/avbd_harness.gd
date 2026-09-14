# Shared harness for the AVBD Godot-side test suites.
#
# A suite subclasses this script and returns its scenarios from `_scenarios()`:
#
#   extends "res://tests/avbd_harness.gd"
#   func _scenarios() -> Array: return [scenario("name", run_thing)]
#
# Each scenario's callable may `await` and returns its evidence (a Dictionary of the
# measurements it made, or null). The harness prints the evidence, counts pass/fail,
# and with `--digest` runs every scenario twice, comparing the evidence text so the
# whole scenario - not just the assertions - is checked for reproducibility.
#
#   Godot --headless --path demo --script res://tests/test_X.gd
#   Godot --headless --path demo --script res://tests/test_X.gd -- --list
#   Godot --headless --path demo --script res://tests/test_X.gd -- --only=my_scenario
#   Godot --headless --path demo --script res://tests/test_X.gd -- --digest --strict
extends SceneTree

const GRAVITY := 9.8


# --- scenario declaration ---------------------------------------------------

func _scenarios() -> Array:
	push_error("a suite must override _scenarios()")
	return []


func scenario(name: String, run: Callable) -> Dictionary:
	return {"name": name, "run": run}


# --- assertions -------------------------------------------------------------

func check(ok: bool, name: String, detail: String) -> void:
	checks += 1
	if not ok:
		failures += 1
	print("%s %-30s %s" % ["  ok " if ok else "FAIL ", name, detail])
	if not ok and strict:
		print("\naborting: --strict stops at the first failure")
		_report_and_quit()


func check_near(value: float, expected: float, tolerance: float, name: String, label := "value") -> void:
	check(absf(value - expected) <= tolerance, name,
			"%s = %.6f (want %.6f +/- %g)" % [label, value, expected, tolerance])


func check_range(value: float, low: float, high: float, name: String, label := "value") -> void:
	check(value >= low and value <= high, name, "%s = %.6f (want %g .. %g)" % [label, value, low, high])


func check_between(value: float, low: float, high: float, name: String, label := "value") -> void:
	check_range(value, low, high, name, label)


# Structural pre-flight: the world creates one solver force per active constraint, so
# a mismatch means a constraint was silently skipped (unresolvable NodePath, a body
# outside the world, ...). The world builds its solver on its first physics tick, so
# this waits for that tick before counting.
#
# Contact manifolds are also solver forces, but their number depends on the scene's
# geometry at that moment, so they are excluded: this compares declared constraints.
func check_forces(world: AVBDWorld3D, expected: int, name := "constraint pre-flight") -> void:
	await steps(world, 1)
	var constraints := world.get_force_count() - world.get_contact_count()
	check(constraints == expected, name,
			"constraint forces = %d, declared = %d (%d contact manifolds ignored)" %
			[constraints, expected, world.get_contact_count()])


# --- simulation helpers -----------------------------------------------------

# Wait for `count` physics ticks.
func frames(count: int) -> void:
	for i in count:
		await physics_frame


# Step the world for exactly `count` solver steps. Waiting on physics frames alone can
# give one world a different number of ticks than another, depending on where in the
# frame it was created, which makes comparisons meaningless.
func steps(world: AVBDWorld3D, count: int) -> void:
	var target := world.get_step_count() + count
	while world.get_step_count() < target:
		await physics_frame


func make_world(parent: Node = null) -> AVBDWorld3D:
	var world := AVBDWorld3D.new()
	if parent == null:
		parent = root
	parent.add_child(world)
	return world


func add_body(parent: Node, body_name: String, size: Vector3, position: Vector3, density := 1.0,
		friction := 0.5, is_static := false) -> AVBDRigidBody3D:
	var body := AVBDRigidBody3D.new()
	body.name = body_name
	body.size = size
	body.density = density
	body.friction = friction
	body.static_body = is_static
	parent.add_child(body)
	body.position = position
	return body


func add_ground(world: AVBDWorld3D, top := 0.0, friction := 0.5, half := 50.0) -> AVBDRigidBody3D:
	# A 1 m thick slab whose top face sits at y = `top`.
	return add_body(world, "Ground", Vector3(half * 2.0, 1.0, half * 2.0), Vector3(0, top - 0.5, 0), 0.0, friction, true)


func add_joint(world: AVBDWorld3D, joint_name: String, a: Node3D, b: Node3D, anchor_a: Vector3,
		anchor_b: Vector3, linear_stiffness := -1.0, angular_stiffness := 0.0,
		fracture_force := -1.0) -> AVBDJoint3D:
	var joint := AVBDJoint3D.new()
	joint.name = joint_name
	joint.anchor_a = anchor_a
	joint.anchor_b = anchor_b
	joint.linear_stiffness = linear_stiffness
	joint.angular_stiffness = angular_stiffness
	joint.fracture_force = fracture_force
	world.add_child(joint)
	if a != null:
		joint.node_a = joint.get_path_to(a)
	joint.node_b = joint.get_path_to(b)
	return joint


# Joint to a fixed point in world space: leave the A side empty.
func add_world_joint(world: AVBDWorld3D, joint_name: String, b: Node3D, world_anchor: Vector3,
		anchor_b := Vector3.ZERO, linear_stiffness := -1.0, angular_stiffness := 0.0) -> AVBDJoint3D:
	return add_joint(world, joint_name, null, b, world_anchor, anchor_b, linear_stiffness, angular_stiffness)


func add_spring(world: AVBDWorld3D, spring_name: String, a: Node3D, b: Node3D, anchor_a: Vector3,
		anchor_b: Vector3, stiffness := 1000.0, rest_length := -1.0) -> AVBDSpring3D:
	var spring := AVBDSpring3D.new()
	spring.name = spring_name
	spring.anchor_a = anchor_a
	spring.anchor_b = anchor_b
	spring.stiffness = stiffness
	spring.rest_length = rest_length
	world.add_child(spring)
	spring.node_a = spring.get_path_to(a)
	spring.node_b = spring.get_path_to(b)
	return spring


func add_ignore(world: AVBDWorld3D, ignore_name: String, a: Node3D, b: Node3D) -> AVBDIgnoreCollision3D:
	var ignore := AVBDIgnoreCollision3D.new()
	ignore.name = ignore_name
	world.add_child(ignore)
	ignore.node_a = ignore.get_path_to(a)
	ignore.node_b = ignore.get_path_to(b)
	return ignore


func add_soft_body(world: AVBDWorld3D, body_name: String, position: Vector3, dimensions := Vector3i(3, 3, 3),
		box_size := Vector3(0.25, 0.25, 0.25)) -> AVBDSoftBody3D:
	var soft := AVBDSoftBody3D.new()
	soft.name = body_name
	soft.dimensions = dimensions
	soft.box_size = box_size
	world.add_child(soft)
	soft.position = position
	return soft


# --- measurement helpers ----------------------------------------------------

func is_finite_body(body: AVBDRigidBody3D) -> bool:
	var p := body.global_position
	var v := body.linear_velocity
	return p.is_finite() and v.is_finite() and body.get_angular_velocity().is_finite()


func all_finite(bodies: Array) -> bool:
	for body in bodies:
		if body is AVBDRigidBody3D and not is_finite_body(body):
			return false
	return true


func max_speed(bodies: Array) -> float:
	var worst := 0.0
	for body in bodies:
		if body is AVBDRigidBody3D:
			worst = maxf(worst, body.linear_velocity.length())
	return worst


func ground_clearance(bodies: Array, floor_y: float) -> float:
	var lowest := INF
	for body in bodies:
		if body is AVBDRigidBody3D:
			lowest = minf(lowest, body.global_position.y)
	return lowest - floor_y


# Pose digest of every simulated body, for comparing two runs.
func digest(world: AVBDWorld3D) -> int:
	var hash := 5381
	for child in world.get_children():
		if child is AVBDRigidBody3D:
			for byte in var_to_bytes(child.global_transform):
				hash = (hash * 33 + byte) & 0xFFFFFFFF
		elif child is AVBDSoftBody3D:
			for i in child.get_cell_count():
				for byte in var_to_bytes(child.to_global(child.get_cell_transform(i).origin)):
					hash = (hash * 33 + byte) & 0xFFFFFFFF
	return hash


# Distance between the two anchors of a joint, in world space.
func joint_error(joint: AVBDJoint3D, a: Node3D, b: Node3D) -> float:
	var point_a: Vector3
	if a == null:
		point_a = joint.global_transform * joint.anchor_a
	else:
		point_a = a.to_global(joint.anchor_a)
	return point_a.distance_to(b.to_global(joint.anchor_b))


# --- runner -----------------------------------------------------------------

var checks := 0
var failures := 0
var strict := false

var _only := ""
var _list := false
var _digest_mode := false
var _evidence := {}


func _initialize() -> void:
	if DisplayServer.get_name() == "headless":
		# Nothing here draws, and the dummy renderer only adds noise.
		pass
	for arg in OS.get_cmdline_user_args():
		if arg == "--list":
			_list = true
		elif arg == "--digest":
			_digest_mode = true
		elif arg == "--strict":
			strict = true
		elif arg.begins_with("--only="):
			_only = arg.split("=")[1]
	_run()


func _run() -> void:
	var scenarios := _scenarios()

	if _list:
		for entry in scenarios:
			print(entry.name)
		quit(0)
		return

	var selected := 0
	for entry in scenarios:
		if not _only.is_empty() and entry.name != _only:
			continue
		selected += 1
		print("\n== %s" % entry.name)
		var before := checks
		var evidence: Variant = await (entry.run as Callable).call()
		_evidence[entry.name] = evidence
		if evidence != null:
			print("     evidence: %s" % _evidence_text(evidence))
		if _digest_mode:
			var repeat: Variant = await (entry.run as Callable).call()
			var same: bool = _evidence_text(repeat) == _evidence_text(evidence)
			check(same, "%s: reproducible" % entry.name,
					"digest %s vs %s" % [_evidence_text(evidence).hash(), _evidence_text(repeat).hash()])
			_evidence["%s#2" % entry.name] = repeat
		if before == checks:
			print("     (this scenario made no assertions)")

	if selected == 0:
		printerr("no scenario matched --only=%s" % _only)

	if _digest_mode:
		print("\nfull-flow digest: %s" % _evidence_text(_evidence).hash())

	_report_and_quit()


func _evidence_text(value: Variant) -> String:
	# Dictionaries print in insertion order, which is deterministic.
	return str(value)


func _report_and_quit() -> void:
	print("\n%d checks, %d failures" % [checks, failures])
	quit(0 if failures == 0 else 1)
