# Shared harness for the AVBD Godot-side test suites.
#
# A suite subclasses this script and returns its scenarios from `_scenarios()`:
#
#   extends "res://tests/avbd_harness.gd"
#   func _scenarios() -> Array: return [scenario("name", run_thing)]
#
# Every suite drives the engine's physics loop with the AVBD physics server
# (physics/3d/physics_engine = "AVBD") through standard Godot nodes: each physics
# frame is exactly one server step at the 60 Hz headless tick, so `steps(n)` waits
# n physics frames. Nothing here touches an AVBD-specific class except the
# registration checks, which is the point: this suite pins the standard-node path.
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


# The engine must have selected AVBD: every scenario below walks the server path.
func check_engine() -> void:
	var engine: String = ProjectSettings.get_setting("physics/3d/physics_engine")
	check(engine == "AVBD", "physics engine is AVBD", "engine=%s" % engine)


# --- simulation helpers -----------------------------------------------------

# Wait for `count` physics ticks; the server steps once per tick (60 Hz headless).
func steps(count: int) -> void:
	for i in count:
		await physics_frame


func frames(count: int) -> void:
	await steps(count)


# A static ground slab: full extents `half * 2` horizontally, 1 m thick, top face at
# y = `top`. Built from standard nodes so the server sees what a game would send.
func add_ground(parent: Node, top := 0.0, friction := 0.5, half := 50.0) -> StaticBody3D:
	var ground := StaticBody3D.new()
	ground.name = "Ground"
	var shape := CollisionShape3D.new()
	shape.name = "Shape"
	var box := BoxShape3D.new()
	box.size = Vector3(half * 2.0, 1.0, half * 2.0)
	shape.shape = box
	ground.add_child(shape)
	# PhysicsMaterial per body (friction pairs through Godot's own node API).
	var material := PhysicsMaterial.new()
	material.friction = friction
	ground.physics_material_override = material
	ground.position = Vector3(0, top - 0.5, 0)
	parent.add_child(ground)
	return ground


# A dynamic box: `size` full extents, `mass` in kg, dropped at `position`.
func add_body(parent: Node, body_name: String, size: Vector3, position: Vector3, mass := 1.0,
		friction := 0.5) -> RigidBody3D:
	var body := RigidBody3D.new()
	body.name = body_name
	var shape := CollisionShape3D.new()
	shape.name = "Shape"
	var box := BoxShape3D.new()
	box.size = size
	shape.shape = box
	body.add_child(shape)
	var material := PhysicsMaterial.new()
	material.friction = friction
	body.physics_material_override = material
	body.mass = mass
	body.position = position
	parent.add_child(body)
	return body


# --- measurement helpers ----------------------------------------------------

func is_finite_body(body: RigidBody3D) -> bool:
	var p := body.global_position
	var v := body.linear_velocity
	return p.is_finite() and v.is_finite() and body.angular_velocity.is_finite()


func all_finite(bodies: Array) -> bool:
	for body in bodies:
		if body is RigidBody3D and not is_finite_body(body):
			return false
	return true


func max_speed(bodies: Array) -> float:
	var worst := 0.0
	for body in bodies:
		if body is RigidBody3D:
			worst = maxf(worst, body.linear_velocity.length())
	return worst


func ground_clearance(bodies: Array, floor_y: float) -> float:
	var lowest := INF
	for body in bodies:
		if body is RigidBody3D:
			lowest = minf(lowest, body.global_position.y)
	return lowest - floor_y


# Pose digest of every given body, for comparing two runs. Transforms are quantized to
# a millimetre so a one-tick measurement offset between two runs does not flip the
# hash, while real divergence (centimetres) does.
func digest_of(bodies: Array) -> int:
	var hash := 5381
	for body in bodies:
		var t: Transform3D = (body as Node3D).global_transform
		t.origin = (t.origin * 1000.0).round() / 1000.0
		for byte in var_to_bytes(t):
			hash = (hash * 33 + byte) & 0xFFFFFFFF
	return hash


# Distance between two bodies' closest anchors (world space), for joint assertions:
# `anchor_a` / `anchor_b` are local offsets on each body.
func anchor_error(a: Node3D, b: Node3D, anchor_a: Vector3, anchor_b: Vector3) -> float:
	return a.to_global(anchor_a).distance_to(b.to_global(anchor_b))


# --- runner -----------------------------------------------------------------

var checks := 0
var failures := 0
var strict := false

var _only := ""
var _list := false
var _digest_mode := false
var _evidence := {}


func _initialize() -> void:
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
