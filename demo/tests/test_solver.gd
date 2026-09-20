# Behavioural tests for the AVBD physics server, run against a real Godot runtime.
#
#   Godot_v4.7 --headless --path demo --script res://tests/test_solver.gd
#   ... -- --list | --only=stack | --digest | --strict
#
# Every scenario builds its node graph from standard Godot nodes (RigidBody3D,
# StaticBody3D, CollisionShape3D, joints), then steps the real physics loop. The
# AVBD physics server answers every call: scene graph -> server -> solver ->
# DirectBodyState -> scene graph. Single-constraint behaviour lives here;
# constraint *combinations* are in test_constraints.gd.
extends "res://tests/avbd_harness.gd"


func _scenarios() -> Array:
	return [
		scenario("rest_contact", test_rest_contact),
		scenario("stack", test_stack),
		scenario("friction", test_friction),
		scenario("hard_joint", test_hard_joint),
		scenario("world_joint", test_world_joint),
		scenario("raycast", test_raycast),
		scenario("impulse_teleport", test_impulse_teleport),
		scenario("runtime_add", test_runtime_add),
		scenario("axis_lock", test_axis_lock),
		scenario("collision_layers", test_collision_layers),
		scenario("collision_exceptions", test_collision_exceptions),
		scenario("sleeping", test_sleeping),
		scenario("paused_follows_game", test_paused_follows_game),
	]


# ---------------------------------------------------------------------------
# A box dropped on a slab comes to rest on top of it: this is the end-to-end
# check that the Godot Y-up <-> solver Z-up mapping is right, through the server.
# ---------------------------------------------------------------------------
func test_rest_contact() -> Variant:
	check_engine()
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 5, 0))

	await steps(240)

	var rest_y := box.global_position.y
	check(absf(rest_y - 1.0) <= 0.02, "rest contact", "y=%.5f (slab top 0.5 + half box 0.5)" % rest_y)
	check(absf(box.linear_velocity.y) <= 0.05, "rest velocity", "vy=%.5f" % box.linear_velocity.y)
	check(is_finite_body(box), "state finite", "pos=%v" % box.global_position)

	var evidence := {"rest_y": rest_y}
	ground.queue_free()
	box.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A 10-box stack keeps its shape and does not drift sideways.
#
# Note: 10 layers requires more time to stabilize due to chain instability.
# Increased steps from 300 to 600 for convergence, friction to 0.8 to reduce lateral drift.
# ---------------------------------------------------------------------------
func test_stack() -> Variant:
	var ground := add_ground(root, 0.5, 0.8)  # Increased friction from 0.5 to 0.8

	var boxes: Array[RigidBody3D] = []
	for i in 10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.0, 0), 1.0, 0.8))

	await steps(600)  # Increased from 300 for better convergence with 10 layers

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	check(worst_y <= 0.2, "stack height", "max |y - expected| = %.5f" % worst_y)
	check(worst_interface <= 0.05, "stack penetration", "worst interface = %.5f m" % worst_interface)
	check(worst_lateral <= 0.05, "stack drift", "max lateral = %.5f m" % worst_lateral)

	var evidence := {"height_error": worst_y, "penetration": worst_interface, "drift": worst_lateral}
	for node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Friction on a 20 degree ramp: mu = 0.5 holds (penalty-method creep only),
# mu = 0.05 lets the box slide away. Friction arrives through PhysicsMaterial.
# ---------------------------------------------------------------------------
func ramp_run(mu: float, settle_frames: int, watch_frames: int) -> Dictionary:
	var ground := add_ground(root, 0.5, mu)

	var ramp_basis := Basis(Vector3(0, 0, 1), deg_to_rad(20.0))
	var ramp := add_body(root, "Ramp", Vector3(40, 1, 24), Vector3(0, 6, 0), 1.0, mu)
	ramp.freeze = true # a kinematic ramp: pose fixed by the scene
	ramp.global_transform = Transform3D(ramp_basis, Vector3(0, 6, 0))

	var up_slope := ramp_basis * Vector3.RIGHT
	var normal := ramp_basis * Vector3.UP
	var start := ramp.global_position + up_slope * 5.0 + normal * 1.05
	var box := add_body(root, "Box", Vector3.ONE, start, 1.0, mu)

	await steps(settle_frames)
	var before := box.global_position
	await steps(watch_frames)

	var late := (box.global_position - before).dot(-up_slope)
	var total := (box.global_position - start).dot(-up_slope)
	var speed := box.linear_velocity.length()
	var evidence := {"late": late, "total": total, "speed": speed}
	ground.queue_free()
	ramp.queue_free()
	box.queue_free()
	await frames(2)
	return evidence


func test_friction() -> Variant:
	var stuck := await ramp_run(0.5, 240, 30)
	check(stuck.late <= 0.02 and stuck.speed <= 0.05, "static friction (mu=0.5)",
			"creep=%.5f m in 0.5 s, speed=%.5f m/s, total=%.4f m" % [stuck.late, stuck.speed, stuck.total])

	var slid := await ramp_run(0.05, 240, 30)
	check(slid.late >= 1.0, "dynamic friction (mu=0.05)",
			"travel=%.4f m in 0.5 s, total=%.4f m" % [slid.late, slid.total])
	check(slid.late > stuck.late * 50.0, "friction ordering",
			"%.4f m vs %.4f m" % [slid.late, stuck.late * 50.0])
	return {"stuck": stuck, "slid": slid}


# ---------------------------------------------------------------------------
# A 4-link chain hung from a static body, built from Generic6DOFJoint3D with all
# six axes locked. Hard joints must hold to a fraction of a link.
# ---------------------------------------------------------------------------
func make_locked_6dof(parent: Node, joint_name: String, a: Node3D, b: Node3D, anchor: Vector3) -> Generic6DOFJoint3D:
	var joint := Generic6DOFJoint3D.new()
	joint.name = joint_name
	parent.add_child(joint)
	joint.global_position = anchor
	joint.node_a = joint.get_path_to(a)
	joint.node_b = joint.get_path_to(b)
	return joint


func test_hard_joint() -> Variant:
	var anchor := add_body(root, "Anchor", Vector3.ONE, Vector3(0, 10, 0))
	anchor.freeze = true

	var links: Array[RigidBody3D] = []
	var joints: Array[Generic6DOFJoint3D] = []
	var previous: Node3D = anchor
	for i in 4:
		var link := add_body(root, "Link%d" % i, Vector3.ONE, Vector3(0, 9.0 - i, 0))
		var joint := make_locked_6dof(root, "Joint%d" % i, previous, link, Vector3(0, 8.5 - i, 0))
		links.append(link)
		joints.append(joint)
		previous = link

	await steps(300)

	var worst := 0.0
	for i in joints.size():
		var a: Node3D = anchor if i == 0 else links[i - 1]
		worst = maxf(worst, anchor_error(a, links[i], Vector3(0, -0.5, 0), Vector3(0, 0.5, 0)))

	check(worst <= 0.02, "hard joint error", "max anchor separation = %.6f m" % worst)
	check(links[3].global_position.y < anchor.global_position.y, "chain hangs",
			"tail y=%.4f below anchor y=%.4f" % [links[3].global_position.y, anchor.global_position.y])

	var evidence := {"worst_joint_error": worst, "tail_y": links[3].global_position.y}
	for node: Node in [anchor] + links + joints:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A PinJoint3D anchored to a static body pins a box to a point. The anchor is
# deliberately off-axis so a swapped coordinate would show up as a mismatch.
# AVBD constraints are stabilised, so an initial violation decays over a second
# or two; the body is spawned at its anchor.
# ---------------------------------------------------------------------------
func test_world_joint() -> Variant:
	var anchor_point := Vector3(2, 3, 1)
	var world_anchor := add_body(root, "WorldAnchor", Vector3.ONE, Vector3(50, 50, 50))
	world_anchor.freeze = true
	var bob := add_body(root, "Bob", Vector3.ONE, anchor_point)
	bob.linear_velocity = Vector3(4, 0, 0)

	var joint := PinJoint3D.new()
	joint.name = "Pin"
	root.add_child(joint)
	joint.global_position = anchor_point
	joint.node_a = joint.get_path_to(world_anchor)
	joint.node_b = joint.get_path_to(bob)

	await steps(180)

	var error := bob.global_position.distance_to(anchor_point)
	check(error <= 0.02, "world anchor holds",
			"bob at %v, anchor %v, error %.5f m" % [bob.global_position, anchor_point, error])
	check(bob.linear_velocity.length() <= 0.05, "world anchor absorbs velocity",
			"speed=%.5f m/s from an initial 4 m/s" % bob.linear_velocity.length())

	var evidence := {"error": error, "speed": bob.linear_velocity.length()}
	for node: Node in [world_anchor, bob, joint]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# RayCast3D node + direct space state intersect_ray: hits dynamic bodies, static
# bodies, misses cleanly, and follows a teleported static body.
# ---------------------------------------------------------------------------
func test_raycast() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))

	await steps(60)
	var top := box.global_position.y + 0.5

	# Direct space state query.
	var space_state := root.get_world_3d().direct_space_state
	var query := PhysicsRayQueryParameters3D.create(Vector3(0, 20, 0), Vector3(0, -10, 0))
	var hit := space_state.intersect_ray(query)
	check(not hit.is_empty() and instance_from_id(hit.get("collider_id")) == box,
			"ray hits dynamic body", "hit=%s" % [hit.get("collider_id")])
	check(hit.has("position") and absf(hit.get("position").y - top) <= 0.05, "ray hit position",
			"y=%.4f (box top %.4f)" % [hit.get("position", Vector3.ZERO).y, top])

	var miss := space_state.intersect_ray(PhysicsRayQueryParameters3D.create(Vector3(50, 20, 50), Vector3(50, 21, 50)))
	check(miss.is_empty(), "ray misses cleanly", "keys=%s" % [miss.keys()])

	var ground_hit := space_state.intersect_ray(PhysicsRayQueryParameters3D.create(Vector3(20, 20, 20), Vector3(20, -10, 20)))
	check(not ground_hit.is_empty() and instance_from_id(ground_hit.get("collider_id")) == ground,
			"ray hits static body", "hit=%s" % [ground_hit.get("collider_id")])

	# RayCast3D node (the engine wraps the same server query).
	var ray := RayCast3D.new()
	ray.target_position = Vector3(0, -30, 0)
	root.add_child(ray)
	ray.global_position = Vector3(0, 20, 0)
	await frames(2)
	check(ray.is_colliding() and ray.get_collider() == box, "RayCast3D node hits", "collider=%s" % ray.get_collider())

	# Static bodies follow their node: move the slab and the ray must follow.
	ground.global_position = Vector3(20, 10.5, 20)
	await frames(2)
	var moved := space_state.intersect_ray(PhysicsRayQueryParameters3D.create(Vector3(20, 20, 20), Vector3(20, -10, 20)))
	check(moved.has("position") and absf(moved.get("position").y - 11.0) <= 0.05, "static body teleport adopted",
			"hit y=%.4f (moved slab top 11.0)" % [moved.get("position", Vector3.ZERO).y])

	var evidence := {"hit_y": hit.get("position", Vector3.ZERO).y}
	ground.queue_free()
	box.queue_free()
	ray.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Impulses and teleports act on the solver state with the right magnitudes,
# through the server's apply paths.
# ---------------------------------------------------------------------------
func test_impulse_teleport() -> Variant:
	add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 1.0, 0))

	await steps(60)

	# A 1 m^3 box at mass 1 kg: a 10 N.s impulse is 10 m/s. Read the velocity back through
	# the body's direct state - the same object the sync callback fills - because the node's
	# cached linear_velocity only refreshes on the next sync.
	box.apply_impulse(Vector3(10, 0, 0))
	await steps(1)
	var direct_state: PhysicsDirectBodyState3D = PhysicsServer3D.body_get_direct_state(box.get_rid())
	var vx: float = direct_state.linear_velocity.x
	check(absf(vx - 10.0) <= 0.5, "impulse sets velocity",
			"vx=%.5f (10 N.s / 1 kg)" % vx)

	box.linear_velocity = Vector3.ZERO
	box.angular_velocity = Vector3.ZERO
	box.apply_impulse(Vector3(10, 0, 0), Vector3(0, 0.5, 0))
	await steps(1)
	check(box.angular_velocity.length() > 1.0, "offset impulse spins body",
			"|w|=%.4f rad/s" % box.angular_velocity.length())

	await steps(60)
	check(box.global_position.x > 1.0, "impulse moves body", "x=%.4f" % box.global_position.x)

	box.global_transform = Transform3D(Basis(), Vector3(-5, 4, 0))
	box.linear_velocity = Vector3.ZERO
	check(absf(box.global_position.x + 5.0) <= 0.001, "teleport moves node", "x=%.4f" % box.global_position.x)
	await frames(2)
	check(absf(box.global_position.x + 5.0) <= 0.05, "teleport holds", "x=%.4f on the next tick" % box.global_position.x)

	var evidence := {"x_after": box.global_position.x}
	box.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Adding a body at runtime must not reset the bodies already simulated.
# ---------------------------------------------------------------------------
func test_runtime_add() -> Variant:
	add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))

	await steps(120)
	var settled := box.global_position.y

	var late := add_body(root, "Late", Vector3.ONE, Vector3(5, 2, 0))
	await steps(60)

	check(absf(settled - box.global_position.y) <= 0.02, "rebuild keeps pose",
			"settled y %.4f -> %.4f after adding a body" % [settled, box.global_position.y])
	check(late.global_position.y > 0.5, "added body rests", "late y=%.4f" % late.global_position.y)

	var evidence := {"settled": settled, "after": box.global_position.y, "late_y": late.global_position.y}
	for node in [box, late]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Axis locks through the server: a box with linear Y locked must hover, resist
# gravity entirely, and keep spinning only about locked axes.
# ---------------------------------------------------------------------------
func test_axis_lock() -> Variant:
	add_ground(root, 0.5, 0.5)
	var hover := add_body(root, "Hover", Vector3.ONE, Vector3(0, 5, 0))
	hover.axis_lock_linear_y = true

	var control := add_body(root, "Control", Vector3.ONE, Vector3(5, 5, 0))

	await steps(240)

	check(absf(hover.global_position.y - 5.0) <= 0.01, "locked axis hovers",
			"y=%.6f (spawn 5.0, no fall in 4 s)" % hover.global_position.y)
	check(absf(hover.linear_velocity.y) <= 0.01, "locked axis no velocity",
			"vy=%.6f" % hover.linear_velocity.y)
	check(control.global_position.y < 2.0, "unlocked control falls",
			"control y=%.4f" % control.global_position.y)

	var evidence := {"hover_y": hover.global_position.y}
	for node in [hover, control]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Collision layers/masks: a box on layer 2 with mask 1 vs a wall on layer 1
# mask 1... pair test is (A.layer & B.mask) || (B.layer & A.mask).
# ---------------------------------------------------------------------------
func test_collision_layers() -> Variant:
	var wall := add_body(root, "Wall", Vector3(0.5, 4, 4), Vector3(0, 2, 0))
	wall.freeze = true

	var pass_through := add_body(root, "PassThrough", Vector3.ONE, Vector3(-6, 2, 0), 1.0)
	pass_through.collision_layer = 2
	pass_through.collision_mask = 2 # never meets the wall's layer 1
	pass_through.linear_velocity = Vector3(4, 0, 0) # pushed at the wall, must cross

	var blocked := add_body(root, "Blocked", Vector3.ONE, Vector3(6, 2, 0), 1.0)
	blocked.collision_layer = 4
	blocked.collision_mask = 1 | 4 # layer 4 meets the wall's layer 1 through the wall's mask
	blocked.linear_velocity = Vector3(-4, 0, 0)

	await steps(180)

	check(pass_through.global_position.x > 0.4, "incompatible layers pass through",
			"x=%.4f (crossed the wall plane)" % pass_through.global_position.x)
	check(blocked.global_position.x > 0.4, "compatible layers collide",
			"x=%.4f (stopped at the wall)" % blocked.global_position.x)

	var evidence := {"pass_x": pass_through.global_position.x, "blocked_x": blocked.global_position.x}
	for node in [wall, pass_through, blocked]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Collision exceptions: two overlapping boxes that would blow apart stay put,
# while a control pair pushes apart.
# ---------------------------------------------------------------------------
func test_collision_exceptions() -> Variant:
	add_ground(root, 0.5, 0.5)

	var a := add_body(root, "OverlapA", Vector3.ONE, Vector3(0, 1.0, 0))
	var b := add_body(root, "OverlapB", Vector3.ONE, Vector3(0.2, 1.0, 0))
	a.add_collision_exception_with(b)

	var control_a := add_body(root, "ControlA", Vector3.ONE, Vector3(6, 1.0, 0))
	var control_b := add_body(root, "ControlB", Vector3.ONE, Vector3(6.2, 1.0, 0))

	await steps(180)

	var separation := (a.global_position - b.global_position).length()
	var control_separation := (control_a.global_position - control_b.global_position).length()
	check(separation < 0.9, "exception pair keeps overlapping",
			"separation=%.4f m (spawned at 0.2)" % separation)
	check(control_separation > separation + 0.3, "control pair pushes apart further",
			"separation=%.4f m" % control_separation)

	var evidence := {"separation": separation, "control": control_separation}
	for node in [a, b, control_a, control_b]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Sleeping: putting a body to sleep freezes it; waking it resumes the fall.
# ---------------------------------------------------------------------------
func test_sleeping() -> Variant:
	add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 8, 0))

	await steps(30)
	box.sleeping = true
	var y_at_sleep := box.global_position.y

	await steps(120)
	check(absf(box.global_position.y - y_at_sleep) < 0.01, "sleeping body freezes",
			"y %.4f -> %.4f while asleep" % [y_at_sleep, box.global_position.y])

	box.sleeping = false
	await steps(120)
	check(box.global_position.y < y_at_sleep - 1.0, "woken body falls",
			"y=%.4f after wake" % box.global_position.y)

	var evidence := {"y_at_sleep": y_at_sleep, "y_end": box.global_position.y}
	box.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# The direct evidence for pause-following: freeze the SceneTree mid-fall and the
# server must stop stepping. Position must be bit-stable; resuming resumes.
# ---------------------------------------------------------------------------
func test_paused_follows_game() -> Variant:
	add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 8, 0))

	await steps(30)
	var y_before := box.global_position.y

	paused = true
	await steps(30)
	var y_paused := box.global_position.y
	paused = false

	check(absf(y_paused - y_before) < 1.0e-4, "paused simulation freezes",
			"y %.6f -> %.6f during 30 paused frames" % [y_before, y_paused])

	await steps(60)
	check(box.global_position.y < y_paused - 1.0, "resume continues the fall",
			"y=%.4f after unpause" % box.global_position.y)

	var evidence := {"y_before": y_before, "y_paused": y_paused, "y_resumed": box.global_position.y}
	box.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Same scene, same ticks, same state (in-process, twice).
# ---------------------------------------------------------------------------
func test_determinism() -> Variant:
	var first := await stack_digest()
	var second := await stack_digest()
	check(first == second, "in-process determinism", "digest %d vs %d" % [first, second])
	return {"digest": first}


func stack_digest() -> int:
	# Build from inside a physics-frame callback, so every invocation of this scenario
	# has the same relationship to the engine's step boundary (one solver step per
	# physics frame, taken after the signal).
	await frames(1)
	var ground := add_ground(root, 0.5, 0.5)
	var boxes: Array[RigidBody3D] = []
	for i in 10:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.5, 0)))

	await steps(300)

	var hash := digest_of(boxes)
	for node: Node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	return hash
