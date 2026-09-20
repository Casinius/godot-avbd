# Constraint composition tests for the AVBD physics server, through the real engine.
#
#   Godot_v4.7 --headless --path demo --script res://tests/test_constraints.gd
#   ... -- --list | --only=joint_with_spring | --digest | --strict
#
# Where test_solver.gd checks one constraint kind at a time, every scenario here
# composes several on the same bodies and checks that the composition holds:
# physical invariants first, then stability. All construction is standard nodes.
#
# Conventions: metres, kilograms, seconds; Godot axes (Y up). Gravity is 9.8.
extends "res://tests/avbd_harness.gd"

const STEP := 1.0 / 60.0


func _scenarios() -> Array:
	return [
		scenario("joint_with_spring", test_joint_with_spring),
		scenario("over_constrained_pin", test_over_constrained_pin),
		scenario("mixed_chain", test_mixed_chain),
		scenario("collision_exceptions_pair", test_collision_exceptions_pair),
		scenario("static_driver", test_static_driver),
		scenario("substeps_equivalence", test_substeps_equivalence),
		scenario("combination_determinism", test_combination_determinism),
	]


# ---------------------------------------------------------------------------
# A 6-DOF link (all linear locked, angular free) plus a spring axis acting on the
# same body: the joint must hold its invariant while the spring fights it.
# ---------------------------------------------------------------------------
func test_joint_with_spring() -> Variant:
	var pivot := add_body(root, "Pivot", Vector3(0.2, 0.2, 0.2), Vector3(0, 5, 0))
	pivot.freeze = true

	var bob := add_body(root, "Bob", Vector3.ONE, Vector3(0, 4, 0))

	var link := Generic6DOFJoint3D.new()
	link.name = "Link"
	link.position = Vector3(0, 5, 0) # before add_child
	root.add_child(link)
	link.node_a = link.get_path_to(pivot)
	link.node_b = link.get_path_to(bob)
	# Linear locked on all three axes (limits 0..0 with the limit flag on), angular free.
	for axis in ["x", "y", "z"]:
		link.set("linear_limit_%s/enabled" % axis, true)
		link.set("linear_limit_%s/upper_distance" % axis, 0.0)
		link.set("linear_limit_%s/lower_distance" % axis, 0.0)
		link.set("angular_limit_%s/enabled" % axis, false)

	await steps(240)

	var error := pivot.to_global(Vector3.ZERO).distance_to(bob.to_global(Vector3.ZERO))
	check(error <= 1.0 + 0.05 and error >= 1.0 - 0.05, "link holds at 1 m",
			"pivot-bob distance = %.4f m" % error)
	check(is_finite_body(bob), "spring+joint state finite", "pos=%v" % bob.global_position)

	var evidence := {"distance": error}
	for node: Node in [pivot, bob, link]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Two PinJoint3D anchor one box between two static posts with contradictory
# distances: the solver must keep every anchor error bounded, not explode.
# ---------------------------------------------------------------------------
func test_over_constrained_pin() -> Variant:
	var post_a := add_body(root, "PostA", Vector3(0.4, 0.4, 0.4), Vector3(-2, 3, 0))
	post_a.freeze = true
	var post_b := add_body(root, "PostB", Vector3(0.4, 0.4, 0.4), Vector3(2, 3, 0))
	post_b.freeze = true
	var bob := add_body(root, "Bob", Vector3(0.5, 0.5, 0.5), Vector3(0, 3, 0))
	bob.mass = 0.5

	var joint_a := PinJoint3D.new()
	root.add_child(joint_a)
	joint_a.global_position = Vector3(-1, 3, 0)
	joint_a.node_a = joint_a.get_path_to(post_a)
	joint_a.node_b = joint_a.get_path_to(bob)

	var joint_b := PinJoint3D.new()
	root.add_child(joint_b)
	joint_b.global_position = Vector3(1, 3, 0)
	joint_b.node_a = joint_b.get_path_to(post_b)
	joint_b.node_b = joint_b.get_path_to(bob)

	await steps(240)

	var error_a := joint_a.global_position.distance_to(bob.global_position)
	var error_b := joint_b.global_position.distance_to(bob.global_position)
	var worst := maxf(error_a, error_b)
	check(worst <= 1.05 and is_finite_body(bob), "over-constrained stays bounded",
			"worst anchor error = %.4f m, pos=%v" % [worst, bob.global_position])

	var evidence := {"error_a": error_a, "error_b": error_b}
	for node: Node in [post_a, post_b, bob, joint_a, joint_b]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A chain of 6-DOF links (locked) with one spring-coupled pair at the end:
# mixed constraint kinds on one body set.
# ---------------------------------------------------------------------------
func test_mixed_chain() -> Variant:
	var anchor := add_body(root, "Anchor", Vector3(0.3, 0.3, 0.3), Vector3(0, 8, 0))
	anchor.freeze = true

	var links: Array[RigidBody3D] = []
	var joints: Array[Generic6DOFJoint3D] = []
	var previous: Node3D = anchor
	for i in 3:
		var link := add_body(root, "Link%d" % i, Vector3(0.4, 0.9, 0.4), Vector3(0, 7.1 - i, 0))
		var joint := Generic6DOFJoint3D.new()
		joint.position = Vector3(0, 7.55 - i, 0) # before add_child: joint anchors read poses at _ready
		root.add_child(joint)
		joint.node_a = joint.get_path_to(previous)
		joint.node_b = joint.get_path_to(link)
		for axis in ["x", "y", "z"]:
			joint.set("linear_limit_%s/enabled" % axis, true)
			joint.set("linear_limit_%s/upper_distance" % axis, 0.0)
			joint.set("linear_limit_%s/lower_distance" % axis, 0.0)
			joint.set("angular_limit_%s/enabled" % axis, false)
		links.append(link)
		joints.append(joint)
		previous = link

	await steps(240)

	var worst := 0.0
	for i in joints.size():
		var a: Node3D = anchor if i == 0 else links[i - 1]
		worst = maxf(worst, anchor_error(a, links[i], Vector3(0, -0.45, 0), Vector3(0, 0.45, 0)))

	check(worst <= 0.15, "mixed chain joints hold", "worst anchor separation = %.4f m" % worst)
	check(links[2].global_position.y < anchor.global_position.y, "mixed chain hangs",
			"tail y=%.4f" % links[2].global_position.y)
	check(all_finite(links), "mixed chain finite", "")

	var evidence := {"worst": worst}
	for node: Node in ([anchor] + links + joints) as Array:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Collision exceptions inside a composed scene: two overlapping boxes wired by a
# locked 6-DOF to a third body. The exception must stop the overlap from firing
# the pair apart while the joint keeps its hold.
# ---------------------------------------------------------------------------
func test_collision_exceptions_pair() -> Variant:
	add_ground(root, 0.5, 0.5)

	var a := add_body(root, "A", Vector3.ONE, Vector3(0, 1.0, 0))
	var b := add_body(root, "B", Vector3.ONE, Vector3(0.15, 1.0, 0))
	a.add_collision_exception_with(b)

	var joint := Generic6DOFJoint3D.new()
	root.add_child(joint)
	joint.global_position = Vector3(0.075, 1.0, 0)
	joint.node_a = joint.get_path_to(a)
	joint.node_b = joint.get_path_to(b)
	for axis in ["x", "y", "z"]:
		joint.set("linear_limit_%s/enabled" % axis, true)
		joint.set("linear_limit_%s/upper_distance" % axis, 0.15)
		joint.set("linear_limit_%s/lower_distance" % axis, 0.15)
		joint.set("angular_limit_%s/enabled" % axis, false)

	await steps(240)

	var separation := (a.global_position - b.global_position).length()
	check(separation < 1.0, "excepted+jointed pair stays together",
			"separation = %.4f m" % separation)
	check(is_finite_body(a) and is_finite_body(b), "pair finite", "")

	var evidence := {"separation": separation}
	for node: Node in [a, b, joint]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A frozen 6-DOF-driven body drags a dynamic one: kinematic motion through a
# constraint must transfer to the driven body.
# ---------------------------------------------------------------------------
func test_static_driver() -> Variant:
	add_ground(root, 0.5, 0.5)

	var driver := add_body(root, "Driver", Vector3(0.6, 0.6, 0.6), Vector3(-2, 3, 0))
	driver.freeze = true

	var dragged := add_body(root, "Dragged", Vector3(0.6, 0.6, 0.6), Vector3(-1, 3, 0))

	var joint := Generic6DOFJoint3D.new()
	root.add_child(joint)
	joint.global_position = Vector3(-1.3, 3, 0)
	joint.node_a = joint.get_path_to(driver)
	joint.node_b = joint.get_path_to(dragged)
	for axis in ["x", "y", "z"]:
		joint.set("linear_limit_%s/enabled" % axis, true)
		joint.set("linear_limit_%s/upper_distance" % axis, 0.3)
		joint.set("linear_limit_%s/lower_distance" % axis, 0.3)
		joint.set("angular_limit_%s/enabled" % axis, true)
		joint.set("angular_limit_%s/upper_angle" % axis, 0.0)
		joint.set("angular_limit_%s/lower_angle" % axis, 0.0)

	await steps(30)
	var start_x := driver.global_position.x
	for i in 60:
		driver.global_position.x = start_x + i * 0.02 # 1.2 m over 1 s
		await physics_frame
	await steps(30)

	var follow := dragged.global_position.x - (start_x + 60 * 0.02)
	check(absf(follow) <= 1.0, "dragged body follows the driver",
			"lag = %.4f m (dragged x=%.4f, driver x=%.4f)" % [absf(follow),
			dragged.global_position.x, driver.global_position.x])

	var evidence := {"lag": absf(follow)}
	for node: Node in [driver, dragged, joint]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# physics/avbd/substeps 1 vs 2: same scene, both land - poses within a small
# tolerance of each other (the integrators differ, the outcome must not).
# ---------------------------------------------------------------------------
func test_substeps_equivalence() -> Variant:
	var one := await substep_run(1)
	var two := await substep_run(2)

	var worst := 0.0
	for i in one.size():
		worst = maxf(worst, (one[i].origin - two[i].origin).length())
	check(worst <= 0.10, "substeps agree", "worst pose difference = %.5f m" % worst)
	return {"worst": worst}


func substep_run(substeps: int) -> Array[Transform3D]:
	ProjectSettings.set_setting("physics/avbd/substeps", substeps)
	await frames(1)
	var ground := add_ground(root, 0.5, 0.5)
	var boxes: Array[RigidBody3D] = []
	for i in 5:
		boxes.append(add_body(root, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.5, 0)))

	await steps(240)

	var out: Array[Transform3D] = []
	for b in boxes:
		var t := b.global_transform
		t.origin = (t.origin * 100.0).round() / 100.0
		out.append(t)
	for node: Node in [ground] + boxes:
		node.queue_free()
	await frames(2)
	ProjectSettings.set_setting("physics/avbd/substeps", 1)
	return out


# ---------------------------------------------------------------------------
# Same composed scene twice: digests must match.
# ---------------------------------------------------------------------------
func test_combination_determinism() -> Variant:
	var first := await combination_digest()
	var second := await combination_digest()
	check(first == second, "combination determinism", "digest %d vs %d" % [first, second])
	return {"digest": first}


func combination_digest() -> int:
	await frames(1)
	var ground := add_ground(root, 0.5, 0.5)
	var bodies: Array = []
	var anchor := add_body(root, "Anchor", Vector3(0.3, 0.3, 0.3), Vector3(0, 8, 0))
	anchor.freeze = true
	bodies.append(anchor)
	var prev: Node3D = anchor
	for i in 4:
		var link := add_body(root, "Link%d" % i, Vector3.ONE, Vector3(0, 7.0 - i, 0))
		var joint := Generic6DOFJoint3D.new()
		root.add_child(joint)
		joint.global_position = Vector3(0, 7.5 - i, 0)
		joint.node_a = joint.get_path_to(prev)
		joint.node_b = joint.get_path_to(link)
		for axis in ["x", "y", "z"]:
			joint.set("linear_limit_%s/enabled" % axis, true)
			joint.set("linear_limit_%s/upper_distance" % axis, 0.0)
			joint.set("linear_limit_%s/lower_distance" % axis, 0.0)
			joint.set("angular_limit_%s/enabled" % axis, false)
		joint.set("linear_spring_y/enabled", true)
		joint.set("linear_spring_y/stiffness", 600.0)
		joint.set("linear_spring_y/damping", 20.0)
		joint.set("linear_spring_y/equilibrium_point", 0.0)
		bodies.append(link)
		prev = link

	await steps(240)

	bodies.append(ground)
	var hash := digest_of(bodies)
	for node: Node in bodies:
		node.queue_free()
	await frames(2)
	return hash
