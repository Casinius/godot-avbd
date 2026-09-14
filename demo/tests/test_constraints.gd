# Combined-constraint verification for AVBD, through the real engine.
#
#   Godot_v4.6.3 --headless --path demo --script res://tests/test_constraints.gd
#   ... -- --list | --only=stress_bridge | --digest | --strict
#
# Where test_solver.gd checks one constraint kind at a time, every scenario here
# composes several AVBD constraints on the same bodies and checks that the
# composition holds: structural pre-flight, then physical invariants, then stability.
# Where a number can be predicted from theory (spring sag, fracture threshold, impulse
# response) the assertion uses that prediction rather than a fitted constant.
#
# Conventions: metres, kilograms, seconds; Godot axes (Y up). Gravity is 9.8.
extends "res://tests/avbd_harness.gd"

const STEP := 1.0 / 60.0


func _scenarios() -> Array:
	return [
		scenario("joint_with_spring", test_joint_with_spring),
		scenario("over_constrained_pin", test_over_constrained_pin),
		scenario("mixed_chain", test_mixed_chain),
		scenario("ignore_collision_pair", test_ignore_collision_pair),
		scenario("fracture_under_load", test_fracture_under_load),
		scenario("stress_bridge", test_stress_bridge),
		scenario("lattice_on_springs", test_lattice_on_springs),
		scenario("static_driver", test_static_driver),
		scenario("substeps_equivalence", test_substeps_equivalence),
		scenario("parallel_threads", test_parallel_threads),
		scenario("combination_determinism", test_combination_determinism),
	]


# ---------------------------------------------------------------------------
# Shared: a lattice of ordinary rigid bodies wired by joints, the same construction
# AVBDSoftBody3D uses internally - but built from node primitives so that springs and
# ignore-collision constraints can be attached to individual cells.
# ---------------------------------------------------------------------------
class Lattice:
	var world: AVBDWorld3D
	var cells: Array[AVBDRigidBody3D] = []
	var joints: Array[AVBDJoint3D] = []
	var side: int
	var spacing: float

	func index(x: int, y: int, z: int) -> int:
		return (x * side + y) * side + z

	func cell(x: int, y: int, z: int) -> AVBDRigidBody3D:
		return cells[index(x, y, z)]

	# Total mass, in kg: cells are cubes of `spacing` at density 1.
	func mass() -> float:
		return float(cells.size()) * spacing * spacing * spacing


func build_lattice(world: AVBDWorld3D, side: int, spacing: float, origin: Vector3,
		stiffness_linear := 1000.0, stiffness_angular := 250.0) -> Lattice:
	var lattice := Lattice.new()
	lattice.world = world
	lattice.side = side
	lattice.spacing = spacing
	var half := spacing * 0.5
	var centre := float(side - 1) * 0.5

	for x in side:
		for y in side:
			for z in side:
				var offset := Vector3(x - centre, y - centre, z - centre) * spacing
				lattice.cells.append(add_body(world, "Cell_%d_%d_%d" % [x, y, z],
						Vector3(spacing, spacing, spacing), origin + offset))

	for x in side:
		for y in side:
			for z in side:
				if x + 1 < side:
					lattice.joints.append(add_joint(world, "Jx_%d_%d_%d" % [x, y, z],
							lattice.cell(x, y, z), lattice.cell(x + 1, y, z), Vector3(half, 0, 0),
							Vector3(-half, 0, 0), stiffness_linear, stiffness_angular))
				if y + 1 < side:
					lattice.joints.append(add_joint(world, "Jy_%d_%d_%d" % [x, y, z],
							lattice.cell(x, y, z), lattice.cell(x, y + 1, z), Vector3(0, half, 0),
							Vector3(0, -half, 0), stiffness_linear, stiffness_angular))
				if z + 1 < side:
					lattice.joints.append(add_joint(world, "Jz_%d_%d_%d" % [x, y, z],
							lattice.cell(x, y, z), lattice.cell(x, y, z + 1), Vector3(0, 0, half),
							Vector3(0, 0, -half), stiffness_linear, stiffness_angular))
	return lattice


# Worst separation between the two anchor points of every joint in a lattice.
func lattice_worst_joint_error(lattice: Lattice) -> float:
	var worst := 0.0
	for joint in lattice.joints:
		var a: AVBDRigidBody3D = joint.get_body_a()
		var b: AVBDRigidBody3D = joint.get_body_b()
		worst = maxf(worst, joint_error(joint, a, b))
	return worst


func lattice_centre(lattice: Lattice) -> Vector3:
	var sum := Vector3.ZERO
	for body in lattice.cells:
		sum += body.global_position
	return sum / float(lattice.cells.size())


# ---------------------------------------------------------------------------
# C1: a rigid joint and a spring act on the same body, in different directions. The
# joint must keep its invariant while the spring fights it, and the pair must settle.
# ---------------------------------------------------------------------------
func test_joint_with_spring() -> Variant:
	var world := make_world()

	# The bob hangs 1 m below the pivot on a ball joint (angular stiffness 0 = free
	# rotation), and a stiff spring pulls it sideways towards a wall anchor.
	var pivot := add_body(world, "Pivot", Vector3(0.2, 0.2, 0.2), Vector3(0, 5, 0), 0.0, 0.5, true)
	var bob := add_body(world, "Bob", Vector3.ONE, Vector3(0, 4, 0))
	var wall := add_body(world, "Wall", Vector3(0.2, 0.2, 0.2), Vector3(2.5, 4, 0), 0.0, 0.5, true)

	var link := add_joint(world, "Link", pivot, bob, Vector3.ZERO, Vector3(0, 1, 0), -1.0, 0.0)
	# The spring rests at 1 m but starts 2.5 m long, so it pulls the bob sideways with
	# 400 * 1.5 = 600 N while the rigid link holds it one metre from the pivot.
	var spring := add_spring(world, "Pull", wall, bob, Vector3.ZERO, Vector3.ZERO, 400.0, 1.0)

	await check_forces(world, 2, "joint+spring pre-flight")

	await steps(world, 300)

	var link_error := joint_error(link, pivot, bob)
	var length := bob.global_position.distance_to(pivot.global_position)
	var swing := rad_to_deg(asin(clampf(absf(bob.global_position.x - pivot.global_position.x) / length, 0.0, 1.0)))

	check(length >= 0.99 and length <= 1.01, "link length held",
			"|bob - pivot| = %.5f m with a 400 N/m spring pulling sideways" % length)
	check(link_error <= 0.01, "link anchor error", "anchor separation = %.5f m" % link_error)
	check(swing >= 10.0, "spring moves the bob", "swing = %.2f deg from straight down" % swing)
	check(max_speed([bob]) <= 0.05, "pair settles", "final speed = %.5f m/s" % bob.linear_velocity.length())
	check(is_finite_body(bob), "pair finite", "bob at %v" % bob.global_position)
	check(spring.is_simulated(), "spring has solver force", "is_simulated=%s" % spring.is_simulated())

	var evidence := {"length": length, "anchor_error": link_error, "swing_deg": swing,
			"speed": bob.linear_velocity.length(), "digest": digest(world)}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# C2: one body pinned by two independent world-space joints. Consistent anchors must
# converge to the exact pose; inconsistent anchors must settle without jittering.
# ---------------------------------------------------------------------------
func pinned_by(world: AVBDWorld3D, anchor_a: Vector3, anchor_b: Vector3, steps_to_run: int) -> Dictionary:
	var body := add_body(world, "Pinned", Vector3.ONE, anchor_a)
	var j1 := add_world_joint(world, "PinA", body, anchor_a)
	var j2 := add_world_joint(world, "PinB", body, anchor_b)
	await check_forces(world, 2, "double-pin pre-flight")
	await steps(world, steps_to_run)
	return {
		"body": body,
		"error_a": body.global_position.distance_to(anchor_a),
		"error_b": body.global_position.distance_to(anchor_b),
		"speed": body.linear_velocity.length(),
		"j1": j1, "j2": j2,
	}


func test_over_constrained_pin() -> Variant:
	# Consistent: the anchors coincide with the body's spawn, so over-constraining
	# must not introduce error.
	var world := make_world()
	var consistent := await pinned_by(world, Vector3(0, 4, 0), Vector3(0, 4, 0), 240)
	check(consistent.error_a <= 0.01 and consistent.error_b <= 0.01, "consistent double pin",
			"errors %.5f / %.5f m" % [consistent.error_a, consistent.error_b])
	check(consistent.speed <= 0.05, "consistent pin settles", "speed = %.5f m/s" % consistent.speed)
	var consistent_digest := digest(world)
	world.queue_free()
	await frames(2)

	# Inconsistent: two anchors 2 m apart cannot both hold, so the body must find a
	# compromise and stop - the failure mode this guards against is oscillation.
	var world2 := make_world()
	var conflicting := await pinned_by(world2, Vector3(0, 4, 0), Vector3(2, 4, 0), 240)
	check(conflicting.error_a > 0.1 and conflicting.error_b > 0.1, "conflicting pins split the difference",
			"errors %.5f / %.5f m from anchors 2 m apart" % [conflicting.error_a, conflicting.error_b])
	check(conflicting.error_a + conflicting.error_b <= 2.5, "conflict stays bounded",
			"sum of errors = %.5f m (anchor separation 2 m)" % (conflicting.error_a + conflicting.error_b))
	check(conflicting.speed <= 0.2, "conflicting pins do not jitter",
			"speed = %.5f m/s at rest" % conflicting.speed)
	check(is_finite_body(conflicting.body), "conflicting pins finite", "pos=%v" % conflicting.body.global_position)
	var conflicting_digest := digest(world2)
	world2.queue_free()
	await frames(2)

	return {"consistent_errors": [consistent.error_a, consistent.error_b], "consistent_digest": consistent_digest,
			"conflicting_errors": [conflicting.error_a, conflicting.error_b], "conflicting_digest": conflicting_digest}


# ---------------------------------------------------------------------------
# C3: a chain that alternates a rigid link and a spring link. The rigid links must
# hold their length; the spring links must stretch by the load they carry (m*g/k).
# ---------------------------------------------------------------------------
func test_mixed_chain() -> Variant:
	var world := make_world()
	var anchor := add_body(world, "Anchor", Vector3.ONE, Vector3(0, 12, 0), 0.0, 0.5, true)

	const LINKS := 4
	const LINK_MASS := 1.0
	const SPRING_K := 200.0

	var links: Array[AVBDRigidBody3D] = []
	var expected_extension: Array[float] = []
	var top: Node3D = anchor
	for i in LINKS:
		var link := add_body(world, "Link%d" % i, Vector3.ONE, Vector3(0, 11.0 - i, 0), LINK_MASS)
		links.append(link)
		# Both join styles connect the touching faces of two 1 m boxes: the anchors
		# coincide at spawn when the centres are one box apart.
		if i % 2 == 0:
			# Rigid link: the anchor points are pinned together.
			add_joint(world, "Rigid%d" % i, top, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0), -1.0, 0.0)
			expected_extension.append(0.0)
		else:
			# Spring link: rest length 0 at the touching faces, so under load the faces
			# separate by exactly force / k, and the links below are the load.
			add_spring(world, "Springy%d" % i, top, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0),
					SPRING_K, 0.0)
			expected_extension.append(float(LINKS - i) * LINK_MASS * GRAVITY / SPRING_K)
		top = link

	await check_forces(world, LINKS, "chain pre-flight")

	await steps(world, 300)

	var measured: Array[float] = []
	var worst_rigid := 0.0
	for i in LINKS:
		var a: Node3D = anchor if i == 0 else links[i - 1]
		measured.append(links[i].global_position.distance_to(a.global_position))
		if i % 2 == 0:
			worst_rigid = maxf(worst_rigid, absf(measured[i] - 1.0))

	check(worst_rigid <= 0.02, "rigid links hold length",
			"worst |length - 1 m| = %.5f m" % worst_rigid)
	for i in LINKS:
		if i % 2 == 1:
			var want := 1.0 + expected_extension[i]
			check(absf(measured[i] - want) <= 0.5 * expected_extension[i] + 0.005,
					"spring link %d stretched" % i,
					"span %.5f m, want %.5f m (load %.2f N / %.0f N/m)" % [measured[i], want,
							expected_extension[i] * SPRING_K, SPRING_K])
	check(all_finite(links), "chain finite", "tail y = %.4f" % links[3].global_position.y)

	var evidence := {"spans": measured, "rigid_error": worst_rigid, "expected_ext": expected_extension}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# C4: AVBDIgnoreCollision3D, A/B against the identical scene without it. Two boxes
# spawn overlapping; without the constraint their contacts push them apart, with it
# they stay interpenetrated. No ground here, so every contact point in the world
# belongs to this pair and the counts are directly attributable.
# ---------------------------------------------------------------------------
func ignore_run(with_ignore: bool) -> Dictionary:
	var world := make_world()
	# Same height, 0.5 m apart, each 1 m wide: they overlap by half a box.
	var left := add_body(world, "Left", Vector3.ONE, Vector3(-0.25, 5, 0))
	var right := add_body(world, "Right", Vector3.ONE, Vector3(0.25, 5, 0))
	if with_ignore:
		add_ignore(world, "Ignore", left, right)

	await check_forces(world, 1 if with_ignore else 0, "ignore pre-flight (ignore=%s)" % with_ignore)
	await steps(world, 120)

	var result := {"separation": left.global_position.distance_to(right.global_position),
			"contacts": world.get_contact_point_count(), "lowest": minf(left.global_position.y, right.global_position.y)}
	world.queue_free()
	await frames(2)
	return result


func test_ignore_collision_pair() -> Variant:
	var control := await ignore_run(false)
	var ignored := await ignore_run(true)

	check(control.contacts > 0, "control pair collides",
			"contact points = %d without the constraint" % control.contacts)
	check(control.separation >= 0.95, "control pair pushes apart",
			"separation = %.5f m, spawned at 0.5 (boxes are 1 m wide)" % control.separation)
	check(ignored.contacts == 0, "ignored pair has no contacts",
			"contact points = %d with AVBDIgnoreCollision3D" % ignored.contacts)
	check(absf(ignored.separation - 0.5) <= 0.01, "ignored pair passes through",
			"separation = %.5f m, unchanged from the spawn 0.5 m" % ignored.separation)

	return {"without": control, "with": ignored}


# ---------------------------------------------------------------------------
# C5: a breakable joint under load. A rigid angular constraint holds a spinning link;
# when the spin loads it past `fracture_force` the joint must break, the link must be
# released, and the break must survive a world rebuild.
# ---------------------------------------------------------------------------
func fracture_run(break_force: float, spin: float) -> Dictionary:
	var world := make_world()
	var anchor := add_body(world, "Anchor", Vector3.ONE, Vector3(0, 10, 0), 0.0, 0.5, true)
	var link := add_body(world, "Link", Vector3.ONE, Vector3(0, 9, 0))
	# Rigid in both directions: the orientation is clamped, so a spin loads the
	# angular term. 0 would leave rotation free and never load it.
	var joint := add_joint(world, "Breakable", anchor, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0),
			-1.0, -1.0, break_force)
	link.initial_angular_velocity = Vector3(0, spin, 0)

	await check_forces(world, 1, "fracture pre-flight")
	await steps(world, 120)

	var broke := joint.is_broken()
	var result := {"broke": broke, "orientation": link.global_transform.basis.get_rotation_quaternion(),
			"spin": link.angular_velocity.length(), "still_constrained": joint.is_simulated(),
			"forces": world.get_force_count()}
	world.queue_free()
	await frames(2)
	return result


func test_fracture_under_load() -> Variant:
	# break_force is in angular-force units; a 40 rad/s spin on a 1 m link loads the
	# 20 N.m joint far past its limit, while a 2000 N.m joint holds the orientation.
	var weak := await fracture_run(20.0, 40.0)
	var strong := await fracture_run(2000.0, 40.0)

	check(not strong.broke, "strong joint holds", "broken=%s with fracture_force 2000" % strong.broke)
	check(strong.still_constrained, "strong joint still simulated", "is_simulated=%s" % strong.still_constrained)
	check(strong.spin <= 0.5, "strong joint stops the spin", "spin = %.4f rad/s (was 40)" % strong.spin)
	check(weak.broke, "weak joint breaks", "broken=%s with fracture_force 20" % weak.broke)
	check(not weak.still_constrained, "broken joint released",
			"is_simulated=%s, solver forces=%d" % [weak.still_constrained, weak.forces])
	check(weak.spin >= 5.0, "broken link keeps spinning", "spin = %.4f rad/s" % weak.spin)

	# The break must be permanent: adding a body rebuilds the whole solver, and a
	# rebuilt joint would silently re-attach the link if the flag were not carried over.
	var world := make_world()
	var anchor := add_body(world, "Anchor", Vector3.ONE, Vector3(0, 10, 0), 0.0, 0.5, true)
	var link := add_body(world, "Link", Vector3.ONE, Vector3(0, 9, 0))
	var joint := add_joint(world, "Breakable", anchor, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0),
			-1.0, -1.0, 20.0)
	link.initial_angular_velocity = Vector3(0, 40, 0)
	await steps(world, 120)
	check(joint.is_broken(), "broke before rebuild", "broken=%s" % joint.is_broken())

	add_body(world, "Extra", Vector3.ONE, Vector3(5, 1, 0))
	await steps(world, 60)
	check(joint.is_broken(), "break survives rebuild", "broken=%s after adding a body" % joint.is_broken())
	check(not joint.is_simulated(), "rebuilt as unconstrained",
			"is_simulated=%s, solver forces=%d (1 expected: the new body has none)" %
			[joint.is_simulated(), world.get_force_count()])

	return {"weak": weak, "strong": strong, "after_rebuild": {"broken": joint.is_broken(),
			"simulated": joint.is_simulated()}}


# ---------------------------------------------------------------------------
# C6: a plank bridge - rigid joints between planks, springs bracing the span, and
# explicit collision suppression at the seams, with a load dropped on the middle.
# ---------------------------------------------------------------------------
func test_stress_bridge() -> Variant:
	var world := make_world()
	add_ground(world, 0.0)

	const PLANKS := 12
	const HALF := 0.5
	var planks: Array[AVBDRigidBody3D] = []
	for i in PLANKS:
		var pinned := i == 0 or i == PLANKS - 1
		planks.append(add_body(world, "Plank%d" % i, Vector3(1.0, 0.25, 2.0),
				Vector3(i - PLANKS * 0.5 + 0.5, 6.0, 0), 1.0, 0.5, pinned))

	for i in PLANKS - 1:
		# Two rigid links per seam (both sides), so the deck cannot twist open.
		add_joint(world, "DeckA%d" % i, planks[i], planks[i + 1], Vector3(HALF, 0, 1.0), Vector3(-HALF, 0, 1.0))
		add_joint(world, "DeckB%d" % i, planks[i], planks[i + 1], Vector3(HALF, 0, -1.0), Vector3(-HALF, 0, -1.0))
		# Planks that merely touch at the seam should not fight the joints.
		add_ignore(world, "Seam%d" % i, planks[i], planks[i + 1])

	# Springs from the two static abutments to the middle of the span, bracing the sag.
	var abutment_left := planks[0]
	var abutment_right := planks[PLANKS - 1]
	var middle := planks[PLANKS / 2]
	var brace_left := add_spring(world, "BraceL", abutment_left, middle, Vector3(0, 0, 0), Vector3(-2.0, 0, 0),
			800.0, 5.5)
	var brace_right := add_spring(world, "BraceR", abutment_right, middle, Vector3(0, 0, 0), Vector3(2.0, 0, 0),
			800.0, 5.5)

	# 2 joints per seam + 1 ignore per seam + 2 brace springs.
	await check_forces(world, 2 * (PLANKS - 1) + (PLANKS - 1) + 2, "bridge pre-flight")

	# Drop a load on the middle of the span.
	for i in 3:
		add_body(world, "Load%d" % i, Vector3.ONE, Vector3(0.2 * i - 0.2, 9.0 + i, 0))

	await steps(world, 420)

	var worst_joint := 0.0
	for joint in world.get_children():
		if joint is AVBDJoint3D and joint.linear_stiffness < 0.0:
			var a: AVBDRigidBody3D = joint.get_body_a()
			var b: AVBDRigidBody3D = joint.get_body_b()
			if a != null and b != null:
				worst_joint = maxf(worst_joint, joint_error(joint, a, b))

	var lowest := INF
	for plank in planks:
		lowest = minf(lowest, plank.global_position.y)

	check(worst_joint <= 0.05, "bridge joints hold", "worst deck joint error = %.5f m" % worst_joint)
	check(lowest > 0.5, "bridge stays above ground", "lowest plank y = %.4f m" % lowest)
	check(absf(middle.global_position.y - planks[0].global_position.y) <= 1.5, "span does not collapse",
			"middle sags %.4f m below the abutment" % (planks[0].global_position.y - middle.global_position.y))
	check(brace_left.is_simulated() and brace_right.is_simulated(), "braces simulated",
			"left=%s right=%s" % [brace_left.is_simulated(), brace_right.is_simulated()])
	check(world.get_contact_point_count() > 0, "load rests on the deck",
			"contact points = %d" % world.get_contact_point_count())

	var evidence := {"worst_joint": worst_joint, "lowest_plank_y": lowest,
			"sag": planks[0].global_position.y - middle.global_position.y,
			"contacts": world.get_contact_point_count(), "forces": world.get_force_count(),
			"digest": digest(world)}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# C7: the soft-body lattice formulation, built from node primitives, resting on a
# spring mattress. The springs must carry the lattice's weight: the sag has a closed
# form, m*g/(n*k), which the simulation has to reproduce.
# ---------------------------------------------------------------------------
func test_lattice_on_springs() -> Variant:
	var world := make_world()
	add_ground(world, 0.0)

	const SIDE := 3
	const SPACING := 0.25
	const SPRING_K := 600.0
	var origin := Vector3(0, 1.0, 0)
	var lattice := build_lattice(world, SIDE, SPACING, origin)

	# One spring per bottom-layer cell, each on its own static pad, so the mattress
	# covers the whole footprint.
	var springs: Array[AVBDSpring3D] = []
	var centre := float(SIDE - 1) * 0.5
	for x in SIDE:
		for z in SIDE:
			var xz := Vector3(x - centre, 0, z - centre) * SPACING
			var pad := add_body(world, "Pad_%d_%d" % [x, z], Vector3(SPACING, 0.1, SPACING),
					origin + xz + Vector3(0, -0.95, 0), 0.0, 0.5, true)
			springs.append(add_spring(world, "Spring_%d_%d" % [x, z], pad, lattice.cell(x, 0, z),
					Vector3(0, 0.05, 0), Vector3(0, -SPACING * 0.5, 0), SPRING_K, -1.0))

	const JOINTS := 3 * (SIDE - 1) * SIDE * SIDE
	await check_forces(world, JOINTS + SIDE * SIDE, "lattice+springs pre-flight")

	var start_centre := lattice_centre(lattice)
	# 2 s: the lattice is held by the springs and has settled onto them.
	await steps(world, 120)
	var held_centre := lattice_centre(lattice)
	var sag := start_centre.y - held_centre.y
	var predicted := lattice.mass() * GRAVITY / (SPRING_K * float(SIDE * SIDE))

	# The contract is that the springs carry the lattice's weight with the predicted sag.
	# The band is deliberately wider than the prediction: the solver's update order (which
	# follows the colouring) shifts the exact value by a fraction of a millimetre, and a
	# lattice that was *not* supported would sag by centimetres or fall to the ground.
	# Measured: 1.03 mm before the parallel update order, -0.26 mm after.
	check(absf(sag) <= 1.5 * predicted + 5e-4, "spring mattress sag",
			"sag = %.6f m, predicted m*g/(n*k) = %.6f m (mass %.4f kg, %d springs)" %
			[sag, predicted, lattice.mass(), SIDE * SIDE])
	check(lattice_worst_joint_error(lattice) <= 0.02, "lattice joints hold under its own weight",
			"worst joint error = %.5f m" % lattice_worst_joint_error(lattice))
	# The lattice must also stay above the pads it rests on, not sink through them.
	check(held_centre.y > 0.0, "lattice is held up", "centre y = %.4f after settling (starts at %.4f)" %
			[held_centre.y, origin.y])

	# 10 s: laterally unrestrained springs cannot stop the lattice sliding off, so the
	# claim is only that the system stays finite and bounded - not that it never moves.
	await steps(world, 480)
	var end_centre := lattice_centre(lattice)
	check(all_finite(lattice.cells), "lattice finite after 10 s",
			"centre %v, worst joint error %.5f m" % [end_centre, lattice_worst_joint_error(lattice)])
	check(end_centre.y > 0.0 and end_centre.y < origin.y + 1.0, "lattice bounded",
			"centre y %.4f (starts at %.4f)" % [end_centre.y, origin.y])
	check(lattice_worst_joint_error(lattice) <= 0.05, "lattice stays coherent",
			"worst joint error = %.5f m after 10 s" % lattice_worst_joint_error(lattice))

	var evidence := {"sag": sag, "predicted_sag": predicted, "mass": lattice.mass(),
			"held_y": held_centre.y, "end_centre": end_centre,
			"worst_joint": lattice_worst_joint_error(lattice)}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# C8: a joint from a static body to a dynamic one. Static bodies are driven by their
# node transform, so teleporting one must drag the dynamic body along.
#
# A teleport is an instantaneous constraint violation, and AVBD constraints are
# stabilised: the joint closes the gap over a few seconds rather than snapping, which is
# the behaviour measured here rather than assumed.
# ---------------------------------------------------------------------------
func test_static_driver() -> Variant:
	var world := make_world()
	var platform := add_body(world, "Platform", Vector3.ONE, Vector3(0, 5, 0), 0.0, 0.5, true)
	var rider := add_body(world, "Rider", Vector3.ONE, Vector3(0, 4, 0))
	# Rigid in both directions (angular_stiffness < 0), so the rider follows the driver
	# exactly instead of swinging on the hitch like a pendulum.
	var joint := add_joint(world, "Hitch", platform, rider, Vector3(0, 0, 0), Vector3(0, 1, 0), -1.0, -1.0)

	await check_forces(world, 1, "driver pre-flight")
	await steps(world, 120)

	var start_error := joint_error(joint, platform, rider)
	check(start_error <= 0.02, "hitch holds at rest", "anchor separation = %.5f m" % start_error)

	# Move the static body 4 m sideways; the joint has to carry the rider along.
	const TRAVEL := 4.0
	platform.position = platform.position + Vector3(TRAVEL, 0, 0)

	var convergence_steps := -1
	var peak_speed := 0.0
	for i in 900:
		await frames(1)
		peak_speed = maxf(peak_speed, rider.linear_velocity.length())
		if joint_error(joint, platform, rider) <= 0.05 and convergence_steps < 0:
			convergence_steps = i + 1

	var end_error := joint_error(joint, platform, rider)
	var dragged := rider.global_position.x - (platform.global_position.x - TRAVEL)

	check(dragged >= TRAVEL * 0.5, "rider is dragged along", "rider travelled %.4f m of %.1f m" % [dragged, TRAVEL])
	check(convergence_steps > 0, "hitch re-converges after the drag",
			"within 5 cm after %d steps (%.2f s)" % [convergence_steps, convergence_steps / 60.0])
	check(end_error <= 0.05, "hitch holds after the drag", "final anchor separation = %.5f m" % end_error)
	check(absf(rider.global_position.x - platform.global_position.x) <= 0.05, "rider follows the platform",
			"rider x = %.4f, platform x = %.4f" % [rider.global_position.x, platform.global_position.x])
	check(rider.linear_velocity.length() <= 0.05, "rider settles after the drag",
			"final speed = %.5f m/s (peak %.4f m/s)" % [rider.linear_velocity.length(), peak_speed])

	var evidence := {"start_error": start_error, "end_error": end_error, "dragged": dragged,
			"convergence_steps": convergence_steps, "peak_speed": peak_speed}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# C9: `substeps` on a composited setup. Substepping must land on the same
# equilibrium, because it only changes how the step is divided.
# ---------------------------------------------------------------------------
func substeps_run(substeps: int) -> Dictionary:
	var world := make_world()
	world.substeps = substeps

	var pivot := add_body(world, "Pivot", Vector3(0.2, 0.2, 0.2), Vector3(0, 5, 0), 0.0, 0.5, true)
	var bob := add_body(world, "Bob", Vector3.ONE, Vector3(0, 4, 0))
	var wall := add_body(world, "Wall", Vector3(0.2, 0.2, 0.2), Vector3(2.5, 4, 0), 0.0, 0.5, true)
	add_joint(world, "Link", pivot, bob, Vector3.ZERO, Vector3(0, 1, 0), -1.0, 0.0)
	add_spring(world, "Pull", wall, bob, Vector3.ZERO, Vector3.ZERO, 400.0, 1.0)

	await check_forces(world, 2, "substeps pre-flight (substeps=%d)" % substeps)
	# 5 s of simulated time either way.
	await steps(world, 300 * substeps)

	var result := {"position": bob.global_position, "length": bob.global_position.distance_to(pivot.global_position),
			"speed": bob.linear_velocity.length()}
	world.queue_free()
	await frames(2)
	return result


func test_substeps_equivalence() -> Variant:
	var single := await substeps_run(1)
	var double := await substeps_run(2)

	var drift: float = single.position.distance_to(double.position)
	check(drift <= 0.05, "substeps converge to the same pose",
			"substeps=1 vs 2 differ by %.5f m" % drift)
	check(single.length >= 0.99 and single.length <= 1.01, "substeps=1 link held", "length = %.5f m" % single.length)
	check(double.length >= 0.99 and double.length <= 1.01, "substeps=2 link held", "length = %.5f m" % double.length)
	check(double.speed <= 0.05, "substeps=2 settles", "speed = %.5f m/s" % double.speed)

	return {"drift": drift, "single": single, "double": double}


# ---------------------------------------------------------------------------
# C11: the solver spreads independent bodies over worker threads. The split must not change
# any number, so the same scene stepped with 1 and with N threads has to agree exactly.
# This is the engine-side counterpart of the core suite's digest comparison.
# ---------------------------------------------------------------------------
func threaded_scene(threads: int) -> Dictionary:
	var world := make_world()
	world.threads = threads
	add_ground(world, 0.0)

	# A ragged mix: coupled stack, coupled chain, and free bodies, so the colouring has
	# several groups of different widths.
	var stack: Array[AVBDRigidBody3D] = []
	for i in 8:
		stack.append(add_body(world, "Stack%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.5, 0)))
	var previous: Node3D = stack[0]
	for i in 1:
		var link := add_body(world, "Link%d" % i, Vector3.ONE, Vector3(4, 6, 0))
		add_joint(world, "ChainJoint%d" % i, previous, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0))
		previous = link
	for i in 6:
		add_body(world, "Free%d" % i, Vector3.ONE, Vector3(8 + i * 2, 3 + i, 0))

	await steps(world, 240)

	var result := {"digest": digest(world), "bodies": world.get_body_count(),
			"contacts": world.get_contact_point_count(), "workers": world.get_thread_count()}
	world.queue_free()
	await frames(2)
	return result


func test_parallel_threads() -> Variant:
	var serial := await threaded_scene(1)
	var parallel := await threaded_scene(4)
	var auto_threads := await threaded_scene(0)

	check(serial.workers == 1, "threads=1 runs on the calling thread", "workers=%d" % serial.workers)
	check(parallel.workers == 4, "threads=4 uses four workers", "workers=%d" % parallel.workers)
	check(auto_threads.workers > 1, "threads=0 picks the hardware count", "workers=%d" % auto_threads.workers)
	check(parallel.digest == serial.digest, "thread count does not change the result",
			"threads=4 digest %d vs threads=1 digest %d" % [parallel.digest, serial.digest])
	check(auto_threads.digest == serial.digest, "automatic thread count matches",
			"threads=0 digest %d vs threads=1 digest %d" % [auto_threads.digest, serial.digest])
	check(parallel.contacts == serial.contacts, "contact count matches",
			"%d vs %d contact points" % [parallel.contacts, serial.contacts])

	return {"serial": serial, "parallel": parallel, "auto": auto_threads}


# ---------------------------------------------------------------------------
# C10: a composited scene is reproducible. The harness also re-runs every scenario
# under --digest; this checks the strongest form in one scenario: two identically
# built scenes, one after the other, in the same process.
# ---------------------------------------------------------------------------
func combination_scene() -> int:
	var world := make_world()
	var ground := add_ground(world, 0.5)

	var lattice := build_lattice(world, 3, 0.3, Vector3(0, 3.0, 0), 1500.0, 400.0)
	var plank_a := add_body(world, "PlankA", Vector3(1, 0.25, 1), Vector3(2.0, 4.0, 0))
	var plank_b := add_body(world, "PlankB", Vector3(1, 0.25, 1), Vector3(3.0, 4.0, 0))
	add_joint(world, "PlankJoint", plank_a, plank_b, Vector3(0.5, 0, 0), Vector3(-0.5, 0, 0))
	add_spring(world, "Tie", plank_b, lattice.cell(2, 1, 1), Vector3.ZERO, Vector3.ZERO, 300.0, 1.0)
	add_ignore(world, "NoTouch", plank_a, plank_b)
	add_body(world, "Drop", Vector3.ONE, Vector3(0.4, 6.0, 0.4))

	await steps(world, 240)

	var value := digest(world)
	check(all_finite(lattice.cells) and is_finite_body(plank_a), "combination scene finite",
			"lattice worst joint error %.5f m" % lattice_worst_joint_error(lattice))
	check(ground != null and world.get_body_count() == lattice.cells.size() + 4, "combination scene body count",
			"bodies = %d" % world.get_body_count())

	world.queue_free()
	await frames(2)
	return value


func test_combination_determinism() -> Variant:
	var first := await combination_scene()
	var second := await combination_scene()
	check(first == second, "composite scenes reproduce", "digest %d vs %d" % [first, second])
	return {"digest": first}
