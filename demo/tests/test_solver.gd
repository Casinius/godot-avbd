# Behavioural tests for the AVBD extension, run against a real Godot runtime.
#
#   Godot_v4.6.3 --headless --path demo --script res://tests/test_solver.gd
#   ... -- --list | --only=stack | --digest | --strict
#
# Every scenario builds its node graph in code, then steps the real physics loop, so
# this exercises the full path: scene graph -> solver -> transforms written back to the
# scene graph. Single-constraint behaviour lives here; constraint *combinations* are in
# test_constraints.gd.
extends "res://tests/avbd_harness.gd"


func _scenarios() -> Array:
	return [
		scenario("rest_contact", test_rest_contact),
		scenario("stack", test_stack),
		scenario("friction", test_friction),
		scenario("joints", test_joints),
		scenario("world_joint", test_world_joint),
		scenario("soft_body", test_soft_body),
		scenario("raycast", test_raycast),
		scenario("runtime_add", test_runtime_add),
		scenario("interaction", test_interaction),
		scenario("determinism", test_determinism),
	]


# ---------------------------------------------------------------------------
# A box dropped on a slab comes to rest on top of it: this is the end-to-end
# check that the Godot Y-up <-> solver Z-up mapping is right.
# ---------------------------------------------------------------------------
func test_rest_contact() -> Variant:
	var world := make_world()
	add_ground(world, 0.5)
	var box := add_body(world, "Box", Vector3.ONE, Vector3(0, 5, 0))

	await check_forces(world, 0, "no constraints")

	await steps(world, 240)

	var rest_y := box.global_position.y
	check(absf(rest_y - 1.0) <= 0.02, "rest contact", "y=%.5f (slab top 0.5 + half box 0.5)" % rest_y)
	check(absf(box.linear_velocity.y) <= 0.05, "rest velocity", "vy=%.5f" % box.linear_velocity.y)
	check(world.get_contact_point_count() > 0, "contact reported",
			"points=%d manifests=%d" % [world.get_contact_point_count(), world.get_contact_count()])
	check(world.get_body_count() == 2, "body count", "bodies=%d" % world.get_body_count())
	check(is_finite_body(box), "state finite", "pos=%v" % box.global_position)

	var evidence := {"rest_y": rest_y, "contacts": world.get_contact_point_count()}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A 10-box stack keeps its shape and does not drift sideways.
# ---------------------------------------------------------------------------
func test_stack() -> Variant:
	var world := make_world()
	add_ground(world, 0.5)

	var boxes: Array[AVBDRigidBody3D] = []
	for i in 10:
		boxes.append(add_body(world, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.5, 0)))

	await steps(world, 300)

	var worst_y := 0.0
	var worst_interface := 0.0
	var worst_lateral := 0.0
	for i in boxes.size():
		var pos := boxes[i].global_position
		worst_y = maxf(worst_y, absf(pos.y - (1.0 + i)))
		# Ideal: the bottom box's centre is 1.0 above the world origin (slab top
		# 0.5 + half box 0.5) and each box above sits one box-height higher.
		var penetration := 1.0 - pos.y if i == 0 else 1.0 - (pos.y - boxes[i - 1].global_position.y)
		worst_interface = maxf(worst_interface, penetration)
		worst_lateral = maxf(worst_lateral, Vector2(pos.x, pos.z).length())

	check(worst_y <= 0.2, "stack height", "max |y - expected| = %.5f" % worst_y)
	check(worst_interface <= 0.05, "stack penetration", "worst interface = %.5f m" % worst_interface)
	check(worst_lateral <= 0.05, "stack drift", "max lateral = %.5f m" % worst_lateral)

	var evidence := {"height_error": worst_y, "penetration": worst_interface, "drift": worst_lateral}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Friction on a 20 degree ramp: mu = 0.5 holds (penalty-method creep only),
# mu = 0.05 lets the box slide away.
# ---------------------------------------------------------------------------
func ramp_run(mu: float, settle_frames: int, watch_frames: int) -> Dictionary:
	var world := make_world()
	add_ground(world, 0.5, mu)

	var ramp_basis := Basis(Vector3(0, 0, 1), deg_to_rad(20.0))
	var ramp := add_body(world, "Ramp", Vector3(40, 1, 24), Vector3(0, 6, 0), 0.0, mu, true)
	ramp.transform = Transform3D(ramp_basis, ramp.position)

	var up_slope := ramp_basis * Vector3.RIGHT
	var normal := ramp_basis * Vector3.UP
	var start := ramp.global_position + up_slope * 5.0 + normal * 1.05
	var box := add_body(world, "Box", Vector3.ONE, start, 1.0, mu)

	await steps(world, settle_frames)
	var before := box.global_position
	await steps(world, watch_frames)

	var late := (box.global_position - before).dot(-up_slope)
	var total := (box.global_position - start).dot(-up_slope)
	var speed := box.linear_velocity.length()
	var evidence := {"late": late, "total": total, "speed": speed}
	world.queue_free()
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
# A 4-link chain hung from a static body: hard joints must hold to a fraction of
# a link, and the chain must hang downwards.
# ---------------------------------------------------------------------------
func test_joints() -> Variant:
	var world := make_world()
	var anchor := add_body(world, "Anchor", Vector3.ONE, Vector3(0, 10, 0), 0.0, 0.5, true)

	var links: Array[AVBDRigidBody3D] = []
	var joints: Array[AVBDJoint3D] = []
	var previous: Node3D = anchor
	for i in 4:
		var link := add_body(world, "Link%d" % i, Vector3.ONE, Vector3(0, 9.0 - i, 0))
		joints.append(add_joint(world, "Joint%d" % i, previous, link, Vector3(0, -0.5, 0), Vector3(0, 0.5, 0)))
		links.append(link)
		previous = link

	await check_forces(world, 4, "joint pre-flight")

	await steps(world, 300)

	var worst := 0.0
	for i in joints.size():
		var a: Node3D = anchor if i == 0 else links[i - 1]
		worst = maxf(worst, joint_error(joints[i], a, links[i]))

	check(worst <= 0.01, "hard joint error", "max anchor separation = %.6f m" % worst)
	check(links[3].global_position.y < anchor.global_position.y, "chain hangs",
			"tail y=%.4f below anchor y=%.4f" % [links[3].global_position.y, anchor.global_position.y])
	check(world.get_force_count() == 4, "joint count", "forces=%d" % world.get_force_count())

	var evidence := {"worst_joint_error": worst, "tail_y": links[3].global_position.y}
	world.queue_free()
	await frames(2)
	return evidence


# A joint anchored to world space (node_a empty) pins a body to a point. The anchor
# is deliberately off-axis so a swapped coordinate would show up as a mismatch.
#
# AVBD constraints are stabilised, so an initial violation decays over a second or
# two instead of snapping; the body is therefore spawned at its anchor.
func test_world_joint() -> Variant:
	var world := make_world()
	var anchor := Vector3(2, 3, 1)
	var bob := add_body(world, "Bob", Vector3.ONE, anchor)
	bob.initial_velocity = Vector3(4, 0, 0)

	var joint := add_world_joint(world, "Pin", bob, anchor)
	await check_forces(world, 1, "joint pre-flight")

	await steps(world, 180)

	var error := bob.global_position.distance_to(anchor)
	check(error <= 0.02, "world anchor holds",
			"bob at %v, anchor %v, error %.5f m" % [bob.global_position, anchor, error])
	check(bob.linear_velocity.length() <= 0.05, "world anchor absorbs velocity",
			"speed=%.5f m/s from an initial 4 m/s" % bob.linear_velocity.length())

	var evidence := {"error": error, "joint_name": joint.name, "speed": bob.linear_velocity.length()}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Lattice soft body: a 3x3x3 blob of jointed boxes dropped on the ground.
# ---------------------------------------------------------------------------
func test_soft_body() -> Variant:
	var world := make_world()
	add_ground(world, 0.0)

	var soft := add_soft_body(world, "Soft", Vector3(0, 4, 0))

	await steps(world, 300)

	var count := soft.get_cell_count()
	var finite := true
	var lowest := INF
	var highest := -INF
	for i in count:
		var t := soft.get_cell_transform(i)
		if not t.origin.is_finite():
			finite = false
		var y := soft.to_global(t.origin).y
		lowest = minf(lowest, y)
		highest = maxf(highest, y)

	check(count == 27, "soft body cells", "cells=%d" % count)
	check(finite, "soft body finite", "all %d cell transforms finite" % count)
	check(lowest > -0.2 and highest < 4.5, "soft body lands", "cell y in [%.3f, %.3f]" % [lowest, highest])
	check(world.get_body_count() == 28, "soft body in solver", "bodies=%d (27 cells + ground)" % world.get_body_count())

	# The drawn instances must match the simulated cells, and the lattice must keep
	# its nominal spacing (0.25 m) after landing. Under --headless the dummy renderer
	# discards MultiMesh buffers (reads come back as identity), so the instance
	# comparison only runs against a real rendering device.
	var mm: MultiMesh = soft.get_multimesh()
	check(mm != null and mm.get_instance_count() == count, "soft body instances",
			"multimesh instances=%d" % [mm.get_instance_count() if mm != null else -1])
	if DisplayServer.get_name() == "headless":
		print("     (instance transforms not comparable under the dummy renderer)")
	else:
		var mismatch := 0.0
		for i in count:
			mismatch = maxf(mismatch, mm.get_instance_transform(i).origin.distance_to(soft.get_cell_transform(i).origin))
		check(mismatch <= 1e-5, "visuals match solver", "max instance/cell mismatch = %.7f m" % mismatch)
	var spacing := soft.get_cell_transform(0).origin.distance_to(soft.get_cell_transform(1).origin)
	check(spacing > 0.2 and spacing < 0.4, "lattice spacing", "cell 0-1 spacing = %.4f m (0.25 nominal)" % spacing)
	# The lattice must fall straight down onto the slab, not drift sideways: the
	# joint anchors are expressed in the solver's Z-up frame, where a sign error
	# would push the lattice along the horizontal axes.
	var worst_lateral := 0.0
	for i in count:
		var cell := soft.to_global(soft.get_cell_transform(i).origin)
		worst_lateral = maxf(worst_lateral, Vector2(cell.x, cell.z).length())
	check(worst_lateral <= 0.6, "lattice stays centred", "max lateral offset = %.4f m" % worst_lateral)

	var evidence := {"cells": count, "lowest_y": lowest, "highest_y": highest, "spacing": spacing,
			"lateral": worst_lateral, "digest": digest(world)}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Raycasts hit dynamic bodies, static bodies and soft bodies.
# ---------------------------------------------------------------------------
func test_raycast() -> Variant:
	var world := make_world()
	var ground := add_ground(world, 0.5)
	var box := add_body(world, "Box", Vector3.ONE, Vector3(0, 3, 0))

	await steps(world, 60)

	var hit := world.raycast(Vector3(0, 20, 0), Vector3(0, -1, 0))
	check(hit.get("body") == box, "raycast hits dynamic body", "body=%s" % [hit.get("body")])
	if hit.has("position"):
		# The box has settled on the slab, so the ray meets its top face.
		var top: float = box.global_position.y + 0.5
		check(absf(hit.position.y - top) <= 0.05, "raycast position", "y=%.4f (box top %.4f)" % [hit.position.y, top])
		check(absf(hit.distance - (20.0 - top)) <= 0.05, "raycast distance",
				"d=%.4f (want %.4f)" % [hit.distance, 20.0 - top])

	var miss := world.raycast(Vector3(50, 20, 50), Vector3(0, 1, 0))
	check(miss.is_empty(), "raycast misses", "keys=%s" % [miss.keys()])

	var ground_hit := world.raycast(Vector3(20, 20, 20), Vector3(0, -1, 0))
	check(ground_hit.get("body") == ground, "raycast hits static body", "body=%s" % [ground_hit.get("body")])

	# Static bodies follow their node: move the slab and the ray must follow.
	ground.position = Vector3(0, 10.5, 0)
	await frames(2)
	var moved := world.raycast(Vector3(20, 20, 20), Vector3(0, -1, 0))
	check(moved.has("position") and absf(moved.position.y - 11.0) <= 0.05, "static body teleport adopted",
			"hit y=%.4f (moved slab top 11.0)" % [moved.get("position", Vector3.ZERO).y])

	var evidence := {"dynamic_hit": hit.get("body").name if hit.get("body") != null else "<none>",
			"distance": hit.get("distance", -1.0)}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Adding a body at runtime must not reset the bodies already simulated.
# ---------------------------------------------------------------------------
func test_runtime_add() -> Variant:
	var world := make_world()
	add_ground(world, 0.5)
	var box := add_body(world, "Box", Vector3.ONE, Vector3(0, 3, 0))

	await steps(world, 120)
	var settled := box.global_position.y

	var late := add_body(world, "Late", Vector3.ONE, Vector3(5, 2, 0))
	await steps(world, 60)

	check(absf(settled - box.global_position.y) <= 0.02, "rebuild keeps pose",
			"settled y %.4f -> %.4f after adding a body" % [settled, box.global_position.y])
	check(world.get_body_count() == 3, "added body simulated", "bodies=%d" % world.get_body_count())
	check(late.global_position.y > 0.5, "added body rests", "late y=%.4f" % late.global_position.y)

	var evidence := {"settled": settled, "after": box.global_position.y, "late_y": late.global_position.y}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Impulses and teleports act on the solver state with the right magnitudes.
# ---------------------------------------------------------------------------
func test_interaction() -> Variant:
	var world := make_world()
	add_ground(world, 0.5)
	var box := add_body(world, "Box", Vector3.ONE, Vector3(0, 1.0, 0))

	await steps(world, 60)

	# A 1 m^3 box at density 1 weighs 1 kg, so a 10 N.s impulse is 10 m/s.
	box.apply_impulse(Vector3(10, 0, 0))
	check(absf(box.linear_velocity.x - 10.0) <= 0.01, "impulse sets velocity",
			"vx=%.5f (10 N.s / 1 kg)" % box.linear_velocity.x)

	# The same impulse at the top face also spins the box: I = m*(1+1)/12 = 1/6,
	# so w = r x J / I = (0,0,-5) / (1/6) = -30 rad/s about Z.
	box.set_linear_velocity(Vector3.ZERO)
	box.set_angular_velocity(Vector3.ZERO)
	box.apply_impulse(Vector3(10, 0, 0), Vector3(0, 0.5, 0))
	var spin := box.angular_velocity.z
	check(absf(spin + 30.0) <= 0.1, "impulse spins body", "wz=%.4f (want -30 rad/s)" % spin)

	await steps(world, 60)
	check(box.global_position.x > 1.0, "impulse moves body", "x=%.4f" % box.global_position.x)

	box.teleport(Vector3(-5, 4, 0))
	check(absf(box.global_position.x + 5.0) <= 0.001, "teleport moves node", "x=%.4f" % box.global_position.x)
	await frames(2)
	check(absf(box.global_position.x + 5.0) <= 0.05, "teleport holds", "x=%.4f on the next tick" % box.global_position.x)

	var evidence := {"impulse_v": 10.0, "spin_w": spin, "teleport_x": box.global_position.x}
	world.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Same scene, same ticks, same state.
# ---------------------------------------------------------------------------
func test_determinism() -> Variant:
	var first := await stack_digest()
	var second := await stack_digest()
	check(first == second, "in-process determinism", "digest %d vs %d" % [first, second])
	return {"digest": first}


func stack_digest() -> int:
	var world := make_world()
	add_ground(world, 0.5)
	var boxes: Array[AVBDRigidBody3D] = []
	for i in 10:
		boxes.append(add_body(world, "Box%d" % i, Vector3.ONE, Vector3(0, 1.0 + i * 1.5, 0)))

	await steps(world, 300)

	var hash := digest(world)
	world.queue_free()
	await frames(2)
	return hash
