# Diagnostic probe: a 4-link chain of fully-locked Generic6DOFJoint3D welds.
# Prints per-tick state for the first ticks so the runaway can be attributed
# (joint force vs. velocity vs. position) instead of guessed at.
#
#   Godot --headless --path demo --script res://tests/probe_chain.gd
extends "res://tests/avbd_harness.gd"

const TICKS := 1300


func _scenarios() -> Array:
	return [scenario("chain", probe_chain)]


func make_locked_6dof(parent: Node, joint_name: String, a: Node3D, b: Node3D, anchor: Vector3) -> Generic6DOFJoint3D:
	var joint := Generic6DOFJoint3D.new()
	joint.name = joint_name
	# Local transform before entering the tree: `root` is the scene root, so local == global
	# here. Assigning `global_position` after add_child() errors out while the root is not
	# yet inside the tree (that is what silently put the old joints at the origin).
	joint.position = anchor
	parent.add_child(joint)
	joint.node_a = joint.get_path_to(a)
	joint.node_b = joint.get_path_to(b)
	return joint


func probe_chain() -> Variant:
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

	print("anchor frozen=%s mass=%.3f  link mass=%.3f" % [anchor.freeze, anchor.mass, links[0].mass])
	print("joint0 global=%v node_a=%s node_b=%s" % [joints[0].global_position,
			joints[0].node_a, joints[0].node_b])

	for tick in TICKS:
		await physics_frame
		if tick < 6 or tick % 100 == 0:
			var line := "t=%4d" % tick
			for i in links.size():
				line += " | L%d y=%9.2f vy=%9.1f s=%d" % [i, links[i].global_position.y,
						links[i].linear_velocity.y, int(links[i].sleeping)]
			print(line)
			print("      anchor sleeping=%d" % int(anchor.sleeping))
			var worst := 0.0
			for i in joints.size():
				var a: Node3D = anchor if i == 0 else links[i - 1]
				worst = maxf(worst, anchor_error(a, links[i], Vector3(0, -0.5, 0), Vector3(0, 0.5, 0)))
			print("      worst anchor separation = %.6f" % worst)

	for link in links:
		link.queue_free()
	for joint in joints:
		joint.queue_free()
	anchor.queue_free()
	await frames(2)
	return {}
