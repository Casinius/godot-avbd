extends "res://tests/avbd_harness.gd"

const SCENE := "res://scenes/pyramid_dbl.tscn"

func _scenarios() -> Array:
	return [scenario("probe", run_probe)]

func _collect(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect(child, out)

func run_probe() -> Variant:
	var scene := (load(SCENE) as PackedScene).instantiate()
	scene.set_script(null)
	root.add_child(scene)
	var bodies: Array[RigidBody3D] = []
	_collect(scene, bodies)
	for body in bodies:
		body.mass = 1.0
		body.sleeping = false
		body.can_sleep = false

	var first_escape := -1
	var escape_info := ""
	for step in 1800:
		await steps(1)
		for i in bodies.size():
			var b := bodies[i]
			var p := b.global_position
			if first_escape < 0 and p.y < -1.0:
				first_escape = step
				escape_info = "step=%d body=%s pos=(%.3f,%.3f,%.3f) vel=(%.3f,%.3f,%.3f)" % [
						step, b.name, p.x, p.y, p.z,
						b.linear_velocity.x, b.linear_velocity.y, b.linear_velocity.z]
		if step % 240 == 239:
			var max_abs_z := 0.0
			var max_abs_x := 0.0
			var max_y := -1.0e30
			var min_y := 1.0e30
			var max_v := 0.0
			var max_v_name := ""
			for b in bodies:
				var p := b.global_position
				max_abs_z = maxf(max_abs_z, absf(p.z))
				max_abs_x = maxf(max_abs_x, absf(p.x))
				max_y = maxf(max_y, p.y)
				min_y = minf(min_y, p.y)
				if b.linear_velocity.length() > max_v:
					max_v = b.linear_velocity.length()
					max_v_name = str(b.name)
			print("t=%4d max|x|=%.2f max|z|=%.2f y=[%.2f,%.2f] max_v=%.2f (%s)" % [
					step + 1, max_abs_x, max_abs_z, min_y, max_y, max_v, max_v_name])
	print("first escape: %s" % (escape_info if first_escape >= 0 else "none"))
	return {}
