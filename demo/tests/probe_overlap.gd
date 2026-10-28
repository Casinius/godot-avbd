extends SceneTree

# Cross-branch probe: measures geometric OBB penetration in the double-pyramid scene
# without relying on any engine-specific diagnostic, so the same numbers can be taken
# on different solver branches.
#
#   Godot --headless --path demo --script res://tests/probe_overlap.gd -- --seconds=10

const SCENE := "res://scenes/pyramid_dbl.tscn"
const HALF := Vector3(0.5, 0.25, 0.5)

var _seconds := 10.0
var _scene := SCENE
var _upper_mass := 1.0

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--seconds="):
			_seconds = float(arg.split("=")[1])
		elif arg.begins_with("--scene="):
			_scene = arg.split("=")[1]
		elif arg.begins_with("--internal="):
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", float(arg.split("=")[1]))
		elif arg.begins_with("--substeps="):
			ProjectSettings.set_setting("physics/avbd/substeps", int(arg.split("=")[1]))
		elif arg.begins_with("--iterations="):
			ProjectSettings.set_setting("physics/avbd/iterations", int(arg.split("=")[1]))
		elif arg.begins_with("--adaptive="):
			ProjectSettings.set_setting("physics/avbd/adaptive_beta", arg.split("=")[1] == "1")
		elif arg.begins_with("--beta="):
			ProjectSettings.set_setting("physics/avbd/beta_linear", float(arg.split("=")[1]))
		elif arg.begins_with("--alpha="):
			ProjectSettings.set_setting("physics/avbd/alpha", float(arg.split("=")[1]))
		elif arg.begins_with("--proj="):
			ProjectSettings.set_setting("physics/avbd/position_projection", arg.split("=")[1] == "1")
		elif arg.begins_with("--proj_iters="):
			ProjectSettings.set_setting("physics/avbd/projection_iterations", int(arg.split("=")[1]))
		elif arg.begins_with("--proj_beta="):
			ProjectSettings.set_setting("physics/avbd/projection_beta", float(arg.split("=")[1]))
		elif arg.begins_with("--proj_cap="):
			ProjectSettings.set_setting("physics/avbd/projection_max_correction", float(arg.split("=")[1]))
		elif arg.begins_with("--proj_slop="):
			ProjectSettings.set_setting("physics/avbd/projection_slop", float(arg.split("=")[1]))
		elif arg.begins_with("--spec="):
			ProjectSettings.set_setting("physics/avbd/speculative_margin", float(arg.split("=")[1]))
		elif arg.begins_with("--upper_mass="):
			_upper_mass = float(arg.split("=")[1])
		elif arg.begins_with("--newton="):
			ProjectSettings.set_setting("physics/avbd/newton_ratio", float(arg.split("=")[1]))
		elif arg.begins_with("--decay="):
			ProjectSettings.set_setting("physics/avbd/stiffness_decay", float(arg.split("=")[1]))
		elif arg.begins_with("--gravity="):
			ProjectSettings.set_setting("physics/3d/default_gravity", float(arg.split("=")[1]))
		elif arg.begins_with("--threads="):
			ProjectSettings.set_setting("physics/avbd/threads", int(arg.split("=")[1]))
	_run()

func _run() -> void:
	var t0 := Time.get_ticks_msec()
	var scene: Node = (load(_scene) as PackedScene).instantiate()
	scene.set_script(null)
	root.add_child(scene)
	var bodies: Array[RigidBody3D] = []
	_collect(scene, bodies)
	for body in bodies:
		body.mass = 1.0 if int(body.name.trim_prefix("Box")) < 210 else _upper_mass
		body.sleeping = false
		body.can_sleep = false

	var hz := float(Engine.physics_ticks_per_second)
	var ticks := int(round(_seconds * hz))
	var samples: Array[float] = []
	var min_y := 1.0e9
	var max_speed := 0.0
	var max_speed_body := ""
	var max_speed_y := 0.0
	var max_speed_x := 0.0
	var max_speed_z := 0.0
	var pairs_over_10cm := 0
	var last_worst := ""
	for tick in ticks:
		await physics_frame
		if tick % 60 != 59:
			continue
		var worst := 0.0
		var worst_pair := ""
		var worst_a: RigidBody3D = null
		var worst_b: RigidBody3D = null
		for i in bodies.size():
			var a := bodies[i]
			min_y = minf(min_y, a.global_position.y)
			if a.linear_velocity.length() > max_speed:
				max_speed = a.linear_velocity.length()
				max_speed_body = a.name
				max_speed_y = a.global_position.y
				max_speed_x = a.global_position.x
				max_speed_z = a.global_position.z
			for j in range(i + 1, bodies.size()):
				var b := bodies[j]
				if a.global_position.distance_squared_to(b.global_position) > 4.0:
					continue
				var d := _sat_depth(a.global_transform, b.global_transform)
				if d > worst:
					worst = d
					worst_pair = "%s/%s" % [a.name, b.name]
					worst_a = a
					worst_b = b
		if worst > 0.2 and worst_a != null:
			print("  DEEP %.2fs %s depth=%.3f A=(%.3f,%.3f,%.3f) v=%.2f B=(%.3f,%.3f,%.3f) v=%.2f" % [
				(tick + 1) / hz, worst_pair, worst,
				worst_a.global_position.x, worst_a.global_position.y, worst_a.global_position.z,
				worst_a.linear_velocity.length(),
				worst_b.global_position.x, worst_b.global_position.y, worst_b.global_position.z,
				worst_b.linear_velocity.length()])
		pairs_over_10cm = pairs_over_10cm if worst <= 0.1 else pairs_over_10cm
		samples.append(worst)
		if worst > 0.1:
			last_worst = worst_pair
		if tick % 30 == 29:
			var watch := ""
			for w in ["Box209", "Box219", "Box419"]:
				var node := scene.get_node_or_null(w)
				if node != null:
					watch += " %s=(%.2f,%.2f,%.2f)" % [w, node.global_position.x, node.global_position.y, node.global_position.z]
			print("t=%.2fs worst=%.3f (%s) min_y=%.3f max_v=%.2f (%s @y=%.2f x=%.2f z=%.2f)%s" % [
				(tick + 1) / hz, worst, worst_pair, min_y, max_speed, max_speed_body, max_speed_y,
				max_speed_x, max_speed_z, watch])
	samples.sort()
	var n := samples.size()
	var p50: float = samples[n / 2] if n > 0 else 0.0
	var p95: float = samples[int(n * 0.95)] if n > 0 else 0.0
	var over := 0
	for s in samples:
		if s > 0.05:
			over += 1
	print("RESULT cfg=%s p50=%.4f p95=%.4f max=%.4f over5cm=%d/%d min_y=%.3f max_v=%.2f worst_pair=%s wall=%.1fs" % [
		OS.get_cmdline_user_args(), p50, p95, samples[n - 1] if n > 0 else 0.0, over, n,
		min_y, max_speed, last_worst, (Time.get_ticks_msec() - t0) / 1000.0])
	quit()

func _collect(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect(child, out)

# Separating-axis depth of two equally sized OBBs (positive = interpenetration).
func _sat_depth(ta: Transform3D, tb: Transform3D) -> float:
	var ax: Array[Vector3] = [ta.basis.x.normalized(), ta.basis.y.normalized(), ta.basis.z.normalized()]
	var bx: Array[Vector3] = [tb.basis.x.normalized(), tb.basis.y.normalized(), tb.basis.z.normalized()]
	var delta := tb.origin - ta.origin
	var best := 1.0e9
	# Face axes.
	for i in 3:
		var ra: float = HALF[i]
		var rb := 0.0
		for j in 3:
			rb += HALF[j] * absf(ax[i].dot(bx[j]))
		var d := absf(delta.dot(ax[i])) - (ra + rb)
		if d > 0.0:
			return 0.0
		best = minf(best, -d)
	for j in 3:
		var rb2: float = HALF[j]
		var ra2 := 0.0
		for i in 3:
			ra2 += HALF[i] * absf(ax[i].dot(bx[j]))
		var d2 := absf(delta.dot(bx[j])) - (ra2 + rb2)
		if d2 > 0.0:
			return 0.0
		best = minf(best, -d2)
	# Edge axes.
	for i in 3:
		for j in 3:
			var axis := ax[i].cross(bx[j])
			var len2 := axis.length_squared()
			if len2 < 1.0e-8:
				continue
			axis /= sqrt(len2)
			var r1 := 0.0
			var r2 := 0.0
			for k in 3:
				r1 += HALF[k] * absf(ax[k].dot(axis))
				r2 += HALF[k] * absf(bx[k].dot(axis))
			var d3 := absf(delta.dot(axis)) - (r1 + r2)
			if d3 > 0.0:
				return 0.0
			best = minf(best, -d3)
	return best
