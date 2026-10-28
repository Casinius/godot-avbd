extends SceneTree

# Loads pyramid_dbl and, when --shift= is given, lowers the upper pyramid so it
# starts in contact instead of falling onto the lower one. Zero-initial-velocity
# settle test: measures whether the pile can *hold* the load, separated from the
# impact transient.

const SCENE := "res://scenes/pyramid_dbl.tscn"
const HALF := Vector3(0.5, 0.25, 0.5)

var _seconds := 5.0
var _shift := 0.0

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--seconds="):
			_seconds = float(arg.split("=")[1])
		elif arg.begins_with("--shift="):
			_shift = float(arg.split("=")[1])
		elif arg.begins_with("--internal="):
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", float(arg.split("=")[1]))
		elif arg.begins_with("--substeps="):
			ProjectSettings.set_setting("physics/avbd/substeps", int(arg.split("=")[1]))
		elif arg.begins_with("--iterations="):
			ProjectSettings.set_setting("physics/avbd/iterations", int(arg.split("=")[1]))
	_run()

func _run() -> void:
	var scene: Node = (load(SCENE) as PackedScene).instantiate()
	scene.set_script(null)
	root.add_child(scene)
	var bodies: Array[RigidBody3D] = []
	_collect(scene, bodies)
	bodies.sort_custom(func(a, b): return _idx(a.name) < _idx(b.name))
	for i in bodies.size():
		var body := bodies[i]
		body.mass = 1.0
		body.sleeping = false
		body.can_sleep = false
		if _idx(body.name) >= 210 and _shift > 0.0:
			body.global_position.y -= _shift

	var hz := float(Engine.physics_ticks_per_second)
	var ticks := int(round(_seconds * hz))
	var samples: Array[float] = []
	var min_y := 1.0e9
	var max_speed := 0.0
	var worst_pair := ""
	for tick in ticks:
		await physics_frame
		if tick % 60 != 59:
			continue
		var worst := 0.0
		var pair := ""
		for i in bodies.size():
			var a := bodies[i]
			min_y = minf(min_y, a.global_position.y)
			max_speed = maxf(max_speed, a.linear_velocity.length())
			for j in range(i + 1, bodies.size()):
				var b := bodies[j]
				if a.global_position.distance_squared_to(b.global_position) > 4.0:
					continue
				var d := _sat_depth(a.global_transform, b.global_transform)
				if d > worst:
					worst = d
					pair = "%s/%s" % [a.name, b.name]
		samples.append(worst)
		worst_pair = pair
		if tick % 240 == 239:
			print("t=%.2fs worst=%.3f (%s) min_y=%.3f max_v=%.2f" % [
				(tick + 1) / hz, worst, pair, min_y, max_speed])
	samples.sort()
	var n := samples.size()
	print("RESULT cfg=%s p50=%.4f p95=%.4f max=%.4f min_y=%.3f max_v=%.2f worst_pair=%s" % [
		OS.get_cmdline_user_args(), samples[n / 2], samples[int(n * 0.95)],
		samples[n - 1], min_y, max_speed, worst_pair])
	quit()

func _idx(name: String) -> int:
	return int(name.trim_prefix("Box"))

func _collect(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect(child, out)

func _sat_depth(ta: Transform3D, tb: Transform3D) -> float:
	var ax: Array[Vector3] = [ta.basis.x.normalized(), ta.basis.y.normalized(), ta.basis.z.normalized()]
	var bx: Array[Vector3] = [tb.basis.x.normalized(), tb.basis.y.normalized(), tb.basis.z.normalized()]
	var delta := tb.origin - ta.origin
	var best := 1.0e9
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
