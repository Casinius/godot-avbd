extends SceneTree

# Controlled impact test: a vertical column of N boxes (2 cm gaps) dropped from
# --height= onto a static ground. Measures the deepest box-box overlap per frame
# so the two solver branches can be compared on a minimal scene.
#
#   Godot --headless --path demo --script res://tests/probe_stack_drop.gd -- --n=5 --height=1
#
# The boxes are spawned already overlapping the drop distance only in the sense
# that the whole column starts higher; their spacing matches the pyramid scenes.

var _n := 5
var _height := 1.0
var _seconds := 4.0

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--n="):
			_n = int(arg.split("=")[1])
		elif arg.begins_with("--height="):
			_height = float(arg.split("=")[1])
		elif arg.begins_with("--seconds="):
			_seconds = float(arg.split("=")[1])
		elif arg.begins_with("--internal="):
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", float(arg.split("=")[1]))
		elif arg.begins_with("--substeps="):
			ProjectSettings.set_setting("physics/avbd/substeps", int(arg.split("=")[1]))
		elif arg.begins_with("--iterations="):
			ProjectSettings.set_setting("physics/avbd/iterations", int(arg.split("=")[1]))
	_run()

func _run() -> void:
	var ground := StaticBody3D.new()
	var gs := CollisionShape3D.new()
	var gshape := BoxShape3D.new()
	gshape.size = Vector3(40, 1, 40)
	gs.shape = gshape
	ground.add_child(gs)
	ground.position = Vector3(0, -0.5, 0)
	root.add_child(ground)

	var boxes: Array[RigidBody3D] = []
	for i in _n:
		var body := RigidBody3D.new()
		body.mass = 1.0
		body.can_sleep = false
		var cs := CollisionShape3D.new()
		var shape := BoxShape3D.new()
		shape.size = Vector3(1, 0.5, 1)
		cs.shape = shape
		body.add_child(cs)
		body.position = Vector3(0, 0.26 + i * 0.52 + _height, 0)
		root.add_child(body)
		boxes.append(body)

	var hz := float(Engine.physics_ticks_per_second)
	var ticks := int(round(_seconds * hz))
	var worst_over := 0.0
	var max_v := 0.0
	for tick in ticks:
		await physics_frame
		if tick % 12 != 11:
			continue
		var deepest := 0.0
		var pair := ""
		for i in boxes.size():
			var a := boxes[i]
			max_v = maxf(max_v, a.linear_velocity.length())
			for j in range(i + 1, boxes.size()):
				var b := boxes[j]
				var d := _sat_depth(a.global_transform, b.global_transform)
				if d > deepest:
					deepest = d
					pair = "%d/%d" % [i, j]
		worst_over = maxf(worst_over, maxf(deepest, maxf(0.0, 0.0)))
		if tick % 60 == 59:
			var ys := PackedFloat32Array()
			for b in boxes:
				ys.append(snappedf(b.global_position.y, 0.0001))
			print("t=%.2fs deepest=%.3f (%s) max_v=%.2f ys=%s" % [
				(tick + 1) / hz, deepest, pair, max_v, ys])
	print("RESULT cfg=%s deepest=%.4f max_v=%.2f" % [OS.get_cmdline_user_args(), worst_over, max_v])
	quit()

func _sat_depth(ta: Transform3D, tb: Transform3D) -> float:
	var half := Vector3(0.5, 0.25, 0.5)
	var ax: Array[Vector3] = [ta.basis.x.normalized(), ta.basis.y.normalized(), ta.basis.z.normalized()]
	var bx: Array[Vector3] = [tb.basis.x.normalized(), tb.basis.y.normalized(), tb.basis.z.normalized()]
	var delta := tb.origin - ta.origin
	var best := 1.0e9
	for i in 3:
		var ra: float = half[i]
		var rb := 0.0
		for j in 3:
			rb += half[j] * absf(ax[i].dot(bx[j]))
		var d := absf(delta.dot(ax[i])) - (ra + rb)
		if d > 0.0:
			return 0.0
		best = minf(best, -d)
	for j in 3:
		var rb2: float = half[j]
		var ra2 := 0.0
		for i in 3:
			ra2 += half[i] * absf(ax[i].dot(bx[j]))
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
				r1 += half[k] * absf(ax[k].dot(axis))
				r2 += half[k] * absf(bx[k].dot(axis))
			var d3 := absf(delta.dot(axis)) - (r1 + r2)
			if d3 > 0.0:
				return 0.0
			best = minf(best, -d3)
	return best
