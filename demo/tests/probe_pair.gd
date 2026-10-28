extends SceneTree

# Minimal impact probe: a hammer box drops onto a resting target box over a static
# ground plane. Logs, per tick, the SAT penetration between the two dynamic boxes and
# the velocities, so the "rigid but penetrates after one collision" report can be
# reproduced on any solver branch with the same numbers.
#
#   Godot --headless --path demo --script res://tests/probe_pair.gd -- --seconds=2

const HALF := Vector3(0.5, 0.25, 0.5)

var _seconds := 2.0
var _drop_y := 3.0
var _hammer_mass := 1.0
var _tick_log := false

func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--seconds="):
			_seconds = float(arg.split("=")[1])
		elif arg.begins_with("--drop_y="):
			_drop_y = float(arg.split("=")[1])
		elif arg.begins_with("--mass="):
			_hammer_mass = float(arg.split("=")[1])
		elif arg.begins_with("--log="):
			_tick_log = arg.split("=")[1] == "1"
	_run()

func _make_box(parent: Node, name: String, half: Vector3, pos: Vector3) -> RigidBody3D:
	var body := RigidBody3D.new()
	body.name = name
	var shape := BoxShape3D.new()
	shape.size = half * 2.0
	var cs := CollisionShape3D.new()
	cs.shape = shape
	body.add_child(cs)
	body.transform = Transform3D(Basis(), pos)
	parent.add_child(body)
	return body

func _run() -> void:
	Engine.physics_ticks_per_second = 240
	ProjectSettings.set_setting("physics/3d/default_gravity", 9.8)
	var scene := Node3D.new()
	scene.name = "Pair"
	root.add_child(scene)

	var ground := _make_box(scene, "Ground", Vector3(4.0, 0.5, 4.0), Vector3(0.0, -0.5, 0.0))
	ground.freeze = true
	var target := _make_box(scene, "Target", HALF, Vector3(0.0, 0.25, 0.0))
	var hammer := _make_box(scene, "Hammer", HALF, Vector3(0.0, _drop_y, 0.0))
	target.can_sleep = false
	hammer.can_sleep = false
	hammer.mass = _hammer_mass

	var ticks := int(round(_seconds * 240.0))
	var max_depth := 0.0
	var max_depth_t := 0.0
	var max_speed := 0.0
	var sum_depth := 0.0
	var depth_ticks := 0
	for tick in ticks:
		await physics_frame
		var depth := _sat_depth(hammer.global_transform, target.global_transform)
		var speed := hammer.linear_velocity.length()
		max_speed = maxf(max_speed, speed)
		if depth > max_depth:
			max_depth = depth
			max_depth_t = float(tick) / 240.0
		if depth > 0.01:
			sum_depth += depth
			depth_ticks += 1
		if _tick_log:
			print("t=%.4f hy=%.4f hvy=%.3f ty=%.4f depth=%.4f speed=%.2f" % [
				float(tick) / 240.0, hammer.global_position.y, hammer.linear_velocity.y,
				target.global_position.y, depth, speed])
	print("RESULT cfg=%s max_depth=%.4f at_t=%.3f depth_ticks=%d mean_depth=%.4f max_speed=%.2f final_hy=%.4f final_ty=%.4f" % [
		OS.get_cmdline_user_args(), max_depth, max_depth_t, depth_ticks,
		sum_depth / float(maxi(1, depth_ticks)), max_speed,
		hammer.global_position.y, target.global_position.y])
	quit()

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
