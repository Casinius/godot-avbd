# Same scene, same frame budget, two engines: time and ground penetration.
#
#   Godot --headless --path demo --fixed-fps 60 --quit-after N --script res://tests/probe_engines.gd
#
# `tools/compare_engines.py` reports wall time only. This adds the number that makes the
# comparison meaningful: how far the settled pyramid sinks into the slab, and how fast the
# pile is still moving at the end. Both engines get the identical scene and frame count; the
# solver budgets are whatever each engine's project settings say (see the printout).
extends SceneTree

const SCENE := "res://scenes/pyramid_dbl.tscn"
# The bottom row's authored centre height is the box half height (1 x 0.5 x 1 box).
const REST_Y := 0.25
const SLAB_HALF_X := 40.0

var _bodies: Array[RigidBody3D] = []


func _collect(node: Node) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			_bodies.append(child)
		_collect(child)


func _initialize() -> void:
	_run()


func _run() -> void:
	await process_frame
	var packed: PackedScene = load(SCENE)
	var scene := packed.instantiate()
	scene.set_script(null)
	root.add_child(scene)
	_collect(scene)
	for body in _bodies:
		body.can_sleep = false
		body.sleeping = false

	var engine: String = ProjectSettings.get_setting("physics/3d/physics_engine")
	var frames := 1200
	var t0 := Time.get_ticks_usec()
	for i in frames:
		await physics_frame
	var wall := float(Time.get_ticks_usec() - t0) / 1.0e6

	var min_y := INF
	var max_speed := 0.0
	var off_slab := 0
	for body in _bodies:
		var p := body.global_position
		if absf(p.x) > SLAB_HALF_X or absf(p.z) > SLAB_HALF_X:
			off_slab += 1
			continue
		min_y = minf(min_y, p.y)
		max_speed = maxf(max_speed, body.linear_velocity.length())

	print("ENGINE=%s bodies=%d frames=%d wall=%.3fs (%.2f ms/frame) sink=%.2f mm max_speed=%.2f m/s off_slab=%d"
			% [engine, _bodies.size(), frames, wall, wall * 1000.0 / float(frames),
			(REST_Y - min_y) * 1000.0, max_speed, off_slab])
	quit(0)
