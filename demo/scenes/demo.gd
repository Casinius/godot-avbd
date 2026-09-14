# Demo driver for the AVBD scenes: prints solver statistics and lets the mouse
# shove bodies around.
#
#   Godot_v4.6.3 --path demo                                   # windowed
#   Godot_v4.6.3 --headless --path demo --quit-after 300       # headless, prints stats
#
# Controls: left click punches along the camera ray, Space pauses, R resets every
# body back to its scene pose.
extends Node3D

@export var impulse_strength := 12.0
@export var stats_interval := 1.0

var world: AVBDWorld3D
var bodies: Array[AVBDRigidBody3D] = []
var initial: Array[Transform3D] = []
var elapsed := 0.0
var steps := 0
var quit_after_steps := 0
var screenshot_path := ""


func _ready() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--steps="):
			quit_after_steps = int(arg.split("=")[1])
		elif arg.begins_with("--screenshot="):
			screenshot_path = arg.split("=")[1]
	world = get_node_or_null("World") as AVBDWorld3D
	if world == null:
		push_warning("demo.gd: no AVBDWorld3D child named 'World'")
		return
	for child in world.get_children():
		if child is AVBDRigidBody3D:
			bodies.append(child)
			initial.append(child.transform)
	print("AVBD demo: %d rigid bodies + soft bodies under %s" % [bodies.size(), world.name])


func _physics_process(_delta: float) -> void:
	if quit_after_steps <= 0:
		return
	steps += 1
	if steps < quit_after_steps:
		return
	print("after %d steps: bodies=%d forces=%d contacts=%d points=%d step=%.3f ms" % [steps,
			world.get_body_count(), world.get_force_count(), world.get_contact_count(),
			world.get_contact_point_count(), world.get_step_time_usec() / 1000.0])
	if not screenshot_path.is_empty():
		await RenderingServer.frame_post_draw
		var image := get_viewport().get_texture().get_image()
		image.save_png(screenshot_path)
		print("saved ", screenshot_path)
	print("DEMO_STEPS_OK")
	get_tree().quit(0)


func _process(delta: float) -> void:
	elapsed += delta
	if world == null or elapsed < stats_interval:
		return
	elapsed = 0.0
	print("bodies=%d forces=%d contacts=%d points=%d step=%.3f ms%s" % [world.get_body_count(),
			world.get_force_count(), world.get_contact_count(), world.get_contact_point_count(),
			world.get_step_time_usec() / 1000.0, "  [paused]" if world.paused else ""])


func _unhandled_input(event: InputEvent) -> void:
	if world == null or not event.is_pressed():
		return
	if event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT:
		_punch(event.position)
	elif event is InputEventKey:
		if event.keycode == KEY_SPACE:
			world.paused = not world.paused
		elif event.keycode == KEY_R:
			_reset()


func _punch(screen_position: Vector2) -> void:
	var camera := get_viewport().get_camera_3d()
	if camera == null:
		return
	var from := camera.project_ray_origin(screen_position)
	var direction := camera.project_ray_normal(screen_position)
	var hit := world.raycast(from, direction, 500.0)
	if hit.is_empty():
		return
	var body := hit.get("body") as AVBDRigidBody3D
	if body == null:
		return
	body.apply_impulse(direction * impulse_strength, hit.get("position") - body.global_position)


func _reset() -> void:
	for i in bodies.size():
		var pose := initial[i]
		bodies[i].teleport(pose.origin, pose.basis.get_rotation_quaternion())
