# Demo driver for the AVBD scenes.
#
# The scene's physics runs through the AVBD physics server (project setting
# physics/3d/physics_engine = "AVBD"): the nodes here are Godot's own RigidBody3D /
# StaticBody3D. The soft-body scene additionally carries an AVBDSoftWorld3D, which
# steps its own solver for the lattices.
#
#   Godot --path demo                                   # windowed
#   Godot --headless --path demo --quit-after 300       # headless, prints stats
#
# Controls: left click punches along the camera ray, Space pauses (the whole
# SceneTree - the physics server follows through set_active), R resets every body
# back to its scene pose.
extends Node3D

@export var impulse_strength := 12.0
@export var stats_interval := 1.0

var bodies: Array[RigidBody3D] = []
var initial: Array[Transform3D] = []
var soft_world: AVBDSoftWorld3D
var elapsed := 0.0
var steps := 0
var quit_after_steps := 0


func _ready() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--steps="):
			quit_after_steps = int(arg.split("=")[1])
	soft_world = get_node_or_null("World") as AVBDSoftWorld3D
	_collect_bodies(self)
	print("AVBD demo: engine=%s, %d rigid bodies%s" % [
			ProjectSettings.get_setting("physics/3d/physics_engine"),
			bodies.size(),
			" + soft lattices" if soft_world != null else ""])


func _collect_bodies(node: Node) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			bodies.append(child)
			initial.append(child.global_transform)
		_collect_bodies(child)


func _physics_process(_delta: float) -> void:
	if quit_after_steps <= 0:
		return
	steps += 1
	if steps < quit_after_steps:
		return
	var server := PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_ACTIVE_OBJECTS)
	var contacts := PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_COLLISION_PAIRS)
	print("after %d steps: active objects=%d collision pairs=%d%s" % [steps, server, contacts,
			" (soft bodies stepped separately)" if soft_world != null else ""])
	print("DEMO_STEPS_OK")
	get_tree().quit(0)


func _process(delta: float) -> void:
	elapsed += delta
	if elapsed < stats_interval:
		return
	elapsed = 0.0
	print("active=%d pairs=%d%s" % [
			PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_ACTIVE_OBJECTS),
			PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_COLLISION_PAIRS),
			"  [paused]" if get_tree().paused else ""])


func _unhandled_input(event: InputEvent) -> void:
	if not event.is_pressed():
		return
	if event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT:
		_punch(event.position)
	elif event is InputEventKey:
		if event.keycode == KEY_SPACE:
			# Pausing the SceneTree is what drives PhysicsServer3D.set_active(false) -
			# the AVBD server must follow and freeze with the game.
			get_tree().paused = not get_tree().paused
		elif event.keycode == KEY_R:
			_reset()


func _punch(screen_position: Vector2) -> void:
	var camera := get_viewport().get_camera_3d()
	if camera == null:
		return
	var from := camera.project_ray_origin(screen_position)
	var to := from + camera.project_ray_normal(screen_position) * 500.0
	# Direct space state query - answered by the AVBD server's AVBDDirectSpaceState3D.
	var query := PhysicsRayQueryParameters3D.create(from, to)
	var hit := get_viewport().get_world_3d().direct_space_state.intersect_ray(query)
	if hit.is_empty():
		return
	var body: RID = hit.get("rid")
	if body == RID():
		return
	var node: Object = instance_from_id(hit.get("collider_id"))
	if node is RigidBody3D:
		node.apply_impulse(camera.project_ray_normal(screen_position) * impulse_strength,
				hit.get("position") - (node as RigidBody3D).global_position)


func _reset() -> void:
	get_tree().paused = false
	for i in bodies.size():
		var body := bodies[i]
		body.global_transform = initial[i]
		body.linear_velocity = Vector3.ZERO
		body.angular_velocity = Vector3.ZERO
		body.sleeping = false
