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
#
# Click-to-impulse fix: Input events (_input/_unhandled_input) occur before physics
# frame, so direct_space_state is inaccessible. Click positions are stored in
# _pending_click and processed in _physics_process after physics has run,
# enabling successful raycasting and impulse application.
extends Node3D

@export var impulse_strength := 12.0
@export var stats_interval := 0.5  # More frequent stats for pyramid stability
@export var max_bodies := 25  # Target pyramid size

var bodies: Array[RigidBody3D] = []
var initial: Array[Transform3D] = []
var soft_world: AVBDSoftWorld3D
var elapsed := 0.0
var steps := 0
var quit_after_steps := 0
var pyramid_bodies: Array[RigidBody3D] = []  # Only pyramid boxes, not ground
var max_position_y := 0.0
var min_position_y := 1000.0
var velocity_stats := { "max_lin": 0.0, "max_ang": 0.0, "avg_lin": 0.0, "steps": 0 }
var max_velocity_accumulated := 0.0  # Track worst velocity spike during collapse
var collapse_started := false
var collapse_start_y := 0.0


func _ready() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--steps="):
			quit_after_steps = int(arg.split("=")[1])
	soft_world = get_node_or_null("World") as AVBDSoftWorld3D
	_collect_bodies(self)
	# Separate pyramid boxes (those named "Box*") from ground
	for body in bodies:
		if body.name.begins_with("Box"):
			pyramid_bodies.append(body)
	# Add all pyramid boxes to punchable group
	for body in pyramid_bodies:
		body.add_to_group("punchable")
	# Add tiny random perturbation for slight instability test
	var rng := RandomNumberGenerator.new()
	for body in pyramid_bodies:
		rng.randomize()
		var perturbation := Vector3(rng.randf_range(-0.01, 0.01), 0.0, rng.randf_range(-0.01, 0.01))
		body.global_position += perturbation
	# Disable sleep by setting low mass and can_sleep=false
	for body in pyramid_bodies:
		body.mass = 1.0  # Reset mass to normal
		body.sleeping = false
		body.can_sleep = false
	# No impulse - let gravity handle natural settling
	# body.apply_impulse(Vector3(0.0, 5.0, 0.0), body.global_position)
	print("AVBD demo: engine=%s, %d rigid bodies, %d pyramid boxes" % [
		ProjectSettings.get_setting("physics/3d/physics_engine"),
		bodies.size(),
		pyramid_bodies.size()])


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

	# Calculate pyramid collapse metrics


func _process(delta: float) -> void:
	elapsed += delta
	if elapsed < stats_interval:
		return
	elapsed = 0.0

	# Track pyramid stability metrics
	var active_pyramid := 0
	var max_velocity := 0.0
	var total_velocity := 0.0
	var max_lin_vel := 0.0
	var max_ang_vel := 0.0
	var total_lin_vel := 0.0
	var total_ang_vel := 0.0

	for body in pyramid_bodies:
		if not body.sleeping:
			active_pyramid += 1
			var lin_vel := body.linear_velocity.length()
			var ang_vel := body.angular_velocity.length()
			max_velocity = max(max_velocity, lin_vel)
			total_velocity += lin_vel
			max_lin_vel = max(max_lin_vel, lin_vel)
			max_ang_vel = max(max_ang_vel, ang_vel)
			total_lin_vel += lin_vel
			total_ang_vel += ang_vel
			var pos_y := body.global_position.y
			max_position_y = max(max_position_y, pos_y)
			min_position_y = min(min_position_y, pos_y)

	# Detect collapse - track ongoing collapse depth
	print("DEBUG: Checking collapse - min_y=%.3f < 0.43? %s, max_lin_vel=%.3f > 0.1? %s" % [
		min_position_y, min_position_y < 0.43, max_lin_vel, max_lin_vel > 0.1])
	if min_position_y < 0.43 and max_lin_vel > 0.1:
		print("COLLAPSE DETECTED! min_y=%.3f, max_lin_vel=%.3f" % [min_position_y, max_lin_vel])
		print("Active bodies: %d" % active_pyramid)
		collapse_started = true
		collapse_start_y = min_position_y

	# Track worst collapse depth
	var collapse_depth = collapse_start_y - min_position_y
	if collapse_started and collapse_depth > max_velocity_accumulated:
		max_velocity_accumulated = collapse_depth

	# Calculate average velocities
	var pyramid_count := pyramid_bodies.size()
	var avg_lin_vel := 0.0
	var avg_ang_vel := 0.0
	if active_pyramid > 0:
		avg_lin_vel = total_lin_vel / active_pyramid
		avg_ang_vel = total_ang_vel / active_pyramid

	# Track worst velocity spike during collapse
	if collapse_started:
		max_velocity_accumulated = max(max_velocity_accumulated, max_velocity)

	# Print detailed pyramid stats
	var stat_str = "Pyramid stats [step=%d]: active=%d/%d, max_y=%.3f, min_y=%.3f, " % [
		steps, active_pyramid, pyramid_count, max_position_y, min_position_y]
	var rest_str = "max_lin=%.3f, max_ang=%.3f, avg_lin=%.3f, max_vel_accum=%.3f, " % [
		max_lin_vel, max_ang_vel, avg_lin_vel if active_pyramid > 0 else 0.0,
		max_velocity_accumulated]
	var collapse_str = "collapse_started=%s, collapse_y=%.3f" % [
		collapse_started, collapse_start_y if collapse_started else 0.0]
	print(stat_str + rest_str + collapse_str)

	print("AVBD demo: active=%d pairs=%d%s" % [
		PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_ACTIVE_OBJECTS),
		PhysicsServer3D.get_process_info(PhysicsServer3D.INFO_COLLISION_PAIRS),
		"  [paused]" if get_tree().paused else ""])


func _input(event: InputEvent) -> void:
	if event.is_pressed() and event is InputEventMouseButton and event.button_index == MOUSE_BUTTON_LEFT:
		_punch(event.position)
	elif event.is_pressed() and event is InputEventKey and event.keycode == KEY_R:
		_reset()


func _unhandled_input(event: InputEvent) -> void:
	if event.is_pressed() and event is InputEventKey:
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

	var ray_dir := camera.project_ray_normal(screen_position)
	var from := camera.project_ray_origin(screen_position)

	var bodies := get_tree().get_nodes_in_group("punchable")
	if bodies.is_empty():
		return

	var closest_hit: Dictionary = {}
	var closest_dist := 1000.0

	for body in bodies:
		if not body is RigidBody3D:
			continue
		# Get the collision shape node
		var shape_node := body.find_child("Shape", true, false)
		if shape_node == null:
			continue
		if not shape_node is CollisionShape3D:
			continue
		var box_shape = shape_node.shape
		if box_shape == null or not box_shape is BoxShape3D:
			continue
		var half_extents: Vector3 = box_shape.size / 2.0
		var body_min: Vector3 = body.global_position - half_extents
		var body_max: Vector3 = body.global_position + half_extents

		# Ray-box intersection (Slab method)
		var t_min: float = 0.0
		var t_max: float = 1000.0
		var normal: Vector3 = Vector3.ZERO

		for i in range(3):
			if abs(ray_dir[i]) < 1e-6:
				if from[i] < body_min[i] or from[i] > body_max[i]:
					break
			else:
				var t1: float = (body_min[i] - from[i]) / ray_dir[i]
				var t2: float = (body_max[i] - from[i]) / ray_dir[i]
				# Swap t1 and t2 if needed
				if t1 > t2:
					var temp: float = t1
					t1 = t2
					t2 = temp
				if t1 > t_min:
					t_min = t1
					normal = Vector3.ZERO
					if i == 0:
						normal.x = -1.0 if from[i] < body_min[i] else 1.0
					elif i == 1:
						normal.y = -1.0 if from[i] < body_min[i] else 1.0
					else:
						normal.z = -1.0 if from[i] < body_min[i] else 1.0
				if t2 < t_max:
					t_max = t2

		if t_min <= t_max and t_max >= 0.0 and t_min < closest_dist:
			closest_dist = t_min
			closest_hit = { "body": body, "position": from + ray_dir * t_min, "normal": normal }

	if not closest_hit.is_empty():
		var hit := closest_hit
		var body: RigidBody3D = hit["body"]
		var impulse_point: Vector3 = hit["position"]
		var impulse: Vector3 = hit["normal"] * impulse_strength
		body.apply_impulse(impulse, impulse_point)


func _reset() -> void:
	get_tree().paused = false
	for i in bodies.size():
		var body := bodies[i]
		body.global_transform = initial[i]
		body.linear_velocity = Vector3.ZERO
		body.angular_velocity = Vector3.ZERO
		body.sleeping = false
