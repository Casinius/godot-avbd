# Single-chassis raycast car demo.
# Launch: godot --path demo res://scenes/car_demo.tscn
# Controls: Up throttle, Down brake, Left/Right steer.
# The four wheels are visual meshes only; the one chassis is the only dynamic body.
extends RigidBody3D

const CHASSIS_SIZE := Vector3(1.8, 0.6, 3.6)
const CHASSIS_MASS := 900.0
const REST_LENGTH := 0.45
const WHEEL_RADIUS := 0.3
const SPRING_STIFFNESS := 22000.0
const DAMPING := 2400.0
const MAX_SUSPENSION_FORCE := 7000.0
const MAX_DRIVE_FORCE := 7000.0
const MAX_BRAKE_FORCE := 10000.0
const MAX_STEER_FORCE := 4000.0
const FRONT_Z := -1.2
const REAR_Z := 1.2
const TRACK_HALF := 0.72
const UP := Vector3.UP

var control_override := false
var override_throttle := 0.0
var override_brake := 0.0
var override_steering := 0.0

var _simulation_seconds := 0.0
var _report_seconds := 0.0
var _rays: Array[RayCast3D] = []
var _wheels: Array[MeshInstance3D] = []
var _wheel_offsets := [
	Vector3(-TRACK_HALF, -0.15, FRONT_Z),
	Vector3(TRACK_HALF, -0.15, FRONT_Z),
	Vector3(-TRACK_HALF, -0.15, REAR_Z),
	Vector3(TRACK_HALF, -0.15, REAR_Z),
]

func _ready() -> void:
	for child in get_children():
		if child is RayCast3D:
			_rays.append(child as RayCast3D)
		elif child is MeshInstance3D and child.name.begins_with("Wheel"):
			_wheels.append(child as MeshInstance3D)
	print("AVBD car demo: engine=%s, one dynamic chassis, controls=arrows (up/down throttle/brake, left/right steer)" % ProjectSettings.get_setting("physics/3d/physics_engine"))


func _integrate_forces(state: PhysicsDirectBodyState3D) -> void:
	# AVBD exposes the direct-state force path as a persistent value. Replace it every
	# callback so the controller output represents exactly one physics frame.
	var total_force := Vector3.ZERO
	var total_torque := Vector3.ZERO
	var throttle := 0.0
	var brake := 0.0
	var steering := 0.0
	if control_override:
		throttle = clampf(override_throttle, -1.0, 1.0)
		brake = clampf(override_brake, -1.0, 1.0)
		steering = clampf(override_steering, -1.0, 1.0)
	else:
		var longitudinal := Input.get_axis("ui_down", "ui_up")
		throttle = clampf(maxf(longitudinal, 0.0), -1.0, 1.0)
		brake = clampf(maxf(-longitudinal, 0.0), -1.0, 1.0)
		steering = clampf(Input.get_axis("ui_left", "ui_right"), -1.0, 1.0)

	var basis := state.transform.basis
	var forward := basis * Vector3(0.0, 0.0, -1.0)
	var right := basis * Vector3.RIGHT
	for ray in _rays:
		if not ray.is_colliding():
			continue
		var distance := ray.global_position.distance_to(ray.get_collision_point())
		var compression := maxf(0.0, REST_LENGTH - (distance - WHEEL_RADIUS))
		var local_offset := ray.position
		var world_offset := basis * local_offset
		var contact_velocity := state.linear_velocity + state.angular_velocity.cross(world_offset)
		var damping_force := DAMPING * contact_velocity.dot(UP)
		var suspension_force := clampf(SPRING_STIFFNESS * compression - damping_force, 0.0, MAX_SUSPENSION_FORCE)
		total_force += UP * suspension_force
		total_torque += world_offset.cross(UP * suspension_force)

	var forward_speed := state.linear_velocity.dot(forward)
	var drive_force := forward * (throttle * MAX_DRIVE_FORCE)
	for side in [-1.0, 1.0]:
		var offset := basis * Vector3(side * TRACK_HALF, 0.0, REAR_Z)
		total_force += drive_force * 0.5
		total_torque += offset.cross(drive_force * 0.5)

	if brake > 0.0 and absf(forward_speed) > 0.001:
		var braking := -forward * signf(forward_speed) * minf(brake * MAX_BRAKE_FORCE, MAX_BRAKE_FORCE)
		for side in [-1.0, 1.0]:
			var offset := basis * Vector3(side * TRACK_HALF, 0.0, REAR_Z)
			total_force += braking * 0.5
			total_torque += offset.cross(braking * 0.5)

	var lateral_force := right * steering * minf(absf(forward_speed) * 1200.0, MAX_STEER_FORCE)
	for side in [-1.0, 1.0]:
		var offset := basis * Vector3(side * TRACK_HALF, 0.0, FRONT_Z)
		total_force += lateral_force * 0.5
		total_torque += offset.cross(lateral_force * 0.5)

	state.set_constant_force(total_force)
	state.set_constant_torque(total_torque)


func _physics_process(delta: float) -> void:
	_simulation_seconds += delta
	_report_seconds += delta
	if _report_seconds < 1.0:
		return
	_report_seconds -= 1.0
	var finite := global_position.is_finite() and linear_velocity.is_finite() and angular_velocity.is_finite()
	var longitudinal := Input.get_axis("ui_down", "ui_up") if not control_override else override_throttle - override_brake
	var throttle := maxf(longitudinal, 0.0)
	var brake := maxf(-longitudinal, 0.0)
	var steering := Input.get_axis("ui_left", "ui_right") if not control_override else override_steering
	print("AVBD car: t=%.2fs speed=%.3f m/s x=%.3f z=%.3f steer=%.2f throttle=%.2f brake=%.2f finite=%s" % [
		_simulation_seconds, linear_velocity.length(), global_position.x, global_position.z,
		steering, throttle, brake, finite])


func _process(_delta: float) -> void:
	for i in mini(_rays.size(), _wheels.size()):
		var ray := _rays[i]
		var wheel := _wheels[i]
		var local_y := -(REST_LENGTH + WHEEL_RADIUS)
		if ray.is_colliding():
			local_y = -0.15 - ray.global_position.distance_to(ray.get_collision_point()) + WHEEL_RADIUS
		wheel.position = Vector3(_wheel_offsets[i].x, local_y, _wheel_offsets[i].z)
		if i < 2:
			wheel.rotation.y = clampf(override_steering if control_override else Input.get_axis("ui_left", "ui_right"), -1.0, 1.0) * deg_to_rad(25.0)
