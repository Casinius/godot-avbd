# Long-run penetration detection for the double-pyramid stress scene.
#
# pyramid_dbl stacks two 20-row pyramids: the lower one settles on the slab and the upper
# one starts ~1.4 m above it, so the pile goes through a *secondary* collision
# (upper-on-lower) after the first ground contacts. Deep shaking/tunnelling would show up
# there first, so the scenario watches penetration for a long wall of simulated time.
#
# Penetration is read from the solver's own contact error through the diagnostic added to
# the AVBD server (PhysicsServer3D.call("get_solver_max_penetration")): a manifold's C0 is
# the geometric overlap captured at the start of the step, so this is the real interpenetration
# and stays correct once the pile rotates (an AABB/transform heuristic would not).
#
# Tracked every physics frame:
#   - the deepest contact anywhere, and the instance-id pair when it is sampled
#   - "secondary collisions": a body that had come to rest is set moving again by a new
#     impact (speed crosses IMPACT_SPEED after REST_FRAMES frames under REST_SPEED)
#   - the deepest penetration seen while such a secondary impact is resolving
#   - sink-through: the lowest body centre against the slab top
#
# The strict scenario (`pyramid_dbl_penetration_strict`) samples the same diagnostic every
# physics frame but keeps the whole distribution instead of only the maximum: percentiles,
# per-threshold frame counts and the longest run above 1 cm distinguish an impact spike
# from a crushed contact. It also measures the lowest body centre, which is an independent
# geometric check of ground penetration (a flat box at rest has its centre at y = 0.25 m).
#
#   Godot --headless --path demo --script res://tests/test_penetration.gd
#   ... -- --only=pyramid_dbl_penetration
#   ... -- --only=pyramid_dbl_penetration_strict
#   ... -- --pen-seconds=180          (longer run)

extends "res://tests/avbd_harness.gd"

const SCENE := "res://scenes/pyramid_dbl.tscn"
const BODY_COUNT := 420
const RUN_SECONDS := 90.0  # logical seconds (60 Hz); override with --pen-seconds=
const REST_FRAMES := 3     # frames under REST_SPEED before a body counts as resting
const REST_SPEED := 0.25   # m/s
const IMPACT_SPEED := 0.5  # m/s that re-wakes a resting body = secondary collision
const SECONDARY_WINDOW := 12
# Bodies above this authored height belong to the dropped (upper) pyramid; the lower
# pyramid's apex is at y ~ 10.1 and the upper base is at y ~ 12.1.
const UPPER_PYRAMID_Y := 11.0

# Strict scenario calibration. The warmup skips the authored 1.4 m drop of the upper
# pyramid, so the assertions describe a settled/persistently-loaded pile, and the settled
# gate (RESTING_PAIR_SPEED below) excludes pairs that are still moving.
#
# Measured on the 420-box double pyramid (both can_sleep settings, 90 s): the deepest pair's
# per-frame maximum sits at p50 ~ 1.4 cm and p99 ~ 2.0 cm. That residual is the augmented
# Lagrangian's equilibrium under a 20-row column - the worst-loaded contact never quite
# reaches the projection slop - not an impact: the failure these checks exist for (loss of
# stiffness under the landing) read 0.5-1.0 m at 13-15 m/s. The bounds therefore sit just
# above the measured resting distribution (3 cm transient, 2.5 cm p99, 2 cm sustained for
# 50 ms) and still fail by an order of magnitude on a crush.
const STRICT_WARMUP_SECONDS := 10.0
const STRICT_WORST_LIMIT := 0.035
const STRICT_P99_LIMIT := 0.025
const STRICT_SUSTAINED_LIMIT := 0.02
const STRICT_SUSTAINED_FRAMES := 12
const GROUND_TOP := 0.0
const FLAT_BOX_HALF_HEIGHT := 0.25
const GROUND_SINK_LIMIT := 0.02
const GROUND_HALF_X := 40.0
# A pair only counts as settled below this speed: the numbers describe the pile's contact
# stiffness at rest, not the overlap of a body still being slammed into the pile. The
# failure these checks guard against (loss of stiffness) reads 0.5-1 m at 13-15 m/s, so a
# 0.5 m/s gate cannot hide it.
const RESTING_PAIR_SPEED := 0.5
# The authored pyramid box is 1 x 0.5 x 1 (Godot frame, Y up).
const BOX_HALF_EXTENTS := Vector3(0.5, 0.25, 0.5)

var _run_seconds := RUN_SECONDS
var _sweep_config := ""
var _ground_z := 20.0
# The scenarios park sleeping by default (can_sleep=false on every box) so a penetration
# number is never measured against a body the solver has stopped solving. `sleep-on` leaves
# the engine default in place to compare against.
var _disable_sleep := true
# pyramid_dbl is two 20-row pyramids, the upper one dropped 1.43 m onto the lower one's
# apex box. `--upper=drop` removes it so the same measurements can be taken on a
# configuration that is physically stable.
var _drop_upper := false


func _initialize() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--pen-seconds="):
			_run_seconds = float(arg.split("=")[1])
		elif arg.begins_with("--sweep-config="):
			_sweep_config = arg.split("=")[1]
		elif arg.begins_with("--ground-z="):
			_ground_z = float(arg.split("=")[1])
		elif arg.begins_with("--upper="):
			_drop_upper = arg.split("=")[1] == "drop"
	super._initialize()


# Solver-budget sweep for the strict scenario. Each case runs in its own process, so this
# applies the settings before the scene (and its bodies) exist; the server re-reads them on
# `settings_changed` and the two awaited frames let that land before the first step.
func _apply_sweep_config() -> void:
	match _sweep_config:
		"", "current":
			return
		"old-budget":
			# The schedule the pyramid regression was authored under: one solver step per
			# 240 Hz engine tick, 12 adaptive iterations (6 Newton rounds).
			ProjectSettings.set_setting("physics/avbd/interleaved_substeps", false)
			ProjectSettings.set_setting("physics/avbd/substeps", 1)
			ProjectSettings.set_setting("physics/avbd/iterations", 12)
			ProjectSettings.set_setting("physics/avbd/internal_ticks_per_second", 0.0)
		"newton-1.0":
			# Current 120 Hz x 2 substeps x 2 rounds, but both rounds are Newton.
			ProjectSettings.set_setting("physics/avbd/newton_ratio", 1.0)
			ProjectSettings.set_setting("physics/avbd/stiffness_decay", 1.0)
		"rounds-8":
			ProjectSettings.set_setting("physics/avbd/rounds_per_substep", 8)
			ProjectSettings.set_setting("physics/avbd/newton_ratio", 1.0)
			ProjectSettings.set_setting("physics/avbd/stiffness_decay", 1.0)
		"adaptive":
			# Explicit budget 0: stepInterleaved falls back to the adaptive target.
			ProjectSettings.set_setting("physics/avbd/rounds_per_substep", 0)
		"sleep-on":
			# Same solver settings as `current`; the only difference is that the boxes keep
			# can_sleep=true, so idle detection, island sleep and manifold freezing are all
			# live. Comparing this against `current` isolates the sleep/wake path.
			_disable_sleep = false
		"no-grace":
			ProjectSettings.set_setting("physics/avbd/contact_grace_ticks", 0)
		"no-spec":
			ProjectSettings.set_setting("physics/avbd/speculative_margin", 0.0)
		"no-proj":
			ProjectSettings.set_setting("physics/avbd/position_projection", false)
		"no-adaptive-beta":
			ProjectSettings.set_setting("physics/avbd/adaptive_beta", false)
		"proj-32":
			ProjectSettings.set_setting("physics/avbd/projection_iterations", 32)
		"proj-64":
			ProjectSettings.set_setting("physics/avbd/projection_iterations", 64)
		"spec-05":
			ProjectSettings.set_setting("physics/avbd/speculative_margin", 0.05)
		"spec-06":
			ProjectSettings.set_setting("physics/avbd/speculative_margin", 0.06)
		"spec-10":
			ProjectSettings.set_setting("physics/avbd/speculative_margin", 0.1)
		_:
			push_error("unknown --sweep-config=%s" % _sweep_config)


func _scenarios() -> Array:
	return [
		scenario("pyramid_dbl_penetration", test_pyramid_dbl_penetration),
		scenario("pyramid_dbl_penetration_strict", test_pyramid_dbl_penetration_strict),
	]


func _collect_rigid_bodies(node: Node, out: Array[RigidBody3D]) -> void:
	for child in node.get_children():
		if child is RigidBody3D:
			out.append(child)
		_collect_rigid_bodies(child, out)


func _describe_instance(instance_id: int) -> String:
	if instance_id == 0:
		return "world"
	var obj := instance_from_id(instance_id)
	return str(obj.name) if obj != null else "body#%d" % instance_id


func _body_speed(instance_id: int) -> float:
	var obj := instance_from_id(instance_id)
	if obj is RigidBody3D:
		return (obj as RigidBody3D).linear_velocity.length()
	return 0.0


# The authored slab is 80 x 1 x 20 (Godot x/y/z), depth 20. Enlarging it lets a
# penetration test ignore boxes that merely slide off the narrow slab and fall.
func _set_ground_depth(scene: Node, depth: float) -> void:
	if is_equal_approx(depth, 20.0):
		return
	var shape_node := scene.get_node_or_null("Ground/Shape") as CollisionShape3D
	if shape_node == null:
		return
	var box := shape_node.shape as BoxShape3D
	if box != null:
		box.size = Vector3(box.size.x, box.size.y, depth)


func test_pyramid_dbl_penetration() -> Variant:
	check_engine()
	check(PhysicsServer3D.has_method("get_solver_max_penetration"),
			"penetration diagnostic present", "PhysicsServer3D.get_solver_max_penetration")

	var packed: PackedScene = load(SCENE)
	var scene := packed.instantiate()
	# The demo driver adds random perturbation and prints stats; a regression run needs the
	# authored poses and no per-frame noise.
	scene.set_script(null)
	root.add_child(scene)

	var bodies: Array[RigidBody3D] = []
	_collect_rigid_bodies(scene, bodies)
	if not _drop_upper:
		check(bodies.size() == BODY_COUNT, "body count",
				"found %d (expected %d)" % [bodies.size(), BODY_COUNT])

	var n := bodies.size()
	var initial_y := PackedFloat32Array()
	initial_y.resize(n)
	for i in n:
		var body := bodies[i]
		body.mass = 1.0
		if _disable_sleep:
			body.sleeping = false
			body.can_sleep = false
		initial_y[i] = body.global_position.y if body.is_inside_tree() else body.position.y
	# Settle the flagged flag array sizing.
	var rest_streak := PackedInt32Array()
	rest_streak.resize(n)
	var secondary_armed := PackedInt32Array()
	secondary_armed.resize(n)
	var secondary_events := PackedInt32Array()
	secondary_events.resize(n)

	var worst_pen := 0.0
	var worst_pen_frame := 0
	var worst_pen_a := 0
	var worst_pen_b := 0
	var pen_at_secondary := 0.0
	var sec_upper := 0
	var sec_lower := 0
	var sec_upper_bodies := 0
	var sec_lower_bodies := 0
	var min_y := 1.0e30
	var max_y := -1.0e30
	var max_xz := 0.0
	var finite := true

	var ticks := int(_run_seconds * float(Engine.physics_ticks_per_second))
	var t0 := Time.get_ticks_msec()
	for frame in ticks:
		await physics_frame

		var frame_pen: float = PhysicsServer3D.call("get_solver_max_penetration")
		if frame_pen > worst_pen:
			worst_pen = frame_pen
			worst_pen_frame = frame
			var deepest: Array = PhysicsServer3D.call("get_solver_deepest_contact")
			if deepest.size() == 3:
				worst_pen_a = int(deepest[0])
				worst_pen_b = int(deepest[1])

		for i in n:
			var body := bodies[i]
			var speed := body.linear_velocity.length()
			if rest_streak[i] >= REST_FRAMES and speed > IMPACT_SPEED:
				secondary_events[i] += 1
				if secondary_events[i] == 1 and initial_y[i] > UPPER_PYRAMID_Y:
					sec_upper_bodies += 1
				elif secondary_events[i] == 1:
					sec_lower_bodies += 1
				if initial_y[i] > UPPER_PYRAMID_Y:
					sec_upper += 1
				else:
					sec_lower += 1
				secondary_armed[i] = SECONDARY_WINDOW
				rest_streak[i] = 0
			elif speed < REST_SPEED:
				rest_streak[i] = mini(rest_streak[i] + 1, REST_FRAMES)
			else:
				rest_streak[i] = 0

			if secondary_armed[i] > 0:
				secondary_armed[i] -= 1
				pen_at_secondary = maxf(pen_at_secondary, frame_pen)

			var pos := body.global_position
			min_y = minf(min_y, pos.y)
			max_y = maxf(max_y, pos.y)
			max_xz = maxf(max_xz, Vector2(pos.x, pos.z).length())
			finite = finite and is_finite_body(body)

	var elapsed := (Time.get_ticks_msec() - t0) / 1000.0
	var logical_time := float(worst_pen_frame) / float(Engine.physics_ticks_per_second)

	check(finite, "pile stays finite", "bodies=%d" % n)
	# A body whose centre drops below the slab surface by more than its own half-height has
	# gone through the ground (or fallen off the slab's edge, which the report makes visible).
	check(min_y >= -0.25, "no deep ground sink", "min centre y=%.4f m" % min_y)
	# Half a box height (0.25 m) is the hard ceiling for a surviving stack; a resting contact
	# should stay around the 1 mm collision margin. 5 cm already means the pile is being crushed.
	check(worst_pen <= 0.05, "penetration bounded",
			"worst=%.4f m at t=%.1f s (%s <-> %s)" % [
				worst_pen, logical_time, _describe_instance(worst_pen_a), _describe_instance(worst_pen_b)])

	var evidence := {
		"seconds": _run_seconds,
		"wall_seconds": elapsed,
		"worst_penetration": worst_pen,
		"worst_penetration_time": logical_time,
		"worst_pair": "%s <-> %s" % [_describe_instance(worst_pen_a), _describe_instance(worst_pen_b)],
		"penetration_at_secondary_collisions": pen_at_secondary,
		"secondary_events_upper": sec_upper,
		"secondary_events_lower": sec_lower,
		"secondary_bodies_upper": sec_upper_bodies,
		"secondary_bodies_lower": sec_lower_bodies,
		"min_y": min_y,
		"max_y": max_y,
		"max_xz": max_xz,
	}
	scene.queue_free()
	await frames(2)
	return evidence


func _percentile(sorted: PackedFloat32Array, q: float) -> float:
	if sorted.is_empty():
		return 0.0
	var idx := clampi(int(round(q * float(sorted.size() - 1))), 0, sorted.size() - 1)
	return sorted[idx]


# SAT penetration depth between two oriented boxes, in metres: positive means overlapping,
# -INF means a separating axis exists. This is the ground truth the solver's reported
# contact error is compared against - the solver's number comes from a cached contact frame
# and anchors, so the two disagree exactly when that cache has gone stale.
func _obb_penetration(a: Transform3D, b: Transform3D, ha: Vector3, hb: Vector3) -> float:
	# Local axes in world space. Basis * Vector3 is the rotation, so this is unambiguous
	# about rows vs columns.
	var ea := [a.basis * Vector3.RIGHT, a.basis * Vector3.UP, a.basis * Vector3.BACK]
	var eb := [b.basis * Vector3.RIGHT, b.basis * Vector3.UP, b.basis * Vector3.BACK]
	var ax: Array[Vector3] = [ea[0], ea[1], ea[2], eb[0], eb[1], eb[2]]
	for i in 3:
		for j in 3:
			ax.append(ea[i].cross(eb[j]))
	var d := b.origin - a.origin
	var best := INF
	for axis in ax:
		var len := axis.length()
		if len < 1.0e-6:
			continue
		var n := axis / len
		var ra := 0.0
		var rb := 0.0
		for i in 3:
			ra += ha[i] * absf(ea[i].dot(n))
			rb += hb[i] * absf(eb[i].dot(n))
		var sep := absf(d.dot(n)) - ra - rb
		if sep > 0.0:
			return -INF
		best = minf(best, -sep)
	return best


func test_pyramid_dbl_penetration_strict() -> Variant:
	check_engine()
	check(PhysicsServer3D.has_method("get_solver_max_penetration"),
			"penetration diagnostic present", "PhysicsServer3D.get_solver_max_penetration")
	_apply_sweep_config()
	# Let the server's settings_changed handler pick the sweep values up before the scene
	# creates its bodies (rebuild_space reads the parameters at the first step).
	await frames(2)

	var packed: PackedScene = load(SCENE)
	var scene := packed.instantiate()
	scene.set_script(null)
	root.add_child(scene)
	_set_ground_depth(scene, _ground_z)

	var bodies: Array[RigidBody3D] = []
	_collect_rigid_bodies(scene, bodies)
	if not _drop_upper:
		check(bodies.size() == BODY_COUNT, "body count",
				"found %d (expected %d)" % [bodies.size(), BODY_COUNT])
	if _drop_upper:
		# The upper pyramid's base row starts at y = 12.06; the lower apex is at 10.13.
		for body in bodies.filter(func(b): return b.global_position.y > 11.0):
			body.free()
		bodies = bodies.filter(func(b): return is_instance_valid(b)) as Array[RigidBody3D]

	var spawn_pos: Array[Vector3] = []
	spawn_pos.resize(bodies.size())
	for i in bodies.size():
		var body := bodies[i]
		body.mass = 1.0
		if _disable_sleep:
			body.sleeping = false
			body.can_sleep = false
		spawn_pos[i] = body.global_position

	var hz := float(Engine.physics_ticks_per_second)
	var warmup_frames := int(STRICT_WARMUP_SECONDS * hz)
	var ticks := int(_run_seconds * hz)

	# Post-warmup distribution. Every engine frame is one sample; with an internal solver
	# rate below the engine rate, a solver tick can be sampled more than once - which is
	# exactly what the sustained-overlap check wants to see.
	var samples := PackedFloat32Array()
	var finite := true
	var min_y := 1.0e30
	var min_y_name := ""
	var min_inside_y := 1.0e30
	var min_inside_name := ""
	var min_inside_y_any := 1.0e30
	var min_inside_name_any := ""
	var max_abs_x := 0.0
	var max_abs_z := 0.0
	var settled_max_speed := 0.0
	var ground_half_z := _ground_z * 0.5

	# Whole-run trackers: the sweep compares configurations by when the collapse starts, not
	# only by the post-warmup distribution.
	var first_above_5cm := -1.0
	var global_min_y := 1.0e30
	var global_max_speed := 0.0

	var worst_warmup := 0.0
	var worst_warmup_frame := 0
	var worst_warmup_a := 0
	var worst_warmup_b := 0

	var worst_settled := 0.0
	var worst_settled_frame := 0
	var worst_settled_a := 0
	var worst_settled_b := 0
	var worst_a_pos := Vector3.ZERO
	var worst_b_pos := Vector3.ZERO
	var worst_pair_boxes := false
	# Post-warmup maximum over every frame, impact transients included: the comparison
	# against the settled distribution shows what the resting filter removed.
	var worst_any_settled := 0.0
	var worst_any_settled_frame := 0
	# What the two reported boxes actually overlapped at that instant (see _obb_penetration).
	var worst_true_depth := -INF

	var above_1mm := 0
	var above_2mm := 0
	var above_5mm := 0
	var above_10mm := 0
	var above_20mm := 0
	var streak_1mm := 0
	var streak_2mm := 0
	var streak_5mm := 0
	var streak_10mm := 0
	var longest_1mm := 0
	var longest_2mm := 0
	var longest_5mm := 0
	var longest_10mm := 0
	var streak_sustained := 0
	var longest_sustained := 0

	# Sleep/wake bookkeeping. Every body here has can_sleep=false, so the solver should never
	# park one: a non-zero `sleeping` count means a body is asleep while still being solved,
	# which skips its velocity update; a non-zero `frozen_manifolds` count means contacts are
	# being left deliberately stale; a rising `rebuilds` count means warm-started contact
	# state is being thrown away mid-run.
	var has_sleep_stats := PhysicsServer3D.has_method("get_solver_sleep_stats")
	var max_sleeping := 0
	var first_sleeping_frame := -1
	var max_frozen := 0
	var min_update_set := 1 << 30
	var max_still_frames := 0
	var sleep_frames_setting := 0
	var sleep_mode_never_min := 1 << 30
	var movable_last := 0
	var rebuilds_start := 0
	var rebuilds_end := 0
	var rebuild_frame := -1
	var rebuild_step := -1
	var rebuild_reason := ""
	var rebuild_reasons := {}
	# 0.25 s timeline: the pile is two pyramids, and the upper one lands on the lower about
	# 0.55 s in. Whether the overlap grows at the impact or afterwards separates "the solver
	# cannot carry the load" from "the impact destroys it".
	var timeline := []
	var timeline_every := maxi(1, int(0.25 * float(Engine.physics_ticks_per_second)))
	if has_sleep_stats:
		var s0: Dictionary = PhysicsServer3D.call("get_solver_sleep_stats")
		rebuilds_start = int(s0.get("rebuilds", 0))

	var t0 := Time.get_ticks_msec()
	for frame in ticks:
		await physics_frame

		# One call, not two: the solver runs on its own thread, so a separate
		# get_solver_max_penetration() and get_solver_deepest_contact() can land on
		# different ticks and report a depth with the wrong pair attached. The deepest
		# contact returns [body_a, body_b, depth] from a single traversal.
		var deep: Array = PhysicsServer3D.call("get_solver_deepest_contact")
		var frame_pen: float = float(deep[2]) if deep.size() == 3 else 0.0
		# Both ends of the deepest pair must be slow for the depth to describe a resting
		# contact; a static body (the slab) counts as zero.
		var pair_speed := 0.0
		if deep.size() == 3:
			pair_speed = maxf(_body_speed(int(deep[0])), _body_speed(int(deep[1])))
		var pair_resting := pair_speed < RESTING_PAIR_SPEED
		if has_sleep_stats and frame > 0:
			# Frame 0 is sampled before the first server step, so the solver's per-step
			# fields (sleepFrames, inUpdateSet) still hold their defaults there.
			var st: Dictionary = PhysicsServer3D.call("get_solver_sleep_stats")
			var asleep := int(st.get("sleeping", 0))
			if asleep > max_sleeping:
				max_sleeping = asleep
			if asleep > 0 and first_sleeping_frame < 0:
				first_sleeping_frame = frame
			max_frozen = maxi(max_frozen, int(st.get("frozen_manifolds", 0)))
			min_update_set = mini(min_update_set, int(st.get("update_set", 0)))
			max_still_frames = maxi(max_still_frames, int(st.get("max_still_frames", 0)))
			sleep_mode_never_min = mini(sleep_mode_never_min, int(st.get("sleep_mode_never", 0)))
			movable_last = int(st.get("movable", 0))
			rebuilds_end = int(st.get("rebuilds", 0))
			if frame == 1:
				# The scene's own setup (mass/can_sleep writes above) marks the space dirty;
				# the first step consumes that and rebuilds once. The invariant below is
				# about rebuilds DURING the run, which would discard warm-started contacts.
				rebuilds_start = rebuilds_end
				# The solver stamps its per-tick fields during the step, so the window is
				# only readable after one has run (frame 0 still shows the defaults).
				sleep_frames_setting = int(st.get("sleep_frames", 0))
			if rebuilds_end > rebuilds_start and rebuild_frame < 0:
				rebuild_frame = frame
				rebuild_step = int(st.get("rebuild_last_step", -1))
				rebuild_reason = String(st.get("rebuild_last_reason", ""))
			rebuild_reasons = st.get("rebuild_reasons", {})
		if frame_pen > 0.05 and first_above_5cm < 0.0:
			first_above_5cm = float(frame) / hz
		if frame < warmup_frames:
			if frame_pen > worst_warmup:
				worst_warmup = frame_pen
				worst_warmup_frame = frame
				worst_warmup_a = int(deep[0])
				worst_warmup_b = int(deep[1])
		else:
			if frame_pen > worst_any_settled:
				worst_any_settled = frame_pen
				worst_any_settled_frame = frame
			# An impact transient on a moving pair is evidence only (worst_any_settled): not a
			# sample and not a streak frame; the warmup already covers the authored landing.
			if pair_resting:
				samples.append(frame_pen)
				if frame_pen > worst_settled:
					worst_settled = frame_pen
					worst_settled_frame = frame
					if deep.size() == 3:
						worst_settled_a = int(deep[0])
						worst_settled_b = int(deep[1])
						var obj_a := instance_from_id(worst_settled_a) as Node3D
						var obj_b := instance_from_id(worst_settled_b) as Node3D
						if obj_a != null:
							worst_a_pos = obj_a.global_position
						if obj_b != null:
							worst_b_pos = obj_b.global_position
						worst_pair_boxes = obj_a is RigidBody3D and obj_b is RigidBody3D
						# The reported depth is -dot(n, xA + rA - xB - rB) over the manifold's
						# cached frame and anchors. Measure what the two boxes actually overlap
						# at the same instant, so a depth that no longer matches the geometry is
						# visible instead of being trusted.
						if worst_pair_boxes:
							worst_true_depth = _obb_penetration(obj_a.global_transform,
									obj_b.global_transform, BOX_HALF_EXTENTS, BOX_HALF_EXTENTS)

				above_1mm += 1 if frame_pen > 0.001 else 0
				above_2mm += 1 if frame_pen > 0.002 else 0
				above_5mm += 1 if frame_pen > 0.005 else 0
				above_10mm += 1 if frame_pen > 0.01 else 0
				above_20mm += 1 if frame_pen > 0.02 else 0

				streak_1mm = streak_1mm + 1 if frame_pen > 0.001 else 0
				streak_2mm = streak_2mm + 1 if frame_pen > 0.002 else 0
				streak_5mm = streak_5mm + 1 if frame_pen > 0.005 else 0
				streak_10mm = streak_10mm + 1 if frame_pen > 0.01 else 0
				longest_1mm = maxi(longest_1mm, streak_1mm)
				longest_2mm = maxi(longest_2mm, streak_2mm)
				longest_5mm = maxi(longest_5mm, streak_5mm)
				longest_10mm = maxi(longest_10mm, streak_10mm)
				streak_sustained = streak_sustained + 1 if frame_pen > STRICT_SUSTAINED_LIMIT else 0
				longest_sustained = maxi(longest_sustained, streak_sustained)

		var frame_max_speed := 0.0
		var frame_max_shift := 0.0
		for i in bodies.size():
			var body := bodies[i]
			var pos := body.global_position
			var speed := body.linear_velocity.length()
			frame_max_speed = maxf(frame_max_speed, speed)
			frame_max_shift = maxf(frame_max_shift, pos.distance_to(spawn_pos[i]))
			global_min_y = minf(global_min_y, pos.y)
			global_max_speed = maxf(global_max_speed, speed)
			if frame < warmup_frames:
				continue
			if pos.y < min_y:
				min_y = pos.y
				min_y_name = str(body.name)
			max_abs_x = maxf(max_abs_x, absf(pos.x))
			max_abs_z = maxf(max_abs_z, absf(pos.z))
			if absf(pos.x) <= GROUND_HALF_X and absf(pos.z) <= ground_half_z:
				if pos.y < min_inside_y_any:
					min_inside_y_any = pos.y
					min_inside_name_any = str(body.name)
				# Same settled gate as the contact depth: a body that is still moving is
				# passing through on an impact, not resting below the slab surface.
				if speed < RESTING_PAIR_SPEED and pos.y < min_inside_y:
					min_inside_y = pos.y
					min_inside_name = str(body.name)
			settled_max_speed = maxf(settled_max_speed, speed)
			finite = finite and is_finite_body(body)

		if frame % timeline_every == 0:
			timeline.append({"t": float(frame) / hz, "pen": frame_pen,
					"max_speed": frame_max_speed, "max_shift": frame_max_shift})

	# Boxes that leave the authored slab footprint are falling, not sinking through it;
	# report them separately so the penetration bound measures what it claims to.
	var bodies_off_slab := 0
	var bodies_below_slab := 0
	for body in bodies:
		var pos := body.global_position
		if absf(pos.x) > GROUND_HALF_X or absf(pos.z) > ground_half_z:
			bodies_off_slab += 1
		if pos.y < -1.5:
			bodies_below_slab += 1

	var elapsed := (Time.get_ticks_msec() - t0) / 1000.0
	var settled_frames := samples.size()
	var sorted := samples.duplicate()
	sorted.sort()
	var mean := 0.0
	for value in samples:
		mean += value
	if settled_frames > 0:
		mean /= float(settled_frames)
	var p50 := _percentile(sorted, 0.50)
	var p95 := _percentile(sorted, 0.95)
	var p99 := _percentile(sorted, 0.99)
	var p999 := _percentile(sorted, 0.999)
	var wanted_center := GROUND_TOP + FLAT_BOX_HALF_HEIGHT
	var ground_sink := wanted_center - min_y
	var ground_sink_on_slab := wanted_center - min_inside_y
	# Independent check of the reported deepest pair: the authored box is 1 x 0.5 x 1.
	# Ignoring rotation, this is the overlap the captured positions can actually have.
	var worst_offset := worst_a_pos - worst_b_pos
	var worst_pair_distance := worst_offset.length()
	var worst_pair_aabb_depth := 0.0
	if worst_pair_boxes:
		var ox := 1.0 - absf(worst_offset.x)
		var oy := 0.5 - absf(worst_offset.y)
		var oz := 1.0 - absf(worst_offset.z)
		if ox > 0.0 and oy > 0.0 and oz > 0.0:
			worst_pair_aabb_depth = minf(ox, minf(oy, oz))

	check(finite, "pile stays finite", "bodies=%d" % bodies.size())
	# Independent of the contact diagnostic: a flat box resting on the slab has its centre
	# at 0.25 m. Only boxes still inside the slab footprint count as sinking; a box that
	# slid past the slab edge (visible in max_abs_x/max_abs_z) is in free fall.
	check(min_inside_y >= wanted_center - GROUND_SINK_LIMIT, "no deep ground sink (on slab)",
			"min centre y=%.4f m (%s) at rest; any-speed min=%.4f m (%s); global=%.1f m, off-slab=%d" % [
				min_inside_y, min_inside_name, min_inside_y_any, min_inside_name_any,
				min_y, bodies_off_slab])
	check(worst_settled <= STRICT_WORST_LIMIT, "settled penetration <= %.0f cm" % (STRICT_WORST_LIMIT * 100.0),
			("worst=%.4f m at rest, t=%.1f s (%s <-> %s); " +
			"any-speed post-warmup max=%.4f m; true SAT overlap of that pair=%.4f m") % [
				worst_settled, float(worst_settled_frame) / hz,
				_describe_instance(worst_settled_a), _describe_instance(worst_settled_b),
				worst_any_settled, worst_true_depth])
	check(p99 <= STRICT_P99_LIMIT, "99th percentile penetration <= %.1f cm" % (STRICT_P99_LIMIT * 100.0),
			"p99=%.4f m (p95=%.4f, p50=%.4f)" % [p99, p95, p50])
	check(longest_sustained < STRICT_SUSTAINED_FRAMES, "no sustained %.0f mm overlap" % (STRICT_SUSTAINED_LIMIT * 1000.0),
			"longest run above %.0f mm = %d engine frames (limit %d); longest above 1 cm = %d" % [
				STRICT_SUSTAINED_LIMIT * 1000.0, longest_sustained,
				STRICT_SUSTAINED_FRAMES, longest_10mm])

	# Sleep/wake invariants. By default every box is created with can_sleep=false, so the
	# solver must see sleep_mode NEVER for all of them, must never park one, and must keep
	# the whole pile in the update set. `sleep-on` flips can_sleep=true instead, so the same
	# checks invert: the pile must actually park once the impact has settled. In both cases
	# the space must not rebuild mid-run (a rebuild destroys every manifold, discarding
	# warm-started contact stiffness).
	var rebuilds_during_run := rebuilds_end - rebuilds_start
	if has_sleep_stats:
		if _disable_sleep:
			check(sleep_mode_never_min == movable_last and movable_last == bodies.size(),
					"can_sleep=false reached the solver",
					"sleep_mode NEVER=%d of movable=%d (bodies=%d)" % [
						sleep_mode_never_min, movable_last, bodies.size()])
			check(max_sleeping == 0, "no body fell asleep",
					"max sleeping=%d (first at frame %d), max stillFrames=%d, sleepFrames=%d" % [
						max_sleeping, first_sleeping_frame, max_still_frames, sleep_frames_setting])
			check(max_frozen == 0, "no manifold frozen",
					"max frozen manifolds=%d" % max_frozen)
			check(min_update_set == bodies.size(), "whole pile in update set",
					"min update set=%d of %d (max sleeping=%d, max frozen=%d)" % [
						min_update_set, bodies.size(), max_sleeping, max_frozen])
		else:
			check(max_sleeping > 0, "settled pile sleeps",
					"max sleeping=%d (first at frame %d), max stillFrames=%d, sleepFrames=%d" % [
						max_sleeping, first_sleeping_frame, max_still_frames, sleep_frames_setting])
			check(sleep_frames_setting > 0, "sleep window exposed",
					"sleepFrames=%d at dt=%.5f s" % [sleep_frames_setting, 1.0 / hz])
		check(rebuilds_during_run == 0, "no space rebuild during the run",
				"rebuilds=%d at frame %d from %s; all reasons=%s" % [
					rebuilds_during_run, rebuild_frame, rebuild_reason, str(rebuild_reasons)])

	var evidence := {
		"config": _sweep_config if not _sweep_config.is_empty() else "current",
		"upper": "drop" if _drop_upper else "keep",
		"body_count": bodies.size(),
		"seconds": _run_seconds,
		"warmup_seconds": STRICT_WARMUP_SECONDS,
		"settled_samples": settled_frames,
		"wall_seconds": elapsed,
		"first_above_5cm": first_above_5cm,
		"global_min_y": global_min_y,
		"global_max_speed": global_max_speed,
		"settled_max_speed": settled_max_speed,
		"worst_warmup": worst_warmup,
		"worst_warmup_time": float(worst_warmup_frame) / hz,
		"worst_warmup_pair": "%s <-> %s" % [_describe_instance(worst_warmup_a), _describe_instance(worst_warmup_b)],
		"worst_settled": worst_settled,
		"worst_settled_time": float(worst_settled_frame) / hz,
		"worst_settled_pair": "%s <-> %s" % [_describe_instance(worst_settled_a), _describe_instance(worst_settled_b)],
		"worst_any_settled": worst_any_settled,
		"worst_any_settled_time": float(worst_any_settled_frame) / hz,
		"resting_pair_speed": RESTING_PAIR_SPEED,
		"worst_pair_boxes": worst_pair_boxes,
		"worst_pair_center_distance": worst_pair_distance,
		"worst_pair_aabb_depth": worst_pair_aabb_depth,
		"worst_pair_true_sat_depth": worst_true_depth,
		"worst_pair_pos_a": worst_a_pos,
		"worst_pair_pos_b": worst_b_pos,
		"mean": mean,
		"p50": p50,
		"p95": p95,
		"p99": p99,
		"p999": p999,
		"frames_above_1mm": above_1mm,
		"frames_above_2mm": above_2mm,
		"frames_above_5mm": above_5mm,
		"frames_above_10mm": above_10mm,
		"frames_above_20mm": above_20mm,
		"longest_run_above_1mm": longest_1mm,
		"longest_run_above_2mm": longest_2mm,
		"longest_run_above_5mm": longest_5mm,
		"longest_run_above_10mm": longest_10mm,
		"longest_run_above_sustained": longest_sustained,
		"sustained_limit": STRICT_SUSTAINED_LIMIT,
		"min_center_y": min_y,
		"min_center_body": min_y_name,
		"ground_sink": ground_sink,
		"ground_z": _ground_z,
		"min_center_y_on_slab": min_inside_y,
		"min_center_body_on_slab": min_inside_name,
		"min_center_y_on_slab_any": min_inside_y_any,
		"min_center_body_on_slab_any": min_inside_name_any,
		"ground_sink_on_slab": ground_sink_on_slab,
		"bodies_off_slab": bodies_off_slab,
		"bodies_below_slab": bodies_below_slab,
		"max_abs_x": max_abs_x,
		"max_abs_z": max_abs_z,
		"sleep_stats": {
			"has_diagnostic": has_sleep_stats,
			"movable": movable_last,
			"sleep_mode_never_min": sleep_mode_never_min,
			"max_sleeping": max_sleeping,
			"first_sleeping_frame": first_sleeping_frame,
			"max_frozen_manifolds": max_frozen,
			"min_update_set": min_update_set if has_sleep_stats else 0,
			"max_still_frames": max_still_frames,
			"sleep_frames": sleep_frames_setting,
			"rebuilds_during_run": rebuilds_end - rebuilds_start,
			"rebuild_frame": rebuild_frame,
			"rebuild_reason": rebuild_reason,
			"rebuild_reasons": rebuild_reasons,
		},
		"timeline": timeline,
	}
	scene.queue_free()
	await frames(2)
	return evidence
