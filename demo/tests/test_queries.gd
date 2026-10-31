# Space-query behaviour for the AVBD physics server.
#
#   Godot --headless --path demo --script res://tests/test_queries.gd
#   ... -- --list | --only=ray_exclusions | --digest | --strict
#
# These scenarios pin the gaps the solver-only suites do not cover: query exclusion
# lists (RayCast3D.add_exception, ShapeCast3D.exclude_parent), Area3D participation,
# capsule query shape mapping, and swept shape casts. They are deliberately written
# against standard Godot nodes and the standard PhysicsDirectSpaceState3D API.
extends "res://tests/avbd_harness.gd"


func _scenarios() -> Array:
	return [
		scenario("ray_exclusions", test_ray_exclusions),
		scenario("shapecast_exclusions", test_shapecast_exclusions),
		scenario("shapecast_capsule", test_shapecast_capsule),
		scenario("area_queries", test_area_queries),
		scenario("direct_cast_motion_sweep", test_direct_cast_motion_sweep),
		scenario("shape_queries", test_shape_queries),
		scenario("motion_query", test_motion_query),
	]


func add_static_box(parent: Node, body_name: String, size: Vector3, position: Vector3) -> StaticBody3D:
	var body := StaticBody3D.new()
	body.name = body_name
	var collision := CollisionShape3D.new()
	collision.name = "Shape"
	var box := BoxShape3D.new()
	box.size = size
	collision.shape = box
	body.add_child(collision)
	body.position = position
	parent.add_child(body)
	return body


func query_params(shape: Shape3D, transform: Transform3D) -> PhysicsShapeQueryParameters3D:
	var params := PhysicsShapeQueryParameters3D.new()
	params.shape = shape
	params.transform = transform
	params.collision_mask = 0xFFFFFFFF
	params.collide_with_bodies = true
	params.collide_with_areas = false
	return params


# ---------------------------------------------------------------------------
# Exclusion lists must be honoured by both the RayCast3D node and the direct
# space-state query. Without it a ray under a character always hits its parent.
# ---------------------------------------------------------------------------
func test_ray_exclusions() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))
	var parent := add_static_box(root, "Parent", Vector3.ONE, Vector3(4, 5, 0))
	await steps(60)

	var space_state := root.get_world_3d().direct_space_state
	var direct_query := PhysicsRayQueryParameters3D.create(
			Vector3(0, 20, 0), Vector3(0, -10, 0))
	direct_query.exclude = [box.get_rid()]
	var direct_hit := space_state.intersect_ray(direct_query)
	check(not direct_hit.is_empty(), "direct ray still hits after exclusion",
			"excluded box at y=3, ground below")
	check(not direct_hit.is_empty() and instance_from_id(direct_hit.get("collider_id")) == ground,
			"direct ray exclusion skips box", "hit=%s" % [direct_hit.get("collider_id")])

	var ray := RayCast3D.new()
	ray.name = "ExcludedRay"
	ray.target_position = Vector3(0, -30, 0)
	ray.add_exception(box)
	root.add_child(ray)
	ray.global_position = Vector3(0, 20, 0)
	ray.force_raycast_update()
	check(ray.get_collider() == ground, "RayCast3D.add_exception skips box",
			"collider=%s" % ray.get_collider())

	var child_ray := RayCast3D.new()
	child_ray.name = "ParentRay"
	child_ray.target_position = Vector3(0, -10, 0)
	child_ray.exclude_parent = true
	parent.add_child(child_ray)
	child_ray.position = Vector3(0, 1.5, 0)
	child_ray.force_raycast_update()
	check(child_ray.get_collider() == null or child_ray.get_collider() != parent,
			"RayCast3D.exclude_parent skips parent", "collider=%s" % child_ray.get_collider())

	var evidence := {
		"direct_collider": direct_hit.get("collider_id"),
		"ray_collider": ray.get_collider(),
		"parent_collider": child_ray.get_collider(),
	}
	for node: Node in [ground, box, parent, ray, child_ray]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# ShapeCast3D uses exclusion lists to report multiple colliders one at a time.
# If the server ignores the list it repeats the same collider until max_results.
# ---------------------------------------------------------------------------
func test_shapecast_exclusions() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))
	await steps(60)

	var cast := ShapeCast3D.new()
	cast.name = "ExcludingCast"
	var sphere := SphereShape3D.new()
	sphere.radius = 0.25
	cast.shape = sphere
	cast.target_position = Vector3(0, -10, 0)
	cast.max_results = 4
	cast.add_exception(box)
	root.add_child(cast)
	cast.global_position = Vector3(0, 10, 0)
	cast.force_shapecast_update()

	check(cast.is_colliding(), "shapecast hits something after excluding box",
			"ground below")
	check(cast.get_collision_count() >= 1, "shapecast reports a collision",
			"count=%d" % cast.get_collision_count())
	var first := cast.get_collider(0) if cast.get_collision_count() > 0 else null
	check(first == ground, "shapecast exclusion skips box", "collider=%s" % first)
	check(cast.get_collision_count() == 1, "shapecast does not repeat the collider",
			"count=%d max_results=%d" % [cast.get_collision_count(), cast.max_results])

	var parent := add_static_box(root, "CastParent", Vector3.ONE, Vector3(5, 5, 0))
	var child_cast := ShapeCast3D.new()
	child_cast.name = "ParentCast"
	var child_sphere := SphereShape3D.new()
	child_sphere.radius = 0.25
	child_cast.shape = child_sphere
	child_cast.target_position = Vector3(0, -1, 0)
	child_cast.exclude_parent = true
	child_cast.max_results = 4
	parent.add_child(child_cast)
	child_cast.position = Vector3(0, 1.5, 0)
	child_cast.force_shapecast_update()
	check(child_cast.get_collider(0) != parent if child_cast.get_collision_count() > 0 else true,
			"ShapeCast3D.exclude_parent skips parent",
			"collider=%s count=%d" % [child_cast.get_collider(0) if child_cast.get_collision_count() > 0 else null,
					child_cast.get_collision_count()])

	var evidence := {
		"count": cast.get_collision_count(),
		"first": first,
		"parent_count": child_cast.get_collision_count(),
	}
	for node: Node in [ground, box, cast, parent, child_cast]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# A CapsuleShape3D is simulated as a cylinder proxy. Queries must use the same
# proxy dimensions; the old default branch treated it as a 1x1x1 box.
#
# Capsule: radius 0.5, full height 2.0 -> bottom is 1.0 below the centre.
# Ground top is 0.5. From y=10, moving -10 m: safe fraction ~ (10 - 1.0 - 0.5) / 10.
# ---------------------------------------------------------------------------
func test_shapecast_capsule() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	await steps(2)

	var space_state := root.get_world_3d().direct_space_state
	var capsule := CapsuleShape3D.new()
	capsule.radius = 0.5
	capsule.height = 2.0
	var params := query_params(capsule, Transform3D(Basis(), Vector3(5, 10, 0)))
	params.motion = Vector3(0, -10, 0)
	var motion := space_state.cast_motion(params)
	check(motion.size() == 2, "capsule cast_motion reports a hit",
			"result=%s" % [motion])
	var safe := motion[0] if motion.size() == 2 else -1.0
	check(safe > 0.80 and safe < 0.89, "capsule cast uses the capsule depth",
			"safe=%.4f (want ~0.85, box proxy would be ~0.90)" % safe)

	var cast := ShapeCast3D.new()
	cast.name = "CapsuleCast"
	cast.shape = capsule
	cast.target_position = Vector3(0, -10, 0)
	cast.max_results = 1
	root.add_child(cast)
	cast.global_position = Vector3(5, 10, 0)
	cast.force_shapecast_update()
	check(cast.is_colliding() and cast.get_collider(0) == ground,
			"ShapeCast3D capsule hits ground", "collider=%s" % cast.get_collider(0))

	var evidence := {"safe": safe, "collider": cast.get_collider(0)}
	ground.queue_free()
	cast.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Area3D must participate in ray/shape/rest queries when collide_with_areas is on.
# ---------------------------------------------------------------------------
func test_area_queries() -> Variant:
	var area := Area3D.new()
	area.name = "QueryArea"
	var collision := CollisionShape3D.new()
	var box := BoxShape3D.new()
	box.size = Vector3(2, 2, 2)
	collision.shape = box
	area.add_child(collision)
	area.collision_layer = 1
	root.add_child(area)
	area.global_position = Vector3(0, 3, 0)
	await frames(2)

	var space_state := root.get_world_3d().direct_space_state
	var ray_query := PhysicsRayQueryParameters3D.create(
			Vector3(0, 20, 0), Vector3(0, -10, 0))
	ray_query.collide_with_bodies = false
	ray_query.collide_with_areas = true
	var ray_hit := space_state.intersect_ray(ray_query)
	check(not ray_hit.is_empty() and ray_hit.get("rid") == area.get_rid(),
			"intersect_ray hits Area3D", "hit=%s" % [ray_hit])

	var shape_params := query_params(SphereShape3D.new(),
			Transform3D(Basis(), Vector3(0, 3, 0)))
	(shape_params.shape as SphereShape3D).radius = 0.2
	shape_params.collide_with_bodies = false
	shape_params.collide_with_areas = true
	var shape_hits := space_state.intersect_shape(shape_params, 8)
	check(shape_hits.size() == 1 and shape_hits[0].get("rid") == area.get_rid(),
			"intersect_shape hits Area3D", "hits=%s" % [shape_hits])

	var cast_params := query_params(SphereShape3D.new(),
			Transform3D(Basis(), Vector3(0, 10, 0)))
	(cast_params.shape as SphereShape3D).radius = 0.2
	cast_params.motion = Vector3(0, -10, 0)
	cast_params.collide_with_bodies = false
	cast_params.collide_with_areas = true
	var cast := space_state.cast_motion(cast_params)
	check(cast.size() == 2, "cast_motion hits Area3D", "result=%s" % [cast])

	var evidence := {"ray": ray_hit, "shape_count": shape_hits.size(), "cast": cast}
	area.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# cast_motion must sweep, not just test the final pose. A fast query crosses the
# box and ends below the ground; the old code only tested the endpoint and missed.
# ---------------------------------------------------------------------------
func test_direct_cast_motion_sweep() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))
	await steps(60)

	var space_state := root.get_world_3d().direct_space_state
	var sphere := SphereShape3D.new()
	sphere.radius = 0.2
	var params := query_params(sphere, Transform3D(Basis(), Vector3(0, 20, 0)))
	params.motion = Vector3(0, -30, 0)
	var motion := space_state.cast_motion(params)
	check(motion.size() == 2, "cast_motion sweeps past the endpoint",
			"result=%s (sphere at y=20 moving to y=-10)" % [motion])
	var safe := motion[0] if motion.size() == 2 else -1.0
	check(safe > 0.0 and safe < 1.0, "cast_motion safe fraction is inside the motion",
			"safe=%.4f" % safe)

	var evidence := {"safe": safe, "motion": motion}
	ground.queue_free()
	box.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# CharacterBody3D.move_and_slide / test_move use PhysicsServer3D.body_test_motion,
# not the direct space state. The server must report the wall's safe fraction and
# collider so the character stops instead of walking through it.
# ---------------------------------------------------------------------------
func test_motion_query() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var wall := add_static_box(root, "Wall", Vector3(1, 3, 4), Vector3(3, 1.5, 0))

	var body := CharacterBody3D.new()
	body.name = "Walker"
	var collision := CollisionShape3D.new()
	var box := BoxShape3D.new()
	box.size = Vector3(0.5, 0.5, 0.5)
	collision.shape = box
	body.add_child(collision)
	body.position = Vector3(1.5, 1.0, 0)
	root.add_child(body)
	await frames(2)

	var blocked := body.test_move(body.global_transform, Vector3(10, 0, 0))
	check(blocked, "test_move sees the wall", "motion=(10,0,0), wall left face x=2.5")

	body.velocity = Vector3(6, 0, 0)
	for i in 60:
		body.move_and_slide()
		await physics_frame
	check(body.global_position.x > 2.0 and body.global_position.x < 2.45,
			"move_and_slide stops at the wall",
			"x=%.4f (wall face 2.5, body half 0.25)" % body.global_position.x)
	check(absf(body.global_position.y - 1.0) < 0.2, "move_and_slide keeps the walker up",
			"y=%.4f" % body.global_position.y)

	var evidence := {"blocked": blocked, "x": body.global_position.x, "y": body.global_position.y}
	for node: Node in [ground, wall, body]:
		node.queue_free()
	await frames(2)
	return evidence


# ---------------------------------------------------------------------------
# Direct shape queries: intersect_point / intersect_shape / collide_shape /
# get_rest_info, including the query exclusion list.
# ---------------------------------------------------------------------------
func test_shape_queries() -> Variant:
	var ground := add_ground(root, 0.5, 0.5)
	var box := add_body(root, "Box", Vector3.ONE, Vector3(0, 3, 0))
	await steps(60)
	var space_state := root.get_world_3d().direct_space_state

	var point_params := PhysicsPointQueryParameters3D.new()
	point_params.position = box.global_position
	point_params.collide_with_bodies = true
	point_params.collide_with_areas = false
	var point_hits := space_state.intersect_point(point_params, 8)
	check(point_hits.size() >= 1 and point_hits[0].get("collider_id") == box.get_instance_id(),
			"intersect_point finds the box", "hits=%s" % [point_hits])

	var sphere := SphereShape3D.new()
	sphere.radius = 0.5
	# Slightly overlapping the box's top face: a realistic rest contact, unlike a shape
	# buried in the box's centre.
	var query_pos := box.global_position + Vector3(0, 0.45, 0)
	var params := query_params(sphere, Transform3D(Basis(), query_pos))
	var shape_hits := space_state.intersect_shape(params, 8)
	check(shape_hits.size() >= 1 and instance_from_id(shape_hits[0].get("collider_id")) == box,
			"intersect_shape finds the box", "hits=%s" % [shape_hits])

	var points := space_state.collide_shape(params, 32)
	check(points.size() >= 2 and points.size() % 2 == 0, "collide_shape returns point pairs",
			"points=%d" % points.size())

	var rest := space_state.get_rest_info(params)
	check(not rest.is_empty() and rest.get("rid") == box.get_rid(),
			"get_rest_info finds the box", "rest=%s" % [rest])
	check(rest.has("normal") and rest.get("normal").length() > 0.9,
			"get_rest_info reports a normal", "normal=%s" % [rest.get("normal")])

	params.exclude = [box.get_rid()]
	var excluded := space_state.intersect_shape(params, 8)
	var saw_box := false
	for hit in excluded:
		if hit.get("rid") == box.get_rid():
			saw_box = true
	check(not saw_box, "intersect_shape honours exclude", "hits=%s" % [excluded])
	var excluded_rest := space_state.get_rest_info(params)
	check(excluded_rest.is_empty() or excluded_rest.get("rid") != box.get_rid(),
			"get_rest_info honours exclude", "rest=%s" % [excluded_rest])

	var evidence := {
		"shape_count": shape_hits.size(),
		"point_count": points.size(),
		"rest_rid": rest.get("rid"),
		"excluded_count": excluded.size(),
	}
	ground.queue_free()
	box.queue_free()
	await frames(2)
	return evidence
