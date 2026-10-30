# Bridge regression: the jointed chain that the plain solver used to lose.
#
# The core scene `Bridge` (tools/core_scenes.cpp) is the same arrangement - 40 planks on a double
# row of ball joints, with boxes dropped on the deck - and `avbd_core_test --scene Bridge` is the
# faster gate for solver work. This script exists for the other half: it runs the *real* scene
# through the GDExtension and the physics server, so it also covers the project settings that
# reach the solver (iterations, relaxations, projection) rather than the core test's env knobs.
#
#   godot --headless --path demo --script res://tests/check_bridge.gd
#   AVBD_BRIDGE_SCENE=res://scenes/bridge2.tscn godot --headless --path demo \
#       --script res://tests/check_bridge.gd
#
# It reports the worst |position| reached by any rigid body, and lists the bodies that left the
# window: a chain that holds keeps its bodies at the geometry it was built from, and a chain that
# is thrown apart puts them in the hundreds or thousands of metres. Override the window with
# AVBD_BRIDGE_LIMIT and the scene with AVBD_BRIDGE_SCENE.
extends SceneTree

const DEFAULT_SCENE := "res://scenes/bridge2.tscn"
const STEPS := 900                 # 15 s at 60 fps, 3.75 s at the demo's 240 Hz tick
# Limit on |position|, in metres from the origin. The demo scenes hang their bridges within a few
# metres; the stress case is the core scene (`avbd_core_test --scene Bridge`, 40 planks spanning
# 40 m), which is why this is an override rather than a constant.
const DEFAULT_LIMIT := 40.0

func _initialize() -> void:
	var path := DEFAULT_SCENE
	var env_scene := OS.get_environment("AVBD_BRIDGE_SCENE")
	if env_scene != "":
		path = env_scene
	var limit := DEFAULT_LIMIT
	var env_limit := OS.get_environment("AVBD_BRIDGE_LIMIT")
	if env_limit != "":
		limit = float(env_limit)

	var packed := load(path)
	if packed == null:
		printerr("BRIDGE_MISSING: ", path)
		quit()
		return
	root.add_child(packed.instantiate())
	await process_frame
	await physics_frame

	var bodies := _bodies()
	print("bridge: scene=", path, " bodies=", bodies.size())

	var worst := 0.0
	for i in STEPS:
		await physics_frame
		var now := _worst_extent(bodies)
		if now > worst:
			worst = now

	# Same measure the core scene is judged by: the chain's own reach is |x| <= 20 plus the deck's
	# width in y, so anything past limit is a body that left the bridge.
	var escaped := 0
	for body in bodies:
		if body is Node3D and body.global_position.length() > limit:
			escaped += 1
			print("  escaped: ", body.name, " at ", body.global_position)

	print("bridge: steps=", STEPS, " worst|p|=", worst, " escaped=", escaped,
			" verdict=", "HELD" if worst <= limit and escaped == 0 else "LOST")
	quit(0 if worst <= limit and escaped == 0 else 1)

func _bodies() -> Array:
	var out := []
	var stack := [root]
	while not stack.is_empty():
		var node = stack.pop_back()
		if node is RigidBody3D:
			out.append(node)
		stack.append_array(node.get_children())
	return out

func _worst_extent(bodies: Array) -> float:
	var worst := 0.0
	for body in bodies:
		if body is Node3D:
			worst = maxf(worst, body.global_position.length())
	return worst
