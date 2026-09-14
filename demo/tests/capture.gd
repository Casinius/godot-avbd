# Renders a scene for a number of physics frames and saves a PNG, so demo scenes can
# be checked visually without a person at the keyboard. Needs a real rendering device
# (do not pass --headless).
#
#   Godot_v4.6.3 --path demo --script res://tests/capture.gd -- \
#       --scene=res://scenes/stack.tscn --frames=180 --out=/tmp/stack.png
extends SceneTree

var scene_path := ""
var out_path := "/tmp/avbd_capture.png"
var frames := 180


func _initialize() -> void:
	_run()


func _run() -> void:
	for arg in OS.get_cmdline_user_args():
		if arg.begins_with("--scene="):
			scene_path = arg.split("=")[1]
		elif arg.begins_with("--out="):
			out_path = arg.split("=")[1]
		elif arg.begins_with("--frames="):
			frames = int(arg.split("=")[1])

	var packed: PackedScene = load(scene_path)
	if packed == null:
		printerr("cannot load scene: ", scene_path)
		quit(1)
		return
	root.add_child(packed.instantiate())

	for i in frames:
		await physics_frame

	await RenderingServer.frame_post_draw
	var image := root.get_texture().get_image()
	var err := image.save_png(out_path)
	if err != OK:
		printerr("cannot save ", out_path, ": ", err)
		quit(1)
		return
	print("CAPTURED ", scene_path, " -> ", out_path, " after ", frames, " frames")
	quit(0)
