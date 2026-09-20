extends Node3D

@export var test_interval := 2.0
var test_timer := 0.0

func _ready() -> void:
	print("TestPunch: Starting click test")

func _process(delta: float) -> void:
	test_timer += delta
	if test_timer >= test_interval:
		test_timer = 0.0
		_send_click_test()

func _send_click_test() -> void:
	print("TestPunch: Simulating click on box at (0, 2, 0)")

	var demo_script = get_parent() as Node
	if demo_script == null or not demo_script.has_method("_punch"):
		print("TestPunch: Demo script not found or missing _punch method")
		return

	# Simulate click at center of screen
	var screen_pos := Vector2(0.5, 0.5)
	demo_script._punch(screen_pos)
