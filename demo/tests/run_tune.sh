#!/bin/bash
cd /home/cyan/Documents/godot-avbd/demo
godot --headless --path . --script res://tests/test_load.gd 2>&1 | tail -15
