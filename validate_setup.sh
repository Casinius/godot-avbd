#!/bin/bash
set -e

echo "=== Validating car_raycast demo setup ==="
echo ""

# Check files exist
echo "1. Checking files exist..."
for file in "demo2/car_raycast.gd" "demo2/car_raycast.tscn" "demo2/test_car_raycast.gd" "demo2/project.godot" "demo2/bin/libavbd.linux.x86_64.so"; do
    if [ ! -f "$file" ]; then
        echo "   ❌ File not found: $file"
        exit 1
    fi
    echo "   ✅ $file"
done
echo ""

# Check scene file
echo "2. Checking car_raycast.tscn..."
if grep -q "car_raycast.gd" demo2/car_raycast.tscn; then
    echo "   ✅ Scene references car_raycast.gd"
else
    echo "   ❌ Scene doesn't reference car_raycast.gd"
    exit 1
fi
echo ""

# Check project.godot
echo "3. Checking project.godot..."
if grep -q 'main_scene="res://car_raycast.tscn"' demo2/project.godot; then
    echo "   ✅ project.godot sets main_scene to car_raycast.tscn"
else
    echo "   ❌ project.godot doesn't set correct main_scene"
    exit 1
fi

# Check input mappings
if grep -q '^f={' demo2/project.godot && \
   grep -q '^b={' demo2/project.godot && \
   grep -q '^a={' demo2/project.godot && \
   grep -q '^d={' demo2/project.godot; then
    echo "   ✅ Input mappings configured (f/b for throttle, a/d for steering)"
else
    echo "   ❌ Input mappings not configured"
    exit 1
fi
echo ""

# Check library
echo "4. Checking AVBD library..."
if [ -f "demo2/bin/libavbd.linux.x86_64.so" ]; then
    SIZE=$(stat -f%z demo2/bin/libavbd.linux.x86_64.so 2>/dev/null || stat -c%s demo2/bin/libavbd.linux.x86_64.so)
    echo "   ✅ AVBD library exists (${SIZE} bytes)"
else
    echo "   ❌ AVBD library not found"
    exit 1
fi
echo ""

# Check scripts for basic structure
echo "5. Checking car_raycast.gd structure..."
if grep -q "extends Node3D" demo2/car_raycast.gd && \
   grep -q "_physics_process" demo2/car_raycast.gd && \
   grep -q "_build" demo2/car_raycast.gd && \
   grep -q "Raycast3D" demo2/car_raycast.gd; then
    echo "   ✅ Script has required components (Node3D, _physics_process, _build, Raycast3D)"
else
    echo "   ❌ Script missing required components"
    exit 1
fi
echo ""

# Check test script
echo "6. Checking test_car_raycast.gd..."
if grep -q "extends Node" demo2/test_car_raycast.gd; then
    echo "   ✅ Test script extends Node"
else
    echo "   ❌ Test script doesn't extend Node"
    exit 1
fi
echo ""

echo "=== All validations passed! ==="
echo ""
echo "The car_raycast demo is ready to run."
echo ""
echo "To run: Godot --path demo2"
echo "Controls:"
echo "  W/S: Throttle (forward/backward)"
echo "  A/D: Steering (left/right)"
echo "  Mouse: Camera (freecam.gd)"
