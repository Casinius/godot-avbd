#!/usr/bin/env sh
# Side-by-side pulse benchmark: identical scene and identical pulse drive, one headless
# run per physics engine, so AVBD can be measured against Jolt (or Godot Physics).
#
#   ./scripts/pulse_bench.sh -- --ticks=60,120,240,480
#   ./scripts/pulse_bench.sh -- --pattern=push --grid=23
#   ENGINES="AVBD,Godot Physics,Jolt Physics" ./scripts/pulse_bench.sh -- --no-probe
#
# Arguments after `--` are forwarded to the bench script unchanged; ENGINES selects which
# engines to run (default: AVBD,Jolt Physics).
#
# The engine is chosen by demo/project.godot's physics/3d/physics_engine. This script
# flips that line for each run and restores the file on exit (including on Ctrl-C), so
# the project is left exactly as it was found.
#
# AVBD needs its GDExtension built first (xmake f -m release && xmake) and the project
# imported once (godot --headless --path demo --import) so Godot sees demo/bin/*.so.
set -eu

demo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
project="$demo/project.godot"
godot=${GODOT:-godot}
engines=${ENGINES:-"AVBD,Jolt Physics"}

backup=$(mktemp)
cp "$project" "$backup"
restore() {
    cp "$backup" "$project"
    rm -f "$backup"
}
trap restore EXIT INT TERM

set_engine() {
    if grep -q '^3d/physics_engine=' "$project"; then
        sed -i "s|^3d/physics_engine=.*|3d/physics_engine=\"$1\"|" "$project"
    else
        sed -i "/^\[physics\]/a 3d/physics_engine=\"$1\"" "$project"
    fi
}

old_ifs=$IFS
IFS=','
for engine in $engines; do
    IFS=$old_ifs
    set_engine "$engine"
    echo "=== $engine ==="
    "$godot" --headless --path "$demo" --script res://tests/test_pulse_bench.gd -- "$@" || true
done
IFS=$old_ifs
