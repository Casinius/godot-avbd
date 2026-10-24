#!/usr/bin/env python3
"""Compare per-frame physics cost between engines on the same demo scenes.

For each scene x engine it runs Godot headless with a fixed frame rate (so every run
executes the identical frame and physics-tick count, independent of machine speed),
times the run, and subtracts the same measurement taken on a generated empty scene as
the control - what remains is the per-frame physics cost of the workload on that
engine. AVBD and the engine(s) named on the command line are compared side by side.

Usage:
  python3 tools/compare_engines.py
  python3 tools/compare_engines.py --scenes pyramid,bridge2 --engines "AVBD,GodotPhysics3D"
  python3 tools/compare_engines.py --frames 600
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

EMPTY_SCENE = """[gd_scene format=3]

[node name="BenchEmpty" type="Node3D"]
"""


def demo_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "demo"))


def write_override(engine):
    path = os.path.join(demo_dir(), "override.cfg")
    with open(path, "w") as f:
        f.write(f'[physics]\n\n3d/physics_engine="{engine}"\n')
    return path


def timed_run(godot, scene_path, frames):
    cmd = [
        godot, "--headless", "--fixed-fps", "60",
        "--quit-after", str(frames),
        "--path", demo_dir(), scene_path,
    ]
    start = time.monotonic()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    wall = time.monotonic() - start
    if result.returncode != 0:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-5:])
        print(f"  warning: nonzero exit, tail:\n{tail}", file=sys.stderr)
    return wall


def measure(godot, scenes, engines, frames, warmup):
    control_path = os.path.join(demo_dir(), "_bench_empty.tscn")
    with open(control_path, "w") as f:
        f.write(EMPTY_SCENE)

    results = {}
    try:
        for engine in engines:
            override = write_override(engine)
            try:
                for scene in scenes:
                    path = control_path if scene == "__control__" else f"res://scenes/{scene}.tscn"
                    if not os.path.exists(control_path) and scene != "__control__":
                        pass
                    if warmup:
                        timed_run(godot, path, frames)
                    wall = timed_run(godot, path, frames)
                    results[(engine, scene)] = wall
                    print(f"  {engine:<16} {scene:<12} {wall:7.2f} s wall")
            finally:
                os.remove(override)
    finally:
        if os.path.exists(control_path):
            os.remove(control_path)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--godot", default=shutil.which("godot") or "godot")
    parser.add_argument("--scenes", default="bridge,bridge2,bridge3,pyramid",
                        help="comma-separated scene names under demo/scenes (without .tscn)")
    parser.add_argument("--engines", default="AVBD,Jolt Physics",
                        help="comma-separated physics/3d/physics_engine values to compare")
    parser.add_argument("--frames", type=int, default=300, help="frames per measurement")
    parser.add_argument("--warmup", action="store_true", default=True)
    parser.add_argument("--no-warmup", dest="warmup", action="store_false")
    args = parser.parse_args()

    scenes = args.scenes.split(",")
    engines = args.engines.split(",")
    results = measure(args.godot, scenes, engines, args.frames, args.warmup)

    print(f"\nper-frame physics cost over {args.frames} fixed frames (wall minus control):")
    header = f"{'scene':<12}" + "".join(f"{e:>18}" for e in engines)
    print(header)
    worst_ratio = None
    for scene in scenes:
        row = f"{scene:<12}"
        costs = {}
        for engine in engines:
            scene_wall = results.get((engine, scene), 0.0)
            control_wall = results.get((engine, "__control__"), 0.0)
            cost_ms = max(0.0, (scene_wall - control_wall)) / args.frames * 1000.0
            costs[engine] = cost_ms
            row += f"{cost_ms:15.2f} ms"
        print(row)
        if len(engines) == 2 and min(costs.values()) > 0:
            ratio = max(costs.values()) / min(costs.values())
            slower = max(costs, key=costs.get)
            worst_ratio = ratio
            print(f"{'':<12}{'':>6}-> {slower} is {ratio:.2f}x slower on this scene")
    return 0


if __name__ == "__main__":
    sys.exit(main())
