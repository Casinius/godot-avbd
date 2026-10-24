#!/usr/bin/env python3
"""Regression gate for the bridge demo scenes.

Runs bridge / bridge2 / bridge3 headless with the AVBD engine for a fixed number of
physics steps, collects the pose digest each scene emits (demo.gd --digest), and
diffs it against a recorded baseline. The simulation is bit-deterministic, so an
equal digest means identical physics; a mismatch is a regression (or an intentional
numerics change - re-record with --update and read the diff first).

Usage:
  python3 tools/bridge_regression.py                    # gate against the baseline
  python3 tools/bridge_regression.py --update           # record a new baseline
  python3 tools/bridge_regression.py --adaptive-beta    # gate with adaptive beta on
  python3 tools/bridge_regression.py --frames 600       # longer settle
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys

SCENES = ["bridge", "bridge2", "bridge3"]
DIGEST_RE = re.compile(r"AVBD_SCENE_DIGEST (\S+)=0x([0-9a-fA-F]+)")


def demo_dir():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "demo"))


def write_override(adaptive_beta):
    # Project settings are read per run; the file wins over project.godot without
    # touching it (the repo's project.godot may deliberately point at Jolt).
    path = os.path.join(demo_dir(), "override.cfg")
    with open(path, "w") as f:
        f.write('[physics]\n\n3d/physics_engine="AVBD"\n')
        if adaptive_beta:
            f.write("avbd/adaptive_beta=true\n")
    return path


def run_scene(godot, scene, frames):
    cmd = [
        godot, "--headless", "--fixed-fps", "60",
        "--quit-after", str(frames + 60),  # safety net; the digest quits earlier
        "--path", demo_dir(), f"res://scenes/{scene}.tscn",
        "++", "--digest", f"--steps={frames}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    match = DIGEST_RE.search(result.stdout)
    if match is None:
        tail = "\n".join((result.stdout + result.stderr).splitlines()[-8:])
        return None, tail
    return f"0x{match.group(2).lower()}", ""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--godot", default=shutil.which("godot") or "godot")
    parser.add_argument("--frames", type=int, default=300, help="physics steps per scene")
    parser.add_argument("--adaptive-beta", action="store_true", help="run with physics/avbd/adaptive_beta on")
    parser.add_argument("--update", action="store_true", help="record the digests as the new baseline")
    parser.add_argument("--baseline", default=os.path.join(os.path.dirname(__file__), "bridge_baseline.json"))
    args = parser.parse_args()

    baseline = {}
    if os.path.exists(args.baseline):
        with open(args.baseline) as f:
            baseline = json.load(f)

    override = write_override(args.adaptive_beta)
    digests = {}
    failures = []
    try:
        for scene in SCENES:
            digest, tail = run_scene(args.godot, scene, args.frames)
            key = scene + ("@adaptive-beta" if args.adaptive_beta else "")
            digests[key] = digest
            label = f"{scene} ({args.frames} steps{' +adaptive-beta' if args.adaptive_beta else ''})"
            if digest is None:
                failures.append(scene)
                print(f"FAIL {label}: no digest emitted\n{tail}")
            elif key in baseline and not args.update and baseline[key] != digest:
                failures.append(scene)
                print(f"FAIL {label}: {digest} != baseline {baseline[key]}")
            else:
                print(f"ok  {label}: {digest}")
    finally:
        os.remove(override)

    if args.update:
        with open(args.baseline, "w") as f:
            json.dump(digests, f, indent=1, sort_keys=True)
            f.write("\n")
        print(f"baseline written: {args.baseline}")
        return 0

    print(f"bridge regression: {'FAIL (' + str(len(failures)) + ')' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
