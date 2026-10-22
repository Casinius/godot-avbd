#!/usr/bin/env python3
"""Per-step physics benchmark: AVBD against Godot's Jolt at the same tick rate.

The measurement problem this solves: `demo/scripts/pulse_bench.sh` runs headless with the
engine's normal frame pacing, and Godot caps catch-up at max_physics_steps_per_frame. At
480 Hz with the cap at 4, the wall time per step is the frame budget (4.17 ms) no matter
how fast the engine is, and realtime_factor cannot exceed 0.5. Here every run gets

    common/max_physics_steps_per_frame=1
    --fixed-fps <tick>

so exactly one physics step happens per main-loop iteration and the loop runs as fast as
the machine allows. The bench script reports steps_per_frame; the driver rejects any
sample where it is not 1, which is what keeps a paced run from being reported as a result.

Everything else is deliberately identical across engines: same generated project, same
bodies, same pulses, same measurement code, same tick. Only physics/3d/physics_engine
changes.

Typical use:

    xmake f -m release && xmake
    python3 tools/pulse_bench.py --ticks 60,120,240,480,960 --grids 11,23,35 --reps 3

    # A/B two AVBD builds (same label, different binary)
    python3 tools/pulse_bench.py --label baseline \
        --avbd-library build/physics_bench/libavbd-baseline.so
    python3 tools/pulse_bench.py --label fixed --profile

Results land in build/physics_bench/ as results.json plus a markdown table.
"""
import argparse
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEMO = REPO / "demo"

# The bench script is self-contained (it builds its scene in code), so the clone needs the
# script, the extension and a generated project - nothing else from demo/.
SCRIPT = "tests/test_pulse_bench.gd"
EXTENSION = "gdavbd/avbd.gdextension"
LIBRARY = "gdavbd/bin/libavbd.linux.x86_64.so"


class BenchError(RuntimeError):
    pass


def project_text(args, ticks, engine):
    lines = [
        "config_version=5",
        "[application]",
        'config/name="AVBD pulse bench"',
        'config/features=PackedStringArray("4.7")',
        "[physics]",
        f'3d/physics_engine="{engine}"',
        "3d/run_on_separate_thread=false",
        "3d/default_gravity=9.8",
        "3d/default_linear_damp=0.0",
        "3d/default_angular_damp=0.0",
        "common/physics_ticks_per_second=%d" % ticks,
        "common/max_physics_steps_per_frame=1",
        "common/physics_interpolation=false",
        "avbd/substeps=1",
        "avbd/threads=%d" % args.threads,
        "avbd/profile=%s" % ("true" if args.profile else "false"),
        "avbd/convergence_threshold=%s" % args.convergence_threshold,
        "avbd/adaptive_iterations=%s" % ("true" if args.adaptive else "false"),
    ]
    if args.iterations > 0:
        lines.append("avbd/iterations=%d" % args.iterations)
    if args.jolt_velocity_steps > 0:
        lines.append("jolt_physics_3d/simulation/velocity_steps=%d" % args.jolt_velocity_steps)
        lines.append("jolt_physics_3d/simulation/position_steps=1")
    lines += [
        "[rendering]",
        'renderer/rendering_method="gl_compatibility"',
        "",
    ]
    return "\n".join(lines)


def prepare_clone(root, args, ticks, engine):
    clone = root / ("%s-%d" % (engine.replace(" ", "_"), ticks))
    (clone / "tests").mkdir(parents=True, exist_ok=True)
    (clone / "gdavbd" / "bin").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(DEMO / SCRIPT, clone / SCRIPT)
    shutil.copyfile(DEMO / EXTENSION, clone / EXTENSION)
    library = Path(args.avbd_library).resolve() if args.avbd_library else DEMO / LIBRARY
    if not library.is_file():
        raise BenchError("AVBD library not found: %s (build first with xmake)" % library)
    shutil.copyfile(library, clone / LIBRARY)
    (clone / "project.godot").write_text(project_text(args, ticks, engine), encoding="utf-8")
    return clone


def run(cmd, log, timeout):
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log.write_text(output, encoding="utf-8")
        raise BenchError("Godot timed out after %ss; log: %s" % (timeout, log)) from error
    except OSError as error:
        raise BenchError("cannot launch Godot: %s" % error) from error
    log.write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise BenchError("Godot exited %d; log: %s\n%s" %
                         (result.returncode, log, result.stdout[-3000:]))
    if "SCRIPT ERROR:" in result.stdout:
        raise BenchError("script error; log: %s\n%s" % (log, result.stdout[-3000:]))
    return result


def import_project(prefix, clone, log, timeout):
    # A headless import is what writes .godot/extension_list.cfg, and without that Godot
    # never loads the .gdextension (and AVBDPhysicsServer3D does not exist in script).
    # Godot 4.7 aborts at editor shutdown after writing the artifacts - with the stock
    # extension too - so success is judged by the artifacts, not the exit code.
    try:
        run(prefix + ["--import", "--quit"], log, timeout)
    except BenchError:
        pass
    if not (clone / ".godot" / "extension_list.cfg").is_file():
        raise BenchError("import did not produce .godot/extension_list.cfg; log: %s" % log)


def measure(args, engine, ticks, work):
    clone = prepare_clone(work, args, ticks, engine)
    prefix = [args.godot, "--headless", "--path", str(clone)]
    import_project(prefix, clone, work / ("%s-%d.import.log" % (engine, ticks)), args.timeout)

    trace = work / ("%s-%d.jsonl" % (engine, ticks))
    if trace.exists():
        trace.unlink()
    command = prefix + [
        "--fixed-fps", str(ticks),
        "--script", "res://" + SCRIPT,
        "--",
        "--ticks=%d" % ticks,
        "--grid=%d" % args.grid,
        "--steps=%d" % args.steps,
        "--warmup=%d" % args.warmup,
        "--pattern=%s" % args.pattern,
        "--pulse-hz=%s" % args.pulse_hz,
        "--reps=%d" % args.reps,
        "--no-probe",
        "--output=%s" % trace.resolve(),
    ]
    if args.threads >= 0:
        command.append("--threads=%d" % args.threads)
    if args.iterations > 0:
        command.append("--iterations=%d" % args.iterations)
    if args.profile:
        command.append("--profile")
    run(command, work / ("%s-%d.run.log" % (engine, ticks)), args.timeout)

    if not trace.is_file():
        raise BenchError("bench produced no trace: %s" % trace)
    samples = [json.loads(line) for line in trace.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not samples:
        raise BenchError("bench produced an empty trace: %s" % trace)
    return samples


def check_sample(sample, args, engine, label):
    where = "%s/%s/tick=%s" % (label, engine, sample.get("ticks"))
    if sample.get("pacing") != "ok":
        raise BenchError("%s: pacing is not one step per frame (steps_per_frame=%s); "
                         "wall time would be the frame budget, not the engine" %
                         (where, sample.get("steps_per_frame")))
    if sample.get("physics_steps") != args.steps:
        raise BenchError("%s: expected %d physics steps, got %s" % (where, args.steps, sample.get("physics_steps")))
    if not sample.get("finite", False):
        raise BenchError("%s: non-finite body state" % where)
    if engine == "AVBD" and args.profile and "solver_rounds" not in sample:
        raise BenchError("%s: profile requested but absent (AVBDPhysicsServer3D.get_instance "
                         "returned nothing)" % where)


def median(samples, key):
    values = [s[key] for s in samples if key in s]
    return statistics.median(values) if values else None


def summarize(args, label, engine, ticks, samples):
    for sample in samples:
        check_sample(sample, args, engine, label)
    row = {
        "label": label,
        "engine": engine,
        "ticks": ticks,
        "bodies": samples[0].get("bodies"),
        "pattern": samples[0].get("pattern"),
        "reps": len(samples),
        "budget_ms": samples[0].get("budget_ms"),
        "wall_ms_per_step": median(samples, "wall_ms_per_step"),
        "engine_ms_per_step": median(samples, "engine_ms_per_step"),
        "physics_ms_per_frame": median(samples, "physics_ms_per_frame"),
        "realtime_factor": median(samples, "realtime_factor"),
        "max_position": median(samples, "max_position"),
        "max_speed": median(samples, "max_speed"),
        "active_objects": samples[0].get("active_objects"),
        "collision_pairs": samples[0].get("collision_pairs"),
    }
    if "solver_ms_per_step" in samples[0]:
        row["solver_ms_per_step"] = median(samples, "solver_ms_per_step")
        for key in ["step_ms", "broad_phase_ms", "warmstart_forces_ms", "colour_graph_ms",
                    "warmstart_bodies_ms", "solve_ms", "finish_ms", "host_pre_ms", "host_sync_ms"]:
            row["solver_" + key] = median(samples, "solver_" + key)
        for key in ["bodies", "forces", "colours", "rounds", "pair_candidates", "substeps"]:
            row["solver_" + key] = samples[0].get("solver_" + key)
    row["spread"] = {}
    for key in ("wall_ms_per_step", "engine_ms_per_step", "solver_ms_per_step"):
        values = [s[key] for s in samples if key in s]
        if len(values) > 1:
            row["spread"][key] = (max(values) - min(values)) / statistics.median(values) if statistics.median(values) else None
    return row


def markdown(rows):
    header = ("| label | engine | tick | bodies | wall ms/step | solver ms/step | "
              "realtime | rounds | pairs | phases (broad/warm/colour/solve/finish) |")
    lines = [header, "|" + "---|" * 10]
    for row in rows:
        phases = ""
        if "solver_broad_phase_ms" in row:
            phases = "%.3f/%.3f/%.3f/%.3f/%.3f" % tuple(
                row.get("solver_" + k, 0.0) for k in
                ("broad_phase_ms", "warmstart_forces_ms", "colour_graph_ms", "solve_ms", "finish_ms"))
        solver = row.get("solver_ms_per_step")
        lines.append("| %s | %s | %d | %s | %.4f | %s | %.3f | %s | %s | %s |" % (
            row["label"], row["engine"], row["ticks"], row["bodies"],
            row["wall_ms_per_step"],
            "%.4f" % solver if solver is not None else "-",
            row["realtime_factor"],
            row.get("solver_rounds", "-"), row.get("solver_pair_candidates", "-"),
            phases or "-"))
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--godot", default="godot")
    parser.add_argument("--avbd-library", help="AVBD .so to install in the clone (A/B builds)")
    parser.add_argument("--engines", default="AVBD,Jolt Physics")
    parser.add_argument("--ticks", default="60,120,240,480,960")
    parser.add_argument("--grids", default="23", help="grid side(s); bodies = grid^2")
    parser.add_argument("--steps", type=int, default=240)
    parser.add_argument("--warmup", type=int, default=60)
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--pattern", choices=("rest", "spin", "push", "both"), default="both")
    parser.add_argument("--pulse-hz", type=float, default=20.0,
                        help="pulses per simulated second; keeps the drive identical across "
                             "tick rates instead of injecting 1/tick impulses")
    parser.add_argument("--threads", type=int, default=-1, help="AVBD worker count; -1 = project default")
    parser.add_argument("--iterations", type=int, default=0, help="0 = each engine's default")
    parser.add_argument("--jolt-velocity-steps", type=int, default=0,
                        help="0 = Jolt's own default; set with --iterations for a pinned run")
    parser.add_argument("--profile", action="store_true", help="collect AVBD phase timings")
    parser.add_argument("--adaptive", action="store_true", help="enable AVBD adaptive iteration count")
    parser.add_argument("--convergence-threshold", default="1e-5")
    parser.add_argument("--label", default="run", help="name for this build/config in the report")
    parser.add_argument("--output-dir", default=str(REPO / "build/physics_bench"))
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--keep-work", action="store_true")
    args = parser.parse_args(argv)

    args.grids = [int(g) for g in args.grids.split(",")]
    ticks_list = [int(t) for t in args.ticks.split(",") if int(t) > 0]
    engines = [e.strip() for e in args.engines.split(",") if e.strip()]
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    work = Path(tempfile.mkdtemp(prefix="avbd-bench-", dir=out))
    try:
        for grid in args.grids:
            args.grid = grid
            for engine in engines:
                for ticks in ticks_list:
                    print("bench %s %d Hz grid=%d ..." % (engine, ticks, grid), flush=True)
                    samples = measure(args, engine, ticks, work)
                    for sample in samples:
                        sample["bodies"] = grid * grid
                    rows.append(summarize(args, args.label, engine, ticks, samples))
                    row = rows[-1]
                    print("  wall=%.4f ms engine=%.4f ms solver=%s realtime=%.3f rounds=%s" % (
                        row["wall_ms_per_step"], row["engine_ms_per_step"],
                        ("%.4f" % row["solver_ms_per_step"]) if "solver_ms_per_step" in row else "-",
                        row["realtime_factor"], row.get("solver_rounds", "-")), flush=True)
    finally:
        if args.keep_work:
            print("work directory kept: %s" % work)
        else:
            shutil.rmtree(work, ignore_errors=True)

    (out / ("results-%s.json" % args.label)).write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    # Accumulate across labels in one directory so an A/B is one file, not a scavenger hunt.
    combined_path = out / "results.json"
    combined = []
    if combined_path.is_file():
        try:
            combined = [r for r in json.loads(combined_path.read_text(encoding="utf-8"))
                        if r.get("label") != args.label]
        except (OSError, json.JSONDecodeError):
            combined = []
    combined.extend(rows)
    combined_path.write_text(json.dumps(combined, indent=2) + "\n", encoding="utf-8")
    report = markdown(combined)
    (out / "results.md").write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BenchError as error:
        print("error: %s" % error, file=sys.stderr)
        sys.exit(1)
