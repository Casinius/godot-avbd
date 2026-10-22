#!/usr/bin/env python3
"""Compare TARS's four real suspension scripts against Jolt, including time rates.

Temporary minimal projects select the backend before startup; TARS is never edited.
Compression change ratio means delta(compression/rest_length)/delta(time), in 1/s.
Strict defaults are diagnostic tolerances, not claims of engine equivalence.
"""
import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

WHEELS = ("fl", "fr", "bl", "br")
VECTOR_FIELDS = {"position": 3, "rotation": 4, "linear_velocity": 3, "angular_velocity": 3}
WHEEL_FIELDS = ("compression", "compression_ratio", "compression_rate", "compression_ratio_rate", "force")


class TraceError(ValueError):
    pass


def finite_number(value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise TraceError(f"{where} must be a finite number")
    return float(value)


def validate_trace(trace, expected_engine=None):
    if not isinstance(trace, dict):
        raise TraceError("trace must be an object")
    if "metadata" in trace:
        if not isinstance(trace["metadata"], dict):
            raise TraceError("metadata must be an object")
        trace = {**trace["metadata"], "samples": trace.get("samples")}
    if trace.get("schema_version") != 1:
        raise TraceError("schema_version must be 1")
    if trace.get("engine") not in ("Jolt", "AVBD"):
        raise TraceError("engine must be Jolt or AVBD")
    if expected_engine and trace["engine"] != expected_engine:
        raise TraceError(f"expected {expected_engine}, got {trace['engine']}")
    hz = trace.get("hz")
    if isinstance(hz, bool) or not isinstance(hz, int) or hz <= 0:
        raise TraceError("hz must be a positive integer")
    if trace.get("scenario") not in ("tars", "symmetric", "com_stress"):
        raise TraceError("unsupported scenario")
    config = trace.get("config")
    if not isinstance(config, dict) or not isinstance(config.get("wheel_lengths"), dict):
        raise TraceError("config.wheel_lengths is required")
    for w in WHEELS:
        if finite_number(config["wheel_lengths"].get(w), f"wheel_lengths.{w}") <= 0:
            raise TraceError("wheel lengths must be positive")
    samples = trace.get("samples")
    if not isinstance(samples, list) or len(samples) < 2:
        raise TraceError("need at least two samples")
    for i, sample in enumerate(samples):
        if not isinstance(sample, dict) or type(sample.get("tick")) is not int or sample["tick"] != i:
            raise TraceError(f"samples[{i}].tick must be contiguous starting at zero")
        time = finite_number(sample.get("time"), "sample time")
        if not math.isclose(time, i / hz, abs_tol=1e-7, rel_tol=1e-7):
            raise TraceError(f"incorrect fixed time grid at sample {i}")
        for field, size in VECTOR_FIELDS.items():
            value = sample.get(field)
            if not isinstance(value, list) or len(value) != size:
                raise TraceError(f"{field} needs {size} components")
            for v in value:
                finite_number(v, field)
        qnorm = math.sqrt(sum(v*v for v in sample["rotation"]))
        if not math.isclose(qnorm, 1, abs_tol=1e-4):
            raise TraceError("rotation quaternion must be unit length")
        if not isinstance(sample.get("wheels"), dict) or set(sample["wheels"]) != set(WHEELS):
            raise TraceError("sample must contain exactly fl/fr/bl/br")
        for w in WHEELS:
            wheel = sample["wheels"][w]
            if not isinstance(wheel, dict):
                raise TraceError("wheel must be an object")
            for f in WHEEL_FIELDS:
                finite_number(wheel.get(f), f"{w}.{f}")
            for f in ("contact", "contact_transition"):
                if type(wheel.get(f)) is not bool:
                    raise TraceError(f"{w}.{f} must be boolean")
            length = config["wheel_lengths"][w]
            if not math.isclose(wheel["compression_ratio"], wheel["compression"] / length, abs_tol=1e-6, rel_tol=1e-6):
                raise TraceError(f"{w} compression ratio has incorrect units")
            old = samples[i-1]["wheels"][w] if i else wheel
            rate = (wheel["compression"] - old["compression"]) * hz
            for f, expected in (("compression_rate", rate), ("compression_ratio_rate", rate / length)):
                # Godot JSON precision can amplify rounding during differentiation.
                if not math.isclose(wheel[f], expected, abs_tol=1e-4, rel_tol=1e-4):
                    raise TraceError(f"{w}.{f} inconsistent at tick {i}")
            if wheel["contact_transition"] != (i > 0 and old["contact"] != wheel["contact"]):
                raise TraceError(f"{w} incorrect contact transition at tick {i}")
    return trace


def compatible(a, b):
    for key in ("hz", "scenario", "config"):
        if a[key] != b[key]:
            raise TraceError(f"incompatible metadata: {key}")
    if a.get("godot_version") != b.get("godot_version"):
        raise TraceError("Godot versions differ")
    for trace in (a, b):
        if trace.get("direct_state_class"):
            actual_avbd = trace["direct_state_class"] == "AVBDDirectBodyState3D"
            if actual_avbd != (trace["engine"] == "AVBD") or (trace["engine"] == "Jolt" and "Jolt" not in trace["direct_state_class"]):
                raise TraceError("backend label does not match the active direct-state class")
    if len(a["samples"]) != len(b["samples"]):
        raise TraceError("sample counts differ")
    for x, y in zip(a["samples"], b["samples"]):
        if x["tick"] != y["tick"] or not math.isclose(x["time"], y["time"], abs_tol=1e-7):
            raise TraceError("sample grids differ")


def rmse(values):
    return math.sqrt(sum(x*x for x in values) / len(values)) if values else None


def summarize(errors, refs, floor=1e-6):
    if not errors:
        return {"count": 0, "rmse": None, "max_abs": None, "normalized_rmse": None}
    rms = rmse(errors)
    return {"count": len(errors), "rmse": rms, "max_abs": max(map(abs, errors)),
            "normalized_rmse": rms / max(rmse(refs), floor)}


def analyze(jolt, avbd, near_zero_floor=1e-3):
    jolt = validate_trace(jolt, "Jolt")
    avbd = validate_trace(avbd, "AVBD")
    compatible(jolt, avbd)
    if near_zero_floor <= 0 or not math.isfinite(near_zero_floor):
        raise TraceError("near-zero floor must be positive")
    js, ats = jolt["samples"], avbd["samples"]
    report = {"schema_version": 1, "scenario": jolt["scenario"], "hz": jolt["hz"],
              "sample_count": len(js), "units": {"compression": "m", "compression_rate": "m/s",
              "compression_ratio": "rest-length fraction", "compression_ratio_rate": "1/s",
              "force": "N", "rotation_error": "rad"}, "state": {}, "derivative": {}, "contact": {}}
    for field in ("position", "linear_velocity", "angular_velocity"):
        errors = [a-b for x,y in zip(js,ats) for a,b in zip(y[field],x[field])]
        refs = [v for x in js for v in x[field]]
        report["state"][field] = summarize(errors, refs)
    report["chassis_rates"] = {}
    for field, units in (("position", "m/s"), ("linear_velocity", "m/s^2"), ("angular_velocity", "rad/s^2")):
        errors, refs = [], []
        for i in range(1, len(js)):
            dt = js[i]["time"] - js[i-1]["time"]
            for component in range(3):
                jr = (js[i][field][component] - js[i-1][field][component]) / dt
                ar = (ats[i][field][component] - ats[i-1][field][component]) / dt
                errors.append(ar - jr)
                refs.append(jr)
        report["chassis_rates"][field] = {"units": units, **summarize(errors, refs)}
    angles = []
    position_errors = []
    for x, y in zip(js, ats):
        dot = abs(sum(a*b for a,b in zip(x["rotation"], y["rotation"])))
        angles.append(2 * math.acos(min(1.0, max(0.0, dot))))
        position_errors.append(math.dist(x["position"], y["position"]))
    report["state"]["rotation_angle"] = summarize(angles, [1]*len(angles))
    first = None
    state_gate, rate_gate, disagreements, coverage = [], [], [], []
    growth_errors = []
    for w in WHEELS:
        length = jolt["config"]["wheel_lengths"][w]
        state = {}
        for field in ("compression", "compression_ratio", "force"):
            refs = [x["wheels"][w][field] for x in js]
            errors = [y["wheels"][w][field]-x["wheels"][w][field] for x,y in zip(js,ats)]
            state[field] = summarize(errors, refs)
        report["state"][w] = state
        state_gate.append(state["compression_ratio"]["rmse"])
        agree = sum(x["wheels"][w]["contact"] == y["wheels"][w]["contact"] for x,y in zip(js,ats))
        ref_contact = sum(x["wheels"][w]["contact"] for x in js)
        report["contact"][w] = {"agreement": agree/len(js), "jolt_contact_samples": ref_contact}
        disagreements.append(1-agree/len(js))
        errors, refs, sustained_errors, sustained_refs, quotients = [], [], [], [], []
        transitions, noncontact, near_zero = 0, 0, 0
        for i in range(1, len(js)):
            x0,x1,y0,y1 = (js[i-1]["wheels"][w],js[i]["wheels"][w],ats[i-1]["wheels"][w],ats[i]["wheels"][w])
            dt = js[i]["time"]-js[i-1]["time"]
            jr = (x1["compression"]-x0["compression"])/dt/length
            ar = (y1["compression"]-y0["compression"])/dt/length
            error = ar-jr
            errors.append(error)
            refs.append(jr)
            growth_errors.append(error)
            transition = x1["contact_transition"] or y1["contact_transition"]
            if transition:
                transitions += 1
            elif not all(v["contact"] for v in (x0,x1,y0,y1)):
                noncontact += 1
            else:
                sustained_errors.append(error)
                sustained_refs.append(jr)
                if abs(jr) < near_zero_floor:
                    near_zero += 1
                else:
                    quotients.append(ar/jr)
        ratio_rate = summarize(errors, refs, near_zero_floor)
        steady_rate = summarize(sustained_errors, sustained_refs, near_zero_floor)
        report["derivative"][w] = {"compression_ratio_rate": ratio_rate,
            "compression_rate": summarize([v*length for v in errors], [v*length for v in refs], near_zero_floor*length),
            "sustained_contact_ratio_rate": steady_rate,
            "avbd_over_jolt_rate_ratio": {"count": len(quotients), "mean": sum(quotients)/len(quotients) if quotients else None,
                "max_deviation_from_one": max(abs(v-1) for v in quotients) if quotients else None},
            "excluded_from_quotient": {"transitions": transitions, "noncontact": noncontact, "near_zero": near_zero}}
        coverage.append(len(sustained_errors))
        rate_gate.append(ratio_rate["normalized_rmse"])
    for i, (x,y) in enumerate(zip(js,ats)):
        if first is None and (position_errors[i] > 1e-6 or any(abs(x["wheels"][w]["compression_ratio"]-y["wheels"][w]["compression_ratio"]) > 1e-6 for w in WHEELS)):
            first = {"tick": x["tick"], "time": x["time"]}
    peak = max(range(len(js)), key=position_errors.__getitem__)
    # Track error growth independently of matching final endpoints.
    growth = [(b-a)*jolt["hz"] for a,b in zip(position_errors, position_errors[1:])]
    report.update({"first_divergence": first, "contact_agreement": 1-max(disagreements),
        "peak_chassis_position_error": {"max_abs": position_errors[peak], "time": js[peak]["time"]},
        "error_growth": {"position_error_rate_m_s": summarize(growth, [1]*len(growth)),
                         "compression_ratio_error_rate_1_s": summarize(growth_errors, [1]*len(growth_errors))},
        "gates": {"max_state_rmse": max(state_gate), "max_derivative_rmse": max(rate_gate),
                  "max_contact_disagreement": max(disagreements), "min_sustained_contact_intervals": min(coverage)}})
    return report


def read_trace(path, engine):
    try:
        return validate_trace(json.loads(Path(path).read_text(encoding="utf-8")), engine)
    except (OSError, json.JSONDecodeError) as error:
        raise TraceError(f"cannot read {path}: {error}") from error


def resource_closure(project):
    """Copy only real fixture dependencies, not unrelated broken scenes/autoloads."""
    classes = {}
    for path in (project/"scenes/vehicle").rglob("*.gd"):
        match = re.search(r"(?m)^class_name\s+(\w+)", path.read_text(encoding="utf-8"))
        if match:
            classes[match.group(1)] = path.relative_to(project).as_posix()
    pending = ["tests/suspension_trace.gd", "tests/zero_tire_model.gd", "def_frontw.tres", "def_rearw.tres"]
    found = set()
    while pending:
        rel = pending.pop()
        if rel in found:
            continue
        path = project/rel
        if not path.is_file():
            raise TraceError(f"fixture dependency is missing: {path}")
        found.add(rel)
        if path.suffix not in (".gd", ".tres", ".tscn"):
            continue
        text = path.read_text(encoding="utf-8")
        pending.extend(re.findall(r'res://([^"\s]+)', text))
        pending.extend(rel for name,rel in classes.items() if re.search(r"\b"+re.escape(name)+r"\b", text))
    return found


def run_command(cmd, log, timeout):
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, timeout=timeout, check=False)
        log.write_text(result.stdout, encoding="utf-8")
    except subprocess.TimeoutExpired as error:
        output = error.stdout or b""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        log.write_text(output, encoding="utf-8")
        raise TraceError(f"Godot timed out after {timeout}s; log: {log}\n{output[-2000:]}") from error
    except OSError as error:
        raise TraceError(f"cannot launch Godot: {error}") from error
    if result.returncode or "SCRIPT ERROR:" in result.stdout or "ERROR:" in result.stdout:
        raise TraceError(f"Godot exited {result.returncode}; log: {log}\n{result.stdout[-4000:]}")


def run_engine(args, engine, output):
    project = Path(args.project).resolve()
    backend = "Jolt Physics" if engine == "Jolt" else "AVBD"
    temp_root = output.parent/".work"
    temp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="tars-suspension-", dir=temp_root) as temp:
        clone = Path(temp)
        for rel in resource_closure(project):
            dest = clone/rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            # Real copies isolate .uid/import writes from the user's subproject.
            dest.write_bytes((project/rel).read_bytes())
        binary = Path(args.avbd_library).resolve() if args.avbd_library else project/"gdavbd/bin/libavbd.linux.x86_64.so"
        if not binary.is_file():
            raise TraceError(f"AVBD library not found: {binary}; build first or use --avbd-library")
        ext = clone/"gdavbd"
        (ext/"bin").mkdir(parents=True)
        (ext/"avbd.gdextension").write_bytes((project/"gdavbd/avbd.gdextension").read_bytes())
        shutil.copyfile(binary, ext/"bin/libavbd.linux.x86_64.so")
        (clone/"project.godot").write_text(f'''config_version=5
[application]
config/name="TARS suspension comparison"
[physics]
3d/physics_engine="{backend}"
3d/run_on_separate_thread=false
common/physics_ticks_per_second={args.hz}
common/max_physics_steps_per_frame=1
3d/default_gravity=9.8
3d/default_linear_damp=0.0
3d/default_angular_damp=0.0
avbd/substeps=1
jolt_physics_3d/simulation/velocity_steps=20
jolt_physics_3d/simulation/position_steps=5
jolt_physics_3d/simulation/allow_sleep=false
[rendering]
renderer/rendering_method="gl_compatibility"
''', encoding="utf-8")
        prefix = [args.godot, "--headless", "--path", str(clone)]
        if output.exists():
            output.unlink()
        run_command(prefix+["--editor", "--import", "--quit"], output.with_suffix(".import.log"), args.timeout)
        run_command(prefix+["--fixed-fps", str(args.hz), "--script", "res://tests/suspension_trace.gd", "--",
                    f"--output={output}", f"--engine={backend}", f"--hz={args.hz}",
                    f"--ticks={args.ticks}", f"--scenario={args.scenario}"], output.with_suffix(".run.log"), args.timeout)
        trace = read_trace(output, engine)
        with binary.open("rb") as stream:
            digest = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
            trace["avbd_library_sha256"] = digest.hexdigest()
        return trace


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--godot", default="godot")
    parser.add_argument("--project", default=str(Path(__file__).resolve().parents[1]/"TARS"))
    parser.add_argument("--avbd-library", help="use a freshly built extension without overwriting TARS's library")
    parser.add_argument("--output-dir", default="build/suspension_comparison")
    parser.add_argument("--hz", type=int, default=500)
    parser.add_argument("--ticks", type=int, default=3000)
    parser.add_argument("--scenario", choices=("tars", "symmetric", "com_stress"), default="tars")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--avbd-trace")
    parser.add_argument("--jolt-trace")
    parser.add_argument("--no-strict", action="store_true")
    parser.add_argument("--max-state-rmse", type=float, default=.05, help="absolute compression/rest-length RMSE")
    parser.add_argument("--max-derivative-rmse", type=float, default=.10, help="rate RMSE / max(Jolt rate RMS, floor)")
    parser.add_argument("--max-contact-disagreement", type=float, default=.05)
    parser.add_argument("--rate-floor", type=float, default=.001, help="near-zero reference threshold, 1/s")
    parser.add_argument("--min-contact-intervals", type=int, default=10)
    args = parser.parse_args(argv)
    if args.hz < 1 or args.ticks < 2 or args.timeout <= 0 or args.min_contact_intervals < 1:
        parser.error("positive hz/timeout/contact coverage and ticks >= 2 required")
    for threshold in (args.max_state_rmse, args.max_derivative_rmse, args.max_contact_disagreement, args.rate_floor):
        if not math.isfinite(threshold) or threshold < 0:
            parser.error("thresholds must be finite and nonnegative")
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    jpath = Path(args.jolt_trace).resolve() if args.jolt_trace else out/"jolt.json"
    apath = Path(args.avbd_trace).resolve() if args.avbd_trace else out/"avbd.json"
    try:
        if args.analyze_only:
            jolt, avbd = read_trace(jpath, "Jolt"), read_trace(apath, "AVBD")
        else:
            jolt = run_engine(args, "Jolt", jpath)
            avbd = run_engine(args, "AVBD", apath)
        report = analyze(jolt, avbd, args.rate_floor)
        report["strict_thresholds"] = {"compression_ratio_rmse": args.max_state_rmse,
            "relative_rate_rmse": args.max_derivative_rmse, "contact_disagreement": args.max_contact_disagreement,
            "min_sustained_contact_intervals": args.min_contact_intervals, "rate_floor_1_s": args.rate_floor}
        gates = report["gates"]
        report["passed"] = (gates["max_state_rmse"] <= args.max_state_rmse and
            gates["max_derivative_rmse"] <= args.max_derivative_rmse and
            gates["max_contact_disagreement"] <= args.max_contact_disagreement and
            gates["min_sustained_contact_intervals"] >= args.min_contact_intervals)
        report["avbd_library_sha256"] = avbd.get("avbd_library_sha256")
        (out/"report.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
        print(f"Jolt reference: {report['sample_count']} samples, {report['hz']} Hz, scenario={report['scenario']}")
        for w in WHEELS:
            state = report["state"][w]["compression_ratio"]["rmse"]
            rate = report["derivative"][w]["compression_ratio_rate"]
            print(f"{w}: compression RMSE={state:.6g} rest length; rate RMSE={rate['rmse']:.6g} 1/s; relative rate error={rate['normalized_rmse']:.2%}")
        print(f"contact agreement={report['contact_agreement']:.2%}; first divergence={report['first_divergence']}; passed={report['passed']}")
        print(f"report: {out/'report.json'}")
        return 0 if report["passed"] or args.no_strict else 1
    except (TraceError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
