import json
import math
import tempfile
import unittest
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parents[1] / "tools"))
import compare_suspension as cs


def trace(engine="Jolt", values=None, contact=None, length=.2):
    values = values if values is not None else [.02, .03, .04, .05]
    contact = contact if contact is not None else [True]*len(values)
    samples = []
    for i, value in enumerate(values):
        rate = (value-values[i-1])*100 if i else 0.0
        wheels = {w: {"compression": value, "compression_ratio": value/length,
            "compression_rate": rate, "compression_ratio_rate": rate/length,
            "force": 100. if contact[i] else 0., "contact": contact[i],
            "contact_transition": i > 0 and contact[i] != contact[i-1]} for w in cs.WHEELS}
        samples.append({"tick": i, "time": i*.01, "position": [i*.1, 0, 0], "rotation": [0, 0, 0, 1],
            "linear_velocity": [1, 0, 0], "angular_velocity": [0, 0, 0], "wheels": wheels})
    return {"schema_version": 1, "engine": engine, "hz": 100, "scenario": "tars",
            "config": {"wheel_lengths": {w: length for w in cs.WHEELS}}, "samples": samples}


class CompareSuspensionTests(unittest.TestCase):
    def test_identical_traces(self):
        report = cs.analyze(trace(), trace("AVBD"))
        self.assertEqual(report["gates"]["max_state_rmse"], 0)
        self.assertEqual(report["gates"]["max_derivative_rmse"], 0)
        self.assertEqual(report["contact_agreement"], 1)
        self.assertAlmostEqual(report["derivative"]["fl"]["avbd_over_jolt_rate_ratio"]["mean"], 1)

    def test_same_endpoint_different_history_has_derivative_error(self):
        a = trace(values=[.02, .03, .04, .05])
        b = trace("AVBD", values=[.02, .01, .02, .05])
        report = cs.analyze(a, b)
        self.assertEqual(a["samples"][-1]["wheels"]["fl"]["compression"], b["samples"][-1]["wheels"]["fl"]["compression"])
        self.assertGreater(report["gates"]["max_derivative_rmse"], 0)
        self.assertEqual(report["first_divergence"]["tick"], 1)
        self.assertGreater(report["error_growth"]["compression_ratio_error_rate_1_s"]["rmse"], 0)

    def test_change_ratio_units_and_relative_rate(self):
        report = cs.analyze(trace(values=[.02, .03, .04, .05]), trace("AVBD", values=[.02, .04, .06, .08]))
        wheel = report["derivative"]["fl"]
        self.assertAlmostEqual(wheel["compression_rate"]["rmse"], 1)
        self.assertAlmostEqual(wheel["compression_ratio_rate"]["rmse"], 5)
        self.assertAlmostEqual(wheel["avbd_over_jolt_rate_ratio"]["mean"], 2)

    def test_zero_reference_uses_floor_without_nan(self):
        report = cs.analyze(trace(values=[.02]*4), trace("AVBD", values=[.02, .03, .04, .05]))
        rate = report["derivative"]["fl"]
        self.assertTrue(math.isfinite(rate["compression_ratio_rate"]["normalized_rmse"]))
        self.assertIsNone(rate["avbd_over_jolt_rate_ratio"]["mean"])
        self.assertEqual(rate["excluded_from_quotient"]["near_zero"], 3)
        json.dumps(report, allow_nan=False)

    def test_contact_transition_exclusion_does_not_hide_state(self):
        contacts = [False, True, True, True]
        report = cs.analyze(trace(contact=contacts), trace("AVBD", contact=contacts))
        wheel = report["derivative"]["fl"]
        self.assertEqual(wheel["compression_ratio_rate"]["count"], 3)
        self.assertEqual(wheel["sustained_contact_ratio_rate"]["count"], 2)
        self.assertEqual(wheel["excluded_from_quotient"]["transitions"], 1)

    def test_all_airborne_has_no_valid_contact_coverage(self):
        a = trace(contact=[False]*4)
        b = trace("AVBD", contact=[False]*4)
        report = cs.analyze(a, b)
        self.assertEqual(report["gates"]["min_sustained_contact_intervals"], 0)
        self.assertIsNone(report["derivative"]["fl"]["sustained_contact_ratio_rate"]["rmse"])

    def test_quaternion_sign_is_same_orientation(self):
        a, b = trace(), trace("AVBD")
        for sample in b["samples"]:
            sample["rotation"] = [0, 0, 0, -1]
        self.assertEqual(cs.analyze(a, b)["state"]["rotation_angle"]["max_abs"], 0)

    def test_invalid_data_and_recorded_rates_rejected(self):
        mutations = [
            lambda t: t["samples"][1]["position"].__setitem__(0, float("nan")),
            lambda t: t["samples"][1].__setitem__("time", .001),
            lambda t: t.__setitem__("engine", "Other"),
            lambda t: t["samples"][1]["wheels"]["fl"].__setitem__("compression_rate", 123),
            lambda t: t["config"]["wheel_lengths"].__setitem__("fl", 0),
            lambda t: t["samples"][1]["wheels"]["fl"].__setitem__("compression_ratio", 123),
            lambda t: t["samples"][1]["wheels"]["fl"].__setitem__("contact", 1),
            lambda t: t["samples"][1]["wheels"]["fl"].__setitem__("contact_transition", True),
        ]
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                value = trace()
                mutate(value)
                with self.assertRaises(cs.TraceError):
                    cs.validate_trace(value)

    def test_metadata_mismatch_and_engine_label_rejected(self):
        b = trace("AVBD")
        b["config"]["mass"] = 123
        with self.assertRaises(cs.TraceError):
            cs.analyze(trace(), b)
        with self.assertRaises(cs.TraceError):
            cs.analyze(trace("AVBD"), trace())
        mislabeled = trace("AVBD")
        mislabeled["direct_state_class"] = "JoltPhysicsDirectBodyState3D"
        with self.assertRaises(cs.TraceError):
            cs.analyze(trace(), mislabeled)

    def test_nested_metadata_and_json_roundtrip(self):
        a = trace()
        samples = a.pop("samples")
        nested = {"metadata": a, "samples": samples}
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"trace.json"
            path.write_text(json.dumps(nested), encoding="utf-8")
            self.assertEqual(cs.read_trace(path, "Jolt")["engine"], "Jolt")
        self.assertEqual(cs.analyze(nested, trace("AVBD"))["gates"]["max_state_rmse"], 0)

    def test_analyze_only_strict_fails_when_airborne(self):
        with tempfile.TemporaryDirectory() as d:
            jolt, avbd = Path(d)/"jolt.json", Path(d)/"avbd.json"
            jolt.write_text(json.dumps(trace(contact=[False]*4)), encoding="utf-8")
            avbd.write_text(json.dumps(trace("AVBD", contact=[False]*4)), encoding="utf-8")
            args = ["--analyze-only", "--jolt-trace", str(jolt), "--avbd-trace", str(avbd), "--output-dir", d]
            with patch("builtins.print"):
                self.assertEqual(cs.main(args), 1)
                self.assertEqual(cs.main(args+["--no-strict"]), 0)

    def test_isolated_runner_backend_and_argument_contract(self):
        project = Path(__file__).parents[1]/"TARS"
        calls = []
        def fake_run(cmd, log, timeout):
            clone = Path(cmd[cmd.index("--path")+1])
            config = (clone/"project.godot").read_text(encoding="utf-8")
            self.assertIn('3d/physics_engine="Jolt Physics"', config)
            self.assertNotIn("[autoload]", config)
            calls.append(cmd)
            if "--script" in cmd:
                path = Path(next(x.split("=", 1)[1] for x in cmd if x.startswith("--output=")))
                path.write_text(json.dumps(trace()), encoding="utf-8")
        with tempfile.TemporaryDirectory(dir=Path(__file__).parents[1]/"build") as d, patch.object(cs, "run_command", side_effect=fake_run):
            binary = Path(d)/"stub.so"
            binary.write_bytes(b"unit-test stub")
            args = SimpleNamespace(project=str(project), godot="fake-godot", hz=100, ticks=3,
                scenario="tars", timeout=1, avbd_library=str(binary))
            original = (project/"project.godot").read_bytes()
            cs.run_engine(args, "Jolt", Path(d)/"jolt.json")
            self.assertEqual((project/"project.godot").read_bytes(), original)
            self.assertIn("--import", calls[0])
            self.assertIn("--engine=Jolt Physics", calls[1])
            self.assertIn("--fixed-fps", calls[1])


if __name__ == "__main__":
    unittest.main()
