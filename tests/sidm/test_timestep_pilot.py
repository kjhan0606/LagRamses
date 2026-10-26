#!/usr/bin/env python3
"""Synthetic tests of conditional pilot admission and matched-rate intervals."""

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from analyze_timestep_diagnostics import analyze_pilot_lines
from validate_pilot_bundle import sha256, validate


GATE_LOG = ["FFTW3 CPU direct Poisson solver enabled", "FFT direct solve DONE",
            "==> Level=11 Step=3 Error=1e-5", "Run completed"]


def pilot_lines():
    lines = GATE_LOG[:-1] + ["Restarting at t=-3 nstep_coarse=10"]
    for step in range(10, 16):
        epoch = .25 + (step - 10) * .001
        lines += [
            f"SIDM_DT: coarse={step} fine={step} level=8 a={epoch} dt_other=.001 "
            "dt_old=.001 pmax=.2 courant=.1 dt_sidm=.0005 dt_local=.0005 "
            "bound_active=T sidm_limited=T",
            f"SIDM_DT_PARENT: coarse={step} fine={step} level=8 a={epoch} "
            "dt_local=.0005 dt_parent=0 dt_selected=.0005 bound_active=F parent_limited=F",
        ]
        if step < 15:
            lines.append("Time elapsed since last coarse step: 10.00 s")
    return lines + ["Main step=15", "Run completed", "Total elapsed time:60.0"]


class PilotRateTests(unittest.TestCase):
    def analyze(self, lines):
        return analyze_pilot_lines(lines, 10, .25, 8, expected_steps=5, warmup_steps=2)

    def test_warmup_is_excluded_and_epochs_are_precise(self):
        result = self.analyze(pilot_lines())
        self.assertEqual(result["measured_start_coarse_step"], 12)
        self.assertEqual(result["measured_end_coarse_step"], 15)
        self.assertEqual(result["measured_coarse_wall_seconds"], 30)
        self.assertAlmostEqual(result["a_start"], .252)
        self.assertAlmostEqual(result["seconds_per_delta_a"], 10000)
        self.assertFalse(result["accepted_as_science_trajectory"])

    def test_final_rounded_epoch_is_not_used_as_precise_endpoint(self):
        lines = [line for line in pilot_lines() if not ("coarse=15" in line)]
        lines.insert(-1, "Fine step=500 a=2.555E-1")
        result = self.analyze(lines)
        self.assertEqual(result["measured_end_coarse_step"], 14)
        self.assertEqual(result["measured_coarse_wall_seconds"], 20)
        self.assertAlmostEqual(result["a_end"], .254)

    def test_missing_timer_or_completion_is_rejected(self):
        for text in ("Time elapsed since last coarse step: 10.00 s", "Run completed"):
            lines = pilot_lines()
            lines.remove(text)
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.analyze(lines)

    def test_wrong_restart_and_stop_are_rejected(self):
        for before, after in (("nstep_coarse=10", "nstep_coarse=9"),
                              ("Main step=15", "Main step=14")):
            with self.subTest(before=before), self.assertRaises(ValueError):
                self.analyze([line.replace(before, after) for line in pilot_lines()])

    def test_failed_gravity_or_banner_only_is_rejected(self):
        for old, new in (("Error=1e-5", "Error=2e-4"),
                         ("FFT direct solve DONE", "unused FFTW banner")):
            with self.subTest(old=old), self.assertRaises(ValueError):
                self.analyze([line.replace(old, new) for line in pilot_lines()])

    def test_wrong_seed_epoch_or_duplicate_base_record_is_rejected(self):
        lines = pilot_lines()
        with self.assertRaises(ValueError):
            analyze_pilot_lines(lines, 10, .3, 8, 5, 2)
        with self.assertRaises(ValueError):
            self.analyze(lines[:4] + lines[4:6] + lines[4:])

    def test_discontinuous_and_reversing_epochs_are_rejected(self):
        for mutated in ([line for line in pilot_lines() if "coarse=11" not in line],
                        [line.replace("a=0.253", "a=0.251") for line in pilot_lines()]):
            with self.subTest(lines=mutated), self.assertRaises(ValueError):
                self.analyze(mutated)

    def test_inconsistent_total_elapsed_is_rejected(self):
        with self.assertRaises(ValueError):
            self.analyze([line.replace("Total elapsed time:60.0", "Total elapsed time:1.0")
                          for line in pilot_lines()])


class BundleGateTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        self.bundle = root / "bundle"
        self.bundle.mkdir()
        self.recovery = root / "recovery"
        self.revision = "eadfc9fdac3e10cac3bf3415c05ab67fb061ae96"
        repository = Path(__file__).resolve().parents[2]
        (self.bundle / "run.nml").write_bytes((repository / "tests/sidm/timestep_pilot_sidm3.nml").read_bytes())
        for name in ("ramses_dt_pilot", "run_l17_timestep_pilot.sbatch", "instrumentation.diff",
                     "analyze_timestep_diagnostics.py", "analyze_mg_convergence_log.py",
                     "audit_poisson_binary.py", "validate_pilot_bundle.py"):
            (self.bundle / name).write_text("fixture, not an executable build\n")
        (self.bundle / "clean_build_binary_gate.json").write_text(json.dumps({
            "fftw_binary_evidence_pass": True,
            "binary_sha256": sha256(self.bundle / "ramses_dt_pilot"),
        }))
        self.manifest = {"purpose": "short_diagnostic_evolution", "accepted_as_science_seed": False,
                         "seed": "/gpfs/kjhan/Hydro/Sidm/validation/z0_production_l17/sidm3/output_00005",
                         "seed_aexp": .250336146250867, "base_source_revision": self.revision}
        self.rehash()
        for label in ("cdm", "sidm3"):
            snapshot = self.recovery / label / "snapshot"
            snapshot.mkdir(parents=True)
            (snapshot / "source_revision.txt").write_text(self.revision)
            for cap in (40, 80):
                case = self.recovery / label / f"cap{cap}"
                case.mkdir()
                (case / "run.log").write_text("\n".join(GATE_LOG))
                (case / "exit_status.txt").write_text(f"cap={cap} runtime_exit=0 analyzer_exit=0")
                (case / "mg_summary.json").write_text(json.dumps({
                    "convergence_gate_pass": True, "fftw_solver_verified": True,
                    "run_completed": True, "maxiter": cap, "epsilon": 1e-4,
                }))
        reader = patch("validate_pilot_bundle.read_seed_info", return_value=
                       "nstep_coarse=838\nncpu=24\nlevelmax=17\naexp=.250336146250867\n")
        reader.start()
        self.addCleanup(reader.stop)

    def rehash(self):
        self.manifest["files_sha256"] = {path.name: sha256(path) for path in self.bundle.iterdir()
                                        if path.name != "manifest.json"}
        (self.bundle / "manifest.json").write_text(json.dumps(self.manifest))

    def test_all_four_outcomes_are_required(self):
        validate(self.bundle, self.recovery)
        path = self.recovery / "sidm3/cap80/mg_summary.json"
        payload = json.loads(path.read_text())
        payload["convergence_gate_pass"] = False
        path.write_text(json.dumps(payload))
        with self.assertRaises(ValueError):
            validate(self.bundle, self.recovery)

    def test_cached_green_summary_cannot_override_bad_raw_log(self):
        (self.recovery / "cdm/cap40/run.log").write_text("\n".join(GATE_LOG).replace("1e-5", "2e-4"))
        with self.assertRaises(ValueError):
            validate(self.bundle, self.recovery)

    def test_changed_binary_or_namelist_is_rejected(self):
        for name in ("ramses_dt_pilot", "run.nml"):
            path = self.bundle / name
            old = path.read_bytes()
            path.write_bytes(old + b"changed")
            with self.subTest(name=name), self.assertRaises(ValueError):
                validate(self.bundle, self.recovery)
            path.write_bytes(old)

    def test_rehashed_unsafe_output_policy_is_rejected(self):
        path = self.bundle / "run.nml"
        path.write_text(path.read_text().replace("noutput=1", "noutput=0"))
        self.rehash()
        with self.assertRaises(ValueError):
            validate(self.bundle, self.recovery)

    def test_source_revision_and_runtime_failure_are_rejected(self):
        path = self.recovery / "cdm/snapshot/source_revision.txt"
        path.write_text("wrong revision")
        with self.assertRaises(ValueError):
            validate(self.bundle, self.recovery)
        path.write_text(self.revision)
        (self.recovery / "sidm3/cap80/exit_status.txt").write_text("runtime_exit=1 analyzer_exit=0")
        with self.assertRaises(ValueError):
            validate(self.bundle, self.recovery)


if __name__ == "__main__":
    unittest.main()
