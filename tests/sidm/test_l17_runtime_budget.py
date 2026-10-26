#!/usr/bin/env python3
"""Small synthetic regression tests for evidence-based runtime projections."""

import tempfile
import unittest
from pathlib import Path

from audit_l17_runtime_budget import audit_interval, forecast, measure_flat_snapshot


class ForecastTests(unittest.TestCase):
    def test_compute_and_calendar_days_are_distinct(self):
        result = forecast(0.2, 0.3, 8640, 1.0, 90, 0.5)
        self.assertAlmostEqual(result["remaining_compute_days_at_fixed_rate"], 0.7)
        self.assertAlmostEqual(result["calendar_days_at_service_fraction"], 1.4)
        self.assertAlmostEqual(result["speedup_required_for_budget_at_service_fraction"], 1.4 / 90)

    def test_already_reached_target_needs_zero_compute(self):
        result = forecast(0.2, 0.3, 8640, 0.3, 90, 1)
        self.assertEqual(result["remaining_compute_days_at_fixed_rate"], 0)

    def test_slow_rate_exceeds_budget_without_a_completion_promise(self):
        result = forecast(0.24, 0.25, 2 * 86400, 1, 90, 1)
        self.assertAlmostEqual(result["remaining_compute_days_at_fixed_rate"], 150)
        self.assertFalse(result["within_budget_under_fixed_rate_assumption"])

    def test_invalid_inputs_fail(self):
        valid = [0.2, 0.3, 8640, 1, 90, 1]
        for index, value in ((0, 0.3), (1, 0.1), (2, 0), (3, 0.25),
                             (4, 0), (5, 0), (5, 1.1), (2, float("nan"))):
            with self.subTest(index=index, value=value), self.assertRaises(ValueError):
                args = valid.copy()
                args[index] = value
                forecast(*args)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.run = Path(self.temporary.name)
        for output, epoch, step in ((2, 0.2, 10), (3, 0.3, 20)):
            directory = self.run / f"output_{output:05d}"
            directory.mkdir()
            (directory / f"info_{output:05d}.txt").write_text(
                f"aexp = {epoch}\nnstep_coarse = {step}\nncpu = 24\nlevelmax = 17\n")
        snapshot = self.run / "snapshot"
        snapshot.mkdir()
        self.namelist = snapshot / "run_123.nml"
        self.namelist.write_text("&RUN_PARAMS\nnrestart=2\n/\n")
        self.runner = snapshot / "run_z0_l17_continue.sbatch"
        self.runner.write_text(
            "#SBATCH --cpus-per-task=2\n#SBATCH --ntasks=24\n"
            "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}\n")
        self.log = self.run / "run_123.log"
        self.log.write_text(
            "Working with nproc = 24 for ndim = 3\n"
            "last commit = " + "a" * 40 + "\n"
            "Restarting at t= -3.0 nstep_coarse= 10\n"
            "FFTW3 CPU direct Poisson solver enabled\nFFT direct solve DONE\n"
            "==> Level=11 Step=3 Error=1.0E-05\n"
            "Main step= 20 mcons=0\nFine step= 100 t=-2 dt=1e-3 a=3.000E-01\n"
            "Run completed\nTotal elapsed time: 8640.0\n")

    def result(self):
        return audit_interval(self.run, 123, 2, 3, 90, 1)

    def test_matched_interval_records_provenance_not_science_acceptance(self):
        result = self.result()
        self.assertAlmostEqual(result["delta_a"], 0.1)
        self.assertEqual(result["configured_omp_threads"], 2)
        self.assertTrue(result["solver_gate"]["fftw_solver_verified"])
        self.assertFalse(result["accepted_as_science_seed"])

    def test_restart_rank_redistribution_is_allowed(self):
        path = self.run / "output_00002" / "info_00002.txt"
        path.write_text(path.read_text().replace("ncpu = 24", "ncpu = 12"))
        self.assertEqual(self.result()["configured_mpi_ranks"], 24)

    def test_failed_physics_still_has_an_operational_rate(self):
        self.log.write_text(self.log.read_text().replace("Error=1.0E-05", "Error=2.0E-04"))
        result = self.result()
        self.assertFalse(result["solver_gate"]["convergence_gate_pass"])
        self.assertFalse(result["accepted_as_science_seed"])
        self.assertGreater(result["seconds_per_delta_a"], 0)

    def test_completed_marker_is_required(self):
        self.log.write_text(self.log.read_text().replace("Run completed", "stopped abruptly"))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            self.result()

    def test_wrong_first_checkpoint_fails(self):
        self.log.write_text(self.log.read_text().replace("nstep_coarse= 10", "nstep_coarse= 11"))
        with self.assertRaisesRegex(ValueError, "coarse-step mismatch"):
            self.result()

    def test_wrong_last_checkpoint_fails(self):
        self.log.write_text(self.log.read_text().replace("Main step= 20", "Main step= 21"))
        with self.assertRaisesRegex(ValueError, "coarse-step mismatch"):
            self.result()

    def test_wrong_epoch_fails(self):
        self.log.write_text(self.log.read_text().replace("a=3.000E-01", "a=4.000E-01"))
        with self.assertRaisesRegex(ValueError, "epoch mismatch"):
            self.result()

    def test_namelist_restart_mismatch_fails(self):
        self.namelist.write_text("nrestart=1\n")
        with self.assertRaisesRegex(ValueError, "namelist/input"):
            self.result()

    def test_changed_solver_settings_are_not_silently_ignored(self):
        self.namelist.write_text(self.namelist.read_text() + "maxiter_fine=40\n")
        with self.assertRaisesRegex(ValueError, "MG settings"):
            self.result()

    def test_missing_thread_contract_fails(self):
        self.runner.write_text("#SBATCH --cpus-per-task=2\n#SBATCH --ntasks=24\n")
        with self.assertRaisesRegex(ValueError, "thread count"):
            self.result()

    def test_flat_snapshot_measurement_is_not_an_integrity_test(self):
        snapshot = self.run / "output_00003"
        (snapshot / "amr_00003.out00001").write_bytes(b"1234")
        (snapshot / "empty").touch()
        result = measure_flat_snapshot(snapshot)
        self.assertEqual(result["files"], 3)
        self.assertEqual(result["rank_payload_file_counts"]["amr"], 1)
        self.assertEqual(result["zero_byte_files"], ["empty"])
        self.assertFalse(result["payload_integrity_verified"])

    def test_nested_payload_requires_a_different_audit(self):
        snapshot = self.run / "output_00003"
        (snapshot / "nested").mkdir()
        with self.assertRaisesRegex(ValueError, "subdirectory"):
            measure_flat_snapshot(snapshot)


if __name__ == "__main__":
    unittest.main()
