#!/usr/bin/env python3
"""Unit tests for the CPU nDGP cap-sweep diagnostic parser."""

from __future__ import annotations

import unittest

import compare_cuda_ndgp_outputs as comparator
import validate_cuda_ndgp_diagnostic as diagnostic


class DiagnosticParserTest(unittest.TestCase):
    def test_capped_and_converged_residuals(self) -> None:
        text = (
            " WARNING: nDGP level  5 NOT converged after 100 iters, "
            "res= 1.721E-01\n"
            " nDGP level  6 converged in  87 iters, res= 9.500E-05\n"
        )
        capped = diagnostic.residuals(text, 5, 100)
        converged = diagnostic.residuals(text, 6, 100)
        self.assertEqual(capped[0]["outcome"], "capped")
        self.assertAlmostEqual(float(capped[0]["residual"]), 0.1721)
        self.assertEqual(converged[0]["outcome"], "converged")

    def test_iteration_beyond_cap_is_rejected(self) -> None:
        with self.assertRaises(diagnostic.DiagnosticError):
            diagnostic.residuals(
                "nDGP level 5 converged in 101 iters, res=1.0E-5", 5, 100
            )

    def test_four_digit_converged_and_capped_rows(self) -> None:
        text = (
            "nDGP level 5 converged in 1047 iters, res=9.9E-5\n"
            "WARNING: nDGP level 5 NOT converged after 1600 iters, res=2.7E-4\n"
        )
        rows = diagnostic.residuals(text, 5, 1600)
        self.assertEqual(rows[0]["iterations"], 1047)
        self.assertEqual(rows[0]["outcome"], "converged")
        self.assertEqual(rows[1]["iterations"], 1600)
        self.assertEqual(rows[1]["outcome"], "capped")

    def test_unknown_warning_spelling_is_rejected(self) -> None:
        with self.assertRaises(diagnostic.DiagnosticError):
            diagnostic.validate_warnings("WARN: unexpected fallback\n", 100)
        with self.assertRaises(diagnostic.DiagnosticError):
            diagnostic.validate_warnings("warning: unexpected fallback\n", 100)

    def test_expected_diagnostic_warnings_are_allowed(self) -> None:
        diagnostic.validate_warnings(
            "WARNING: IC header carries no omega_b (legacy grafic); "
            "using namelist omega_b\n"
            "WARNING: nDGP level  5 NOT converged after 100 iters, "
            "res= 1.721E-01\n",
            100,
        )

    def test_cpu_diagnostic_rejects_cuda_dispatch(self) -> None:
        for marker in (
            "[CUDA_MG] B=64 C=8 gs=1",
            "[CUDA_NGR] B=64 C=8 uploads=1",
            "[CUDA_PM_GATHER] B=64 C=8 gather=1",
            "CUDA pool: MPI local rank 0 -> GPU 0",
            "Adaptive loop: CUDA pool initialized",
        ):
            with self.subTest(marker=marker):
                with self.assertRaises(diagnostic.DiagnosticError):
                    diagnostic.validate_cpu_only_markers(marker)

    def test_smooth_residual_acceptance(self) -> None:
        rows = {
            "5": [
                {"outcome": "converged", "iterations": 300, "residual": 9.0e-5},
                {"outcome": "converged", "iterations": 250, "residual": 8.0e-5},
            ],
            "6": [
                {"outcome": "converged", "iterations": 700, "residual": 9.9e-5},
                {"outcome": "converged", "iterations": 600, "residual": 7.0e-5},
            ],
        }
        self.assertEqual(diagnostic.validate_smooth_residuals(rows, 1.0e-4), 700)

    def test_smooth_rejects_capped_or_slow_or_missing_solves(self) -> None:
        good = {"outcome": "converged", "iterations": 300, "residual": 9.0e-5}
        cases = (
            {"5": [good, good], "6": [good, {**good, "outcome": "capped"}]},
            {"5": [good, good], "6": [good, {**good, "iterations": 801}]},
            {"5": [good], "6": [good, good]},
            {"5": [good, good], "6": [good, {**good, "residual": 1.0e-4}]},
        )
        for rows in cases:
            with self.subTest(rows=rows):
                with self.assertRaises(diagnostic.DiagnosticError):
                    diagnostic.validate_smooth_residuals(rows, 1.0e-4)

    def test_l5_control_separates_diagnostic_from_scientific_status(self) -> None:
        passed = [
            {"outcome": "converged", "iterations": 700, "residual": 9.0e-5},
            {"outcome": "converged", "iterations": 400, "residual": 8.0e-5},
        ]
        failed = [
            {"outcome": "capped", "iterations": 1600, "residual": 2.7e-4},
            {"outcome": "converged", "iterations": 100, "residual": 8.0e-5},
        ]
        self.assertEqual(
            diagnostic.assess_l5_control_residuals(passed, 1.0e-4)["status"],
            "BASE_SOLVER_PASS",
        )
        assessment = diagnostic.assess_l5_control_residuals(failed, 1.0e-4)
        self.assertEqual(assessment["status"], "BASE_SOLVER_FAIL")
        self.assertEqual(assessment["max_iteration"], 1600)

    def test_l5_control_requires_exactly_two_rows(self) -> None:
        with self.assertRaises(diagnostic.DiagnosticError):
            diagnostic.assess_l5_control_residuals(
                [{"outcome": "converged", "iterations": 1, "residual": 0.0}],
                1.0e-4,
            )

    def test_l5_control_requires_both_rank_owners_and_no_l6(self) -> None:
        valid = [
            {"source_cpu": 1, "owner_level_counts": {"5": 2048}},
            {"source_cpu": 2, "owner_level_counts": {"5": 2048, "6": 0}},
            {"kind": "info"},
        ]
        diagnostic.validate_l5_owner_inputs(valid)
        invalid_cases = (
            [
                {"source_cpu": 1, "owner_level_counts": {"5": 4096}},
                {"source_cpu": 2, "owner_level_counts": {"5": 0}},
            ],
            [
                {"source_cpu": 1, "owner_level_counts": {"5": 2048, "6": 1}},
                {"source_cpu": 2, "owner_level_counts": {"5": 2048}},
            ],
            [
                {"source_cpu": 1, "owner_level_counts": {"5": 2048}},
                {"source_cpu": 1, "owner_level_counts": {"5": 2048}},
            ],
        )
        for owners in invalid_cases:
            with self.subTest(owners=owners):
                with self.assertRaises(diagnostic.DiagnosticError):
                    diagnostic.validate_l5_owner_inputs(owners)

    def test_pinned_particle_set_is_noncontiguous_and_exact_size(self) -> None:
        identities = comparator.pinned_particle_ids()
        self.assertEqual(len(identities), 61440)
        self.assertNotIn(1 + 8 + 32 * (8 + 32 * 8), identities)
        self.assertIn(32769, identities)
        self.assertIn(65536, identities)


if __name__ == "__main__":
    unittest.main()
