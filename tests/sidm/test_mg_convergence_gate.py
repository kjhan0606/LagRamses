#!/usr/bin/env python3
"""Regression tests for rejecting unsafe scientific runs from log evidence."""

import unittest

from analyze_mg_convergence_log import analyze_lines


class ConvergenceGateTests(unittest.TestCase):
    def summarize(self, *lines):
        return analyze_lines(list(lines), 1e-4, 10)

    def test_empty_log_has_no_convergence_evidence(self):
        self.assertFalse(self.summarize()["convergence_gate_pass"])

    def test_last_iteration_can_converge(self):
        result = self.summarize("==> Level=11 Step=10 Error=8.0D-05")
        self.assertTrue(result["convergence_gate_pass"])
        self.assertEqual(result["true_nonconvergences"], 0)

    def test_legacy_warning_still_fails_strict_gate(self):
        result = self.summarize("==> Level=11 Step=10 Error=8.0E-05",
                                "WARN: Fine multigrid Poisson failed to converge...")
        self.assertEqual(result["true_nonconvergences"], 0)
        self.assertFalse(result["convergence_gate_pass"])

    def test_real_startup_failure(self):
        result = self.summarize("==> Level=11 Step=10 Error=1.706E-04")
        self.assertEqual(result["true_nonconvergences"], 1)
        self.assertFalse(result["convergence_gate_pass"])

    def test_ambiguous_high_residual_without_iteration_cap_fails(self):
        self.assertFalse(self.summarize(
            "==> Level=11 Step=4 Error=2.0E-04")["convergence_gate_pass"])

    def test_nonfinite_residual_is_not_serialized_as_nan(self):
        for value in ("NaN", "Infinity", "-Inf"):
            with self.subTest(value=value):
                result = self.summarize(f"==> Level=11 Step=4 Error={value}")
                self.assertFalse(result["convergence_gate_pass"])
                self.assertEqual(result["nonfinite_residuals"], 1)
                self.assertIsNone(result["failures"][0]["error"])

    def test_fftw_request_without_compiled_support_fails(self):
        result = self.summarize("WARNING: use_fftw ignored (not compiled with USE_FFTW)")
        self.assertFalse(result["fftw_base_solver_enabled"])
        self.assertTrue(result["fftw_requested_but_ignored"])

    def test_fftw_enabled_requires_banner(self):
        self.assertTrue(self.summarize("FFTW3 CPU direct Poisson solver enabled (use_fftw=T)")
                        ["fftw_base_solver_enabled"])

    def test_fftw_banner_alone_does_not_prove_execution(self):
        self.assertFalse(self.summarize("FFTW3 CPU direct Poisson solver enabled (use_fftw=T)")
                         ["fftw_solver_verified"])
        self.assertTrue(self.summarize("FFTW3 CPU direct Poisson solver enabled (use_fftw=T)",
                                     "==> Level=8 FFT direct solve DONE (fast path)")
                        ["fftw_solver_verified"])

    def test_invalid_thresholds_are_rejected(self):
        for epsilon, maxiter in ((0, 10), (float("nan"), 10), (1e-4, 0)):
            with self.assertRaises(ValueError):
                analyze_lines([], epsilon, maxiter)


if __name__ == "__main__":
    unittest.main()
