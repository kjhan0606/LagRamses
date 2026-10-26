#!/usr/bin/env python3
"""Distinguish RAMSES Poisson support from linked-library/banner evidence."""

import unittest

from audit_poisson_binary import classify_symbols


class BinaryGateTests(unittest.TestCase):
    def test_dynamic_fftw_build_has_linked_code_evidence(self):
        result = classify_symbols(
            "0001 T fftw_poisson_solve_uniform_\n"
            " U fftw_plan_dft_r2c_3d@FFTW_3\n U fftw_execute_dft_r2c@FFTW_3\n")
        self.assertTrue(result["fftw_binary_evidence_pass"])
        self.assertFalse(result["solver_dispatch_verified"])
        self.assertFalse(result["physics_convergence_verified"])

    def test_static_fftw_build_is_supported(self):
        self.assertTrue(classify_symbols(
            "0001 T fftw_poisson_solve_uniform_\n"
            "0002 T fftw_plan_dft_r2c_3d\n0003 T fftw_execute_dft_r2c\n"
        )["fftw_binary_evidence_pass"])

    def test_mpi_plans_and_execute_are_supported(self):
        self.assertTrue(classify_symbols(
            "0001 T fftw_poisson_solve_uniform_\n"
            " U fftw_mpi_plan_dft_r2c_3d\n U fftw_mpi_execute_dft_r2c\n"
        )["fftw_binary_evidence_pass"])

    def test_library_without_poisson_routine_fails(self):
        result = classify_symbols(
            "0001 T multigrid_fine_\n0002 T fftw_plan_dft_r2c_3d\n"
            "0003 T fftw_execute_dft_r2c\n0004 B amr_parameters_mp_use_fftw_\n")
        self.assertTrue(result["regular_code_symbols_present"])
        self.assertFalse(result["fftw_binary_evidence_pass"])
        self.assertEqual(result["fftw_poisson_routines"], [])

    def test_undefined_poisson_routine_is_not_a_definition(self):
        self.assertFalse(classify_symbols(
            " U fftw_poisson_solve_uniform_\n"
            " U fftw_plan_dft_r2c_3d\n U fftw_execute_dft_r2c\n"
        )["fftw_binary_evidence_pass"])

    def test_entry_point_alone_is_insufficient(self):
        self.assertFalse(classify_symbols("0001 T fftw_poisson_solve_uniform_\n")
                         ["fftw_binary_evidence_pass"])

    def test_stripped_or_empty_symbol_table_is_insufficient(self):
        for text in ("", "nm: binary: no symbols", " U fftw_execute\n"):
            with self.subTest(text=text):
                self.assertFalse(classify_symbols(text)["fftw_binary_evidence_pass"])

    def test_feature_variable_or_banner_is_not_a_solver(self):
        self.assertFalse(classify_symbols(
            "0001 B amr_parameters_mp_use_fftw_\n"
            "FFTW3 CPU direct Poisson solver enabled (use_fftw=T)\n"
        )["fftw_binary_evidence_pass"])


if __name__ == "__main__":
    unittest.main()
