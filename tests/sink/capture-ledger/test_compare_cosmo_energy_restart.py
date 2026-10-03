#!/usr/bin/env python3
"""Regression tests for the cosmological restart roundoff discriminator."""

import unittest

import numpy as np

from compare_cosmo_energy_restart import compare_array


class RoundoffGateTest(unittest.TestCase):
    def test_machine_roundoff_passes(self):
        left = np.array([1.0, 0.0, -0.5])
        right = np.array([1.0 + 3.0e-15, 0.0, -0.5])
        self.assertTrue(compare_array(left, right, True)["pass"])

    def test_physical_change_fails(self):
        left = np.array([1.0, 0.0, -0.5])
        right = np.array([1.0 + 1.0e-9, 0.0, -0.5])
        self.assertFalse(compare_array(left, right, True)["pass"])

    def test_exact_fields_reject_roundoff(self):
        left = np.array([1.0])
        right = np.array([1.0 + 3.0e-15])
        self.assertFalse(compare_array(left, right, False)["pass"])

    def test_nonfinite_values_fail(self):
        left = np.array([1.0, np.nan])
        self.assertFalse(compare_array(left, left, True)["pass"])


if __name__ == "__main__":
    unittest.main()
