#!/usr/bin/env python3
"""Check closed-form angular moments against independent adaptive quadrature."""

import math
import unittest

from scipy.integrate import quad

from test_sidm_angular import EPS2, ISOTROPIC_MOMENTS, matched_isotropic_normalizations, \
    mean_identical_particle_transfer, mean_momentum_transfer, mean_viscosity_transfer, \
    yukawa_model_contract


class AngularIntegralTests(unittest.TestCase):
    def test_normalization_and_transport_weights(self):
        for regulator in (*EPS2, 19.999, 20.001, 200.0, 20000.0, 2e8):
            norm = regulator * (regulator + 2.0) / 2.0
            expectations = (
                ("normalization", lambda mu: 1.0, 1.0),
                ("momentum", lambda mu: 1.0 - mu, mean_momentum_transfer(regulator)),
                ("viscosity", lambda mu: 1.0 - mu**2, mean_viscosity_transfer(regulator)),
                ("identical particle", lambda mu: 1.0 - abs(mu),
                 mean_identical_particle_transfer(regulator)),
            )
            for name, weight, analytic in expectations:
                with self.subTest(eps2=regulator, weight=name):
                    measured, error = quad(
                        lambda mu: norm * weight(mu) / (1.0 - mu + regulator)**2,
                        -1.0, 1.0, points=[0.0], epsabs=1e-11, epsrel=1e-11,
                    )
                    self.assertLess(error, 1e-10)
                    self.assertLess(abs(measured - analytic), 3e-10)

    def test_isotropic_reference_weights(self):
        for name, weight in (("momentum", lambda mu: 1.0 - mu),
                             ("viscosity", lambda mu: 1.0 - mu**2),
                             ("identical_particle", lambda mu: 1.0 - abs(mu))):
            integral, _ = quad(lambda mu: weight(mu) / 2.0, -1.0, 1.0, points=[0.0])
            self.assertAlmostEqual(integral, ISOTROPIC_MOMENTS[name], places=14)

    def test_matched_cross_sections_against_quadrature(self):
        weights = {"momentum": lambda mu: 1.0 - mu,
                   "viscosity": lambda mu: 1.0 - mu**2,
                   "identical_particle": lambda mu: 1.0 - abs(mu)}
        for regulator in EPS2:
            normalizations = matched_isotropic_normalizations(regulator)
            for name, weight in weights.items():
                with self.subTest(b=regulator, weight=name):
                    norm = regulator * (regulator + 2.0) / 2.0
                    moment, _ = quad(lambda mu: norm * weight(mu) /
                                     (1.0 - mu + regulator)**2, -1.0, 1.0,
                                     points=[0.0], epsabs=1e-11, epsrel=1e-11)
                    iso_total = 3.0
                    matching = normalizations[name]["rutherford_total_over_isotropic_total"]
                    self.assertAlmostEqual(iso_total * matching * moment,
                                           iso_total * ISOTROPIC_MOMENTS[name], places=10)

    def test_default_matching_factors(self):
        ratios = matched_isotropic_normalizations(0.02)
        expected = {"momentum": 13.656457044536154,
                    "viscosity": 6.094965806548731,
                    "identical_particle": 7.619590874470207}
        for name, value in expected.items():
            self.assertAlmostEqual(ratios[name]["rutherford_total_over_isotropic_total"],
                                   value, places=10)
        # Doubling the identical-particle weight in BOTH models cancels.
        identical = ratios["identical_particle"]
        self.assertEqual((2 * identical["isotropic_weight_mean"]) /
                         (2 * identical["rutherford_weight_mean"]),
                         identical["rutherford_total_over_isotropic_total"])
        viscosity = ratios["viscosity"]
        self.assertAlmostEqual(viscosity["isotropically_normalized_weight_mean"],
                               1.5 * viscosity["rutherford_weight_mean"], places=14)
        self.assertAlmostEqual(1.0 / viscosity["isotropically_normalized_weight_mean"],
                               viscosity["rutherford_total_over_isotropic_total"], places=14)

    def test_isotropic_limit(self):
        for regulator in (2e8, 2e20, 2e100):
            for result in matched_isotropic_normalizations(regulator).values():
                self.assertLess(abs(result["rutherford_total_over_isotropic_total"] - 1.0), 1e-8)

    def test_invalid_regulators(self):
        for regulator in (0.0, -1.0, math.nan, math.inf):
            for function in (mean_momentum_transfer, mean_viscosity_transfer,
                             mean_identical_particle_transfer):
                with self.subTest(b=regulator, function=function.__name__):
                    with self.assertRaises(ValueError):
                        function(regulator)

    def test_yukawa_born_total_integral_and_code_difference(self):
        contract = yukawa_model_contract()
        v0 = contract["v0_kms"]
        for row in contract["cases"]:
            velocity = row["velocity_kms"]
            r2 = (velocity / v0)**2
            # Robertson et al. angular differential law, normalized by sigma0.
            born, _ = quad(lambda mu: 0.5 / (1.0 + r2 * (1.0 - mu) / 2.0)**2,
                           -1.0, 1.0, epsabs=1e-12, epsrel=1e-12)
            self.assertAlmostEqual(born, row["born_total_over_sigma0_without_floor"], places=12)
            self.assertAlmostEqual(row["implemented_total_over_sigma0"] / born,
                                   row["implemented_over_born_total"], places=12)
        self.assertTrue(contract["cases"][0]["angular_floor_active"])
        self.assertFalse(contract["cases"][1]["angular_floor_active"])
        self.assertEqual(contract["cases"][0]["eps2_with_code_floor"], 20000.0)
        self.assertEqual(contract["cases"][1]["eps2_with_code_floor"], 20000.0)


if __name__ == "__main__":
    unittest.main()
