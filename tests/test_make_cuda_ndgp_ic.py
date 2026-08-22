#!/usr/bin/env python3
"""Unit checks for the CUDA/nDGP deterministic IC profiles."""

from __future__ import annotations

import math
import pathlib
import tempfile
import unittest

import make_cuda_ndgp_ic as fixture


class SmoothProfileTest(unittest.TestCase):
    @staticmethod
    def manifest_entries(name: str) -> dict[str, str]:
        path = pathlib.Path(__file__).resolve().parent / name
        return {
            relative: digest
            for digest, relative in (
                line.split(maxsplit=1)
                for line in path.read_text(encoding="ascii").splitlines()
            )
        }

    def test_global_cell_centres_and_level_shifts(self) -> None:
        self.assertEqual(fixture.smooth_coordinate(0, 1.0, 0.0), 0.5)
        self.assertEqual(fixture.smooth_coordinate(0, 0.5, 8.0), 8.25)
        cases = (
            (1.0, 0.0, (1.0 / 16.0) * math.cos(math.pi / 16.0)),
            (0.5, 8.0, (1.0 / 8.0) * math.cos(math.pi / 32.0)),
        )
        for dx, offset, expected_max in cases:
            shifts = [
                fixture.smooth_displacement_base_cells(i, dx, offset) / dx
                for i in range(fixture.N)
            ]
            self.assertTrue(any(value > 0.0 for value in shifts))
            self.assertTrue(any(value < 0.0 for value in shifts))
            self.assertAlmostEqual(
                max(abs(value) for value in shifts), expected_max, places=15
            )
            self.assertAlmostEqual(math.fsum(shifts), 0.0, places=14)
            self.assertTrue(all(abs(value) < 1.0 for value in shifts))

    def test_physical_displacement_is_level_independent(self) -> None:
        q = 8.25
        base_value = fixture.SMOOTH_AMPLITUDE * (fixture.H0 / 100.0) * math.sin(
            2.0 * math.pi * fixture.SMOOTH_MODE * q / fixture.N
        )
        self.assertEqual(fixture.smooth_position_value(0, 0.5, 8.0), base_value)
        self.assertLessEqual(abs(base_value), 0.04375)

    def test_cic_density_is_finite_signed_and_mass_conserving(self) -> None:
        for dx, offset in ((1.0, 0.0), (0.5, 8.0)):
            weights = fixture.smooth_cic_axis_weights(dx, offset)
            self.assertAlmostEqual(math.fsum(weights), fixture.N, places=14)
            values = [
                fixture.smooth_density_contrast(ix, iy, iz, weights)
                for iz in range(fixture.N)
                for iy in range(fixture.N)
                for ix in range(fixture.N)
            ]
            self.assertTrue(all(math.isfinite(value) for value in values))
            self.assertLess(min(values), 0.0)
            self.assertGreater(max(values), 0.0)
            self.assertGreater(math.fsum(value * value for value in values), 0.0)
            self.assertAlmostEqual(math.fsum(values) / len(values), 0.0, places=14)

    def test_harsh_profile_oracle_is_unchanged(self) -> None:
        histogram: dict[int, int] = {}
        for iz in range(fixture.N):
            for iy in range(fixture.N):
                for ix in range(fixture.N):
                    count = int(fixture.density_contrast(ix, iy, iz) + 1)
                    histogram[count] = histogram.get(count, 0) + 1
        self.assertEqual(histogram, {0: 18944, 1: 4096, 2: 6144, 4: 3072, 8: 512})

    def test_both_pinned_manifests_regenerate_exactly(self) -> None:
        cases = (
            (fixture.HARSH_PROFILE, "cuda_ndgp_ic.sha256"),
            (fixture.SMOOTH_PROFILE, "cuda_ndgp_smooth_ic.sha256"),
            (fixture.MATCHED_PROFILE, "cuda_ndgp_matched_ic.sha256"),
        )
        with tempfile.TemporaryDirectory() as temporary:
            for profile, manifest_name in cases:
                with self.subTest(profile=profile):
                    actual = fixture.generate(pathlib.Path(temporary) / profile, profile)
                    expected = (
                        pathlib.Path(__file__).resolve().parent / manifest_name
                    ).read_text(encoding="ascii")
                    self.assertEqual(actual, expected)

    def test_ids_refmap_and_velocities_are_byte_identical_between_profiles(self) -> None:
        harsh = self.manifest_entries("cuda_ndgp_ic.sha256")
        smooth = self.manifest_entries("cuda_ndgp_smooth_ic.sha256")
        matched = self.manifest_entries("cuda_ndgp_matched_ic.sha256")
        common = [
            path
            for path in harsh
            if any(
                path.endswith(suffix)
                for suffix in (
                    "ic_particle_ids",
                    "ic_refmap",
                    "ic_velcx",
                    "ic_velcy",
                    "ic_velcz",
                )
            )
        ]
        self.assertEqual(len(common), 10)
        for path in common:
            self.assertEqual(harsh[path], smooth[path], path)
            self.assertEqual(harsh[path], matched[path], path)

    def test_matched_children_reconstruct_the_coarse_source(self) -> None:
        for parent_local in range(16):
            parent = 8 + parent_local
            base_position = fixture.decoded_position_base_cells(
                parent,
                1.0,
                0.0,
                fixture.smooth_position_value(parent, 1.0, 0.0),
            )
            children = [
                fixture.decoded_position_base_cells(
                    2 * parent_local + bit,
                    0.5,
                    8.0,
                    fixture.matched_position_value(2 * parent_local + bit),
                )
                for bit in (0, 1)
            ]
            self.assertLess(max(abs(value - base_position) for value in children), 1.0e-7)
        metrics = fixture.matched_coarse_diagnostics()
        self.assertEqual(metrics["uniform_integral"], 1.0)
        self.assertEqual(metrics["matched_integral"], 1.0)
        self.assertLess(metrics["linf"], fixture.MATCHED_COARSE_CEILING)


if __name__ == "__main__":
    unittest.main()
