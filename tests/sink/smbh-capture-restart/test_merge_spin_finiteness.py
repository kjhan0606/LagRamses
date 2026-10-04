#!/usr/bin/env python3
"""Guard the finite remnant-spin lifecycle for degenerate encounters."""

import math
import unittest
from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[3]
    / "patch"
    / "lagRamses"
    / "sink_particle.kjhan.f90"
).read_text(encoding="utf-8")


def remnant_vector(q, a1, a2, spin1, spin2, orbital, lmodana, maxspin=0.998):
    """Small numerical model of the guarded Fortran vector construction."""
    lmod = math.sqrt(sum(component * component for component in orbital))
    scale = 1.0 / (1.0 + q) ** 2
    raw = [scale * (a1 * spin1[axis] + a2 * spin2[axis] * q**2)
           for axis in range(3)]
    if lmod > 0.0:
        raw = [raw[axis] + scale * orbital[axis] / lmod * lmodana * q
               for axis in range(3)]
    magnitude = math.sqrt(sum(component * component for component in raw))
    if magnitude > 0.0:
        return min(magnitude, maxspin), tuple(component / magnitude for component in raw)
    return 0.0, (0.0, 0.0, 0.0)


class MergeSpinFinitenessTest(unittest.TestCase):
    def setUp(self) -> None:
        start = SOURCE.index("Lmod =SQRT(Lx**2+Ly**2+Lz**2)")
        end = SOURCE.index("!YDspin", start)
        self.block = SOURCE[start:end]

    def test_orbital_direction_is_used_only_for_nonzero_angular_momentum(self) -> None:
        guard = self.block.index("if(Lmod>0d0)then")
        fallback = self.block.index("else", guard)
        finish = self.block.index("endif", fallback)
        orbital = self.block[guard:fallback]
        zero_orbit = self.block[fallback:finish]
        for component in ("Lx/Lmod", "Ly/Lmod", "Lz/Lmod"):
            self.assertIn(component, orbital)
            self.assertNotIn(component, zero_orbit)
        self.assertIn("a1*ax1+a2*ax2*q**2", zero_orbit)
        self.assertIn("a1*ay1+a2*ay2*q**2", zero_orbit)
        self.assertIn("a1*az1+a2*az2*q**2", zero_orbit)

        # Protect every term of the established nonzero-orbit construction.
        self.assertIn("a1*ax1+a2*ax2*q**2+(Lx/Lmod)*Lmodana*q", orbital)
        self.assertIn("a1*ay1+a2*ay2*q**2+(Ly/Lmod)*Lmodana*q", orbital)
        self.assertIn("a1*az1+a2*az2*q**2+(Lz/Lmod)*Lmodana*q", orbital)

    def test_zero_spin_vector_is_not_normalized(self) -> None:
        norm = self.block.index("a1mod=SQRT(ax1**2+ay1**2+az1**2)")
        guard = self.block.index("if(a1mod>0d0)then", norm)
        fallback = self.block.index("else", guard)
        guarded = self.block[guard:fallback]
        zero_spin = self.block[fallback:self.block.index("endif", fallback)]
        for component in ("ax1/a1mod", "ay1/a1mod", "az1/a1mod"):
            self.assertIn(component, guarded)
            self.assertNotIn(component, zero_spin)
        self.assertIn("spinmag_new(igrp)=0d0", zero_spin)
        self.assertIn("bhspin_new(igrp,1:ndim)=0d0", zero_spin)

    def test_nonzero_remnant_spin_retains_magnitude_cap(self) -> None:
        self.assertIn("spinmag_new(igrp)=MIN(a1mod,maxspin)", self.block)

    def test_zero_orbit_and_zero_progenitor_spins_are_exactly_zero(self) -> None:
        magnitude, direction = remnant_vector(
            1.0, 0.0, 0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0), 3.0,
        )
        self.assertEqual(magnitude, 0.0)
        self.assertEqual(direction, (0.0, 0.0, 0.0))

    def test_zero_orbit_retains_finite_progenitor_spin_vector(self) -> None:
        magnitude, direction = remnant_vector(
            0.5, 0.4, 0.2, (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0), 3.0,
        )
        raw = (0.4 / 2.25, 0.2 * 0.25 / 2.25, 0.0)
        expected_magnitude = math.hypot(raw[0], raw[1])
        self.assertTrue(math.isfinite(magnitude))
        self.assertAlmostEqual(magnitude, expected_magnitude)
        self.assertEqual(direction[2], 0.0)
        self.assertAlmostEqual(direction[0], raw[0] / expected_magnitude)
        self.assertAlmostEqual(direction[1], raw[1] / expected_magnitude)

    def test_nonzero_orbit_matches_full_reference_expression(self) -> None:
        q, a1, a2, lmodana = 0.5, 0.4, 0.2, 2.5
        spin1, spin2, orbital = (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 2.0)
        magnitude, direction = remnant_vector(q, a1, a2, spin1, spin2,
                                               orbital, lmodana)
        scale = 1.0 / (1.0 + q) ** 2
        expected = (
            scale * (a1 * spin1[0] + a2 * spin2[0] * q**2
                     + orbital[0] / 2.0 * lmodana * q),
            scale * (a1 * spin1[1] + a2 * spin2[1] * q**2
                     + orbital[1] / 2.0 * lmodana * q),
            scale * (a1 * spin1[2] + a2 * spin2[2] * q**2
                     + orbital[2] / 2.0 * lmodana * q),
        )
        expected_magnitude = math.sqrt(sum(component**2 for component in expected))
        self.assertAlmostEqual(magnitude, min(expected_magnitude, 0.998))
        for actual, component in zip(direction, expected):
            self.assertAlmostEqual(actual, component / expected_magnitude)


if __name__ == "__main__":
    unittest.main()
