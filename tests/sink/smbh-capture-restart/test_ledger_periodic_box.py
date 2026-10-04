#!/usr/bin/env python3
"""Guard the per-axis periodic-box contract in capture-ledger events."""

import math
import unittest
from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[3]
    / "patch"
    / "lagRamses"
    / "sink_particle.kjhan.f90"
).read_text(encoding="utf-8")
START = SOURCE.index("\nsubroutine write_smbh_capture_ledger(") + 1
END = SOURCE.index("end subroutine write_smbh_capture_ledger", START)
WRITER = SOURCE[START:END]


def periodic_extents(scale: float, xbound: tuple[float, float, float]) -> tuple[float, ...]:
    extents = tuple(scale * value for value in xbound)
    if any(not math.isfinite(value) or value <= 0.0 for value in extents):
        raise ValueError("non-finite or non-positive periodic box extent")
    return extents


class LedgerPeriodicBoxTest(unittest.TestCase):
    def test_anisotropic_extents_are_scaled_per_axis(self) -> None:
        self.assertEqual(periodic_extents(0.25, (4.0, 8.0, 12.0)), (1.0, 2.0, 3.0))

    def test_invalid_extents_fail_closed(self) -> None:
        for xbound in ((1.0, 0.0, 1.0), (1.0, -1.0, 1.0),
                       (1.0, math.inf, 1.0), (1.0, math.nan, 1.0)):
            with self.subTest(xbound=xbound), self.assertRaises(ValueError):
                periodic_extents(1.0, xbound)

    def test_writer_validates_extents_before_opening_ledger(self) -> None:
        initialization = WRITER.index("periodic_box_size(idim)=scale*xbound(idim)")
        finite_guard = WRITER.index("ieee_is_finite(periodic_box_size(idim))", initialization)
        positive_guard = WRITER.index("periodic_box_size(idim) <= 0d0", finite_guard)
        fatal = WRITER.index("non-finite or non-positive periodic box extent", positive_guard)
        ledger_open = WRITER.index("open(newunit=ledger_unit", fatal)
        self.assertLess(initialization, finite_guard)
        self.assertLess(finite_guard, positive_guard)
        self.assertLess(positive_guard, fatal)
        self.assertLess(fatal, ledger_open)

    def test_event_writes_all_three_axis_extents(self) -> None:
        field = WRITER.index("periodic_box_size_code")
        event_text = WRITER[field:WRITER.index("omega_m", field)]
        for axis in (1, 2, 3):
            self.assertIn(f"json_real(periodic_box_size({axis}))", event_text)

    def test_minimum_image_uses_the_same_validated_extents(self) -> None:
        self.assertEqual(WRITER.count("box_size=periodic_box_size(idim)"), 4)
        self.assertEqual(WRITER.count("box_size=scale*xbound(idim)"), 0)


if __name__ == "__main__":
    unittest.main()
