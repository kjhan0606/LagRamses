#!/usr/bin/env python3
"""Source-level guards for failures exercised by the two-phase smoke test."""

import re
import unittest
from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[3]
    / "patch"
    / "lagRamses"
    / "sink_particle.kjhan.f90"
).read_text(encoding="utf-8")


class MergeCloudRegressionTest(unittest.TestCase):
    def test_cloud_particles_are_linked_to_actual_parent_grids(self) -> None:
        self.assertIn(
            "ind_cloud_grid(j)=ind_grid(ind_grid_part(j))",
            SOURCE,
        )
        self.assertIn(
            "call add_list(ind_cloud,ind_cloud_grid,ok_true,np)",
            SOURCE,
        )
        self.assertNotIn(
            "call add_list(ind_cloud,ind_grid_part,ok_true,np)",
            SOURCE,
        )

    def test_zero_angular_momentum_spin_update_has_no_zero_norm_division(self) -> None:
        merger = re.search(
            r"subroutine merge_sink\(ilevel\)(.*?)end subroutine merge_sink",
            SOURCE,
            flags=re.DOTALL | re.IGNORECASE,
        )
        self.assertIsNotNone(merger)
        body = merger.group(1).lower()
        self.assertRegex(body, r"if\s*\(lmod>0d0\)\s*then")
        self.assertRegex(body, r"if\s*\(a1mod>0d0\)\s*then")
        self.assertIn("bhspin_new(igrp,1:3)=0d0", body)


if __name__ == "__main__":
    unittest.main()
