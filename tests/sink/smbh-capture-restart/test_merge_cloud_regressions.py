#!/usr/bin/env python3
"""Source-level guards for failures exercised by the two-phase smoke test."""

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

    def test_sink_and_cloud_removal_use_actual_parent_grids(self) -> None:
        self.assertEqual(
            SOURCE.count("ind_parent_grid(j)=ind_grid(ind_grid_part(j))"),
            2,
        )
        self.assertEqual(
            SOURCE.count("call remove_list(ind_part,ind_parent_grid,ok,np)"),
            2,
        )
        self.assertNotIn(
            "call remove_list(ind_part,ind_grid_part,ok,np)",
            SOURCE,
        )

if __name__ == "__main__":
    unittest.main()
