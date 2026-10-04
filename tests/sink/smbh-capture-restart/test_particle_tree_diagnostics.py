#!/usr/bin/env python3
"""Guard the bounded level-1 particle-tree failure diagnostic."""

import unittest
from pathlib import Path


SOURCE = (
    Path(__file__).resolve().parents[3]
    / "patch"
    / "lagRamses"
    / "particle_tree.kjhan.f90"
).read_text(encoding="utf-8")


class ParticleTreeDiagnosticTest(unittest.TestCase):
    def test_make_tree_distinguishes_source_from_staging_failure(self) -> None:
        self.assertIn("MAKE_TREE_SOURCE_GRID_INVALID", SOURCE)
        self.assertIn("MAKE_TREE_STAGE_GRID_MISMATCH", SOURCE)
        self.assertIn("if(igrid<1.or.igrid>ngridmax)then", SOURCE)
        self.assertEqual(SOURCE.count("if(ind_grid(1)/=igrid)then"), 2)
        self.assertIn("' staged_grid=',ind_grid(1)", SOURCE)

    def test_level1_parent_failure_is_bounded_and_fail_closed(self) -> None:
        marker = SOURCE.index("CHECK_TREE_LEVEL1_PARENT_INVALID")
        diagnostic = SOURCE.rindex("tree_parent_invalid=.false.", 0, marker)
        neighbor_lookup = SOURCE.index("call get3cubefather", diagnostic)
        block = SOURCE[diagnostic:neighbor_lookup]

        self.assertIn("igrid_diag<1.or.igrid_diag>ngridmax", block)
        self.assertIn("ifather_diag<1.or.ifather_diag>ncoarse", block)
        self.assertIn("min(5,max(0,numbp(igrid_diag)))", block)
        self.assertIn("idp(ipart_diag)", block)
        self.assertIn("xp(ipart_diag,1:ndim)", block)
        self.assertIn("levelp(ipart_diag)", block)
        self.assertIn("if(tree_parent_invalid)call clean_stop", block)


if __name__ == "__main__":
    unittest.main()
