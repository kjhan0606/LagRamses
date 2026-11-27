"""A module dependency must not become the default production target."""
from pathlib import Path
import subprocess
import unittest


class MakeDefaultTargetTest(unittest.TestCase):
    def test_default_is_ramses_and_sidm_geometry_is_linked(self):
        root = Path(__file__).resolve().parents[2]
        # -q never runs a recipe; an out-of-date target has exit status 1.
        result = subprocess.run(['make', '-C', str(root / 'bin'), '-qp'],
                                text=True, capture_output=True)
        self.assertIn(result.returncode, (0, 1), result.stderr)
        self.assertIn('.DEFAULT_GOAL := ramses\n', result.stdout)
        self.assertIn('SIDM_OBJ = sidm_scatter.o sidm_cell_geometry.o ', result.stdout)
        self.assertIn('sidm_cell_geometry.o: pm_commons.o amr_commons.kjhan.o', result.stdout)


if __name__ == '__main__':
    unittest.main()
