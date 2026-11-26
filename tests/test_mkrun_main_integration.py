"""Regression for the text/template entry point retained when merging main."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'patch' / 'cuRamses' / 'aux'))
import ramses_nml_generator as rng


class MainMergeTest(unittest.TestCase):
    def run_generator(self, *args):
        return subprocess.run([sys.executable, str(ROOT / 'mkrun.py'), *args],
                              cwd=ROOT, text=True, capture_output=True, check=True)

    def test_template_preserves_schedule_and_explicit_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self.run_generator('--mode', 'text', '--outdir', str(tmp / 'base'),
                               '--name', 'base', '--ic', 'none', '--boxlen', '50',
                               '--set', 'courant_factor=0.7', '--redshifts', '1,0')
            base = tmp / 'base' / 'base.nml'
            self.run_generator('--mode', 'text', '--outdir', str(tmp / 'derived'),
                               '--name', 'derived', '--ic', 'none', '--template',
                               str(base), '--set', 'courant_factor=0.6')
            original, _ = rng.parse_namelist(base.read_text())
            derived, _ = rng.parse_namelist((tmp / 'derived' / 'derived.nml').read_text())
            self.assertEqual(original['noutput'], '2')
            self.assertEqual(derived['noutput'], '2')
            self.assertEqual(derived['aout'], original['aout'])
            self.assertEqual(float(derived['boxlen']), 50.0)
            self.assertEqual(float(derived['courant_factor']), 0.6)

    def test_text_help_is_reachable_without_interaction(self):
        result = self.run_generator('--mode', 'text', '--help')
        self.assertIn('--template', result.stdout)
        self.assertIn('--set', result.stdout)


if __name__ == '__main__':
    unittest.main()
