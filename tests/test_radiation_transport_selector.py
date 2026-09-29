#!/usr/bin/env python3
"""Small regressions for the single radiation_transport namelist control."""
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'patch' / 'cuRamses' / 'aux'))

import mkrun  # noqa: E402
import ramses_nml_generator as rng  # noqa: E402


class RadiationTransportSelectorTests(unittest.TestCase):
    def test_generator_emits_one_transport_selector_and_separate_physics(self):
        modes = ('none', 'snrt_sn', 'snrt_mn', 'ramses_rt', 'aton')
        for mode in modes:
            values = {'radiation_transport': mode}
            if mode == 'snrt_mn':
                values['snrt_moment_order'] = 5
            rendered = rng.format_namelist(values)
            self.assertIn("radiation_transport='{}'".format(mode), rendered)
            self.assertIsNone(re.search(
                r'(?im)^\s*(?:rt|aton|snrt_transport_model)\s*=', rendered
            ))

        physics = rng.format_namelist({
            'radiation_transport': 'none',
            'haardt_madau': True,
            'neq_chem': True,
            'dust_mass_enabled': True,
        })
        for entry in ("radiation_transport='none'", 'haardt_madau=.true.',
                      'neq_chem=.true.', 'dust_mass_enabled=.true.'):
            self.assertIn(entry, physics)

    def test_comparison_choice_overrides_template_and_drops_irrelevant_order(self):
        template = (
            "&RUN_PARAMS\n"
            "  radiation_transport='snrt_sn'\n"
            "  snrt_moment_order=3\n"
            "/\n"
        )
        mn = {'radiation_transport': 'snrt_mn', 'snrt_moment_order': 5}
        rendered = mkrun.set_group_values(template, 'RUN_PARAMS', mn)
        rendered = mkrun.merge_into_group(rendered, 'RUN_PARAMS', list(mn), mn)
        self.assertIn("radiation_transport='snrt_mn'", rendered)
        self.assertIn('snrt_moment_order=5', rendered)

        sn = {'radiation_transport': 'snrt_sn'}
        rendered = mkrun.set_group_values(
            rendered, 'RUN_PARAMS', sn, remove=('snrt_moment_order',)
        )
        rendered = mkrun.merge_into_group(rendered, 'RUN_PARAMS', list(sn), sn)
        self.assertIn("radiation_transport='snrt_sn'", rendered)
        self.assertNotIn('snrt_moment_order', rendered)
        self.assertEqual(rendered.count('radiation_transport='), 1)

    def test_duplicate_selector_is_rejected(self):
        duplicate = "&RUN_PARAMS\nradiation_transport='none'\nradiation_transport='snrt_sn'\n/\n"
        with self.assertRaisesRegex(ValueError, 'duplicate radiation_transport'):
            mkrun.set_group_values(
                duplicate, 'RUN_PARAMS', {'radiation_transport': 'snrt_sn'}
            )


if __name__ == '__main__':
    unittest.main()
