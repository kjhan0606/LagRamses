#!/usr/bin/env python3
"""Regressions for concise physics switches expanded into namelist fields."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'patch' / 'cuRamses' / 'aux'))

import ramses_nml_generator as rng  # noqa: E402


class SimplePhysicsOptionTests(unittest.TestCase):
    def test_dmo_disables_hydro_physics_as_a_complete_mode(self):
        values = rng.apply_simple_physics_options({
            'hydro': True,
            'sink': True,
            'create_sinks': True,
            'agn': True,
            'sink_AGN': True,
            'random_jet': True,
            'haardt_madau': True,
            'radiation_transport': 'snrt_mn',
            'dust_mass_enabled': True,
            'ir_feedback': True,
            'mhd_enabled': True,
        }, dmo=True)

        for key in ('cosmo', 'pic', 'poisson'):
            self.assertTrue(values[key], key)
        self.assertFalse(values['hydro'])
        for key in ('sink', 'create_sinks', 'agn', 'sink_AGN', 'random_jet', 'haardt_madau',
                    'dust_mass_enabled', 'ir_feedback', 'mhd_enabled'):
            self.assertFalse(values[key], key)
        self.assertEqual(values['radiation_transport'], 'none')

    def test_hydro_and_independent_off_switches_preserve_uvb_by_default(self):
        values = rng.apply_simple_physics_options({
            'haardt_madau': True,
            'radiation_transport': 'snrt_mn',
            'snrt_moment_order': 5,
            'dust_mass_enabled': True,
            'dust_growth': True,
            'ir_feedback': True,
            'sink': True,
            'create_sinks': True,
            'sink_AGN': True,
            'bondi': True,
        }, hydro=True, no_dust=True, no_ray=True, no_agn=True)

        for key in ('cosmo', 'pic', 'poisson', 'hydro'):
            self.assertTrue(values[key], key)
        self.assertEqual(values['radiation_transport'], 'none')
        self.assertEqual(values['snrt_moment_order'], 3)
        self.assertTrue(values['haardt_madau'])
        for key in ('dust_mass_enabled', 'dust_growth', 'ir_feedback',
                    'sink', 'create_sinks', 'sink_AGN', 'bondi'):
            self.assertFalse(values[key], key)

        no_uvb = rng.apply_simple_physics_options(values, no_uvb=True)
        self.assertFalse(no_uvb['haardt_madau'])

        # Isolate --no-agn so dust disabling cannot mask its IR-feedback effect.
        agn_off = rng.apply_simple_physics_options({
            'hydro': True,
            'radiation_transport': 'snrt_mn',
            'haardt_madau': True,
            'dust_mass_enabled': True,
            'sink': True,
            'create_sinks': True,
            'agn': True,
            'sink_AGN': True,
            'bondi': True,
            'drag': True,
            'spin_bh': True,
            'mad_jet': True,
            'random_jet': True,
            'ir_feedback': True,
        }, no_agn=True)
        for key in ('sink', 'create_sinks', 'agn', 'sink_AGN', 'bondi',
                    'drag', 'spin_bh', 'mad_jet', 'random_jet', 'ir_feedback'):
            self.assertFalse(agn_off[key], key)
        self.assertEqual(agn_off['radiation_transport'], 'snrt_mn')
        self.assertTrue(agn_off['haardt_madau'])
        self.assertTrue(agn_off['dust_mass_enabled'])
        rendered = rng.format_namelist(agn_off, emit_defaults=True).lower()
        self.assertRegex(rendered, r'\bagn\s*=\s*\.false\.')
        self.assertRegex(rendered, r'\bir_feedback\s*=\s*\.false\.')

    def test_basic_species_selects_existing_grackle_neq_route(self):
        table = '/fixture/CloudyData_noUVB.h5'
        values = rng.apply_simple_physics_options(
            {'haardt_madau': False, 'radiation_transport': 'snrt_sn',
             'dust_mass_enabled': True},
            nspecies='basic', grackle_data_file=table)

        self.assertTrue(values['hydro'])
        self.assertTrue(values['cooling'])
        self.assertTrue(values['metal'])
        self.assertEqual(values['cooling_method'], 'grackle_neq')
        self.assertFalse(values['neq_chem'])
        self.assertEqual(values['radiation_transport'], 'none')
        self.assertFalse(values['dust_mass_enabled'])
        self.assertEqual(values['grackle_neq_data_file'], table)
        self.assertFalse(values['haardt_madau'])

    def test_conflicting_or_incomplete_shortcuts_are_rejected(self):
        invalid_calls = (
            {'dmo': True, 'hydro': True},
            {'dmo': True, 'nspecies': 'basic'},
            {'grackle_data_file': '/fixture/rates.h5'},
            {'nspecies': 'extended'},
        )
        for options in invalid_calls:
            with self.subTest(options=options), self.assertRaises(ValueError):
                rng.apply_simple_physics_options({}, **options)


if __name__ == '__main__':
    unittest.main()
