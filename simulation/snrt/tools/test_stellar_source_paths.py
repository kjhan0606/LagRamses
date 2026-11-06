#!/usr/bin/env python3
"""Focused converter routing/projection tests. Execute on Grammar Slurm.

With --source-root and the five source arguments, also exercise the exact
RAR-local raw-source readers. No generated source package is imported.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import build_parsec_mixed_lowmass as mixed
import build_sukhbold_native_lowmass as sn
import read_karakas_lugaro2016 as agb
from stellar_source_paths import SourcePathError, contract_root, source_child, source_path


SOURCE_INPUTS = None


class SourceRoutingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'rar'
        self.root.mkdir()
        self.assets = self.root / 'agb'
        self.assets.mkdir()
        self.matrix = self.root / 'matrix.json'
        self.payload = dict(source_root=str(self.root), candidates=[dict(
            candidate_id='karakas_lugaro2016_agb', source_asset_path='agb',
            approval_id=None, selected_gross_normalization_policy=agb.GROSS_NORMALIZATION_POLICY)])
        self.write_matrix()

    def write_matrix(self):
        self.matrix.write_text(json.dumps(self.payload))

    def test_paths_reject_absolute_parent_and_symlink_escapes(self):
        local = self.assets / 'local'
        local.write_text('local')
        self.assertEqual(source_path('agb/local', root=self.root), local)
        view = self.root / 'view'
        view.mkdir()
        (view/'raw').symlink_to(local)
        self.assertEqual(source_child(view, 'raw', root=self.root), local)
        outside = self.root.parent / 'other-project'
        outside.write_text('external')
        (self.assets / 'link').symlink_to(outside)
        for name in ('../other-project', outside, self.assets / 'link'):
            with self.subTest(name=name), self.assertRaises(SourcePathError):
                source_path(name, root=self.root)
        for name in ('../matrix.json', str(local)):
            for root in (None, self.root):
                with self.subTest(name=name, root=root), self.assertRaises(SourcePathError):
                    source_child(self.assets, name, root=root)
        self.assertEqual(source_child(self.assets, 'link'), outside)
        self.assertEqual(source_child(self.assets, 'link').read_text(), 'external')
        with self.assertRaises(SourcePathError):
            source_child(self.assets, 'link', root=self.root)
        with self.assertRaises(FileNotFoundError):
            source_path('missing', root=self.root).read_bytes()

    def test_contract_root_cannot_be_omitted_or_disagree(self):
        with self.assertRaises(SourcePathError):
            contract_root(self.matrix, {}, default_path=agb.DEFAULT_MATRIX)
        with self.assertRaises(SourcePathError):
            contract_root(self.matrix, self.payload, root=self.assets)
        self.assertIsNone(contract_root(agb.DEFAULT_MATRIX, {}, default_path=agb.DEFAULT_MATRIX))

    def test_selected_external_asset_is_rejected_before_read(self):
        self.payload['candidates'][0]['source_asset_path'] = str(self.root.parent / 'LRD_JWST')
        self.write_matrix()
        with self.assertRaises(SourcePathError):
            agb.read_karakas_lugaro2016(self.matrix)

    def test_scoped_lifetime_has_no_legacy_fallback(self):
        with self.assertRaisesRegex(ValueError, 'explicit lifetime_path'):
            agb.read_karakas_lugaro2016(self.matrix, include_lifetimes=True)
        with self.assertRaisesRegex(ValueError, 'explicit lifetime_path'):
            agb.attach_kl16_lifetimes([], source_root=self.root)
        with self.assertRaises(SourcePathError):
            agb.attach_kl16_lifetimes([], agb.LIFETIME_PATH, source_root=self.root)

    def test_fishlock_routes_explicit_bytes_but_preserves_content_pin(self):
        target = self.root / 'fishlock.txt'
        target.write_text('not the pinned Fishlock table')
        with self.assertRaisesRegex(ValueError, 'fingerprint mismatch'):
            agb.read_fishlock2014(lifetime_model='raiteri96_padova', matrix_path=self.matrix,
                                  yield_path=target)
        with self.assertRaises(FileNotFoundError):
            agb.read_fishlock2014(lifetime_model='raiteri96_padova', matrix_path=self.matrix)

    def test_mixed_rejects_incomplete_local_source_set(self):
        with self.assertRaisesRegex(ValueError, 'KL16 lifetimes and Sukhbold'):
            mixed.terminal_sources(2, matrix_path=self.matrix)

    def test_manifest_hashes_cannot_conflict(self):
        self.assertEqual(mixed.merged_source_hashes({'a': 'one'}, {'a': 'one', 'b': 'two'}),
                         {'a': 'one', 'b': 'two'})
        with self.assertRaisesRegex(ValueError, 'conflicting source fingerprints'):
            mixed.merged_source_hashes({'a': 'one'}, {'a': 'two'})

    def test_sukhbold_contract_rejects_escaping_release_or_files(self):
        contract = json.loads(sn.source.DEFAULT_CONTRACT.read_text())
        path = self.root / 'sukhbold.json'
        contract['source']['release_root_relative_path'] = '../external'
        path.write_text(json.dumps(contract))
        with self.assertRaises(SourcePathError):
            sn.source.audit_sukhbold2016_candidate(root=self.root, contract_path=path)
        contract['source']['release_root_relative_path'] = 'release'
        contract['source']['files'] = {'../outside': {'bytes': 0, 'sha256': ''},
                                     **contract['source']['files']}
        path.write_text(json.dumps(contract))
        with self.assertRaises(SourcePathError):
            sn.source.audit_sukhbold2016_candidate(root=self.root, contract_path=path)
        del contract['source']['files']['INDEX.html']
        path.write_text(json.dumps(contract))
        with self.assertRaisesRegex(ValueError, 'omits a consumed source file'):
            sn.source.audit_sukhbold2016_candidate(root=self.root, contract_path=path)

    def test_sukhbold_relocation_preserves_physical_identity(self):
        contract = json.loads(sn.source.DEFAULT_CONTRACT.read_text())
        contract['source']['release_root_relative_path'] = 'relocated'
        path = self.root / 'sukhbold.json'
        path.write_text(json.dumps(contract))
        release = self.root / 'relocated'
        release.mkdir()
        for name in contract['source']['files']:
            (release / name).write_text('routing fixture; source audit is mocked')

        def record(mass):
            return dict(stable_wind_sum_msun=1., stable_ejecta_sum_msun=mass-2.,
                        baryonic_mass_cut_after_fallback_msun=1., final_kinetic_energy_erg=1e51,
                        stable_wind_by_tracked_element_msun={e: 1. if e == 'H' else 0. for e in sn.ELEMENTS},
                        stable_ejecta_by_tracked_element_msun={e: mass-2. if e == 'H' else 0. for e in sn.ELEMENTS},
                        selected_radioactive_inventory={})

        report = {'z96_grid': {'models': {str(9+i*.25): record(9+i*.25) for i in range(13)}}}
        kwargs = dict(solar_coordinate=.02, wind_speed_km_s=1000., lifetime_model='raiteri96_padova',
                      remnant_policy='stable_segment_residual', max_remnant_difference_msun=.1)
        with patch.object(sn.source, 'audit_sukhbold2016_candidate', return_value=report) as audit, \
             patch.object(sn.tarfile, 'open') as archive, \
             patch.object(sn.source, '_extract_text', return_value='fixture'), \
             patch.object(sn.source, '_parse_z96_results', return_value={m: {} for m in (12.25,12.5,12.75,13.)}), \
             patch.object(sn.source, '_parse_yield_table', side_effect=lambda text, mass, **kw: record(mass)):
            default = sn.build_lowmass(**kwargs)
            provenance = {}
            relocated = sn.build_lowmass(**kwargs, root=self.root, contract_path=path,
                                         source_root=self.root, provenance=provenance)
            self.assertEqual(default, relocated)
            audit.assert_called_with(root=self.root, contract_path=path, source_root=self.root)
            self.assertEqual([call.args[0] for call in archive.call_args_list[-2:]],
                             [release/'explosion_results_PHOTB.tar.gz', release/'nucleosynthesis_yields.tar.gz'])
            self.assertEqual(provenance['release_root'], str(release))
            self.assertEqual(provenance['contract']['path'], str(path))


class ProjectionTests(unittest.TestCase):
    def test_all_element_normalization_preserves_raw_rows(self):
        row = dict(commented_out=False, coordinate={'initial_mass_msun': 2.},
                   selected_mass_expelled_msun=1.5, selected_final_mass_msun=.5,
                   selected_mass_source='yield_table_header',
                   elements_by_atomic_number={z: {'gross_mass_msun': float(z)} for z in range(1,79)})
        original = copy.deepcopy(row)
        result = agb.normalize_selected_ejecta(row)
        self.assertEqual(row, original)
        self.assertEqual(result['returned_mass_msun'], 1.5)
        self.assertEqual(result['remnant_mass_msun'], .5)
        self.assertAlmostEqual(sum(result['gross_mass_msun_by_atomic_number'].values()), 1.5)
        self.assertAlmostEqual(result['gross_mass_msun_by_atomic_number'][26] /
                               result['gross_mass_msun_by_atomic_number'][1], 26.)
        self.assertGreater(result['untracked_ejecta_msun'], 0.)

    def test_existing_endpoint_channels_ties_and_envelope_policy(self):
        row = dict(coordinate=dict(initial_mass_msun=7., metallicity_mass_fraction=.001),
                   evolution=dict(stellar_lifetime_yr=1e8, lifetime_source='fixture', core_kind='ONe'),
                   selected_ejecta=dict(tracked_ejected_mass_msun=[1.]*11,
                                       returned_mass_msun=6., remnant_mass_msun=1.))
        later = copy.deepcopy(row)
        later['evolution']['stellar_lifetime_yr'] = 2e8
        endpoint = mixed.endpoint(8., .001, [row, later], {})
        self.assertEqual(endpoint[3:8], (1e8, 2, 3, 0, 0.))
        self.assertAlmostEqual(endpoint[1]+endpoint[2], 8.)
        self.assertTrue(endpoint[8]['envelope_proxy'])
        snmodel = dict(mass=10., wind_elements=[.1]*11, terminal_elements=[.2]*11,
                       wind=2., terminal=6., remnant=2., age=2e7, terminal_energy=1e51)
        endpoint = mixed.endpoint(10., .001, [], {'models': [snmodel]})
        self.assertEqual(endpoint[1:8], (8., 2., 2e7, 3, 4, 1, 1e51))
        self.assertTrue(all(abs(v-.3)<1e-15 for v in endpoint[0]))

    def test_rar_local_pinned_source_projection(self):
        if SOURCE_INPUTS is None:
            self.skipTest('supply explicit RAR source arguments for the raw-source regression')
        provenance = {}
        rows, selection = mixed.terminal_sources(2, provenance=provenance, **SOURCE_INPUTS)
        kl = agb.read_karakas_lugaro2016(SOURCE_INPUTS['matrix_path'], include_lifetimes=True,
            source_root=SOURCE_INPUTS['source_root'], lifetime_path=SOURCE_INPUTS['lifetime_path'])
        self.assertEqual(len(kl['records']), 62)
        self.assertEqual(len(kl['excluded_commented_records']), 2)
        sample = next(r for r in kl['records'] if r['coordinate']['initial_mass_msun']==3.5
                      and r['coordinate']['metallicity_mass_fraction']==.03)
        self.assertEqual(sample['selected_ejecta']['returned_mass_msun'], 2.773)
        self.assertEqual(sample['selected_ejecta']['remnant_mass_msun'], .727)
        overshoot = next(r for r in kl['records'] if r['coordinate']['initial_mass_msun']==1.75
                         and r['coordinate']['metallicity_mass_fraction']==.014)
        self.assertEqual(overshoot['evolution']['stellar_lifetime_yr'], 1755e6)
        self.assertEqual(len(selection['models']), 17)
        self.assertEqual(sum(r.get('overshoot_label')=='Fishlock2014_original' for r in rows), 16)
        self.assertEqual(provenance['fishlock']['sha256'], agb.FISHLOCK_SHA256)
        self.assertEqual(provenance['kl16']['lifetimes']['sha256'], agb.LIFETIME_SHA256)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('source-root', 'agb-matrix', 'kl16-lifetimes', 'fishlock-yields',
                   'sukhbold-root', 'sukhbold-contract'):
        parser.add_argument('--'+option, type=Path)
    args, remaining = parser.parse_known_args()
    if any(value is not None for value in vars(args).values()):
        if any(value is None for value in vars(args).values()):
            parser.error('raw-source regression requires all six source arguments')
        SOURCE_INPUTS = dict(source_root=args.source_root, matrix_path=args.agb_matrix,
            lifetime_path=args.kl16_lifetimes, fishlock_path=args.fishlock_yields,
            sukhbold_root=args.sukhbold_root, sukhbold_contract=args.sukhbold_contract)
    unittest.main(argv=[__file__, *remaining])
