"""Focused converter tests using pinned real sources; never launch RAMSES."""
import argparse
import io
import json
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from build_parsec_pair_feedback import (BISTABILITY_MODEL, BISTABILITY_SELECTOR,
    MSUN, PARSEC_LSUN, Z_GRIDS, bistability_wind_speed, phase_wind_speed,
    source_name, wind_history, read_ejecta, project)


def elementary():
    temperatures = [8000, 10000, 10001, 12500, 22500, 25000, 27500, 50000,
                    25000, 25000]
    track = np.zeros((len(temperatures), 38))
    track[:, 0] = 40
    track[:, 3] = 5
    track[:, 4] = np.log10(temperatures)
    track[:, 29] = .7
    track[-2:, 29] = [.399, .4]
    track[:, 5] = np.sqrt(1e5*PARSEC_LSUN /
                          (4*np.pi*5.670374419e-5*(10**track[:, 4])**4))
    old, _, _ = phase_wind_speed(track, .014)
    new, branches, _ = bistability_wind_speed(track, .014)
    np.testing.assert_allclose(new/old, [1, 1, .5, .5, .5, .75, 1, 1, 1, .75],
                               rtol=5e-15, atol=0)
    # log10 roundoff at either ramp endpoint can affect its diagnostic label,
    # but not its continuous velocity. All five branches must be represented.
    assert set(branches) == set(range(5))
    assert branches[8] == 1 and branches[9] == 4
    for z in [.008, .03]:
        bistability_wind_speed(track, z)
    for kind in ['low_z', 'high_z', 'nan_z', 'radius', 'gamma', 'hydrogen', 'nan']:
        bad, z = track.copy(), .014
        if kind == 'low_z': z = .007999
        if kind == 'high_z': z = .030001
        if kind == 'nan_z': z = np.nan
        if kind == 'radius': bad[3, 5] *= 2
        if kind == 'gamma':
            bad[3, 3] += 2
            bad[3, 5] *= 10
        if kind == 'hydrogen': bad[3, 29] = 1.1
        if kind == 'nan': bad[3, 0] = np.nan
        try:
            bistability_wind_speed(bad, z)
        except ValueError:
            pass
        else:
            raise AssertionError(kind)
    print('ELEMENTARY_LAW_BOUNDARIES_AND_ADMISSION_PASS')


def check(source, package, baseline):
    manifest = json.loads((package/'manifest.json').read_text())
    old_manifest = json.loads((baseline/'manifest.json').read_text())
    assert manifest['model'] == BISTABILITY_MODEL
    assert manifest['wind_model']['id'] == BISTABILITY_SELECTOR
    assert manifest['wind_model']['transition_temperature_K'] == [22500, 27500]
    history = (package/'history.nml').read_text()
    for token in ['version=4', BISTABILITY_MODEL, "timing_policy='wind_linear_terminal_step'",
                  'net_channel_available=T,F,T,F,T', 'input_imf_max=600', 'mass_domain_min=14']:
        assert token in history, token
    old_nodes = {(n['mass'], n['z']): n for n in old_manifest['nodes']}
    rows, old_rows = np.loadtxt(package/'yields.dat'), np.loadtxt(baseline/'yields.dat')
    # Terminal histories must be bitwise identical; only channel1 energy and
    # energy-selected compression knots may change.
    np.testing.assert_array_equal(rows[rows[:, 0] != 1], old_rows[old_rows[:, 0] != 1])
    max_error, changed, affected_mass = 0., 0, 0.
    branches_seen = np.zeros(5)
    grid = manifest.get('metallicity_grid', 'solar_pair')
    assert len(manifest['nodes']) == 45*len(Z_GRIDS[grid])
    with zipfile.ZipFile(source/'all_ejecta.zip') as ejecta:
        for zs, ys in Z_GRIDS[grid]:
            initial, wind = read_ejecta(ejecta, f'ejecta/Z{zs}_Y{ys}_winds_ejecta.dat')
            with zipfile.ZipFile(source/source_name(zs, ys, 'tracks')) as archive:
                for node in manifest['nodes']:
                    if node['z'] != float(zs): continue
                    mass = node['mass']
                    track = np.loadtxt(io.BytesIO(archive.read(source_name(zs, ys, 'tracks', mass))), skiprows=3)
                    m, t, x = track[:, 0], 10**track[:, 4], track[:, 29]
                    lum = 10**track[:, 3]*PARSEC_LSUN
                    radius = np.sqrt(lum/(4*np.pi*5.670374419e-5*t**4))
                    gamma = .2*(1+x)*lum/(4*np.pi*2.99792458e10*6.67430e-8*m*MSUN)
                    # Independent scalar piecewise definition, not the helper.
                    d = np.array([1.3 if v <= 22500 else 2.6 if v >= 27500 else
                                  1.3 + (v-22500)*1.3/5000 for v in t])
                    d[x < .4] = 1.6
                    v2 = d**2*2*6.67430e-8*m*MSUN/radius*(1-gamma)*(float(zs)/.02)**.26
                    v2[t <= 10000] = 1e12
                    dm = -np.diff(m)*node['wind_mass']/(m[0]-m[-1])
                    energy = .25*MSUN*np.sum(dm*(v2[1:]+v2[:-1]))
                    max_error = max(max_error, abs(energy/node['wind_energy_erg']-1))
                    original = old_nodes[mass, node['z']]
                    for key in ['wind_mass', 'terminal_mass', 'baryonic_remnant',
                                'terminal_energy', 'age_yr', 'fate', 'he_core']:
                        assert node[key] == original[key], (mass, node['z'], key)
                    assert node['wind_energy_erg'] <= original['wind_energy_erg']*(1+2e-14)
                    changed += node['wind_energy_erg'] < original['wind_energy_erg']*(1-1e-12)
                    fractions = np.array(node['closure_branch_mass_fraction'])
                    assert len(fractions) == 5 and abs(fractions.sum()-1) < 2e-13
                    assert abs(sum(node['closure_branch_energy_fraction'])-1) < 2e-13
                    branches_seen += fractions
                    affected_mass += node['wind_mass']*sum(fractions[3:])
                    assert node['energy_in_compression_bound']
                    assert node['max_cumulative_error_per_final_component'] <= manifest['wind_tolerance_per_final_component']
                    wind_rows = rows[(rows[:, 0] == 1) & (rows[:, 1] == mass) & (rows[:, 2] == node['z'])]
                    old_wind = old_rows[(old_rows[:, 0] == 1) & (old_rows[:, 1] == mass) & (old_rows[:, 2] == node['z'])]
                    np.testing.assert_array_equal(np.delete(wind_rows[-1], 6), np.delete(old_wind[-1], 6))
                    assert wind_rows[-1, 6] == node['wind_energy_erg']
                    # All knots refer to real source times (plus origin).
                    assert np.isin(wind_rows[:, 3], np.r_[0., track[:, 1]]).all()
                    cumulative_energy = np.r_[0., np.cumsum(.25*MSUN*dm*(v2[1:]+v2[:-1]))]
                    unique = np.r_[np.diff(track[:, 1]) > 0, True]
                    reconstructed = np.interp(track[unique, 1], wind_rows[:, 3], wind_rows[:, 6])
                    assert np.max(np.abs(reconstructed-cumulative_energy[unique]))/energy <= \
                        manifest['wind_tolerance_per_final_component']+1e-12
                    if mass == 40 and node['z'] == .014:
                        _, constant, _ = wind_history(track, mass, node['wind_mass'],
                            project(wind[mass]), initial, 1e-4, np.full(len(track), 1e8))
                        np.testing.assert_allclose(constant[:, 12], .5*MSUN*constant[:, 0]*1e16,
                                                   rtol=2e-13, atol=0)
    assert max_error < 1e-12 and changed > 0 and np.all(branches_seen > 0)
    print('ACTUAL_SOURCE_CONVERSION_PASS nodes=', len(manifest['nodes']),
          'energy_changed_nodes=', changed, 'max_relative_energy_error=', max_error,
          'affected_unweighted_wind_mass_Msun=', affected_mass)
    print('MATERIAL_ENDPOINTS_TERMINAL_HISTORIES_AND_NATIVE_V4_FLAGS_PASS')


def legacy_bytes(regenerated, retained):
    # Manifest is included, not just the native payload.
    for name in ['yields.dat', 'history.nml', 'manifest.json']:
        assert (regenerated/name).read_bytes() == (retained/name).read_bytes(), name
    print('LEGACY_ALL_FILES_BYTE_IDENTICAL', regenerated)


def cli_rejections(source, package):
    converter = Path(__file__).resolve().parents[1]/'tools/build_parsec_pair_feedback.py'
    output = package.parent/'must-not-be-created'
    assert not output.exists()
    for extra, message in [
            (['--wind-model', BISTABILITY_SELECTOR, '--wind-km-s', '1000'], 'do not supply'),
            (['--wind-model', BISTABILITY_SELECTOR, '--metallicity-grid', 'precision_eleven'], 'no extreme-Z'),
            (['--wind-model', 'not_a_model'], 'invalid choice')]:
        result = subprocess.run([sys.executable, '-B', str(converter), '--source-dir', str(source),
                                 '--output', str(output), *extra], capture_output=True, text=True)
        assert result.returncode != 0 and message in result.stderr, result.stderr
        assert not output.exists()
    print('CLI_CONFLICT_EXTREME_Z_UNKNOWN_SELECTOR_REJECTIONS_PASS')


def sed_interface(sed, baseline_sed):
    lines = (sed/'nodes.dat').read_text().splitlines()
    old = (baseline_sed/'nodes.dat').read_text().splitlines()
    assert lines[0] == 'SNRT_PARSEC_QE_INTERVAL_V1'
    assert lines[2] == BISTABILITY_MODEL and len(lines[2]) <= 64
    assert lines[:2]+lines[3:] == old[:2]+old[3:]
    nml = (sed/'source.nml').read_text()
    for token in ['version=4', "population_binding='match_feedback_high_mass_only'",
                  "energy_semantics='photon_number_and_energy_v1'", 'imf_max=600',
                  "interpolation='linear_cumulative_age_linear_Z'"]:
        assert token in nml, token
    print('MATCHED_SED_UNCHANGED_EXCEPT_MODEL_NATIVE_INTERFACE_FLAGS_PASS')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--package', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--legacy-fixed', type=Path, nargs=2)
    p.add_argument('--legacy-phase', type=Path, nargs=2)
    p.add_argument('--sed', type=Path, nargs=2, help='new and retained phase SED directories')
    a = p.parse_args()
    elementary()
    check(a.source, a.package, a.baseline)
    for pair in [a.legacy_fixed, a.legacy_phase]:
        if pair: legacy_bytes(*pair)
    cli_rejections(a.source, a.package)
    if a.sed: sed_interface(*a.sed)
