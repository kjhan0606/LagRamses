#!/usr/bin/env python3
"""Compare compiled fixed-source QUMOND with held-out analytic NFW curves."""
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from ngc2403_cpu512_radial_holdout import (
    C200_GRID, H0, V200_GRID, fit_nfw, nfw_curve, velocity_summary,
)
from ngc2403_lagr_rar import (
    CM_PER_KPC, M_PER_KPC, N_AZIMUTH, empirical_rar,
    read_uniform, trilinear,
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def check_inputs(run_root, gate_root):
    if not (gate_root / 'COMPLETE').is_file():
        raise ValueError(f'independent force analysis is not complete: {gate_root}')
    gate_path = gate_root / 'analysis.json'
    gate = json.loads(gate_path.read_text())
    if gate.get('status') != 'PERIODIC_SOURCE_AND_FORCE_GATE_PASS':
        raise ValueError(f'unexpected force-gate status: {gate.get("status")}')
    if int(gate.get('grid', {}).get('nmesh', -1)) != 512:
        raise ValueError('the compiled run is not the expected 512^3 grid')
    if not math.isclose(float(gate['grid']['box_kpc']), 64.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError('unexpected periodic box size')
    if set(gate.get('cases', {})) != {'newtonian', 'qumond'}:
        raise ValueError('missing matched force cases')
    for model, metrics in gate['cases'].items():
        path = run_root/model/'output_00001/data_00001.h5'
        if Path(metrics['snapshot']).resolve() != path.resolve() or not path.is_file():
            raise ValueError('force analysis refers to another snapshot')
        if not re.fullmatch(r'[0-9a-f]{64}', metrics.get('snapshot_sha256', '')):
            raise ValueError('force analysis lacks a valid snapshot digest')
        names = ['median_relative_force_error', 'p90_relative_force_error']
        if model == 'newtonian':
            names.append('source_rms_relative_error')
        if any(not 0 <= float(metrics.get(name, math.nan)) <= 1e-6 for name in names):
            raise ValueError('force/source analysis threshold not met')
        log = (run_root/model/'run.log').read_text()
        if 'Run completed' not in log or re.search(r'failed to converge|NaN/Inf detected|ERROR:|Aborting', log, re.I):
            raise ValueError('underlying fixed-source simulation did not finish cleanly')

    policy_path = run_root / 'output_policy.json'
    policy = json.loads(policy_path.read_text())
    if policy.get('cases') != 2 or policy.get('nstepmax') != 1:
        raise ValueError('the run does not match the two-case one-step force diagnostic')
    if policy.get('source') != 'same THINGS-like baryon source as compiled 256^3 comparison':
        raise ValueError('unexpected baryon-source contract')

    contract_path = run_root / 'source_contract.json'
    contract = json.loads(contract_path.read_text())
    if contract.get('scenario') != 'things':
        raise ValueError('expected the THINGS-like source contract')
    if not math.isclose(float(contract.get('stellar_ml', math.nan)), 0.5,
                        rel_tol=0.0, abs_tol=1e-12):
        raise ValueError('unexpected stellar mass-to-light ratio')
    if set(gate.get('source_files_sha256', {})) != {'star_profile_code.dat', 'hi_profile_code.dat', 'source_contract.json'}:
        raise ValueError('force analysis source digests are incomplete')
    for name, expected in gate['source_files_sha256'].items():
        candidate = run_root / name
        if not candidate.is_file() or sha256(candidate) != expected:
            raise ValueError(f'force-gate source hash mismatch: {candidate}')
    return gate, policy, contract


def sample_radial_acceleration(force, radius_kpc, meta):
    azimuth = np.arange(N_AZIMUTH, dtype=np.float64) * (2.0 * math.pi / N_AZIMUTH)
    radius_code = radius_kpc * CM_PER_KPC / meta['unit_length_cm']
    if radius_code >= 0.5 * meta['boxlen_code']:
        raise ValueError(f'radius {radius_kpc} kpc reaches the periodic half-box')
    points = np.column_stack((
        0.5 + radius_code * np.cos(azimuth),
        0.5 + radius_code * np.sin(azimuth),
        np.full(N_AZIMUTH, 0.5),
    ))
    unit_radial = np.column_stack((np.cos(azimuth), np.sin(azimuth)))
    acceleration_si = meta['unit_length_cm'] / meta['unit_time_s']**2 / 100.0
    sampled = -np.sum(trilinear(force, points, meta['boxlen_code'])[:2].T * unit_radial,
                      axis=1) * acceleration_si
    if np.any(~np.isfinite(sampled)) or np.median(sampled) <= 0.0:
        raise ValueError(f'invalid median inward acceleration at {radius_kpc} kpc')
    return {
        'median': float(np.median(sampled)),
        'p16': float(np.quantile(sampled, 0.16)),
        'p84': float(np.quantile(sampled, 0.84)),
    }


def save_plot(out, observed, radius, gbar, gq, gnfw, gobs):
    x_obs = np.asarray([float(row['gbar_m_s2']) for row in observed])
    gobs = np.asarray(gobs)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.5, 4.4), layout='constrained')
    for axis in (ax, bx):
        axis.tick_params(labelsize=9)

    xgrid = np.geomspace(min(np.min(x_obs), np.min(gbar)) * 0.8,
                         max(np.max(x_obs), np.max(gbar)) * 1.2, 256)
    ax.scatter(x_obs, gobs, s=16, color='black', label='Observed SPARC')
    ax.scatter(gbar, gq, s=17, marker='s', color='#247f7f',
               label='Compiled QUMOND')
    ax.scatter(gbar, gnfw, s=17, marker='^', color='#d45b25',
               label='NFW halo, held out')
    ax.plot(xgrid, empirical_rar(xgrid), color='0.5', ls='--', lw=1.0,
            label='Empirical RAR reference')
    ax.plot(xgrid, xgrid, color='0.7', ls=':', lw=0.9,
            label='Newtonian baryons')
    ax.set(xscale='log', yscale='log', xlabel=r'$g_{\rm bar}$ (m s$^{-2}$)',
           ylabel=r'$g_R$ (m s$^{-2}$)')
    ax.text(0.03, 0.97, r'NGC 2403 · fixed source · $512^3$',
            transform=ax.transAxes, ha='left', va='top', fontsize=10)
    ax.legend(loc='lower right', frameon=False, fontsize=10)

    q_residual = np.log10(np.asarray(gq) / gobs)
    nfw_residual = np.log10(np.asarray(gnfw) / gobs)
    bx.axhline(0.0, color='0.6', lw=0.8)
    bx.plot(radius, q_residual, 's-', ms=3, lw=0.8, color='#247f7f',
            label='QUMOND / observed')
    bx.plot(radius, nfw_residual, '^-', ms=3, lw=0.8, color='#d45b25',
            label='NFW / observed')
    bx.set(xlabel='Radius (kpc)', ylabel=r'$\log_{10}(g_{\rm model}/g_{\rm obs})$ (dex)')
    bx.text(0.97, 0.97, 'Four contiguous radial holdouts',
            transform=bx.transAxes, ha='right', va='top', fontsize=10)
    bx.legend(loc='upper right', bbox_to_anchor=(0.98, 0.91), frameon=False, fontsize=10)
    fig.savefig(out / 'rar_comparison.pdf')
    fig.savefig(out / 'rar_comparison.png', dpi=180)
    plt.close(fig)


def main():
    if len(sys.argv) not in (4, 5):
        raise SystemExit(
            'usage: script OUTPUT_DIRECTORY LAGRAMSES_RUN OBSERVED_CSV [COMPLETED_FORCE_ANALYSIS]')
    out, run_root, observed_path = map(Path, sys.argv[1:4])
    gate_root = Path(sys.argv[4]) if len(sys.argv) == 5 else run_root
    observed_output = out / 'ngc2403_observed.csv'
    if (not out.is_dir() or not observed_output.is_file() or
            {path.name for path in out.iterdir()} != {observed_output.name} or
            observed_path.resolve() != observed_output.resolve()):
        raise SystemExit(
            f'output must contain only the staged observed CSV: {out}')
    gate, policy, contract = check_inputs(run_root, gate_root)

    observed_sha256 = sha256(observed_path)
    with observed_path.open(newline='') as stream:
        all_rows = list(csv.DictReader(stream))
    observed = [row for row in all_rows
                if row['passes_exploratory_point_cut'] == '1' and
                float(row['radius_kpc']) >= 1.0]
    observed.sort(key=lambda row: float(row['radius_kpc']))
    radius = np.asarray([float(row['radius_kpc']) for row in observed])
    vobs = np.asarray([float(row['vobs_km_s']) for row in observed])
    error = np.asarray([float(row['vobs_error_km_s']) for row in observed])
    gobs = np.asarray([float(row['gobs_m_s2']) for row in observed])
    gbar_sparc = np.asarray([float(row['gbar_m_s2']) for row in observed])
    vgas = np.asarray([float(row['vgas_km_s']) for row in observed])
    vdisk = np.asarray([float(row['vdisk_km_s']) for row in observed])
    vbul = np.asarray([float(row['vbul_km_s']) for row in observed])
    if len(radius) != 64 or np.any(error <= 0.0):
        raise ValueError('expected the documented 64-point sample with positive velocity errors')
    calc_gobs = (vobs * 1000.0)**2 / (radius * M_PER_KPC)
    if not np.allclose(calc_gobs, gobs, rtol=2e-12, atol=0.0):
        raise ValueError('observed gobs does not match Vobs^2/R in the staged CSV')
    calc_gbar = (vgas * np.abs(vgas) + 0.5 * vdisk * np.abs(vdisk) +
                 0.7 * vbul * np.abs(vbul)) * 1.0e6 / (radius * M_PER_KPC)
    if not np.allclose(calc_gbar, gbar_sparc, rtol=2e-12, atol=0.0):
        raise ValueError('observed gbar does not match the signed SPARC mass-model convention')

    metadata = {}
    gbar = np.empty(len(radius))
    gq = np.empty(len(radius))
    gbar_azimuth = []
    gq_azimuth = []
    for model in ('newtonian', 'qumond'):
        force, current = read_uniform(run_root / model)
        if current['snapshot_sha256'] != gate['cases'][model]['snapshot_sha256']:
            raise ValueError('stored force changed since the independent analysis')
        metadata[model] = current
        if current['nmesh'] != 512:
            raise ValueError(f'{model} snapshot has unexpected resolution')
        if model == 'qumond':
            for key in ('boxlen_code', 'unit_length_cm', 'unit_time_s'):
                if not math.isclose(metadata['newtonian'][key], current[key], rel_tol=1e-12):
                    raise ValueError(f'Newtonian and QUMOND snapshots disagree in {key}')
        target = gbar if model == 'newtonian' else gq
        azimuth_records = gbar_azimuth if model == 'newtonian' else gq_azimuth
        for index, rad in enumerate(radius):
            sample = sample_radial_acceleration(force, rad, current)
            target[index] = sample['median']
            azimuth_records.append(sample)
        del force

    vbar2 = gbar * radius * M_PER_KPC / 1.0e6
    q_velocity = np.sqrt(gq * radius * M_PER_KPC) / 1000.0
    nfw_velocity = np.full(len(radius), np.nan)
    fold_results = []
    records = []
    folds = [(16 * index, 16 * (index + 1)) for index in range(4)]
    for fold, (start, stop) in enumerate(folds):
        test_mask = np.zeros(len(radius), dtype=bool)
        test_mask[start:stop] = True
        train_mask = ~test_mask
        v200, c200, train_chi2 = fit_nfw(
            radius[train_mask], vbar2[train_mask], vobs[train_mask], error[train_mask])
        nfw_velocity_all = nfw_curve(radius, vbar2, v200, c200)
        nfw_velocity[test_mask] = nfw_velocity_all[test_mask]
        q_stats = velocity_summary(q_velocity[test_mask], vobs[test_mask], error[test_mask])
        nfw_stats = velocity_summary(nfw_velocity[test_mask], vobs[test_mask], error[test_mask])
        fold_results.append({
            'fold': fold,
            'heldout_radius_kpc': [float(radius[start]), float(radius[stop - 1])],
            'training_points': int(train_mask.sum()),
            'heldout_points': int(test_mask.sum()),
            'nfw_fit': {
                'v200_km_s': v200,
                'c200': c200,
                'training_diagonal_chi2': train_chi2,
                'grid_edge': bool(v200 in (V200_GRID[0], V200_GRID[-1]) or
                                  c200 in (C200_GRID[0], C200_GRID[-1])),
            },
            'heldout_fixed_qumond': q_stats,
            'heldout_nfw': nfw_stats,
        })
        for index in range(start, stop):
            gnfw = float(nfw_velocity[index]**2 * 1.0e6 /
                         (radius[index] * M_PER_KPC))
            records.append({
                'fold': fold,
                'radius_kpc': float(radius[index]),
                'vobs_km_s': float(vobs[index]),
                'vobs_error_km_s': float(error[index]),
                'gobs_m_s2': float(gobs[index]),
                'gbar_sparc_m_s2': float(gbar_sparc[index]),
                'gbar_compiled_newtonian_m_s2': float(gbar[index]),
                'gqumond_compiled_m_s2': float(gq[index]),
                'gnfw_heldout_m_s2': gnfw,
                'vqumond_compiled_km_s': float(q_velocity[index]),
                'vnfw_heldout_km_s': float(nfw_velocity[index]),
                'vbar_compiled_km_s': float(math.sqrt(vbar2[index])),
                'nfw_v200_km_s': v200,
                'nfw_c200': c200,
                'log10_gbar_compiled_over_sparc_dex': float(
                    np.log10(gbar[index] / gbar_sparc[index])),
                'log10_gqumond_over_observed_dex': float(np.log10(gq[index] / gobs[index])),
                'log10_gnfw_over_observed_dex': float(np.log10(gnfw / gobs[index])),
                'gbar_azimuth_p16_m_s2': gbar_azimuth[index]['p16'],
                'gbar_azimuth_p84_m_s2': gbar_azimuth[index]['p84'],
                'qumond_azimuth_p16_m_s2': gq_azimuth[index]['p16'],
                'qumond_azimuth_p84_m_s2': gq_azimuth[index]['p84'],
            })

    if np.any(~np.isfinite(nfw_velocity)):
        raise ValueError('an out-of-fold NFW prediction is missing')
    q_residual = np.log10(gq / gobs)
    nfw_g = nfw_velocity**2 * 1.0e6 / (radius * M_PER_KPC)
    nfw_residual = np.log10(nfw_g / gobs)
    with (out / 'holdout.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    summary = {
        'status': 'COMPILED_512_FIXED_SOURCE_QUMOND_VS_NFW_HOLDOUT',
        'scope': 'Compiled lagRamses uniform 512^3 periodic force fields sampled at 64 SPARC radii. Four contiguous radial blocks are held out. NFW parameters are fit on the other 48 radii at fixed stellar M/L=0.5. QUMOND remains fixed.',
        'run_root': str(run_root.resolve()),
        'force_analysis_root': str(gate_root.resolve()),
        'grid': gate['grid'],
        'output_policy': policy,
        'source_contract': contract,
        'sample': {
            'points': len(radius),
            'radius_kpc': [float(radius.min()), float(radius.max())],
            'selection': 'SPARC exploratory point cut and radius >= 1 kpc',
            'force_sample': f'median inward radial acceleration over {N_AZIMUTH} midplane azimuths per radius',
        },
        'nfw_fit_grid': {
            'h0_km_s_mpc': H0,
            'v200_km_s': [float(V200_GRID[0]), float(V200_GRID[-1]),
                          float(V200_GRID[1] - V200_GRID[0])],
            'c200': [float(C200_GRID[0]), float(C200_GRID[-1]),
                     float(C200_GRID[1] - C200_GRID[0])],
        },
        'out_of_fold_velocity': {
            'fixed_qumond': velocity_summary(q_velocity, vobs, error),
            'analytic_nfw': velocity_summary(nfw_velocity, vobs, error),
        },
        'out_of_fold_rar_residual': {
            'median_log10_gqumond_over_observed_dex': float(np.median(q_residual)),
            'rms_log10_gqumond_over_observed_dex': float(np.sqrt(np.mean(q_residual**2))),
            'median_log10_gnfw_over_observed_dex': float(np.median(nfw_residual)),
            'rms_log10_gnfw_over_observed_dex': float(np.sqrt(np.mean(nfw_residual**2))),
            'median_log10_compiled_gbar_over_sparc_dex': float(
                np.median(np.log10(gbar / gbar_sparc))),
        },
        'folds': fold_results,
        'inputs_sha256': {
            str(gate_root / 'analysis.json'): sha256(gate_root / 'analysis.json'),
            str(run_root / 'output_policy.json'): sha256(run_root / 'output_policy.json'),
            str(run_root / 'source_contract.json'): sha256(run_root / 'source_contract.json'),
            str(observed_path.resolve()): observed_sha256,
            'newtonian_snapshot': metadata['newtonian']['snapshot_sha256'],
            'qumond_snapshot': metadata['qumond']['snapshot_sha256'],
        },
        'limitations': [
            'The NFW component is a two-parameter analytic halo fitted on training radii and evaluated on excluded radii. It is not a live CDM halo and the calculation is not a cosmological Lambda CDM simulation.',
            'The spherical isolated NFW acceleration is added to the periodic baryonic force; this is not a full periodic total-matter Poisson solve.',
            'Each NFW fold fits two parameters on 48 radii while QUMOND remains fixed at a0 and stellar M/L=0.5. The predictive comparison has unequal fitting freedom.',
            'The four folds share most training radii. Fold scores are correlated and do not provide independent model-selection evidence.',
            'Pointwise velocity errors provide diagonal weights. Radial covariance, distance, inclination, stellar-population uncertainty, gas reconstruction, beam effects, and non-circular motions are not propagated.',
            'The SPARC baryonic accelerations use a different gas-source reconstruction from the prescribed lagRamses source. The measured source mismatch remains explicit.',
            'The force fields use a periodic uniform box. The comparison does not establish AMR-boundary accuracy, a continuum limit, galaxy evolution, or cosmological initial-condition consistency.',
            'The Newtonian baryon field is used only to construct a matched-source analytic NFW prediction. The Newtonian baryon-only curve is not identified as the Lambda CDM model.',
        ],
    }
    if sha256(observed_path) != observed_sha256:
        raise ValueError('staged observed CSV changed during analysis')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    save_plot(out, observed, radius, gbar, gq, nfw_g, gobs)
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
