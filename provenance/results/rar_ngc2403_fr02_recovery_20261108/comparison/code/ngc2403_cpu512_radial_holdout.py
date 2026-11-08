#!/usr/bin/env python3
"""Contiguous radial holdout for NGC 2403 CPU reference force curves."""
import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

H0 = 70.0
M_PER_KPC = 3.0856775814913673e19
G_DAGGER = 1.2e-10
V200_GRID = np.linspace(30.0, 250.0, 881)
C200_GRID = np.linspace(2.0, 30.0, 561)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def nfw_f(x):
    return np.log1p(x) - x / (1.0 + x)


def empirical_rar(gbar):
    return gbar / (-np.expm1(-np.sqrt(gbar / G_DAGGER)))


def nfw_curve(radius, vbar2, v200, c200):
    x = radius * H0 / (100.0 * v200)
    vhalo2 = v200**2 * nfw_f(c200 * x) / (x * float(nfw_f(c200)))
    return np.sqrt(vbar2 + vhalo2)


def velocity_summary(predicted, observed, error):
    residual = predicted - observed
    return {
        'rms_km_s': float(np.sqrt(np.mean(residual**2))),
        'median_absolute_km_s': float(np.median(np.abs(residual))),
        'mean_residual_km_s': float(np.mean(residual)),
        'diagonal_chi2': float(np.sum((residual / error)**2)),
        'diagonal_chi2_per_point': float(np.mean((residual / error)**2)),
    }


def fit_nfw(radius, vbar2, observed, error, *, v200_grid=None, c200_grid=None):
    v200_grid = V200_GRID if v200_grid is None else np.asarray(v200_grid)
    c200_grid = C200_GRID if c200_grid is None else np.asarray(c200_grid)
    best_score = math.inf
    best_v200 = None
    best_c200 = None
    c = c200_grid[None, :, None]
    fc = nfw_f(c200_grid)[None, :, None]
    for start in range(0, len(v200_grid), 24):
        v200 = v200_grid[start:start + 24]
        x = radius[None, None, :] * H0 / (100.0 * v200[:, None, None])
        vhalo2 = v200[:, None, None]**2 * nfw_f(c * x) / (x * fc)
        model = np.sqrt(vbar2[None, None, :] + vhalo2)
        score = np.sum(((model - observed[None, None, :]) /
                        error[None, None, :])**2, axis=2)
        iv, ic = np.unravel_index(int(np.argmin(score)), score.shape)
        candidate = float(score[iv, ic])
        if candidate < best_score:
            best_score = candidate
            best_v200 = float(v200[iv])
            best_c200 = float(c200_grid[ic])

    return best_v200, best_c200, best_score


def main():
    if len(sys.argv) != 5:
        raise SystemExit('usage: script OUTPUT CPU512_RADII BASELINE_RAR SOURCE_CONTRACT')
    out, reference_path, baseline_path, contract_path = map(Path, sys.argv[1:])
    if not out.is_dir() or any(out.iterdir()):
        raise SystemExit(f'output must exist and be empty: {out}')

    with reference_path.open(newline='') as stream:
        reference = list(csv.DictReader(stream))
    with baseline_path.open(newline='') as stream:
        baseline = [row for row in csv.DictReader(stream) if row['scenario'] == 'THINGS']
    if len(reference) != 64 or len(baseline) != 64:
        raise ValueError(f'expected 64 matched radii, got {len(reference)} and {len(baseline)}')
    reference.sort(key=lambda row: float(row['radius_kpc']))
    baseline.sort(key=lambda row: float(row['radius_kpc']))

    radius = np.asarray([float(row['radius_kpc']) for row in reference])
    radius_base = np.asarray([float(row['radius_kpc']) for row in baseline])
    if not np.array_equal(radius, radius_base):
        raise ValueError('CPU reference and compiled baseline radii differ')
    observed = np.asarray([float(row['vobs_km_s']) for row in baseline])
    error = np.asarray([float(row['vobs_error_km_s']) for row in baseline])
    gobs = np.asarray([float(row['gobs_m_s2']) for row in baseline])
    if not np.allclose(np.asarray([float(row['vobs_km_s']) for row in reference]),
                       observed, rtol=0.0, atol=1e-13):
        raise ValueError('observed velocities differ between source files')
    if not np.allclose(np.asarray([float(row['gobs_m_s2']) for row in reference]),
                       gobs, rtol=1e-12, atol=1e-22):
        raise ValueError('observed accelerations differ between source files')
    if np.any(error <= 0.0):
        raise ValueError('observed velocity errors must be positive')

    folds = [(16 * k, 16 * (k + 1)) for k in range(4)]
    outputs = {}
    records = []
    for nmesh in (256, 512):
        gbar = np.asarray([float(row[f'gbar_{nmesh}_m_s2']) for row in reference])
        gq = np.asarray([float(row[f'gq_{nmesh}_m_s2']) for row in reference])
        vq = np.asarray([float(row[f'vq_{nmesh}_km_s']) for row in reference])
        vbar2 = gbar * radius * M_PER_KPC / 1.0e6
        nfw_prediction = np.full(64, np.nan)
        fold_results = []
        for fold, (start, stop) in enumerate(folds):
            test_mask = np.zeros(64, dtype=bool)
            test_mask[start:stop] = True
            train_mask = ~test_mask
            v200, c200, train_chi2 = fit_nfw(
                radius[train_mask], vbar2[train_mask],
                observed[train_mask], error[train_mask])
            model_all = nfw_curve(radius, vbar2, v200, c200)
            nfw_prediction[test_mask] = model_all[test_mask]
            holdout_q = velocity_summary(vq[test_mask], observed[test_mask], error[test_mask])
            holdout_nfw = velocity_summary(nfw_prediction[test_mask], observed[test_mask], error[test_mask])
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
                'heldout_fixed_qumond': holdout_q,
                'heldout_nfw': holdout_nfw,
            })
            for index in range(start, stop):
                vnfw = float(model_all[index])
                gnfw = vnfw**2 * 1.0e6 / (radius[index] * M_PER_KPC)
                records.append({
                    'nmesh': nmesh,
                    'fold': fold,
                    'radius_kpc': float(radius[index]),
                    'vobs_km_s': float(observed[index]),
                    'vobs_error_km_s': float(error[index]),
                    'vqumond_km_s': float(vq[index]),
                    'vnfw_km_s': vnfw,
                    'v200_km_s': v200,
                    'c200': c200,
                    'gbar_m_s2': float(gbar[index]),
                    'gqumond_m_s2': float(gq[index]),
                    'gnfw_m_s2': float(gnfw),
                    'log10_gqumond_over_gobs_dex': float(np.log10(gq[index] / gobs[index])),
                    'log10_gnfw_over_gobs_dex': float(np.log10(gnfw / gobs[index])),
                })

        aggregate_q = velocity_summary(vq, observed, error)
        aggregate_nfw = velocity_summary(nfw_prediction, observed, error)
        outputs[str(nmesh)] = {
            'grid': {'nmesh': nmesh, 'box_kpc': 64.0, 'cell_width_kpc': 64.0 / nmesh},
            'folds': fold_results,
            'all_radii_out_of_fold_predictions': {
                'fixed_qumond': aggregate_q,
                'nfw': aggregate_nfw,
                'median_log10_gqumond_over_gobs_dex': float(np.median(np.log10(gq / gobs))),
                'median_log10_gnfw_over_gobs_dex': float(
                    np.median(np.log10(np.asarray([r['gnfw_m_s2'] for r in records if r['nmesh'] == nmesh]) / gobs))),
            },
        }

    summary = {
        'status': 'INDEPENDENT_CPU_PERIODIC_NGC2403_RADIAL_HOLDOUT',
        'scope': 'Four contiguous 16-radius holdouts. NFW parameters are fitted on the other 48 radii at fixed stellar M/L=0.5. QUMOND curves remain fixed.',
        'radii_kpc': [float(radius.min()), float(radius.max())],
        'nfw_fit': {'h0_km_s_mpc': H0, 'v200_grid': [30.0, 250.0, 0.25], 'c200_grid': [2.0, 30.0, 0.05]},
        'resolutions': outputs,
        'limitations': [
            'NFW predictions use two parameters fitted on 48 radii, while QUMOND remains fixed at M/L=0.5. The comparison has unequal fitting freedom.',
            'Pointwise velocity errors provide diagonal training weights. Radial covariance, distance, inclination, gas-profile, beam, and non-circular-motion uncertainties are not propagated.',
            'The four folds share training data. Fold scores are correlated and do not form independent evidence or a model-selection statistic.',
            'The NFW halo is an analytic post-processing model, not a live CDM halo or a cosmological simulation.',
            'The 512^3 QUMOND curve is an independent CPU reconstruction, not a compiled lagRamses result. Periodic-boundary and continuum uncertainties remain.',
        ],
    }
    with (out / 'holdout.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    for source in (reference_path, baseline_path, contract_path):
        shutil.copy2(source, out / source.name)
    summary['input_sha256'] = {str(path): sha256(path)
                               for path in (reference_path, baseline_path, contract_path)}
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
