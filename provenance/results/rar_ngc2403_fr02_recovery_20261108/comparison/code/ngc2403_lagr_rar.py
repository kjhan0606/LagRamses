#!/usr/bin/env python3
"""Compare static lagRamses QUMOND forces with SPARC RAR points."""
import csv
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE = Path('/gpfs/kjhan/RAR_QUMOND')
SPARC = BASE / 'analysis/ngc2403_input_448195/ngc2403_observed.csv'
SCENARIOS = {
    'THINGS': Path(os.environ.get(
        'RAR_THINGS_ROOT', str(BASE / 'runs/ngc2403_matched_gate_things_448337'))),
    'HALOGAS': Path(os.environ.get(
        'RAR_HALOGAS_ROOT', str(BASE / 'runs/ngc2403_matched_gate_halogas_448338'))),
}
CM_PER_KPC = 3.0856775814913673e21
M_PER_KPC = CM_PER_KPC / 100.0
G_DAGGER = 1.2e-10
N_AZIMUTH = 64


def info_units(path):
    text = path.read_text()
    result = {}
    for key in ('unit_l', 'unit_t'):
        match = re.search(r'(?im)^\s*' + key + r'\s*=\s*([^\s]+)', text)
        if not match:
            raise ValueError(f'missing {key} in {path}')
        result[key] = float(match.group(1))
    return result


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_uniform(case_root):
    snapshot = case_root / 'output_00001/data_00001.h5'
    info = case_root / 'output_00001/info_00001.txt'
    units = info_units(info)
    with h5py.File(snapshot, 'r') as handle:
        header = handle['header'].attrs
        levelmin = int(np.asarray(header['levelmin']).item())
        levelmax = int(np.asarray(header['nlevelmax']).item())
        boxlen = float(np.asarray(header['boxlen']).item())
        if levelmin != levelmax:
            raise ValueError(f'nonuniform reference grid: {levelmin}:{levelmax}')
        level = levelmin
        nmesh = 2**level
        dx = boxlen / nmesh
        amr = handle[f'amr/level_{level}']
        gravity = handle[f'gravity/level_{level}']
        xg = np.column_stack([np.asarray(amr[f'xg_{axis}']) for axis in (1, 2, 3)])
        ngrid = len(xg)
        if 8 * ngrid != nmesh**3:
            raise ValueError(f'incomplete uniform mesh: {ngrid} grids at {nmesh}^3')
        if np.any(np.asarray(amr['son_flag']) != 0):
            raise ValueError('uniform-grid reference unexpectedly contains refined cells')
        child = np.arange(8)
        child_xyz = np.column_stack((child % 2, (child // 2) % 2, child // 4))
        positions = (xg[:, None, :] + (child_xyz[None, :, :] - 0.5) * dx) % boxlen
        ijk = np.rint(positions / dx - 0.5).astype(np.int64) % nmesh
        flat = np.ravel_multi_index(ijk.reshape(-1, 3).T, (nmesh,) * 3)
        if not np.all(np.bincount(flat, minlength=nmesh**3) == 1):
            raise ValueError('AMR coordinates do not map one-to-one to the uniform grid')
        forces = np.empty((3, nmesh, nmesh, nmesh), dtype=np.float64)
        for axis in range(3):
            values = np.asarray(gravity[f'f_{axis + 1}'], dtype=np.float64).reshape(ngrid, 8)
            forces[axis, ijk[..., 0], ijk[..., 1], ijk[..., 2]] = values
        del positions, ijk, flat, xg, values
    if not np.all(np.isfinite(forces)):
        raise ValueError(f'nonfinite forces in {snapshot}')
    snapshot_sha256 = sha256_file(snapshot)
    gate_summary = case_root.parent / 'analysis.json'
    if gate_summary.is_file():
        gate = json.loads(gate_summary.read_text())
        expected = gate['cases'][case_root.name]['snapshot_sha256']
        if snapshot_sha256 != expected:
            raise ValueError(f'snapshot hash differs from the force-gate record: {snapshot}')
    return forces, {'snapshot': str(snapshot), 'snapshot_sha256': snapshot_sha256,
                    'level': level, 'nmesh': nmesh, 'boxlen_code': boxlen, 'dx_code': dx,
                    'dx_kpc': dx * units['unit_l'] / CM_PER_KPC,
                    'unit_length_cm': units['unit_l'], 'unit_time_s': units['unit_t']}


def trilinear(field, points, boxlen):
    nmesh = field.shape[1]
    dx = boxlen / nmesh
    index = (np.mod(points, boxlen) / dx) - 0.5
    lower = np.floor(index).astype(np.int64)
    weight = index - lower
    sampled = np.zeros((3, len(points)), dtype=np.float64)
    for ox in (0, 1):
        wx = (1.0 - weight[:, 0]) if ox == 0 else weight[:, 0]
        ix = (lower[:, 0] + ox) % nmesh
        for oy in (0, 1):
            wy = (1.0 - weight[:, 1]) if oy == 0 else weight[:, 1]
            iy = (lower[:, 1] + oy) % nmesh
            for oz in (0, 1):
                wz = (1.0 - weight[:, 2]) if oz == 0 else weight[:, 2]
                iz = (lower[:, 2] + oz) % nmesh
                sampled += field[:, ix, iy, iz] * (wx * wy * wz)[None, :]
    return sampled


def load_sparc():
    with SPARC.open(newline='') as stream:
        all_rows = list(csv.DictReader(stream))
    rows = [r for r in all_rows
            if r['passes_exploratory_point_cut'] == '1' and float(r['radius_kpc']) >= 1.0]
    if len(rows) != 64:
        raise ValueError(f'expected the documented 64-point NGC 2403 selection, got {len(rows)}')
    return rows


def empirical_rar(gbar):
    return gbar / (-np.expm1(-np.sqrt(gbar / G_DAGGER)))


def summarize(values):
    values = np.asarray(values, dtype=np.float64)
    return {'median': float(np.median(values)), 'p16': float(np.quantile(values, 0.16)),
            'p84': float(np.quantile(values, 0.84)), 'rms': float(np.sqrt(np.mean(values**2)))}


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: ngc2403_lagr_rar.py OUTPUT_DIRECTORY')
    out = Path(sys.argv[1])
    if out.exists():
        raise SystemExit(f'output already exists: {out}')
    observed = load_sparc()
    radius_kpc = np.array([float(r['radius_kpc']) for r in observed])
    vobs = np.array([float(r['vobs_km_s']) for r in observed])
    everr = np.array([float(r['vobs_error_km_s']) for r in observed])
    vgas = np.array([float(r['vgas_km_s']) for r in observed])
    vdisk = np.array([float(r['vdisk_km_s']) for r in observed])
    vbul = np.array([float(r['vbul_km_s']) for r in observed])
    gbar_sparc = np.array([float(r['gbar_m_s2']) for r in observed])
    gobs = np.array([float(r['gobs_m_s2']) for r in observed])
    calc_gobs = (vobs * 1000.0)**2 / (radius_kpc * M_PER_KPC)
    if not np.allclose(calc_gobs, gobs, rtol=2e-12, atol=0.0):
        raise ValueError('SPARC gobs columns do not match Vobs^2/R')
    calc_gbar = (vgas * np.abs(vgas) + 0.5 * vdisk * np.abs(vdisk) +
                 0.7 * vbul * np.abs(vbul)) * 1.0e6 / (radius_kpc * M_PER_KPC)
    if not np.allclose(calc_gbar, gbar_sparc, rtol=2e-12, atol=0.0):
        raise ValueError('SPARC gbar columns do not match the signed mass-model convention')

    azimuth = np.arange(N_AZIMUTH, dtype=np.float64) * (2.0 * math.pi / N_AZIMUTH)
    empirical_sparc = empirical_rar(gbar_sparc)
    rows = []
    per_scenario = {}
    source_hashes = {str(SPARC): hashlib.sha256(SPARC.read_bytes()).hexdigest()}
    for scenario, root in SCENARIOS.items():
        fields = {}
        meta = {}
        contract_path = root / 'source_contract.json'
        contract = json.loads(contract_path.read_text())
        if not math.isclose(float(contract['unit_length_cm']),
                            float(info_units((root / 'qumond/output_00001/info_00001.txt'))['unit_l']),
                            rel_tol=1e-12):
            raise ValueError(f'{scenario}: source contract and output length units differ')
        for source in (contract_path, root / 'star_profile_code.dat', root / 'hi_profile_code.dat'):
            source_hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        gate_summary_path = root / 'analysis.json'
        if gate_summary_path.is_file():
            gate_sources = json.loads(gate_summary_path.read_text())['source_files_sha256']
            for source in (contract_path, root / 'star_profile_code.dat', root / 'hi_profile_code.dat'):
                expected = gate_sources[source.name]
                if source_hashes[str(source)] != expected:
                    raise ValueError(f'{scenario}: input source hash differs from the force-gate record: {source}')
        for model in ('newtonian', 'qumond'):
            case_root = root / model
            field, current = read_uniform(case_root)
            fields[model] = field
            meta[model] = current
        if meta['newtonian']['nmesh'] != meta['qumond']['nmesh'] or meta['newtonian']['level'] != meta['qumond']['level']:
            raise ValueError(f'{scenario}: Newtonian and QUMOND grids differ')
        for key in ('unit_length_cm', 'unit_time_s', 'boxlen_code'):
            if not math.isclose(meta['newtonian'][key], meta['qumond'][key], rel_tol=1e-12):
                raise ValueError(f'{scenario}: mismatched {key}')
        model_rows = []
        for i, r in enumerate(observed):
            rad_code = radius_kpc[i] * CM_PER_KPC / meta['qumond']['unit_length_cm']
            if rad_code >= 0.5 * meta['qumond']['boxlen_code']:
                raise ValueError('observed radius reaches periodic half-box')
            points = np.column_stack((
                0.5 + rad_code * np.cos(azimuth),
                0.5 + rad_code * np.sin(azimuth),
                np.full(N_AZIMUTH, 0.5),
            ))
            unit_radial = np.column_stack((np.cos(azimuth), np.sin(azimuth)))
            # unit_l is in cm and unit_t in seconds, so the code
            # acceleration is first in cm/s^2; convert to SI by dividing by 100.
            acceleration = meta['qumond']['unit_length_cm'] / meta['qumond']['unit_time_s']**2 / 100.0
            sampled = {}
            for model in ('newtonian', 'qumond'):
                force = trilinear(fields[model], points, meta[model]['boxlen_code'])
                inward_code = -np.sum(force[:2].T * unit_radial, axis=1)
                inward_si = inward_code * acceleration
                if np.any(~np.isfinite(inward_si)) or np.median(inward_si) <= 0:
                    raise ValueError(f'{scenario} {model}: invalid radial acceleration at R={radius_kpc[i]}')
                sampled[model] = inward_si
            gbar_model = float(np.median(sampled['newtonian']))
            gq = float(np.median(sampled['qumond']))
            gn = gbar_model
            vq = math.sqrt(gq * radius_kpc[i] * M_PER_KPC) / 1000.0
            gemp_model = float(empirical_rar(np.array([gbar_model]))[0])
            record = {
                'scenario': scenario, 'radius_kpc': radius_kpc[i], 'vobs_km_s': vobs[i],
                'vobs_error_km_s': everr[i], 'gbar_sparc_m_s2': gbar_sparc[i],
                'gobs_m_s2': gobs[i], 'gbar_lagramses_newton_m_s2': gbar_model,
                'gnewton_baryons_m_s2': gn, 'gqumond_lagramses_m_s2': gq,
                'vqumond_km_s': vq,
                'gbar_model_over_sparc': gbar_model / gbar_sparc[i],
                'log10_gbar_model_over_sparc': math.log10(gbar_model / gbar_sparc[i]),
                'log10_gobs_over_empirical_sparc_rar': math.log10(gobs[i] / empirical_sparc[i]),
                'log10_gqumond_over_observed': math.log10(gq / gobs[i]),
                'log10_gqumond_over_empirical_model_rar': math.log10(gq / gemp_model),
                'gqumond_azimuth_p16_m_s2': float(np.quantile(sampled['qumond'], 0.16)),
                'gqumond_azimuth_p84_m_s2': float(np.quantile(sampled['qumond'], 0.84)),
            }
            rows.append(record)
            model_rows.append(record)
        per_scenario[scenario] = {
            'snapshot_and_units': meta,
            'median_log10_gqumond_over_observed_dex': summarize([r['log10_gqumond_over_observed'] for r in model_rows]),
            'median_log10_gqumond_over_empirical_model_rar_dex': summarize([r['log10_gqumond_over_empirical_model_rar'] for r in model_rows]),
            'median_log10_gbar_source_ratio_dex': summarize([r['log10_gbar_model_over_sparc'] for r in model_rows]),
            'median_abs_velocity_residual_km_s': float(np.median([abs(r['vqumond_km_s']-r['vobs_km_s']) for r in model_rows])),
        }
        for model in ('newtonian', 'qumond'):
            path = Path(meta[model]['snapshot'])
            source_hashes[str(path)] = meta[model]['snapshot_sha256']

    grid_summary = {}
    resolution_specs = set()
    for scenario, scenario_data in per_scenario.items():
        grid_summary[scenario] = {}
        for model, meta in scenario_data['snapshot_and_units'].items():
            dx_kpc = float(meta['dx_kpc'])
            box_kpc = float(meta['boxlen_code']) * float(meta['unit_length_cm']) / CM_PER_KPC
            grid_summary[scenario][model] = {
                'nmesh': int(meta['nmesh']), 'dx_kpc': dx_kpc, 'box_kpc': box_kpc,
            }
            resolution_specs.add((int(meta['nmesh']), round(dx_kpc, 6)))
    resolution_text = ', '.join(
        f'{nmesh}^3 (dx={dx_kpc:.3f} kpc)' for nmesh, dx_kpc in sorted(resolution_specs)
    )

    out.mkdir(parents=True)
    csv_path = out / 'ngc2403_lagr_rar.csv'
    with csv_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    # For each observed point, plot one QUMOND prediction per gas-map scenario.
    lookup = {name: [r for r in rows if r['scenario'] == name] for name in SCENARIOS}
    empirical_x = np.geomspace(min(gbar_sparc.min(), min(r['gbar_lagramses_newton_m_s2'] for r in rows))*0.8,
                               max(gbar_sparc.max(), max(r['gbar_lagramses_newton_m_s2'] for r in rows))*1.2, 256)
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10.4, 4.4), layout='constrained')
    for axis in (ax, bx):
        axis.tick_params(labelsize=9)
    ax.scatter(gbar_sparc, gobs, s=18, color='black', label='Observed SPARC', zorder=3)
    for name, color, marker in (('THINGS', '#d45b25', 's'), ('HALOGAS', '#247f7f', '^')):
        model = lookup[name]
        gmod = np.array([r['gbar_lagramses_newton_m_s2'] for r in model])
        gq = np.array([r['gqumond_lagramses_m_s2'] for r in model])
        ax.scatter(gmod, gq, s=18, marker=marker, color=color, label=f'lagRamses QUMOND · {name}')
    ax.plot(empirical_x, empirical_rar(empirical_x), color='0.45', ls='--', lw=1.1,
            label='Empirical RAR reference')
    ax.plot(empirical_x, empirical_x, color='0.65', ls=':', lw=1.0,
            label='Baryon-only Newtonian')
    ax.set(xscale='log', yscale='log', xlabel=r'$g_{\rm bar}$ (m s$^{-2}$)',
           ylabel=r'$g_R$ (m s$^{-2}$)')
    ax.text(0.03, 0.97, 'NGC 2403 · static source', transform=ax.transAxes,
            ha='left', va='top', fontsize=10)
    ax.legend(loc='lower right', frameon=False, fontsize=8)

    bx.axhline(0.0, color='0.65', lw=0.8)
    bx.plot(radius_kpc, [r['log10_gobs_over_empirical_sparc_rar'] for r in observed_rows(observed, gobs, empirical_sparc)],
            color='black', marker='o', ms=3, lw=0.8, label='Observed / empirical RAR')
    for name, color, marker in (('THINGS', '#d45b25', 's'), ('HALOGAS', '#247f7f', '^')):
        model = lookup[name]
        bx.plot(radius_kpc, [r['log10_gqumond_over_observed'] for r in model],
                color=color, marker=marker, ms=3, lw=0.8, label=f'{name} QUMOND / observed')
    bx.set(xlabel='Radius (kpc)', ylabel=r'$\log_{10}(g_1/g_2)$ (dex)')
    bx.text(0.03, 0.97, f'Fixed baryons · {resolution_text}', transform=bx.transAxes,
            ha='left', va='top', fontsize=10)
    bx.legend(loc='lower right', frameon=False, fontsize=8)
    fig.savefig(out / 'ngc2403_lagr_rar.pdf')
    fig.savefig(out / 'ngc2403_lagr_rar.png', dpi=180)
    plt.close(fig)

    summary = {
        'status': 'STATIC_LAGR_RAR_AUDIT',
        'points': len(observed), 'radius_kpc': [float(radius_kpc.min()), float(radius_kpc.max())],
        'selection': 'SPARC exploratory point cut and R >= 1 kpc; 64 points',
        'force_sample': f'trilinear midplane sample at {N_AZIMUTH} azimuths per observed radius; median inward radial force',
        'units': 'lagRamses acceleration converted with unit_l/unit_t^2 from each info file; gbar/gobs use SPARC signed mass model and observed velocity',
        'grid_resolution': grid_summary,
        'empirical_reference': 'gbar / (1 - exp(-sqrt(gbar / 1.2e-10 m s^-2))); descriptive reference only',
        'scenarios': per_scenario,
        'gas_map_sensitivity_median_log10_gq_halogas_over_things_dex': summarize([
            math.log10(lookup['HALOGAS'][i]['gqumond_lagramses_m_s2'] / lookup['THINGS'][i]['gqumond_lagramses_m_s2'])
            for i in range(len(observed))]),
        'inputs_sha256': source_hashes,
        'limitations': [
            f'Static, prescribed baryons in periodic uniform lagRamses runs at {resolution_text}; no galaxy evolution or cosmological initial conditions.',
            'The horizontal SPARC gbar and simulated Newtonian gbar come from different gas-source reconstructions; the ratio is retained as a source mismatch, not an error bar.',
            'Published velocity errors are not a full covariance; distance, inclination, stellar M/L and gas-map uncertainties are not propagated.',
            'The empirical RAR reference is not independent of the SPARC data and is not a theory fit.',
            'Baryon-only Newtonian is a solver control, not Lambda CDM; no halo is included in these force outputs.',
            'Periodic boundary conditions remain an explicit model systematic.',
        ],
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def observed_rows(rows, gobs, empirical):
    return [{'log10_gobs_over_empirical_sparc_rar': math.log10(g / e)} for g, e in zip(gobs, empirical)]


if __name__ == '__main__':
    main()
