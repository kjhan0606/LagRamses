#!/usr/bin/env python3
"""Compare two compiled 512-cubed gas-source prescriptions at matched radii."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from ngc2403_lagr_rar import M_PER_KPC, empirical_rar, read_uniform
from ngc2403_compiled512_cdm_nfw_holdout import sample_radial_acceleration
from ngc2403_cpu512_radial_holdout import (
    C200_GRID, V200_GRID, H0, fit_nfw, nfw_curve, velocity_summary,
)

# Exact checksum of the immutable 1279415 CSV, not the normalized tracked CSV.
BASELINE_SHA = '702aa0a4c46f78d61f8c0971983163073bc99d3fc4106d68622965f3307544db'
SOURCE_SHA = {
    'star_profile_code.dat': '1cf451b8f9ff983ec083d332f7b4a32ea40ff803752e48e5d756328b297d89d3',
    'hi_profile_code.dat': 'e041eb1e22db9c1b39e58a7d08259f5f9bb4c5dc138ff314c590c03010a1dd2d',
    'source_contract.json': '007e745db540161b701e924f3d66475721953a11010cfcd7c658e5a1e7eeb591',
}
OBS_SHA = '852eb3bcbf10e0e088f70598a25075b991867eb64a0c7ca0cd46d6ff7eb6fb95'
BINARY_SHA = '7baf42cc884f33e9babb09b21aca5c46c380b3c32186cd562d20971c6414061c'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def summary(values):
    a = np.asarray(values)
    if not np.all(np.isfinite(a)):
        raise ValueError('nonfinite summary')
    return dict(mean=float(a.mean()), median=float(np.median(a)),
                rms=float(np.sqrt(np.mean(a*a))), min=float(a.min()), max=float(a.max()))


def plot(out, r, vobs, err, oldq, q, vnfw, gobs):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.5), layout='constrained')
    ax.errorbar(r, vobs, err, fmt='.', color='black', ms=3, lw=.6, label='SPARC observed')
    ax.plot(r, oldq, color='#247f7f', lw=1.2, label='QUMOND: THINGS-like gas')
    ax.plot(r, q, color='#7951a0', lw=1.2, label='QUMOND: FR02 gas')
    ax.plot(r, vnfw, color='#d45b25', ls='--', lw=1, label='NFW: FR02 gas, held out')
    ax.set(xlabel='Radius (kpc)', ylabel=r'$V_c$ (km s$^{-1}$)', ylim=(40, 165))
    ax.text(.03, .97, r'NGC 2403 · matched $512^3$ fields', transform=ax.transAxes,
            va='top', fontsize=10)
    ax.legend(loc='lower right', frameon=False, fontsize=10)
    for v, c, label, ls in ((oldq, '#247f7f', 'QUMOND: THINGS-like gas', '-'),
                            (q, '#7951a0', 'QUMOND: FR02 gas', '-'),
                            (vnfw, '#d45b25', 'NFW: FR02 gas, held out', '--')):
        bx.plot(r, np.log10(v*v*1e6/(r*M_PER_KPC)/gobs), color=c, label=label, ls=ls, lw=1.2)
    bx.axhline(0, color='.65', lw=.8)
    bx.set(xlabel='Radius (kpc)', ylabel=r'$\log_{10}(g_{\rm model}/g_{\rm obs})$', ylim=(-.25, .17))
    bx.text(.03, .97, 'Residuals at the same radii', transform=bx.transAxes, va='top', fontsize=10)
    bx.legend(loc='lower left', frameon=False, fontsize=10)
    for a in (ax, bx):
        a.set_xlim(.5, 21.4)
        a.tick_params(labelsize=10)
    fig.savefig(out/'source_comparison.pdf')
    fig.savefig(out/'source_comparison.png', dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    for name in ('run', 'baseline', 'baseline_summary', 'observed', 'out'):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    run, out = args.run, args.out
    if not out.is_dir() or any(out.iterdir()):
        raise ValueError('output directory must exist and be empty')
    if not (run/'COMPLETE').is_file():
        raise ValueError('compiled force run incomplete')
    gate = json.loads((run/'analysis.json').read_text())
    contract = json.loads((run/'source_contract.json').read_text())
    policy = json.loads((run/'output_policy.json').read_text())
    old_summary = json.loads(args.baseline_summary.read_text())
    if gate['status'] != 'PERIODIC_SOURCE_AND_FORCE_GATE_PASS' or gate['grid']['nmesh'] != 512:
        raise ValueError('expected passed compiled 512 field check')
    if gate['grid']['box_kpc'] != 64.0 or gate['grid']['cell_width_kpc'] != .125:
        raise ValueError('unexpected geometry')
    if policy['binary_sha256'] != BINARY_SHA or old_summary['output_policy']['binary_sha256'] != BINARY_SHA:
        raise ValueError('source comparison requires identical compiled binaries')
    if old_summary['grid'] != gate['grid']:
        raise ValueError('source comparison requires matched physical grids')
    old_contract = old_summary['source_contract']
    for key in ('unit_length_cm', 'star_exp_scaleheight_code', 'gas_exp_scaleheight_code'):
        if not math.isclose(contract[key], old_contract[key], rel_tol=1e-13):
            raise ValueError(f'nonradial-source contract changed: {key}')
    for name, expected in SOURCE_SHA.items():
        if sha(run/name) != expected or gate['source_files_sha256'][name] != expected:
            raise ValueError(f'unexpected FR02 source: {name}')
    if sha(args.baseline) != BASELINE_SHA or sha(args.observed) != OBS_SHA:
        raise ValueError('baseline/observed input changed')
    for label, metrics in gate['cases'].items():
        keys = ['median_relative_force_error', 'p90_relative_force_error']
        if label == 'newtonian': keys.append('source_rms_relative_error')
        if any(not 0 <= metrics[k] <= 1e-6 for k in keys):
            raise ValueError('force gate failed')
        if Path(metrics['snapshot']).resolve() != (run/label/'output_00001/data_00001.h5').resolve():
            raise ValueError('force gate refers to a different snapshot')
    old = load_csv(args.baseline)
    obs = [x for x in load_csv(args.observed) if x['passes_exploratory_point_cut'] == '1' and float(x['radius_kpc']) >= 1]
    old.sort(key=lambda x: float(x['radius_kpc']))
    obs.sort(key=lambda x: float(x['radius_kpc']))
    if len(old) != 64 or len(obs) != 64:
        raise ValueError('expected unchanged 64-radius sample')
    def arr(rows, key): return np.array([float(x[key]) for x in rows])
    r, vobs, err, gobs, gs = (arr(obs,k) for k in ('radius_kpc','vobs_km_s','vobs_error_km_s','gobs_m_s2','gbar_m_s2'))
    for k, a in [('radius_kpc',r),('vobs_km_s',vobs),('vobs_error_km_s',err),('gobs_m_s2',gobs),('gbar_sparc_m_s2',gs)]:
        if not np.array_equal(arr(old,k), a): raise ValueError(f'baseline observation mismatch: {k}')
    if np.any(np.diff(r) <= 0) or np.any(err <= 0): raise ValueError('invalid radius/error')
    if not np.allclose(vobs**2*1e6/(r*M_PER_KPC),gobs,rtol=2e-12,atol=0): raise ValueError('gobs convention mismatch')
    calc_gs = sum(f*arr(obs,k)*np.abs(arr(obs,k)) for f,k in ((1.,'vgas_km_s'),(.5,'vdisk_km_s'),(.7,'vbul_km_s')))*1e6/(r*M_PER_KPC)
    if not np.allclose(calc_gs,gs,rtol=2e-12,atol=0): raise ValueError('signed gbar convention mismatch')
    sample, metadata = {}, {}
    for label in ('newtonian','qumond'):
        force, meta = read_uniform(run/label)
        if meta['nmesh'] != 512 or not math.isclose(meta['dx_kpc'], .125, abs_tol=1e-12):
            raise ValueError('snapshot resolution mismatch')
        if metadata:
            for k in ('unit_length_cm','unit_time_s','boxlen_code'):
                if meta[k] != metadata['newtonian'][k]: raise ValueError(f'paired {k} mismatch')
        metadata[label] = meta
        sample[label] = [sample_radial_acceleration(force, rad, meta) for rad in r]
        del force
    gb, gq = (np.array([x['median'] for x in sample[k]]) for k in ('newtonian','qumond'))
    vq = np.sqrt(gq*r*M_PER_KPC)/1000
    vnfw = np.empty(64)
    folds = []
    for fold in range(4):
        held = np.arange(64)//16 == fold
        v200, c200, chi = fit_nfw(r[~held], (gb*r*M_PER_KPC/1e6)[~held], vobs[~held], err[~held])
        vnfw[held] = nfw_curve(r[held], (gb*r*M_PER_KPC/1e6)[held], v200, c200)
        folds.append(dict(fold=fold, training_points=48, heldout_points=16,
                          radius_kpc=[float(r[held][0]),float(r[held][-1])],v200_km_s=v200,c200=c200,
                          training_diagonal_chi2=chi, grid_edge=bool(v200 in (V200_GRID[0],V200_GRID[-1]) or c200 in (C200_GRID[0],C200_GRID[-1]))))
    gnfw = vnfw**2*1e6/(r*M_PER_KPC)
    oldq = arr(old,'gqumond_compiled_m_s2')
    oldb = arr(old,'gbar_compiled_newtonian_m_s2')
    measures = {
        'log10_qumond_fr02_over_things_dex':np.log10(gq/oldq),
        'log10_gbar_fr02_over_things_dex':np.log10(gb/oldb),
        'qumond_fr02_minus_observed_dex':np.log10(gq/gobs),
        'qumond_things_minus_observed_dex':np.log10(oldq/gobs),
        'nfw_fr02_minus_observed_dex':np.log10(gnfw/gobs),
        'gbar_fr02_minus_sparc_dex':np.log10(gb/gs),
        'gbar_things_minus_sparc_dex':np.log10(oldb/gs),
        'qumond_fr02_minus_reference_dex':np.log10(gq/empirical_rar(gb)),
        'reference_baryon_input_shift_dex':np.log10(empirical_rar(gb)/empirical_rar(gs)),
        'reference_sparc_minus_observed_dex':np.log10(empirical_rar(gs)/gobs),
    }
    identity = measures['qumond_fr02_minus_reference_dex']+measures['reference_baryon_input_shift_dex']+measures['reference_sparc_minus_observed_dex']-measures['qumond_fr02_minus_observed_dex']
    if np.max(np.abs(identity)) > 1e-12: raise ValueError('residual decomposition fails')
    rows = []
    for i, radius in enumerate(r):
        row = dict(radius_kpc=float(radius),fold=i//16,vobs_km_s=vobs[i],vobs_error_km_s=err[i],gobs_m_s2=gobs[i],gbar_sparc_m_s2=gs[i],
                   gbar_fr02_m_s2=gb[i],gqumond_fr02_m_s2=gq[i],gnfw_fr02_heldout_m_s2=gnfw[i],
                   gbar_things_m_s2=oldb[i],gqumond_things_m_s2=oldq[i],vqumond_fr02_km_s=vq[i],vnfw_fr02_km_s=vnfw[i])
        row.update({k:float(v[i]) for k,v in measures.items()})
        for label in sample:
            for quantile in ('p16','p84'): row[f'{label}_azimuth_{quantile}_m_s2']=sample[label][i][quantile]
        rows.append(row)
    with (out/'comparison.csv').open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    result = dict(status='MATCHED_COMPILED_512_GAS_SOURCE_COMPARISON',run_root=str(run),
                  grid=gate['grid'],source_contract=contract,baseline_contract=old_contract,
                  sample=dict(points=64, radius_kpc=[float(r[0]),float(r[-1])],azimuths_per_radius=64),
                  all_radii={k:summary(v) for k,v in measures.items()},
                  outer_radii=dict(radius_kpc=[float(r[48]),float(r[-1])], statistics={k:summary(v[48:]) for k,v in measures.items()}),
                  velocity=dict(things_qumond=velocity_summary(arr(old,'vqumond_compiled_km_s'),vobs,err),
                                fr02_qumond=velocity_summary(vq,vobs,err),fr02_nfw=velocity_summary(vnfw,vobs,err)),
                  nfw_folds=folds, nfw_fit_grid=dict(H0=H0,V200=[float(V200_GRID[0]),float(V200_GRID[-1])],c200=[float(C200_GRID[0]),float(C200_GRID[-1])]),
                  identity_max_abs_dex=float(np.max(np.abs(identity))), force_metadata=metadata,
                  inputs_sha256={str(path):sha(path) for path in (args.baseline,args.baseline_summary,args.observed,run/'analysis.json',run/'source_contract.json')},
                  limits=[
                      'Same binary, grid, stars, heights, a0 and sampling. The radial gas source prescription changes, including profile and taper; this is not a pure density-normalization effect.',
                      'The FR02 graph profile is a candidate and not established as the exact SPARC gas input. Its uncertainty is not quantified.',
                      'Both sources use a periodic 64-kpc box and finite 0.125-kpc cells. Operator verification does not prove continuum or isolated-boundary accuracy.',
                      'QUMOND is a fixed full field solve; NFW fits two parameters on 48 radii per fold. Folds share training data, and diagonal errors omit observational covariance and nuisance parameters.',
                      'NFW is an isolated analytic contribution added to the periodic baryon force, not a live Lambda-CDM or cosmological simulation.',
                      'The reference decomposition is algebraic accounting, not an error budget or causal separation of physics. One galaxy does not measure intrinsic population RAR scatter.',
                  ])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(out,r,vobs,err,arr(old,'vqumond_compiled_km_s'),vq,vnfw,gobs)
    print(json.dumps({k:result[k] for k in ('status','velocity','all_radii','outer_radii','nfw_folds')},indent=2))


if __name__ == '__main__': main()
