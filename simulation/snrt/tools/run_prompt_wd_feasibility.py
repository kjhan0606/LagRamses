#!/usr/bin/env python3
"""Bounded OFFLINE COSMIC comparison, never a native yield/DTD exporter.

Requires isolated cosmic-popsynth==4.2.0 and its pinned example Params.ini.
Grid weights describe an explicitly assumed population, NOT an observed rate.
WD inventories include ineligible WDs and cannot certify N100 event capacity.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import cosmic
from cosmic import utils
from cosmic.evolve import Evolve
from cosmic.sample.initialbinarytable import InitialBinaryTable
import numpy as np
import pandas as pd
from helium_retention_reference import helium_regime, MODEL_ID

PARAMS_SHA256 = '68b2c0a7a90cd179e4f98c27428935472713e161ea0a19e5c8bf7eef8d3d04e9'
BCM_COLUMNS = ['tphys', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2',
               'RRLO_1', 'RRLO_2', 'deltam_1', 'deltam_2',
               'porb', 'sep', 'SN_1', 'SN_2', 'bin_state']  # bin_num is appended by COSMIC


def imf_integral(lo, hi, moment=0):
    """Unnormalised Kroupa-2001 two-slope IMF, .08--120 Msun."""
    value = 0.
    for a, b, slope, factor in ((.08, .5, 1.3, 1.), (.5, 120., 2.3, .5)):
        l, h = max(lo, a), min(hi, b)
        if h > l:
            q = 1 + moment - slope
            value += factor * (h**q-l**q)/q
    return value


def summarize(initial, bpp, bcm):
    births, disappearances = [], []
    for bid, rows in bpp.groupby('bin_num', sort=False):
        # Keep solver event order, including multiple transitions at one time.
        if np.any(np.diff(rows.tphys) < -1e-8):
            raise ValueError('nonmonotone event history')
        meta = initial.loc[int(bid)]
        for component in (1, 2):
            k, m = f'kstar_{component}', f'mass_{component}'
            wd = rows[rows[k] == 11]
            if not wd.empty:
                first = wd.iloc[0]
                births.append(dict(bin_num=int(bid), component=component,
                                   age_myr=float(first.tphys), mass_msun=float(first[m]),
                                   weight_per_initial_msun=float(meta.weight_per_initial_msun)))
            previous = rows[k].shift()
            for _, row in rows[(previous == 11) & (rows[k] == 15)].iterrows():
                disappearances.append(dict(bin_num=int(bid), component=component,
                                            age_myr=float(row.tphys), evol_type=int(row.evol_type)))
    snapshots = []
    weights = initial.weight_per_initial_msun
    for age in (40., 45., 50., 69., 100.):
        rows = bcm[np.isclose(bcm.tphys, age, rtol=0., atol=1e-7)]
        # Merged/disrupted systems can terminate their regular output early.
        # Do not extrapolate missing survivors into a fabricated population.
        if rows.bin_num.duplicated().any():
            raise ValueError('ambiguous snapshot rows')
        w = rows.bin_num.astype(int).map(weights).to_numpy()
        inventory = sum(np.sum(w*np.where(rows[f'kstar_{j}'] == 11,
                                         rows[f'mass_{j}'], 0.)) for j in (1, 2))
        snapshots.append(dict(age_myr=age, recorded_systems=len(rows),
                              recorded_CO_mass_per_initial_msun=float(inventory)))
    return dict(first_CO_age_myr=min((r['age_myr'] for r in births), default=None),
                CO_births_by_40=sum(r['age_myr'] <= 40 for r in births),
                CO_birth_records=births,
                CO_disappearances_not_classified_as_SNIa=disappearances,
                snapshots_not_complete_population=snapshots)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--params', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--followup-grid', type=Path,
                   help='Replay three Z=.01 representative histories and a wide-binary control')
    p.add_argument('--retention-screen', type=Path,
                   help='Classify existing fine histories; no counterfactual binary evolution')
    args = p.parse_args()
    if args.followup_grid is not None and args.retention_screen is not None:
        p.error('select only one follow-up mode')
    if cosmic.__version__ != '4.2.0':
        p.error('requires pinned COSMIC 4.2.0')
    if hashlib.sha256(args.params.read_bytes()).hexdigest() != PARAMS_SHA256:
        p.error('requires unchanged v4.2.0 example Params.ini')
    if args.output_dir.exists():
        p.error('refusing existing output directory')
    bse, sse, _, _, _, _ = utils.parse_inifile(str(args.params))
    # All physical flags are from the pinned file. Only the initial population
    # and output sampling below replace that file's illustrative sampler.
    args.output_dir.mkdir(parents=True, exist_ok=False)
    if args.retention_screen is not None:
        summary = dict(model_id=MODEL_ID, status='local_reference_screen_not_binary_reevolution',
                       parent_selection_sha256=hashlib.sha256((args.retention_screen/'selection.json').read_bytes()).hexdigest(),
                       kernel_sha256=hashlib.sha256(Path(__file__).with_name('helium_retention_reference.py').read_bytes()).hexdigest(),
                       assumptions=['W17 solar-composition hot nonrotating WD reference applied diagnostically.',
                                    'Actual parent Z=.01; no source metallicity recalibration.',
                                    'Both rates are proxies, NOT measured incident He transfer.',
                                    'No continuation beyond original WD destruction or donor collapse.',
                                    'No automatic expansion wind cap, flash efficiency extrapolation, or N100 events.'],
                       cases={})
        for label in ('early_CO','early_disappearance','donor_growth','wide_control'):
            path = args.retention_screen/f'{label}-bcm.csv.gz'
            bcm = pd.read_csv(path,float_precision='round_trip')
            records = []
            for j in (1,2):
                k = 3-j
                selected = bcm[(bcm[f'kstar_{j}']==11)&bcm[f'kstar_{k}'].between(7,9)
                               &(bcm[f'RRLO_{k}']>=1)]
                for _, row in selected.iterrows():
                    for proxy,rate in [('donor_total_loss',-row[f'deltam_{k}']),
                                       ('WD_net_gain',row[f'deltam_{j}'])]:
                        if rate <= 0:
                            continue
                        records.append(dict(age_myr=float(row.tphys),component=j,
                                            mass_msun=float(row[f'mass_{j}']),rate_proxy=proxy,
                                            rate_msun_yr=float(rate),**helium_regime(row[f'mass_{j}'],rate)))
            frame = pd.DataFrame(records)
            # JSON preserves explicit null for unsupported efficiencies.
            (args.output_dir/f'{label}.json').write_text(json.dumps(records,indent=2)+'\n')
            summary['cases'][label] = dict(input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                           classified_proxy_samples=len(records),
                                           regimes=frame.regime.value_counts().to_dict() if records else {})
        (args.output_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print('RETENTION_REFERENCE_SCREEN_COMPLETE_NO_N100_ADMISSION',flush=True)
        return
    if args.followup_grid is not None:
        grid = pd.read_csv(args.followup_grid/'z0.01-grid.csv', float_precision='round_trip').set_index('bin_num')
        original = pd.read_csv(args.followup_grid/'z0.01-initC.csv.gz', float_precision='round_trip').set_index('bin_num')
        selected = []
        for label, bid in [('early_CO',283), ('early_disappearance',289),
                           ('donor_growth',164), ('wide_control',283)]:
            row = grid.loc[bid]
            period = 1e9 if label == 'wide_control' else row.period
            selected.append(dict(label=label, original_bin_num=bid, m1=row.m1, m2=row.m2,
                                 period_days=period, randomseed=int(original.loc[bid,'randomseed'])))
        selection = dict(status='followup_comparison_not_source_admission',
                         cosmic_version=cosmic.__version__, params_sha256=PARAMS_SHA256,
                         BSEDict=bse, SSEDict=sse, metallicity=.01, bcm_spacing_myr=.005,
                         runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                         parent_selection_sha256=hashlib.sha256((args.followup_grid/'selection.json').read_bytes()).hexdigest(),
                         cases=selected)
        (args.output_dir/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
        for case in selected:
            init = InitialBinaryTable.InitialBinaries(m1=case['m1'],m2=case['m2'],
                                                      porb=case['period_days'],ecc=0.,tphysf=100.,
                                                      kstar1=1,kstar2=1,metallicity=.01)
            bpp,bcm,initc,kicks = Evolve.evolve(init,BSEDict=bse.copy(),SSEDict=sse.copy(),
                                              randomseed=case['randomseed'],nproc=1,dtp=.005,
                                              bcm_columns=BCM_COLUMNS)
            for name,table in [('bpp',bpp),('bcm',bcm),('initC',initc),('kicks',kicks)]:
                table.to_csv(args.output_dir/f'{case["label"]}-{name}.csv.gz',index=False)
            print('FOLLOWUP_DONE',case['label'],len(bpp),len(bcm),flush=True)
        return
    binary_fraction = .5  # fraction of SYSTEMS in this comparison
    ni, mi = imf_integral(.08, 120), imf_integral(.08, 120, 1)
    # Uniform secondary MASS on [.08, m1] => mean m2=(m1+.08)/2.
    mass_normalization = mi + binary_fraction*(mi+.08*ni)/2
    manifest = dict(status='offline_comparison_only_not_source_admission',
                    cosmic_version=cosmic.__version__, params_sha256=PARAMS_SHA256,
                    runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    BSEDict=bse, SSEDict=sse, random_seed=260907,
                    initial_population=dict(primary_imf='Kroupa2001', primary_support=[.08,120],
                                            binary_fraction_by_system=binary_fraction,
                                            secondary_mass='uniform [.08,m1]',
                                            log10_period_days='uniform [0,4.5]', eccentricity=0.,
                                            grid_primary_bin_edges=list(range(4,14)),
                                            donor_midpoint_fractions=[.1,.3,.5,.7,.9],
                                            period_midpoints=list(np.arange(.25,4.5,.5)),
                                            metallicities=[.007,.01,.014,.03]),
                    evolution_end_myr=100., bcm_spacing_myr=1.,
                    limitations=['Targeted midpoint quadrature; not a converged population rate.',
                                 'Primary masses outside 4--13 and singles are not evolved.',
                                 'No reweighting to observed DTD; denominator includes all initial systems.',
                                 'CO WD formation/inventory is not eligible N100 fuel or an explosion count.',
                                 'CO disappearance can be a merger or sub-Chandrasekhar event.',
                                 'No KL16 mass/chemistry splice or production-model changes.'])
    (args.output_dir/'selection.json').write_text(json.dumps(manifest, indent=2)+'\n')
    np.random.seed(260907)
    summary = {}
    for z in manifest['initial_population']['metallicities']:
        inputs = []
        for low, u, logp in itertools.product(range(4,13), [.1,.3,.5,.7,.9], np.arange(.25,4.5,.5)):
            m1 = low+.5
            inputs.append(dict(m1=m1, m2=.08+u*(m1-.08), period=10.**logp,
                               weight_per_initial_msun=binary_fraction*imf_integral(low,low+1)
                               /mass_normalization/5/9))
        grid = pd.DataFrame(inputs)
        grid.index.name = 'bin_num'
        init = InitialBinaryTable.InitialBinaries(m1=grid.m1.to_numpy(), m2=grid.m2.to_numpy(),
                                                  porb=grid.period.to_numpy(), ecc=np.zeros(len(grid)),
                                                  tphysf=np.full(len(grid),100.),
                                                  kstar1=np.ones(len(grid)), kstar2=np.ones(len(grid)),
                                                  metallicity=np.full(len(grid),z))
        print(f'RUN Z={z} systems={len(grid)}', flush=True)
        bpp, bcm, initc, kicks = Evolve.evolve(init, BSEDict=bse.copy(), SSEDict=sse.copy(),
                                             nproc=2, dtp=1., bcm_columns=BCM_COLUMNS)
        if bpp.bin_num.nunique() != len(grid) or not np.all(np.isfinite(bpp[['tphys','mass_1','mass_2']])):
            raise ValueError('missing/nonfinite evolution')
        prefix = args.output_dir/f'z{z:g}'
        grid.to_csv(str(prefix)+'-grid.csv')
        for name, table in [('bpp',bpp), ('bcm',bcm), ('initC',initc), ('kicks',kicks)]:
            table.to_csv(str(prefix)+f'-{name}.csv.gz', index=False)
        result = summarize(grid, bpp, bcm)
        summary[str(z)] = result
        print(f'DONE Z={z} first_CO={result["first_CO_age_myr"]} CO_by40={result["CO_births_by_40"]}', flush=True)
    (args.output_dir/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('OFFLINE_FEASIBILITY_COMPLETE_NOT_N100_QUALIFICATION', flush=True)


if __name__ == '__main__':
    main()
