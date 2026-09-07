#!/usr/bin/env python3
"""Offline real KL16+LC18 input, with explicit wind/release approximations.

SNIa mode binds the existing Kroupa effective-binary population; it does NOT
change its 40 Myr DTD or guarantee sufficient WDs. Native causality may reject
that combination. No raw source, historical approval or runtime default edits.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from build_lc18_native_wind import build_native_wind, SOLAR_MASS_CGS
from read_karakas_lugaro2016 import (read_karakas_lugaro2016, read_fishlock2014, DEFAULT_MATRIX,
                                    LIFETIME_SHA256, FISHLOCK_SHA256)


def build_combined(*, rotation: int, massive_wind_speed_km_s: float, agb_wind_speed_km_s: float,
                   agb_release: str, agb_energy: str, population: str, imf_id: int = 2,
                   low_z_agb: str = 'none'):
    if low_z_agb not in ('none', 'fishlock2014_raiteri96'):
        raise ValueError('unsupported low-Z source/timing choice')
    if population not in ('single', 'snia_baseline') or (population == 'snia_baseline' and imf_id != 1):
        raise ValueError('SNIa baseline requires explicit Kroupa imf_id=1; no implicit IMF conversion')
    if agb_release != 'terminal_envelope' or agb_energy != 'isotropic_thermalized':
        raise ValueError('unsupported AGB release/energy selection')
    if not math.isfinite(agb_wind_speed_km_s) or not 0 < agb_wind_speed_km_s < 29979.2458:
        raise ValueError('supply positive nonrelativistic AGB wind speed; no default')
    table, history = build_native_wind(rotation=rotation, wind_speed_km_s=massive_wind_speed_km_s,
                                      timing='uniform_until_terminal', composition='as_tabulated_mean',
                                      energy='isotropic_thermalized', imf_id=imf_id)
    report = read_karakas_lugaro2016(include_lifetimes=True)
    # Common ordinary CO-AGB support is 1--6 Msun. Do not turn ONe or hybrid
    # CO(Ne) cores into an N100 WD supplier or resurrect commented 8 Msun nodes.
    co = [r for r in report['records'] if r['evolution']['core_kind'] == 'CO']
    zs = sorted({r['coordinate']['metallicity_mass_fraction'] for r in co})
    upper = min(max(r['coordinate']['initial_mass_msun'] for r in co
                    if r['coordinate']['metallicity_mass_fraction'] == z) for z in zs)
    nodes = [r for r in co if r['coordinate']['initial_mass_msun'] <= upper]
    if upper != 6 or len(nodes) != 58:
        raise ValueError('unexpected ordinary CO-AGB source support')
    candidate = next(c for c in json.loads(DEFAULT_MATRIX.read_text())['candidates']
                     if c['candidate_id'] == 'karakas_lugaro2016_agb')
    selection = dict(version=1, population=population, imf_id=imf_id, rotation=rotation,
                     agb_release=agb_release, agb_energy=agb_energy, agb_wind_speed_km_s=agb_wind_speed_km_s,
                     massive_wind_speed_km_s=massive_wind_speed_km_s,
                     massive_table_sha256=hashlib.sha256(table.encode()).hexdigest(),
                     kl16_files=candidate['source_file_sha256'], lifetime_sha256=LIFETIME_SHA256,
                     normalization=report['normalization_policy'], agb_mass_support=[1,upper],
                     core_selection='ordinary_CO_only_common_support',
                     snia_policy='unchanged_approved_DTD_strict_causal_WD_budget',
                     agb_decay='source_fully_decayed_elements', net='unavailable_diagnostic_zero',
                     kl16_model_coordinates=[dict(**r['coordinate'], overshoot=r['overshoot_label']) for r in nodes])
    if low_z_agb != 'none':
        low = read_fishlock2014(lifetime_model='raiteri96_padova')
        nodes = low + nodes
        selection['low_z_agb'] = dict(model=low_z_agb,source_sha256=FISHLOCK_SHA256,
            mass_authority='sum_all_source_gross_elements_no_normalization',
            lifetime='Raiteri96_Padova_fit; not source-matched Monash lifetimes',
            lifetime_coefficients=[[10.13,.07547,-.008084],[-4.424,-.7939,-.1187],[1.262,.3385,.05417]],
            core_selection='Fishlock2014_CO_1_to_6; exclude_ONe_7',
            z=.001,inter_source_mixture='linear_Z_cumulative_mixture',
            common_z_range=[.001,.01345],snia='no_new_microscopic_WD_supplier_claim')
    identity = hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest()
    rows = []
    for node in nodes:
        c = node['coordinate']; e = node['selected_ejecta']
        age = node['evolution']['stellar_lifetime_yr']
        energy = .5*e['returned_mass_msun']*SOLAR_MASS_CGS*(agb_wind_speed_km_s*1e5)**2
        for time in (0., age, 2e10):
            f = float(time >= age)
            row = [2., c['initial_mass_msun'], c['metallicity_mass_fraction'], time,
                   f*e['returned_mass_msun'], f*e['remnant_mass_msun'], f*energy,
                   0.,0.,0., *(f*v for v in e['tracked_ejected_mass_msun']), *([0.]*11)]
            rows.append(' '.join(format(v,'.17g') for v in row))
    history = re.sub(r"model_id='[^']*'", f"model_id='kl16_lc18_user_v1:{identity}'", history)
    history = re.sub(r"model_coordinates='[^']*'", "model_coordinates='KL16 normalized; actual M/Z/Y/overshoot lifetimes; explicit isotropic winds; ordinary CO envelope/WD step'", history)
    if low_z_agb != 'none':
        history = re.sub(r"model_coordinates='[^']*'",
            "model_coordinates='KL16 source ages + Fishlock Z=.001/Raiteri96 Padova ages; gross-sum mass; CO<=6; linear Z source mixture'", history)
    history = history.replace("net_yield_policy='unavailable_diagnostic_zero',",
                              "net_yield_policy='unavailable_diagnostic_zero',agb_release_policy='terminal_step',")
    if population == 'snia_baseline':
        history = history.replace('input_population_id=0,input_binary_fraction=0,',
                                  'input_population_id=1,input_binary_fraction=.5,')
    header = '# KL16+LC18 real yields with explicit release/energy/population approximations; not all-channel approval.\n'
    header += '# '+json.dumps(selection,sort_keys=True)+'\n'
    table = header+table+'\n'.join(rows)+'\n'
    return table, history


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--rotation', type=int, choices=(0,150,300), required=True)
    p.add_argument('--massive-wind-speed-km-s', type=float, required=True)
    p.add_argument('--agb-wind-speed-km-s', type=float, required=True)
    p.add_argument('--agb-release', choices=('terminal_envelope',), required=True)
    p.add_argument('--agb-energy', choices=('isotropic_thermalized',), required=True)
    p.add_argument('--population', choices=('single','snia_baseline'), required=True)
    p.add_argument('--imf-id', type=int, choices=(0,1,2,4), default=2)
    p.add_argument('--low-z-agb', choices=('none','fishlock2014_raiteri96'), default='none',
                   help='Optional explicit source/lifetime mixture; never extrapolates yields')
    args = p.parse_args()
    if args.output_dir.exists():
        p.error('existing output directory refused')
    table, history = build_combined(rotation=args.rotation, massive_wind_speed_km_s=args.massive_wind_speed_km_s,
                                   agb_wind_speed_km_s=args.agb_wind_speed_km_s, agb_release=args.agb_release,
                                   agb_energy=args.agb_energy, population=args.population, imf_id=args.imf_id,
                                   low_z_agb=args.low_z_agb)
    args.output_dir.mkdir(parents=True,exist_ok=False)
    for name, content in [('yields.dat',table),('history.nml',history)]:
        with (args.output_dir/name).open('x') as f:
            f.write(content)
        print(name,hashlib.sha256(content.encode()).hexdigest())


if __name__ == '__main__':
    main()
