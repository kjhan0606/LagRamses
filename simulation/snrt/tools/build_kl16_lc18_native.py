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
from build_lc18_native_wind import build_native_wind, SOLAR_MASS_CGS, DECAY_MODELS
from read_karakas_lugaro2016 import (read_karakas_lugaro2016, read_fishlock2014, DEFAULT_MATRIX,
                                    LIFETIME_SHA256, FISHLOCK_SHA256, attach_fishlock_pulses)


def build_combined(*, rotation: int, massive_wind_speed_km_s: float, agb_wind_speed_km_s: float,
                   agb_release: str, agb_energy: str, population: str, imf_id: int = 2,
                   low_z_agb: str = 'none', massive_source: str = 'wind_only',
                   snii_energy_erg: float | None = None, massive_wind_timing: str = 'uniform_until_terminal',
                   agb_source_scope: str = 'co_common', agb_net: str = 'unavailable',
                   massive_decay: str = 'as_tabulated_no_decay'):
    if agb_net not in ('unavailable','normalized_initial_MZY'):
        raise ValueError('unsupported AGB net-yield model')
    if agb_source_scope not in ('co_common','kl16_envelopes_to7'):
        raise ValueError('unsupported AGB source scope')
    if low_z_agb not in ('none', 'fishlock2014_raiteri96'):
        raise ValueError('unsupported low-Z source/timing choice')
    if population not in ('single', 'snia_baseline') or (population == 'snia_baseline' and imf_id != 1):
        raise ValueError('SNIa baseline requires explicit Kroupa imf_id=1; no implicit IMF conversion')
    if agb_release not in ('terminal_envelope','fishlock_tp_mass_loss') or agb_energy != 'isotropic_thermalized':
        raise ValueError('unsupported AGB release/energy selection')
    if agb_release=='fishlock_tp_mass_loss' and low_z_agb=='none':
        raise ValueError('Fishlock thermal pulses require the explicit low-Z source')
    if not math.isfinite(agb_wind_speed_km_s) or not 0 < agb_wind_speed_km_s < 29979.2458:
        raise ValueError('supply positive nonrelativistic AGB wind speed; no default')
    table, history = build_native_wind(rotation=rotation, wind_speed_km_s=massive_wind_speed_km_s,
                                      timing=massive_wind_timing, composition='as_tabulated_mean',
                                      energy='isotropic_thermalized', imf_id=imf_id,
                                      massive_source=massive_source, snii_energy_erg=snii_energy_erg,
                                      decay=massive_decay)
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
    non_co = []
    if agb_source_scope == 'kl16_envelopes_to7':
        upper = 7
        nodes = [r for r in report['records'] if r['coordinate']['initial_mass_msun'] <= upper]
        if len(nodes) != 61:
            raise ValueError('unexpected KL16 1--7 envelope grid')
        non_co = [r for r in nodes if r['evolution']['core_kind'] != 'CO']
        if [(r['coordinate']['initial_mass_msun'],r['coordinate']['metallicity_mass_fraction'],
             r['evolution']['core_kind']) for r in non_co] != [(7.,.007,'CO_Ne')]:
            raise ValueError('non-CO source map changed')
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
    if massive_decay != DECAY_MODELS[0]:
        selection['massive_decay'] = massive_decay
    if agb_source_scope != 'co_common':
        selection.update(agb_source_scope=agb_source_scope,
            core_selection='all_source_envelopes_1_to7; hybrid_CO_Ne_excluded_from_strict_Ia_CO_inventory',
            non_co_nodes=[dict(**r['coordinate'],core_kind=r['evolution']['core_kind']) for r in non_co])
    if low_z_agb != 'none':
        low = read_fishlock2014(lifetime_model='raiteri96_padova',include_one=upper==7,
                               net_yields=agb_net!='unavailable')
        nodes = low + nodes
        selection['low_z_agb'] = dict(model=low_z_agb,source_sha256=FISHLOCK_SHA256,
            mass_authority='sum_all_source_gross_elements_no_normalization',
            lifetime='Raiteri96_Padova_fit; not source-matched Monash lifetimes',
            lifetime_coefficients=[[10.13,.07547,-.008084],[-4.424,-.7939,-.1187],[1.262,.3385,.05417]],
            core_selection='Fishlock2014_CO_1_to_6; exclude_ONe_7',
            z=.001,inter_source_mixture='linear_Z_cumulative_mixture',
            common_z_range=[.001,.01345],snia='no_new_microscopic_WD_supplier_claim')
        if upper==7:
            non_co = [r for r in nodes if r['evolution']['core_kind']!='CO']
            selection['core_selection']='all_selected_envelopes_1_to7; CO_Ne_and_ONe_excluded_from_strict_Ia_CO_inventory'
            selection['low_z_agb']['core_selection']='Fishlock2014_1_to7_envelopes; ONe_7_excluded_from_strict_CO_inventory'
            selection['non_co_nodes']=[dict(**r['coordinate'],core_kind=r['evolution']['core_kind']) for r in non_co]
    if agb_net != 'unavailable':
        initial_models = []
        for n in nodes:
            if 'initial_net_model' in n:
                initial_models.append(n['initial_net_model'])
                continue
            # Initial abundances precede overshoot and the later AGB mixing
            # mass prescription. Match actual initial M,Z,Y; permit only
            # duplicate entries with exactly identical abundance vectors.
            keys = ('initial_mass_msun','metallicity_mass_fraction','initial_helium_label')
            matches = [b for b in report['initial_composition_records']
                       if all(b['coordinate'][k]==n['coordinate'][k] for k in keys)]
            vectors = [{k:v['mass_fraction'] for k,v in b['elements_by_atomic_number'].items()} for b in matches]
            if not vectors or any(v!=vectors[0] for v in vectors):
                raise ValueError('missing or ambiguous initial M/Z/Y abundance vector')
            x = vectors[0]; total = math.fsum(x.values())
            if not math.isfinite(total) or total<=0:
                raise ValueError('invalid initial composition sum')
            x = {k:v/total for k,v in x.items()}
            n['selected_net'] = [n['selected_ejecta']['gross_mass_msun_by_atomic_number'][k]-
                x[k]*n['selected_ejecta']['returned_mass_msun'] for k in (1,2,6,7,8,10,12,14,16,20,26)]
            initial_models.append(dict(coordinate=n['coordinate'],source_lines=[b['source_line'] for b in matches],
                                       initial_sum_before_normalization=total,initial_fractions=x))
        selection.update(net='AGB_only_gross_minus_normalized_source_X0_times_returned',
                         initial_composition_matching='initial_M_Z_Y; identical duplicates only',
                         initial_composition_normalization='sum_all_source_element_mass_fractions',
                         initial_models=initial_models)
        if low_z_agb!='none':
            selection['initial_composition_matching']='KL16_initial_M_Z_Y; Fishlock_same_model_X0_column'
    if agb_release=='fishlock_tp_mass_loss':
        selection['agb_timing']=attach_fishlock_pulses(nodes)
    identity = hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest()
    rows = []
    for node in nodes:
        c = node['coordinate']; e = node['selected_ejecta']
        age = node['evolution']['stellar_lifetime_yr']
        energy = .5*e['returned_mass_msun']*SOLAR_MASS_CGS*(agb_wind_speed_km_s*1e5)**2
        for time,f in node.get('wind_knots',[(0.,0.),(age,1.),(2e10,1.)]):
            row = [2., c['initial_mass_msun'], c['metallicity_mass_fraction'], time,
                   f*e['returned_mass_msun'], float(time>=age)*e['remnant_mass_msun'], f*energy,
                   0.,0.,0., *(f*v for v in e['tracked_ejected_mass_msun']),
                   *(f*v for v in node.get('selected_net',[0.]*11))]
            rows.append(' '.join(format(v,'.17g') for v in row))
    history = re.sub(r"model_id='[^']*'", f"model_id='kl16_lc18_user_v1:{identity}'", history)
    history = re.sub(r"model_coordinates='[^']*'", "model_coordinates='KL16 normalized; actual M/Z/Y/overshoot lifetimes; explicit isotropic winds; ordinary CO envelope/WD step'", history)
    if low_z_agb != 'none':
        history = re.sub(r"model_coordinates='[^']*'",
            "model_coordinates='KL16 source ages + Fishlock Z=.001/Raiteri96 Padova ages; gross-sum mass; CO<=6; linear Z source mixture'", history)
    history = history.replace("net_yield_policy='unavailable_diagnostic_zero',",
                              "net_yield_policy='unavailable_diagnostic_zero',agb_release_policy='terminal_step',")
    if non_co:
        history = history.replace("agb_release_policy='terminal_step',",
            "agb_release_policy='terminal_step',agb_non_co_count=1,\n"
            " agb_non_co_mass=7,agb_non_co_z=.007,agb_non_co_kind=1,")
        history = re.sub(r"model_coordinates='[^']*'",
            "model_coordinates='KL16 1--7 source ages/envelopes; CO(Ne) at 7/Z=.007 excluded from strict Ia inventory; LC18 source wind/SN'",history)
        if low_z_agb!='none':
            history = history.replace('agb_non_co_count=1,','agb_non_co_count=2,').replace(
                'agb_non_co_mass=7,agb_non_co_z=.007,agb_non_co_kind=1,',
                'agb_non_co_mass=7,7,agb_non_co_z=.001,.007,agb_non_co_kind=2,1,')
            history = re.sub(r"model_coordinates='[^']*'",
                "model_coordinates='KL16+Fishlock 1--7 envelopes; Fishlock Raiteri96 ages; ONe/.001 and CO(Ne)/.007 excluded from strict Ia CO supply'",history)
    if agb_net != 'unavailable':
        history = history.replace("net_yield_policy='unavailable_diagnostic_zero',",
            "net_yield_policy='channel_mask',net_channel_available=.false.,.true.,.false.,.false.,.false.,")
    if population == 'snia_baseline':
        history = history.replace('input_population_id=0,input_binary_fraction=0,',
                                  'input_population_id=1,input_binary_fraction=.5,')
    if agb_release=='fishlock_tp_mass_loss':
        selected=[n for n in nodes if 'terminal_jump_fraction' in n]
        masses=','.join(format(n['coordinate']['initial_mass_msun'],'.17g') for n in selected)
        fractions=','.join(format(n['terminal_jump_fraction'],'.17g') for n in selected)
        history=history.replace("agb_release_policy='terminal_step',",
            f"agb_release_policy='wind_history_terminal_remnant',agb_wind_jump_count={len(selected)},\n"
            f" agb_wind_jump_mass={masses},agb_wind_jump_z={len(selected)}*.001,\n"
            f" agb_wind_jump_fraction={fractions},")
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
    p.add_argument('--agb-release', choices=('terminal_envelope','fishlock_tp_mass_loss'), required=True)
    p.add_argument('--agb-energy', choices=('isotropic_thermalized',), required=True)
    p.add_argument('--population', choices=('single','snia_baseline'), required=True)
    p.add_argument('--imf-id', type=int, choices=(0,1,2,4), default=2)
    p.add_argument('--low-z-agb', choices=('none','fishlock2014_raiteri96'), default='none',
                   help='Optional explicit source/lifetime mixture; never extrapolates yields')
    p.add_argument('--massive-source', choices=('wind_only','lc18_set_r'), default='wind_only')
    p.add_argument('--snii-energy-erg', type=float, help='Required explicit comparison parameter for lc18_set_r')
    p.add_argument('--massive-wind-timing', choices=('uniform_until_terminal','phase_mass_loss_or_uniform'),
                   default='uniform_until_terminal')
    p.add_argument('--agb-source-scope', choices=('co_common','kl16_envelopes_to7'),default='co_common')
    p.add_argument('--agb-net',choices=('unavailable','normalized_initial_MZY'),default='unavailable')
    p.add_argument('--massive-decay',choices=DECAY_MODELS,default=DECAY_MODELS[0],
                   help='LC18-only isotope projection; fully decayed AGB sources stay unchanged')
    args = p.parse_args()
    if args.output_dir.exists():
        p.error('existing output directory refused')
    table, history = build_combined(rotation=args.rotation, massive_wind_speed_km_s=args.massive_wind_speed_km_s,
                                   agb_wind_speed_km_s=args.agb_wind_speed_km_s, agb_release=args.agb_release,
                                   agb_energy=args.agb_energy, population=args.population, imf_id=args.imf_id,
                                   low_z_agb=args.low_z_agb, massive_source=args.massive_source,
                                   snii_energy_erg=args.snii_energy_erg, massive_wind_timing=args.massive_wind_timing,
                                   agb_source_scope=args.agb_source_scope, agb_net=args.agb_net,
                                   massive_decay=args.massive_decay)
    args.output_dir.mkdir(parents=True,exist_ok=False)
    for name, content in [('yields.dat',table),('history.nml',history)]:
        with (args.output_dir/name).open('x') as f:
            f.write(content)
        print(name,hashlib.sha256(content.encode()).hexdigest())


if __name__ == '__main__':
    main()
