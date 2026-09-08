#!/usr/bin/env python3
"""LC18 -> native user-selected wind-only model, NOT production approval.

Reuse the checksum-verified reader and phase lifetime aggregation. Offline
input adapter only: no Python inside RAMSES. Raw sources and historical
approval contracts remain unchanged. Derived files are local-use artifacts;
redistribution permission is not inferred from public access.
"""
from __future__ import annotations
import argparse
from bisect import bisect_right
import hashlib
import json
import math
from pathlib import Path
import re
from adapt_g2_candidate_sources import adapt_candidate, LIMONGI_ID
from fp1_limongi_phase_history import build_phase_histories

ELEMENTS = ('H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe')
SOLAR_MASS_CGS = 1.98847e33


def build_native_wind(*, rotation: int, wind_speed_km_s: float,
                      timing: str, composition: str, energy: str,
                      imf_id: int = 2, massive_source: str = 'wind_only',
                      snii_energy_erg: float | None = None) -> tuple[str, str]:
    if massive_source not in ('wind_only', 'lc18_set_r'):
        raise ValueError('unsupported massive source')
    full = massive_source == 'lc18_set_r'
    if full:
        if snii_energy_erg is None or not math.isfinite(snii_energy_erg) or snii_energy_erg <= 0:
            raise ValueError('Set R requires explicit positive SN energy; LC18 binding energy is NOT injection energy')
    elif snii_energy_erg is not None:
        raise ValueError('SN energy is only meaningful for lc18_set_r')
    if rotation not in (0, 150, 300) or imf_id not in (0, 1, 2, 4):
        raise ValueError('unsupported rotation/IMF selection')
    if not math.isfinite(wind_speed_km_s) or not 0 < wind_speed_km_s < 29979.2458:
        raise ValueError('supply a positive nonrelativistic (<0.1c) effective wind speed; no default is inferred')
    if timing not in ('uniform_until_terminal', 'phase_mass_loss_or_uniform') or \
            (composition, energy) != ('as_tabulated_mean', 'isotropic_thermalized'):
        raise ValueError('unsupported release/composition/energy model')
    report = adapt_candidate(LIMONGI_ID, include_records=True)
    components = report['source_components']
    histories, _ = build_phase_histories(components['evolutionary_properties']['records'],
                                       ['MS', 'H', 'He', 'C', 'Ne', 'O', 'Si', 'PSN'])
    zmap = report['source_axes']['metallicity_mass_fraction_from_source_article']
    def coordinate(record):
        c = record['source_model_coordinate']
        return (c['rotation_velocity_km_s'], c['metallicity_feh'], c['initial_mass_msun'])
    winds = {coordinate(r): r['source_reported_isotopic_yields']
             for r in components['wind_yields']['records']}
    def project(isotopes):
        grouped = {e: [] for e in ELEMENTS}
        for isotope, value in isotopes.items():
            match = re.fullmatch(r'([A-Z][a-z]?)(\d*)', isotope)
            if match is None or not math.isfinite(value) or value < 0:
                raise ValueError(f'invalid source isotope/value: {isotope}')
            if match[1] in grouped:
                grouped[match[1]].append(value)
        return math.fsum(isotopes.values()), [math.fsum(grouped[e]) for e in ELEMENTS]
    nodes = []
    wind_shapes = {}
    unresolved_phase_nodes = []
    for record in components['recommended_yields']['records']:
        c = record['source_model_coordinate']
        m = c['initial_mass_msun']; feh = c['metallicity_feh']
        if c['rotation_velocity_km_s'] != rotation or not (13 if full else 40) <= m <= 120:
            continue
        # Set R table8 is wind-only above 25 Msun. Do not subtract table9
        # (available only for lower masses) or merge BR26 terminal ejecta.
        recommended = record['source_reported_isotopic_yields']
        isotopes = winds[coordinate(record)] if full and m <= 25 else recommended
        total, ejecta = project(isotopes)
        terminal_total, terminal_ejecta = 0., [0.]*len(ELEMENTS)
        if full and m <= 25:
            if recommended.keys() != isotopes.keys():
                raise ValueError('Set R and wind isotope sets differ')
            # Table 8 includes winds. Subtract once, isotope by isotope;
            # reject negative differences, never clip or double-count wind.
            terminal_total, terminal_ejecta = project({k: recommended[k]-isotopes[k] for k in isotopes})
            if terminal_total <= 0 or total + terminal_total > m:
                raise ValueError('non-closing LC18 terminal endpoint')
        age = histories[(rotation, feh, m)]['terminal_age_yr']
        if not math.isfinite(age) or age <= 0 or not 0 < total <= m or sum(ejecta) > total*(1+1e-12):
            raise ValueError(f'invalid wind budget/lifetime at {(rotation, feh, m)}')
        # Integrated isotope sum is the declared mass authority. Rounded
        # phase endpoint masses are NOT forced to agree with that total.
        kinetic = .5*total*SOLAR_MASS_CGS*(wind_speed_km_s*1e5)**2
        if not math.isfinite(kinetic):
            raise ValueError('wind energy overflow')
        node = (float(zmap[str(feh)]), m, age, total, ejecta, kinetic)
        if full:
            node += (terminal_total, terminal_ejecta, snii_energy_erg if m <= 25 else 0.)
        nodes.append(node)
        if timing == 'phase_mass_loss_or_uniform':
            phase = histories[(rotation, feh, m)]
            endpoint = phase['terminal_cumulative_wind_mass_msun']
            if endpoint == 0:
                # Printed phase masses cannot resolve these small winds.
                # The explicitly named option retains a uniform fallback;
                # never invent a phase or erase nonzero isotopic wind mass.
                unresolved_phase_nodes.append(dict(mass=m, feh=feh, wind_mass_msun=total))
                shape = [(0., 0.), (age, 1.)]
            else:
                shape = [(0., 0.)] + [(p['cumulative_age_yr'], p['cumulative_wind_mass_msun']/endpoint)
                                      for p in phase['nodes']]
                if any(b[0] <= a[0] or b[1] < a[1] for a,b in zip(shape,shape[1:])):
                    raise ValueError('invalid phase mass-loss shape')
            wind_shapes[(node[0],m)] = shape
    nodes.sort(key=lambda n: (n[0], n[1]))
    masses = [13,15,20,25,30,40,60,80,120] if full else [40,60,80,120]
    if len(nodes) != 4*len(masses) or any([n[1] for n in nodes if n[0] == z] != masses for z in set(n[0] for n in nodes)):
        raise ValueError('incomplete LC18 mass/Z branch')
    source_files = report['verified_acquisition']['verified_files']
    selection = dict(source_files=source_files, rotation_km_s=rotation, wind_speed_km_s=wind_speed_km_s,
                     timing=timing, composition=composition, energy=energy, imf_id=imf_id,
                     decay='as_tabulated_no_decay', net='unavailable_diagnostic_zero',
                     mass_authority='sum_all_table8_isotopes', version=1)
    if full:
        selection.update(version=2, massive_source=massive_source, snii_energy_erg=snii_energy_erg,
            terminal='table8_minus_table9_at_13_15_20_25; wind_only_at_30_and_above',
            injected_energy='explicit_constant_comparison_parameter_not_LC18_binding_energy',
            remnant='baryonic_residual_initial_minus_wind_minus_terminal_not_gravitational_mass',
            mass_cells='nearest_source_node_fraction_scaling; separate_40_seam; 25/30_boundary=27.5',
            unsupported_mass_msun=[8,13])
    if timing == 'phase_mass_loss_or_uniform':
        selection.update(wind_timing='piecewise_linear_table5_cumulative_loss_fraction',
            wind_endpoint='isotope_sum_unchanged; phase_masses_set_shape_only',
            unresolved_phase_uniform_nodes=unresolved_phase_nodes)
    identity = hashlib.sha256(json.dumps(selection,sort_keys=True).encode()).hexdigest()
    header = ['# LC18 Set R high-mass wind; local user-selected approximation, NOT publication/production approval.',
              '# Citation: Limongi & Chieffi 2018, ApJS 237:13; VizieR/CDS J/ApJS/237/13.',
              '# No isotope decay, net-yield inference, radial momentum or terminal explosion.',
              '# Uniform release and fixed mean composition; isotropic kinetic energy thermalized.',
              '# net columns are unavailable diagnostic placeholders, NOT physical zero net production.',
              '# ' + json.dumps(selection,sort_keys=True)]
    if full:
        header[0] = '# LC18 Set R wind + CCSN, 13--120 Msun; explicit comparison model, not publication approval.'
        header[2] = '# No isotope decay, net-yield inference or directed momentum; terminal ejecta at source lifetimes.'
    if wind_shapes:
        header[3] = '# Table5 mass-loss timing; explicit uniform fallback for unresolved phase mass losses; fixed composition/speed.'
    # A common union of actual lifetimes satisfies the existing rectangular
    # age grid. No dense time sampling or large synthetic atlas is needed.
    ages = sorted({0., 2e10, *(n[2] for n in nodes)})
    rows = []
    for node in nodes:
        z,m,life,total,ejecta,kinetic = node[:6]
        terminal_total, terminal_ejecta, terminal_energy = node[6:] if full else (0., [0.]*11, 0.)
        for ch in (1,3):
            channel_ages = ages
            if wind_shapes:
                # Native source-node interpolation consumes each star's own
                # knots. Do not multiply every phase by all other M/Z nodes.
                channel_ages = sorted({0.,life,2e10,*(t for t,_ in wind_shapes[(z,m)])}) if ch==1 else [0.,life,2e10]
            for age in channel_ages:
                f = min(age/life,1.) if ch == 1 else 0.
                if ch == 1 and wind_shapes:
                    shape = wind_shapes[(z,m)]
                    index = bisect_right([t for t,_ in shape], age)-1
                    if index == len(shape)-1:
                        f = 1.
                    else:
                        t0,f0 = shape[index]; t1,f1 = shape[index+1]
                        f = f0+(f1-f0)*(age-t0)/(t1-t0)
                terminal_f = float(ch == 3 and age >= life)
                remnant = m-total-terminal_total if terminal_f else 0.
                row = [ch,m,z,age,f*total+terminal_f*terminal_total,remnant,
                       f*kinetic+terminal_f*terminal_energy,0.,0.,0.,
                       *(f*v+terminal_f*t for v,t in zip(ejecta,terminal_ejecta)),*([0.]*11)]
                rows.append(' '.join(str(v) if isinstance(v,int) else format(v,'.17g') for v in row))
    def array(index):
        return ','.join(format(n[index],'.17g') for n in nodes)
    table8 = next(f['sha256'] for f in source_files if f['relative_path'].endswith('/table8.dat'))
    history = f'''! Explicit LC18 local-use approximation. No physical-package approval claimed.
&stellar_high_mass_history
 version=1,node_count={len(nodes)},
 model_id='lc18_user_wind_v1:{identity}',
 wind_source_id='LC18:table8:{table8}',terminal_source_id='user_wind_only_collapse',
 model_coordinates='rotation={rotation}; vwind={wind_speed_km_s:g}km/s; uniform; mean-composition; no-decay; thermalized',
 timing_policy='wind_linear_terminal_step',metallicity_policy='linear_Z_cumulative_mixture',
 net_yield_policy='unavailable_diagnostic_zero',
 input_imf_id={imf_id},input_population_id=0,input_binary_fraction=0,
 input_imf_min=.08,input_imf_max=120,
 mass_msun={array(1)},
 metallicity={array(0)},
 terminal_age_yr={array(2)},
 terminal_outcome=16*0
/
'''
    if full:
        history = history.replace('version=1,', 'version=2,')
        history = history.replace('lc18_user_wind_v1:', 'lc18_set_r_v2:')
        history = history.replace(f"wind_source_id='LC18:table8:{table8}',terminal_source_id='user_wind_only_collapse'",
            f"wind_source_id='LC18:SetR:{identity}',terminal_source_id='LC18:SetR:{identity}'")
        history = history.replace('terminal_outcome=16*0',
            'terminal_outcome='+','.join(str(int(n[6]>0)) for n in nodes))
    if wind_shapes:
        history = history.replace('; uniform;', '; table5-shape/fallback;')
    return '\n'.join(header+rows)+'\n', history


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True, help='NEW local directory; existing paths refused')
    p.add_argument('--rotation', type=int, choices=(0,150,300), required=True)
    p.add_argument('--wind-speed-km-s', type=float, required=True)
    p.add_argument('--timing', choices=('uniform_until_terminal','phase_mass_loss_or_uniform'), required=True)
    p.add_argument('--composition', choices=('as_tabulated_mean',), required=True)
    p.add_argument('--energy', choices=('isotropic_thermalized',), required=True)
    p.add_argument('--imf-id', type=int, choices=(0,1,2,4), default=2)
    p.add_argument('--massive-source', choices=('wind_only','lc18_set_r'), default='wind_only')
    p.add_argument('--snii-energy-erg', type=float, help='Required explicit comparison parameter for lc18_set_r')
    args = p.parse_args()
    if args.output_dir.exists():
        p.error('output directory already exists; source artifacts must not be overwritten')
    table, history = build_native_wind(rotation=args.rotation, wind_speed_km_s=args.wind_speed_km_s,
                                     timing=args.timing, composition=args.composition, energy=args.energy, imf_id=args.imf_id,
                                     massive_source=args.massive_source, snii_energy_erg=args.snii_energy_erg)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for name, content in [('yields.dat',table),('history.nml',history)]:
        with (args.output_dir/name).open('x') as f:
            f.write(content)
        print(name, hashlib.sha256(content.encode()).hexdigest())


if __name__ == '__main__':
    main()
