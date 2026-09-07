#!/usr/bin/env python3
"""LC18 -> native user-selected wind-only model, NOT production approval.

Reuse the checksum-verified reader and phase lifetime aggregation. Offline
input adapter only: no Python inside RAMSES. Raw sources and historical
approval contracts remain unchanged. Derived files are local-use artifacts;
redistribution permission is not inferred from public access.
"""
from __future__ import annotations
import argparse
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
                      imf_id: int = 2) -> tuple[str, str]:
    if rotation not in (0, 150, 300) or imf_id not in (0, 1, 2, 4):
        raise ValueError('unsupported rotation/IMF selection')
    if not math.isfinite(wind_speed_km_s) or not 0 < wind_speed_km_s < 29979.2458:
        raise ValueError('supply a positive nonrelativistic (<0.1c) effective wind speed; no default is inferred')
    if (timing, composition, energy) != ('uniform_until_terminal', 'as_tabulated_mean', 'isotropic_thermalized'):
        raise ValueError('unsupported release/composition/energy model')
    report = adapt_candidate(LIMONGI_ID, include_records=True)
    components = report['source_components']
    histories, _ = build_phase_histories(components['evolutionary_properties']['records'],
                                       ['MS', 'H', 'He', 'C', 'Ne', 'O', 'Si', 'PSN'])
    zmap = report['source_axes']['metallicity_mass_fraction_from_source_article']
    nodes = []
    for record in components['recommended_yields']['records']:
        c = record['source_model_coordinate']
        m = c['initial_mass_msun']; feh = c['metallicity_feh']
        if c['rotation_velocity_km_s'] != rotation or not 40 <= m <= 120:
            continue
        # Set R table8 is wind-only above 25 Msun. Do not subtract table9
        # (available only for lower masses) or merge BR26 terminal ejecta.
        isotopes = record['source_reported_isotopic_yields']
        grouped = {e: [] for e in ELEMENTS}
        for isotope, value in isotopes.items():
            match = re.fullmatch(r'([A-Z][a-z]?)(\d*)', isotope)
            if match is None or not math.isfinite(value) or value < 0:
                raise ValueError(f'invalid source isotope/value: {isotope}')
            if match[1] in grouped:
                grouped[match[1]].append(value)
        ejecta = [math.fsum(grouped[e]) for e in ELEMENTS]
        total = math.fsum(isotopes.values())
        age = histories[(rotation, feh, m)]['terminal_age_yr']
        if not math.isfinite(age) or age <= 0 or not 0 < total <= m or sum(ejecta) > total*(1+1e-12):
            raise ValueError(f'invalid wind budget/lifetime at {(rotation, feh, m)}')
        # Integrated isotope sum is the declared mass authority. Rounded
        # phase endpoint masses are NOT forced to agree with that total.
        kinetic = .5*total*SOLAR_MASS_CGS*(wind_speed_km_s*1e5)**2
        if not math.isfinite(kinetic):
            raise ValueError('wind energy overflow')
        nodes.append((float(zmap[str(feh)]), m, age, total, ejecta, kinetic))
    nodes.sort(key=lambda n: (n[0], n[1]))
    if len(nodes) != 16 or any([n[1] for n in nodes if n[0] == z] != [40,60,80,120] for z in set(n[0] for n in nodes)):
        raise ValueError('expected complete 4-Z x 4-mass LC18 high-mass branch')
    source_files = report['verified_acquisition']['verified_files']
    selection = dict(source_files=source_files, rotation_km_s=rotation, wind_speed_km_s=wind_speed_km_s,
                     timing=timing, composition=composition, energy=energy, imf_id=imf_id,
                     decay='as_tabulated_no_decay', net='unavailable_diagnostic_zero',
                     mass_authority='sum_all_table8_isotopes', version=1)
    identity = hashlib.sha256(json.dumps(selection,sort_keys=True).encode()).hexdigest()
    header = ['# LC18 Set R high-mass wind; local user-selected approximation, NOT publication/production approval.',
              '# Citation: Limongi & Chieffi 2018, ApJS 237:13; VizieR/CDS J/ApJS/237/13.',
              '# No isotope decay, net-yield inference, radial momentum or terminal explosion.',
              '# Uniform release and fixed mean composition; isotropic kinetic energy thermalized.',
              '# net columns are unavailable diagnostic placeholders, NOT physical zero net production.',
              '# ' + json.dumps(selection,sort_keys=True)]
    # A common union of actual lifetimes satisfies the existing rectangular
    # age grid. No dense time sampling or large synthetic atlas is needed.
    ages = sorted({0., 2e10, *(n[2] for n in nodes)})
    rows = []
    for z,m,life,total,ejecta,kinetic in nodes:
        for ch in (1,3):
            for age in ages:
                f = min(age/life,1.) if ch == 1 else 0.
                remnant = m-total if ch == 3 and age >= life else 0.
                row = [ch,m,z,age,f*total,remnant,f*kinetic,0.,0.,0.,
                       *(f*v for v in ejecta),*([0.]*11)]
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
    return '\n'.join(header+rows)+'\n', history


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True, help='NEW local directory; existing paths refused')
    p.add_argument('--rotation', type=int, choices=(0,150,300), required=True)
    p.add_argument('--wind-speed-km-s', type=float, required=True)
    p.add_argument('--timing', choices=('uniform_until_terminal',), required=True)
    p.add_argument('--composition', choices=('as_tabulated_mean',), required=True)
    p.add_argument('--energy', choices=('isotropic_thermalized',), required=True)
    p.add_argument('--imf-id', type=int, choices=(0,1,2,4), default=2)
    args = p.parse_args()
    if args.output_dir.exists():
        p.error('output directory already exists; source artifacts must not be overwritten')
    table, history = build_native_wind(rotation=args.rotation, wind_speed_km_s=args.wind_speed_km_s,
                                     timing=args.timing, composition=args.composition, energy=args.energy, imf_id=args.imf_id)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for name, content in [('yields.dat',table),('history.nml',history)]:
        with (args.output_dir/name).open('x') as f:
            f.write(content)
        print(name, hashlib.sha256(content.encode()).hexdigest())


if __name__ == '__main__':
    main()
