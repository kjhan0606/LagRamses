#!/usr/bin/env python3
"""Matched two-carrier LC18 comparison; offline only, no simulation launch.

Consume an EXISTING prompt-projected yields.dat/history.nml pair. Verify its
LC18 rows against the existing converter and pinned isotopic sources, then
write a NEW companion without changing either input. This is abundance-only:
Al26 -> Mg26, Fe60 -> Ni60 (Co60 delay collapsed), transparent MeV energy,
no dust/CHIMES or full nuclide-network claim. LC18 integrated mean isotope
composition follows its selected wind-mass knots, not measured isotope timing.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

from adapt_g2_candidate_sources import adapt_candidate, LIMONGI_ID
from build_lc18_native_wind import build_native_wind

MODEL = 'lc18_al26_fe60_transparent_v1'
PROMPT = 'prompt_t12_le_100yr_baryonic_v1'
NUCLEAR_SHA = '1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859'
HALF_LIFE_S = (717e3*31557600., 2.62e6*31557600.)
ROOT = Path(__file__).resolve().parents[3]


def native_rows(text):
    rows = [[float(x) for x in line.split()] for line in text.splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    if not rows or any(len(r) != 32 or not all(math.isfinite(x) for x in r) for r in rows):
        raise ValueError('invalid native row format')
    return rows


def history_string(text, key):
    matches = re.findall(r'\b'+key+r"\s*=\s*'([^']*)'", text)
    if len(matches) != 1:
        raise ValueError('missing/duplicate history identity: '+key)
    return matches[0]


def history_array(text, key):
    match = re.search(r'\b'+key+r'\s*=\s*(.*?)(?=\b\w+\s*=|/)', text, re.S)
    if not match:
        raise ValueError('missing history coordinate: '+key)
    values = []
    for token in match[1].strip(' \n,').split(','):
        token = token.strip().lower().replace('d', 'e')
        if '*' in token:
            count, value = token.split('*')
            values.extend([float(value)]*int(count))
        else:
            values.append(float(token))
    return values


def build_companion(table_text, history_text):
    metadata = [json.loads(line[2:]) for line in table_text.splitlines() if line.startswith('# {')]
    selections = [m for m in metadata if m.get('massive_source') == 'lc18_set_r' and 'decay' in m]
    if len(selections) != 1 or selections[0]['decay'] != PROMPT:
        raise ValueError('requires matched prompt-projected LC18 Set R, not no-decay/PARSEC/decayed AGB')
    selection = selections[0]
    reference, ref_history = build_native_wind(rotation=selection['rotation_km_s'],
        wind_speed_km_s=selection['wind_speed_km_s'], timing=selection['timing'],
        composition=selection['composition'], energy=selection['energy'], imf_id=selection['imf_id'],
        massive_source='lc18_set_r', snii_energy_erg=selection['snii_energy_erg'], decay=PROMPT)
    rows, ref_rows = native_rows(table_text), native_rows(reference)
    selected_rows = [r for r in rows if r[0] in (1., 3.)]
    if selected_rows != ref_rows:
        raise ValueError('LC18 channel/M/Z/age/material rows differ from matched converter')
    identity_keys = ('model_id', 'wind_source_id', 'terminal_source_id', 'model_coordinates')
    identities = [history_string(history_text, k) for k in identity_keys]
    for key in ('wind_source_id', 'terminal_source_id', 'timing_policy', 'metallicity_policy'):
        if history_string(history_text, key) != history_string(ref_history, key):
            raise ValueError('source/history mismatch: '+key)
    for key in ('mass_msun', 'metallicity', 'terminal_age_yr', 'terminal_outcome', 'input_imf_id'):
        if history_array(history_text, key) != history_array(ref_history, key):
            raise ValueError('source/history coordinate mismatch: '+key)
    # The combined adapter puts its own selection first and changes model_id;
    # the high-mass source IDs remain the independently regenerated LC18 IDs.
    first = metadata[0]
    prefix = 'kl16_lc18_user_v1:' if 'kl16_files' in first else 'lc18_set_r_v2:'
    if identities[0] != prefix+hashlib.sha256(json.dumps(first, sort_keys=True).encode()).hexdigest():
        raise ValueError('table/history model identity mismatch')
    nuclear_path = ROOT/'external/g2_candidates/nuclear_decay/nubase_4.mas20'
    if hashlib.sha256(nuclear_path.read_bytes()).hexdigest() != NUCLEAR_SHA:
        raise ValueError('pinned NUBASE2020 bytes differ')
    report = adapt_candidate(LIMONGI_ID, include_records=True)
    components = report['source_components']
    def coordinate(record):
        c = record['source_model_coordinate']
        return c['rotation_velocity_km_s'], c['metallicity_feh'], c['initial_mass_msun']
    wind = {coordinate(r): r['source_reported_isotopic_yields'] for r in components['wind_yields']['records']}
    zmap = report['source_axes']['metallicity_mass_fraction_from_source_article']
    endpoints = {}
    for record in components['recommended_yields']['records']:
        rotation, feh, mass = coordinate(record)
        if rotation != selection['rotation_km_s']:
            continue
        total = record['source_reported_isotopic_yields']
        released_wind = wind[coordinate(record)] if mass <= 25 else total
        z = float(zmap[str(feh)])
        endpoints[(1., mass, z)] = [released_wind[n] for n in ('Al26', 'Fe60')]
        endpoints[(3., mass, z)] = [total[n]-released_wind[n] if mass <= 25 else 0.
                                   for n in ('Al26', 'Fe60')]
    final_return = {key: max(r[4] for r in ref_rows if tuple(r[:3]) == key) for key in endpoints}
    inventory_rows = []
    for row in rows:
        key = tuple(row[:3])
        parents = [0., 0.]
        if row[0] in (1., 3.):
            total = final_return[key]
            fraction = row[4]/total if total > 0 else 0.
            parents = [fraction*p for p in endpoints[key]]
        if min(parents) < 0 or not all(math.isfinite(p) for p in parents):
            raise ValueError('invalid source parent yield')
        tolerance = 64*math.ulp(1.)*max(row[4], 1e-300)
        if parents[0] > row[4]-math.fsum(row[10:21])+tolerance or parents[1] > row[20]+tolerance:
            raise ValueError('radioactive parent is not a subset of its host metal')
        inventory_rows.append(' '.join(format(x, '.17g') for x in row+parents))
    payload = '\n'.join(inventory_rows)+'\n'
    header = (f"&stellar_radioactive_companion\n version=1,model_id='{MODEL}',\n"
        f" source_projection='{PROMPT}',row_count={len(rows)},\n"
        ' source_identity='+','.join("'"+s+"'" for s in identities)+',\n'
        f" source_yields_sha256='{hashlib.sha256(table_text.encode()).hexdigest()}',\n"
        f" source_history_sha256='{hashlib.sha256(history_text.encode()).hexdigest()}',\n"
        f" nuclear_sha256='{NUCLEAR_SHA}',\n"
        f" inventory_sha256='{hashlib.sha256(payload.encode()).hexdigest()}',\n"
        ' half_life_s='+','.join(format(x, '.17g') for x in HALF_LIFE_S)+'\n/\n')
    return header+payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--yield-table', required=True, type=Path)
    parser.add_argument('--history', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path, help='NEW companion file; existing paths refused')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists')
    text = build_companion(args.yield_table.read_bytes().decode(), args.history.read_bytes().decode())
    with args.output.open('x') as stream:
        stream.write(text)
    print(MODEL, args.output, hashlib.sha256(text.encode()).hexdigest())


if __name__ == '__main__':
    main()
