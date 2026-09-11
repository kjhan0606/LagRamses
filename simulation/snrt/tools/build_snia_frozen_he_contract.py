#!/usr/bin/env python3
"""One native effective-hybrid sidecar from retained COSMIC system169.

No binary evolution is run. Contact and incident He transfer are prescribed,
NOT inferred from deltam/RRLO. stdout is the complete ordered native handoff.
The event replaces empirical Ia counts; ordinary SSP returns are untouched.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re

from helium_retention_reference import helium_regime

MODEL = 'tabulated_frozen_he_contact_v1'
CLOSURE = 'imposed_contact_stable_eta1_frozen135_v1'
APPROVAL = 'SNIA-FROZEN-HE-HYBRID-2026-09-11'
COSMIC_COMMIT = 'f9b90f451bca014e9e0adb3a410bee8752e30e53'
PARAMS_SHA = '68b2c0a7a90cd179e4f98c27428935472713e161ea0a19e5c8bf7eef8d3d04e9'
RATE = 2e-6
THRESHOLD = 1.4004633930489443
SYSTEM = 169


def csv_rows(path):
    raw = path.read_bytes()
    text = gzip.decompress(raw).decode() if path.suffix == '.gz' else raw.decode()
    return list(csv.DictReader(io.StringIO(text)))


def build(project: Path, history: Path) -> str:
    grid_path = history/'grid-v2/z0.01-grid.csv'
    selection_path = history/'grid-v2/selection.json'
    bpp_path = history/'grid-v2/z0.01-bpp.csv.gz'
    bcm_path = history/'grid-v2/z0.01-bcm.csv.gz'
    params_path = history/'Params-v4.2.0.ini'
    baseline_path = project/'simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml'
    n100_path = project/'simulation/snrt/data/fp2_snia_hesma_n100_approved_event_source_v1.json'
    selection = json.loads(selection_path.read_text())
    if (selection['cosmic_version'] != '4.2.0' or selection['params_sha256'] != PARAMS_SHA
            or hashlib.sha256(params_path.read_bytes()).hexdigest() != PARAMS_SHA):
        raise ValueError('requires retained pinned COSMIC4.2.0 comparison')
    grid = [r for r in csv_rows(grid_path) if int(r['bin_num']) == SYSTEM]
    if len(grid) != 1:
        raise ValueError('missing/duplicate selected system')
    grid = grid[0]
    primary, secondary = float(grid['m1']), float(grid['m2'])
    weight = float(grid['weight_per_initial_msun'])
    if (primary, secondary) != (7.5, 5.274) or not 0 < weight*(primary+secondary) < 1:
        raise ValueError('selected birth population changed')
    rows = [r for r in csv_rows(bpp_path) if int(float(r['bin_num'])) == SYSTEM]
    if not rows or float(rows[0]['tphys']) != 0:
        raise ValueError('original ZAMS clock missing')
    if (float(rows[0]['mass_1']), float(rows[0]['mass_2'])) != (primary, secondary):
        raise ValueError('component birth identity changed')
    if any(float(b['tphys']) < float(a['tphys']) for a, b in zip(rows, rows[1:])):
        raise ValueError('unordered provider phase records')
    contact = next(r for r in rows if float(r['kstar_1']) == 11
                   and 7 <= float(r['kstar_2']) <= 9)
    contact_age = float(contact['tphys'])*1e6
    wd, donor = float(contact['mass_1']), float(contact['mass_2'])
    if not float(contact['RRLO_2']) < 1:
        raise ValueError('this named comparison assumes detached post-CE seed, then imposes contact')
    event_age = contact_age+(THRESHOLD-wd)/RATE
    end_age = float(rows[-1]['tphys'])*1e6
    if not 0 < contact_age < event_age <= end_age:
        raise ValueError('threshold event is outside retained history')
    # BPP seed plus regular BCM knots bracket the imposed event. Only donor
    # background losses are retained; COSMIC WD accretion is NOT added again.
    background = [(contact_age, donor)]
    for r in csv_rows(bcm_path):
        if int(float(r['bin_num'])) != SYSTEM:
            continue
        age = float(r['tphys'])*1e6
        if age <= contact_age:
            continue
        if float(r['kstar_1']) != 11 or not 7 <= float(r['kstar_2']) <= 9:
            raise ValueError('retained pair changes channel before the imposed event bracket')
        background.append((age, float(r['mass_2'])))
        if age >= event_age:
            break
    if background[-1][0] < event_age:
        raise ValueError('donor history does not bracket event')
    for a, b in zip(background, background[1:]):
        if b[0] <= a[0] or b[1] > a[1]+1e-12:
            raise ValueError('nonmonotone donor background')
    a, b = background[-2:]
    donor_at_event = a[1]+(b[1]-a[1])*(event_age-a[0])/(b[0]-a[0])
    loss = donor-donor_at_event
    survivor = donor_at_event-(THRESHOLD-wd)
    if loss < 0 or survivor < 0:
        raise ValueError('prescribed transfer overdraws donor')
    for mass in (wd, 1.35):
        ref = helium_regime(mass, RATE)
        if ref['regime'] != 'stable_central_growth_candidate_only' or ref['eta_conditional'] != 1:
            raise ValueError('prescribed rate is outside adopted stable/central reference')

    n100 = json.loads(n100_path.read_text())
    event = n100['event']
    if (event['returned_mass_msun_per_event'] != THRESHOLD
            or event['wd_debit_msun_per_event'] != THRESHOLD
            or event['terminal_remnant_msun_per_event'] != 0):
        raise ValueError('N100 zero-remnant event mass changed')
    # Copy the existing physical/coupling groups verbatim. The separate new
    # comparison approval does not pretend to replace the physical-source ID.
    baseline = baseline_path.read_text()
    split = baseline.index('&snia_physical_contract')
    population, physical = baseline[:split], baseline[split:]
    for key in ('returned_mass_per_event', 'wd_debit_per_event'):
        match = re.search(rf'\b{key}\s*=\s*([^\s,]+)', physical)
        if match is None or float(match[1].replace('d', 'e')) != THRESHOLD:
            raise ValueError('runtime N100 template does not match threshold')
    population = population[population.index('&snia_population_realization'):]
    population = population.replace('&snia_population_realization',
                                    f"&snia_population_realization\n event_model='{MODEL}'", 1)
    population = re.sub(r'population_source_id\s*=\s*"[^"]*"',
                        'population_source_id="cosmic:system169:frozen_he_contact:effective_hybrid"', population)
    population = re.sub(r'events_per_initial_msun\s*=\s*[^\s,]+',
                        f'events_per_initial_msun={weight:.17e}', population)
    parameters = dict(birth_z=.01, age_yr=event_age, weight=weight, contact_age_yr=contact_age,
                      seed_wd=wd, seed_donor=donor, donor_background_loss=loss,
                      transfer_rate=RATE, threshold=THRESHOLD, history_end_age=end_age,
                      initial_binary_mass=primary+secondary, retention_reference_z=.02,
                      retention_mass_ceiling=1.35)
    if not all(math.isfinite(v) for v in parameters.values()):
        raise ValueError('nonfinite event parameter')
    paths = [bpp_path, bcm_path, grid_path, selection_path, params_path, Path(__file__),
             Path(__file__).with_name('helium_retention_reference.py'), n100_path, baseline_path]
    lines = ['&snia_frozen_he_event']
    lines += [f' row%{key}={value:.17e}' for key, value in parameters.items()]
    lines += [f' row%history_bin_num={SYSTEM}', f" row%closure_id='{CLOSURE}'",
              f" row%comparison_approval='{APPROVAL}'", f" row%source_commit='{COSMIC_COMMIT}'"]
    lines += [f" row%source_sha256({i})='{hashlib.sha256(path.read_bytes()).hexdigest()}'"
              for i, path in enumerate(paths, 1)]
    lines += ['/', '']
    return ('! Effective hybrid, NOT a self-consistent WD/donor/ordinary-SSP population.\n'
            '! Contact and 2e-6Msun/yr incident He feed are prescribed, NOT inferred.\n'
            '! Original-ZAMS age; full-initial-SSP weight, no extra binary/IMF scaling.\n'
            f'! Imposed-event donor survivor (comparison only): {survivor:.17g} Msun.\n'
            + population+'\n'.join(lines)+'\n'+physical)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument('--history-dir', type=Path, required=True)
    args = parser.parse_args()
    print(build(args.project, args.history_dir), end='')


if __name__ == '__main__':
    main()
