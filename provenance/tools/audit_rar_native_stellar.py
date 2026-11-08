#!/usr/bin/env python3
"""Independent endpoint source, birth, KDK and saved-field checks for the bounded fixture."""
import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
from rar_native_particles import epoch
from rar_evolving_thermal import thermal_accounting

p = argparse.ArgumentParser()
p.add_argument('run', type=Path)
p.add_argument('out', type=Path)
p.add_argument('--sinks', action='store_true')
a = p.parse_args()
text = (a.run/'run.log').read_text()
if not (a.run/'COMPLETE').is_file() or 'Run completed' not in text:
    raise ValueError('incomplete native stellar execution')


def records(content, marker, count, with_step=False):
    found = []
    step = -1
    lines = iter(content.splitlines())
    for line in lines:
        if 'RAR_QUMOND_EVOLVE_END_V1' in line:
            step = int(line.split('RAR_QUMOND_EVOLVE_END_V1', 1)[1].split()[0])
        if marker not in line:
            continue
        tokens = line.split(marker, 1)[1].split()
        while len(tokens) < count:
            tokens.extend(next(lines).split())
        if len(tokens) != count:
            raise ValueError('malformed record: '+marker)
        # ES output may omit the E for a three-digit exponent (e.g. 1.2-102).
        # Preserve that Fortran value rather than rejecting a tiny receipt.
        row = [float(re.sub(r'^([+-]?\d*\.?\d+)([+-]\d+)$', r'\1E\2', v.replace('D', 'E')))
               for v in tokens]
        if with_step:
            row = [step, *row]
        found.append(row)
    result = np.asarray(found).reshape(-1, count+int(with_step))
    if not np.isfinite(result).all():
        raise ValueError('nonfinite record: '+marker)
    return result


def close(got, expected, label, atol=1e-12, rtol=1e-10):
    if not np.allclose(got, expected, rtol=rtol, atol=atol):
        error = float(np.max(abs(np.asarray(got)-np.asarray(expected)), initial=0.))
        raise ValueError(f'{label}: maximum residual {error}')


def rows(path, columns):
    lines = path.read_text().splitlines()
    result = np.loadtxt(path, ndmin=2) if any(v.strip() and not v.startswith('#') for v in lines) else np.empty((0, columns))
    if result.shape[1] != columns or not np.isfinite(result).all():
        raise ValueError('invalid export '+str(path))
    return result[np.argsort(result[:, 1])]


starts = records(text, 'RAR_QUMOND_EVOLVE_BEGIN_V1', 5)
ends = records(text, 'RAR_QUMOND_EVOLVE_END_V1', 5)
steps = len(ends)
if steps != 4 or len(starts) != steps+1:
    raise ValueError('unexpected completed-step coverage')
dt = ends[:, 3]
times = np.r_[0., np.cumsum(dt)]
close(ends[:, 0], np.arange(steps), 'step sequence', atol=0, rtol=0)
close(ends[:, 1:3], np.column_stack((times[:-1], times[1:])), 'endpoint clocks', atol=1e-16, rtol=0)
close(starts[:steps, 2], dt, 'common timestep', atol=0, rtol=0)
close(starts[:, 1], times, 'begin time', atol=1e-16, rtol=0)
close(starts[:, 3], times, 'begin field source time', atol=1e-16, rtol=0)
generation_rows = []
lines = iter(text.splitlines())
for line in lines:
    if 'RAR_QUMOND_SOURCE_GENERATION_V1' not in line:
        continue
    tokens = line.split('RAR_QUMOND_SOURCE_GENERATION_V1', 1)[1].split()
    while len(tokens) < 4:
        tokens.extend(next(lines).split())
    if len(tokens) != 4 or tokens[3] not in ('T', 'F'):
        raise ValueError('invalid source generation receipt')
    generation_rows.append((int(tokens[0]), float(tokens[1]), int(tokens[2]), tokens[3]=='T'))
generation_expected = [(0, 0., 0, False)]
for k in range(steps):
    generation_expected.extend([(k, times[k+1], k, False), (k, times[k+1], k+1, True)])
if generation_rows != generation_expected:
    raise ValueError('same-time source generation was skipped or repeated')
paths = sorted(a.run.glob('composite_epoch_*_leaf.tsv'))
if len(paths) != 2*steps+1:
    raise ValueError('expected pre-source and post-source fields at each endpoint')
leaves = [np.loadtxt(path) for path in paths]
reference = leaves[0]
volume = reference[:, 5]**3
lookup = {int(v): i for i, v in enumerate(reference[:, 0])}
particles, gas, particle_reports, sinks, sink_reports = [], [], [], [], []
for k, (path, state) in enumerate(zip(paths, leaves)):
    if not np.array_equal(state[:, :6], reference[:, :6]):
        raise ValueError('fixed topology changed')
    time = float(path.read_text().splitlines()[1].split()[3])
    close(time, times[(k+1)//2], 'field clock', atol=1e-16, rtol=0)
    part, density, report = epoch(a.run/f'particles_epoch_{k:05d}.tsv', state, role='baryon')
    if a.sinks:
        body, body_density, body_report = epoch(a.run/f'sinks_epoch_{k:05d}.tsv', state, role='physical_sink')
        sinks.append(body); sink_reports.append(body_report)
        density += body_density
    particles.append(part)
    gas.append(state[:, 6]-density)
    particle_reports.append(report)
    snap = k in (4, 8)
    args = [sys.executable, str(a.out/'verify.py'), str(a.run.parent), str(a.out/f'field_{k:02d}.json'),
            '--case', a.run.name, '--evolving', '--particles', '--particle-role', 'baryon',
            '--leaf-file', path.name, '--snapshots', '1' if snap else '0',
            '--snapshot-glob', f'output_{k//4:05d}/data_*.h5' if snap else 'no_snapshot_for_this_epoch']
    if a.sinks:
        args.append('--sinks')
    subprocess.run(args, check=True)
if len(particles[0]) != 0:
    raise ValueError('fixture must begin with no stellar particles')
total_mass = np.array([state[:, 6] @ volume for state in leaves])
if not a.sinks:
    close(total_mass, total_mass[0], 'physical gas plus stellar mass', atol=1e-12, rtol=0)

births = records(text, 'RAR_STELLAR_BIRTH_V1', 13)
returns = records(text, 'RAR_STELLAR_COMMIT_V1', 28, with_step=True)
if not len(births) or not len(returns) or set(births[:, 1]) != set(reference[:, 1]):
    raise ValueError('birth/return must actually occur on both refinement levels')
if set(returns[:, 2]) != set(reference[:, 1]):
    raise ValueError('stellar return did not exercise both refinement levels')
if len(np.unique(births[:, 3])) != len(births):
    raise ValueError('duplicate newborn identity')
birth_by_id = {int(row[3]): row for row in births}
for birth in births:
    step, level, cell, identity = map(int, birth[:4])
    when, vol, mass, before, after = birth[4:9]
    close(when, times[step+1], 'birth clock', atol=1e-16, rtol=0)
    close(vol, volume[lookup[cell]], 'donor volume', atol=0, rtol=0)
    close(mass, (before-after)*vol, 'accepted birth debit/credit', atol=2e-18, rtol=1e-12)
    if level != int(reference[lookup[cell], 1]) or mass <= 0 or after < .5*before-1e-14:
        raise ValueError('invalid birth cap/donor')
    close(mass/1e-5, round(mass/1e-5), 'integer accepted star count', atol=1e-10, rtol=0)
cap_fraction = births[:, 6]/(births[:, 7]*births[:, 5])
if np.max(cap_fraction) < .49 or 'Lost=' not in text:
    raise ValueError('synthetic case did not exercise high depletion')

inventories = {}
inventory_text = re.sub(r'(RAR_GLOBAL_INVENTORY_V1)\s*', r'\1', text)
for label in ('before_native_return', 'after_native_return', 'before_native_birth', 'after_native_birth'):
    for rec in records(inventory_text, 'RAR_GLOBAL_INVENTORY_V1'+label, 33):
        # 3 header numbers + NVAR=18 gas totals + 12 census values.
        step, background, nv = rec[:3]
        if nv != 18 or background != 1:
            raise ValueError('inventory units/shape mismatch')
        key = (int(step), label)
        if key in inventories:
            raise ValueError('duplicate source inventory')
        inventories[key] = rec[3:]
expected = {(k, label) for k in range(steps) for label in
            ('before_native_return', 'after_native_return', 'before_native_birth', 'after_native_birth')}
if set(inventories) != expected:
    raise ValueError('source inventory coverage mismatch')
budgets = []
max_kdk = 0.
max_gas_kick = 0.
for k in range(steps):
    previous, pre, post = particles[2*k:2*k+3]
    moved = rows(a.run/f'particle_step_{k:05d}.tsv', 26)
    if moved.shape[0] != len(previous) or not np.array_equal(moved[:, :4], previous[:, :4]):
        raise ValueError('newborn entered an earlier KDK roster')
    if not np.array_equal(pre[:, :4], previous[:, :4]):
        raise ValueError('source mass changed inside KDK')
    close(moved[:, 4], dt[k], 'particle step length', atol=0, rtol=0)
    x0, v0, f0, vh, x1, f1, v1 = [moved[:, j:j+3] for j in (5, 8, 11, 14, 17, 20, 23)]
    comparisons = [(x0, previous[:, 4:7]), (v0, previous[:, 7:10]), (f0, previous[:, 10:13]),
                   (f1, pre[:, 10:13]), (vh, v0+.5*dt[k]*f0), (x1, np.mod(x0+dt[k]*vh, 1.)),
                   (v1, vh+.5*dt[k]*f1), (pre[:, 4:7], x1), (pre[:, 7:10], vh)]
    for got, wanted in comparisons:
        close(got, wanted, 'old-roster KDK', atol=1e-13, rtol=0)
        max_kdk = max(max_kdk, float(np.max(abs(got-wanted), initial=0.)))
    after = {int(row[1]): row for row in post}
    old_ids = set(map(int, previous[:, 1]))
    for i, old in enumerate(previous):
        now = after[int(old[1])]
        close(now[4:7], x1[i], 'survivor position', atol=0, rtol=0)
        close(now[7:10], v1[i], 'isotropic return velocity', atol=0, rtol=0)
        if now[3] > old[3] or now[3] < 0:
            raise ValueError('invalid surviving stellar mass')
    new_ids = set(after)-old_ids
    batch_births = births[births[:, 0] == k]
    if new_ids != set(map(int, batch_births[:, 3])):
        raise ValueError('unrecorded birth or deleted old body')
    for birth in batch_births:
        now = after[int(birth[3])]
        close(now[3], birth[6], 'newborn mass', atol=0, rtol=0)
        close(now[4:7], reference[lookup[int(birth[2])], 2:5], 'newborn not drifted', atol=0, rtol=0)
        close(now[7:10], birth[9:12], 'newborn synchronized gas velocity', atol=0, rtol=0)
    back = returns[returns[:, 0] == k, 1:]
    sums = back[:, 10:].sum(axis=0)
    before = inventories[(k, 'before_native_return')]
    returned = inventories[(k, 'after_native_return')]
    born0 = inventories[(k, 'before_native_birth')]
    born1 = inventories[(k, 'after_native_birth')]
    for inv in (before, returned, born0, born1):
        if inv[26] != 0 or inv[27] != 0 or inv[18] != len(reference):
            raise ValueError('invalid/duplicated list or leaf census')
    close(returned[:18]-before[:18], sums, 'whole return gas receipt', atol=2e-13, rtol=1e-8)
    close(before[22]-returned[22], np.sum(back[:, 6]), 'whole stellar debit', atol=2e-13, rtol=1e-8)
    close(born0[0]-born1[0], np.sum(batch_births[:, 6]), 'whole birth gas debit', atol=2e-13, rtol=1e-10)
    close(born1[22]-born0[22], np.sum(batch_births[:, 6]), 'whole newborn credit', atol=2e-13, rtol=1e-10)
    birth_dp = np.sum(batch_births[:, 6, None]*batch_births[:, 9:12], axis=0)
    birth_de = np.sum(batch_births[:, 6]*(batch_births[:, 12]+.5*np.sum(batch_births[:, 9:12]**2, axis=1)))
    close(born0[1:4]-born1[1:4], birth_dp, 'birth momentum exchange', atol=2e-13, rtol=1e-10)
    close(born0[4]-born1[4], birth_de, 'birth gas internal and kinetic energy removal', atol=2e-13, rtol=1e-10)
    old_by_slot = {int(row[0]): (i, row) for i, row in enumerate(previous)}
    for receipt in back:
        i, old = old_by_slot[int(receipt[2])]
        identity = int(old[1]);birth = birth_by_id[identity]
        close(receipt[4], times[k+1]-birth[4], 'SSP age progress', atol=1e-16, rtol=0)
        close(receipt[11:14], receipt[6]*v1[i], 'isotropic returned momentum', atol=1e-17, rtol=1e-8)
        close(old[3]-after[identity][3], receipt[6], 'individual stellar mass debit', atol=1e-17, rtol=1e-10)
    kicks = np.loadtxt(a.run/f'gravity_kicks_{k:05d}.tsv')
    if kicks.shape != (2*len(reference), 17) or set(kicks[:, 0]) != {0, 1}:
        raise ValueError('extra or missing gas gravity kick')
    for phase in (0, 1):
        row = kicks[kicks[:, 0] == phase]
        if len(row) != len(reference) or len(np.unique(row[:, 1])) != len(reference):
            raise ValueError('gravity kick does not cover each physical leaf once')
        row = row[np.argsort([lookup[int(v)] for v in row[:, 1]])]
        source = 2*k+phase
        close(row[:, 3], gas[source], 'gas kick source density', atol=2e-14, rtol=1e-13)
        close(row[:, 4], gas[2*k+1], 'gas kick receiver density', atol=2e-14, rtol=1e-13)
        close(row[:, 5], .5*dt[k], 'gas half-step', atol=0, rtol=0)
        close(row[:, 6:9], leaves[source][:, 10:13], 'gas force epoch before stellar source', atol=0, rtol=0)
        dp = row[:, 13:16]-row[:, 9:12]
        wanted_dp = row[:, 3, None]*row[:, 5, None]*row[:, 6:9]
        de = row[:, 16]-row[:, 12]
        wanted_de = (np.sum(row[:, 13:16]**2, axis=1)-np.sum(row[:, 9:12]**2, axis=1))/(2*row[:, 4])
        close(dp, wanted_dp, 'gas gravity momentum source', atol=2e-14, rtol=0)
        close(de, wanted_de, 'gas gravity work', atol=2e-14, rtol=0)
        max_gas_kick = max(max_gas_kick, float(np.max(abs(dp-wanted_dp))), float(np.max(abs(de-wanted_de))))
    budgets.append({'step': k, 'new_stars': len(batch_births), 'birth_mass': float(batch_births[:, 6].sum()),
                    'return_events': len(back), 'return_mass': float(back[:, 6].sum()),
                    'maximum_return_inventory_residual': float(np.max(abs(returned[:18]-before[:18]-sums)))})

nml = (a.run/'effective.nml').read_text()
unit_t = float(re.search(r'(?im)^\s*units_time\s*=\s*(\S+)', nml)[1].replace('d', 'e'))
thermal = thermal_accounting(text, reference, steps, dt, unit_t, records)['statistics']
if thermal['cooling_records'] == 0 or thermal['heating_records'] == 0:
    raise ValueError('cooling/heating coverage incomplete')
for output, k in ((1, 2), (2, 4)):
    with h5py.File(a.run/f'output_{output:05d}/data_{output:05d}.h5') as h:
        g = h['particles'];ix = np.argsort(g['identity'][:]);last = particles[2*k]
        if a.sinks:
            ix = ix[g['ptypep'][:][ix] == 1]
        close(g['identity'][:][ix], last[:, 1], 'saved identities', atol=0, rtol=0)
        close(g['mass'][:][ix], last[:, 3], 'saved stellar masses', atol=0, rtol=0)
        for j in range(3):
            close(g[f'x_{j+1}'][:][ix], last[:, 4+j], 'saved endpoint positions', atol=0, rtol=0)
            close(g[f'v_{j+1}'][:][ix], last[:, 7+j], 'saved endpoint velocities', atol=0, rtol=0)
        for j, identity in enumerate(last[:, 1].astype(int)):
            birth = birth_by_id[identity]
            close(g['tpp'][:][ix][j], birth[4], 'saved birth time', atol=0, rtol=0)
            close(g['mp0'][:][ix][j], birth[6], 'saved birth mass', atol=0, rtol=0)
            close(g['indtab'][:][ix][j], times[k]-birth[4], 'saved committed SSP age', atol=1e-16, rtol=0)
subprocess.run([sys.executable, str(a.out/'stellar.py'), str(a.run/'run.log'), str(a.out/'elements.json'), '--require-elements'], check=True)
sink_result = None
if a.sinks:
    from rar_native_sinks import audit
    sink_result = audit(a.run, leaves, sinks, dt, times, total_mass, records, close)
    sink_result['source_epochs'] = sink_reports
result = {'schema': 'rar_native_stellar_endpoint_audit_v1', 'steps': steps, 'end_time': float(times[-1]),
          'timestep_sequence': dt.tolist(), 'leaves': len(reference), 'births': len(births),
          'max_accepted_depletion_fraction': float(max(cap_fraction)), 'budgets': budgets,
          'maximum_total_mass_difference': float(max(abs(total_mass-total_mass[0]))),
          'maximum_kdk_residual': max_kdk, 'particle_epochs': particle_reports, 'thermal': thermal,
          'maximum_gas_gravity_transaction_residual': max_gas_kick,
          'physical_sinks': sink_result,
          'source_revision': (a.run/'source_commit.txt').read_text().strip(),
          'inputs_sha256': {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in
                            [a.run/'run.log', a.run/'effective.nml', a.run/'inputs.sha256', a.run/'outputs.sha256']},
          'limits': 'Synthetic depletion stress on a fixed serial two-level periodic mesh. Isothermal birth removes gas internal energy. Physical sinks/AGN are checked only when requested. No calibrated SF, MPI/regrid/restart, cosmological IC, galaxy-RAR, or zero AMR self-force claim.'}
(a.out/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
print(json.dumps({k:result[k] for k in ('steps','end_time','leaves','births','maximum_total_mass_difference','maximum_kdk_residual')},indent=2))
