#!/usr/bin/env python3
"""Compare transfer tables on common log-k points; not a realized-IC test."""
import json
import os
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

base = Path('/gpfs/kjhan/RAR_QUMOND')
old = np.loadtxt(base / 'initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084/transfer.dat')
paths = [base / f'transfer_probes/lcdm_g205_z99_highk_1290351_{i}/transfer_z99.dat' for i in (0, 1)]
new = [np.loadtxt(p) for p in paths]
out = Path(os.environ['RAR_COMPARE_OUTPUT'])
names = {'cdm_density': 1, 'baryon_density': 2, 'cdm_velocity': 10, 'baryon_velocity': 11}
bands = [(0.25, 1), (1, 4), (4, 16), (16, 32), (32, 64), (64, 128), (128, 155), (155, 257), (257, 445)]

def interpolate(table, col, kval, allow_extrapolation):
    x, y = table[:, 0], table[:, col]
    if not np.isfinite(table).all() or np.any(np.diff(x) <= 0):
        raise ValueError('invalid transfer data')
    if np.min(kval) < x[0] or (not allow_extrapolation and np.max(kval) > x[-1]):
        raise ValueError('evaluation beyond permitted table range')
    lx, lk = np.log10(x), np.log10(kval)
    positive = bool(np.all(y > 0))
    values = np.log10(y) if positive else y
    result = CubicSpline(lx, values, bc_type='natural', extrapolate=False)(np.minimum(lk, lx[-1]))
    right = lk > lx[-1]
    result[right] = values[-1] + (values[-1]-values[-2])/(lx[-1]-lx[-2])*(lk[right]-lx[-1])
    return np.power(10.0, result) if positive else result

report = {'scope': 'GR CAMB transfer comparison on equally spaced log-k points, not an IC power-spectrum or solver-accuracy convergence test',
          'reference': str(paths[1]), 'sampling_100': str(paths[0]),
          'original_table': 'lcdm_g205_multifluid_l11_z99_diag_1290084/transfer.dat',
          'interpolation': 'natural cubic in log(k) and log(T) for strictly positive columns, otherwise signed T; high-k original uses last-two-node extrapolation',
          'bands': []}
for low, high in bands:
    kval = np.geomspace(low, high, 2048)
    entry = {'k_h_Mpc': [low, high], 'fields': {}}
    for name, col in names.items():
        ref = interpolate(new[1], col, kval, False)
        a = interpolate(new[0], col, kval, False)
        b = interpolate(old, col, kval, True)
        norm = np.linalg.norm(ref)
        if norm == 0 or not np.isfinite(norm):
            raise ValueError('zero or invalid reference norm')
        entry['fields'][name] = {
            'sample100_vs_200_relative_l2': float(np.linalg.norm(a-ref)/norm),
            'original_vs_200_relative_l2': float(np.linalg.norm(b-ref)/norm),
            'sample100_sign_disagreement_fraction': float(np.mean(np.sign(a) != np.sign(ref))),
            'original_sign_disagreement_fraction': float(np.mean(np.sign(b) != np.sign(ref)))}
    report['bands'].append(entry)
(out / 'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
