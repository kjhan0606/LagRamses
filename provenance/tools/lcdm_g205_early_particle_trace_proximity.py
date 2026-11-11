#!/usr/bin/env python3
"""Measure periodic proximity of zoom particles to parent trace-cell centres."""
import json
import os
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

snapshot = Path(os.environ['RAR_SNAPSHOT'])
trace = Path(os.environ['RAR_PARENT_TRACE'])
out = Path(os.environ['RAR_ANALYSIS_OUTPUT'])
with np.load(trace) as archive:
    cells = np.asarray(archive['g205_r3_cell_xyz'], dtype=np.int64)
if cells.shape != (4536, 3) or len(np.unique(cells, axis=0)) != 4536:
    raise ValueError('unexpected Group 205 trace cells')
shift = np.array([46, -54, 56])
centres = ((cells + shift) % 256 + 0.5) / 256.0
tree = cKDTree(centres, boxsize=1.0)
thresholds = np.array([0.5, 1.0, 2.0, 4.0])
expected_mass = (1.0 - 0.0493017/0.315)*2.0**-33
hist = {str(level): {'total': 0, 'within': np.zeros(len(thresholds), dtype=np.int64),
                     'd_min': np.inf, 'd_max': 0.0} for level in (10, 11)}
with h5py.File(snapshot, 'r') as h:
    a = float(np.asarray(h['header'].attrs['aexp']).item())
    if abs(a-0.0101) > 1e-6:
        raise ValueError('wrong snapshot epoch')
    p = h['particles']
    n = len(p['mass'])
    for start in range(0, n, 250_000):
        end = min(n, start+250_000)
        mass = np.asarray(p['mass'][start:end])
        lev = np.asarray(p['levelp'][start:end])
        fine = np.isclose(mass, expected_mass, rtol=1e-5, atol=0)
        if not fine.any():
            continue
        xyz = np.column_stack([p[f'x_{axis}'][start:end] for axis in (1,2,3)])
        if not np.isfinite(xyz).all() or np.any((xyz < 0) | (xyz >= 1)):
            raise ValueError('invalid periodic particle coordinates')
        for level in (10, 11):
            mask = fine & (lev == level)
            if not mask.any():
                continue
            distance = tree.query(xyz[mask], workers=4)[0] * 256.0
            group = hist[str(level)]
            group['total'] += len(distance)
            group['within'] += np.sum(distance[:, None] <= thresholds, axis=0)
            group['d_min'] = min(group['d_min'], float(distance.min()))
            group['d_max'] = max(group['d_max'], float(distance.max()))
expected = {'10': 42860, '11': 6229836}
if any(hist[k]['total'] != v for k,v in expected.items()):
    raise ValueError('particle counts differ from independent snapshot inventory')
result = {'scale_factor': a, 'parent_trace': 'Group 205 union within 3 R200c and HOP members',
          'shift_base_cells_xyz': shift.tolist(), 'trace_cells': len(cells),
          'distance_unit': 'base-cell width = 25/256 comoving Mpc/h',
          'thresholds_base_cells': thresholds.tolist(), 'levels': {},
          'limits': 'Distances compare early evolved zoom positions to shifted parent Lagrangian cell centres. They are an Eulerian proximity proxy, not particle-ID membership, halo contamination, AMR-leaf coverage, or a late-time resolution measurement.'}
for level, group in hist.items():
    result['levels'][level] = {'total_finest_mass_particles': group['total'],
        'within_counts': group['within'].tolist(),
        'within_fractions': (group['within']/group['total']).tolist(),
        'minimum_distance_base_cells': group['d_min'], 'maximum_distance_base_cells': group['d_max']}
(out/'proximity.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
