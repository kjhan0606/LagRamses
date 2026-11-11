#!/usr/bin/env python3
"""Check finest-mass particles against actual finest-level AMR grid boxes."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

snapshot = Path(os.environ['RAR_SNAPSHOT'])
out = Path(os.environ['RAR_ANALYSIS_OUTPUT'])
mass_fine = (1.0 - 0.0493017/0.315)*2.0**-33
nside = 1024  # level-11 grids are boxes of two level-11 cells per axis

def keys(xyz):
    return (xyz[:, 0]*nside + xyz[:, 1])*nside + xyz[:, 2]

with h5py.File(snapshot, 'r') as h:
    a = float(np.asarray(h['header'].attrs['aexp']).item())
    if abs(a-0.0101) > 1e-6:
        raise ValueError('wrong snapshot epoch')
    g = h['amr/level_11']
    centres = np.column_stack([g[f'xg_{i}'][:] for i in (1,2,3)])
    lattice = centres*nside - 0.5
    if np.max(np.abs(lattice-np.rint(lattice))) > 1e-7:
        raise ValueError('unexpected AMR grid-centre convention')
    grid_xyz = np.rint(lattice).astype(np.int64) % nside
    grid_keys = np.sort(keys(grid_xyz))
    if len(np.unique(grid_keys)) != len(grid_keys):
        raise ValueError('duplicate level-11 grid boxes')
    p = h['particles']
    total = len(p['mass'])
    counts = {'level10_inside': 0, 'level10_outside': 0,
              'level11_inside': 0, 'level11_outside': 0, 'other': 0}
    for lo in range(0, total, 500_000):
        hi = min(lo+500_000, total)
        mass = np.asarray(p['mass'][lo:hi])
        level = np.asarray(p['levelp'][lo:hi])
        fine = np.isclose(mass, mass_fine, rtol=1e-5, atol=0)
        if not fine.any():
            continue
        xyz = np.column_stack([p[f'x_{i}'][lo:hi] for i in (1,2,3)])[fine]
        if not np.isfinite(xyz).all() or np.any((xyz < 0) | (xyz >= 1)):
            raise ValueError('invalid particle coordinates')
        level = level[fine]
        coord = np.floor(xyz*nside).astype(np.int64)
        pkey = keys(coord)
        idx = np.searchsorted(grid_keys, pkey)
        valid = idx < len(grid_keys)
        inside = np.zeros(len(idx), dtype=bool)
        inside[valid] = grid_keys[idx[valid]] == pkey[valid]
        for lev in (10, 11):
            counts[f'level{lev}_inside'] += int(np.count_nonzero((level==lev)&inside))
            counts[f'level{lev}_outside'] += int(np.count_nonzero((level==lev)&~inside))
        counts['other'] += int(np.count_nonzero((level!=10)&(level!=11)))
if sum(counts.values()) != 6272696:
    raise ValueError('finest-mass count mismatch')
summary = {'scale_factor': a, 'finest_grid_boxes': len(grid_keys),
           'finest_mass_particles': sum(counts.values()), 'counts': counts,
           'inside_fraction': (counts['level10_inside']+counts['level11_inside'])/sum(counts.values()),
           'method': 'independently map periodic positions into level-11 grid boxes from AMR grid-centre coordinates; level-11 is maximum level',
           'limits': 'One early output only. Grid-box containment is geometric, not CIC-stencil or force accuracy, future halo membership, tracer coverage, or convergence.'}
(out/'leaf_coverage.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
