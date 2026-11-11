#!/usr/bin/env python3
"""Bound distances of level-11 particles outside stored level-11 boxes."""
import json
import os
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial import cKDTree

snapshot = Path(os.environ['RAR_SNAPSHOT'])
out = Path(os.environ['RAR_ANALYSIS_OUTPUT'])
mass_fine = (1.0 - 0.0493017/0.315)*2.0**-33
nside = 1024

def keys(x):
    return (x[:, 0]*nside+x[:, 1])*nside+x[:, 2]

with h5py.File(snapshot, 'r') as h:
    a = float(np.asarray(h['header'].attrs['aexp']).item())
    if abs(a-0.0101) > 1e-6:
        raise ValueError('unexpected epoch')
    g = h['amr/level_11']
    centers = np.column_stack([g[f'xg_{i}'][:] for i in (1,2,3)])
    lattice = centers*nside-0.5
    if np.max(np.abs(lattice-np.rint(lattice))) > 1e-7:
        raise ValueError('unexpected grid coordinate convention')
    grid_keys = np.sort(keys(np.rint(lattice).astype(np.int64) % nside))
    if len(np.unique(grid_keys)) != len(grid_keys):
        raise ValueError('duplicate grid boxes')
    p = h['particles']
    outside = []
    for lo in range(0,len(p['mass']),500_000):
        hi = min(lo+500_000,len(p['mass']))
        mass = p['mass'][lo:hi]
        lev = p['levelp'][lo:hi]
        selected = np.isclose(mass,mass_fine,rtol=1e-5,atol=0) & (lev==11)
        if not selected.any():
            continue
        xyz = np.column_stack([p[f'x_{i}'][lo:hi] for i in (1,2,3)])[selected]
        if not np.isfinite(xyz).all() or np.any((xyz<0)|(xyz>=1)):
            raise ValueError('invalid periodic position')
        particle_keys = keys(np.floor(xyz*nside).astype(np.int64))
        ix = np.searchsorted(grid_keys,particle_keys)
        valid = ix < len(grid_keys)
        inside = np.zeros(len(ix), dtype=bool)
        inside[valid] = grid_keys[ix[valid]] == particle_keys[valid]
        outside.append(xyz[~inside])
outside = np.concatenate(outside, axis=0)
if len(outside) != 679:
    raise ValueError(f'expected 679 mismatches, found {len(outside)}')
tree = cKDTree(centers, boxsize=1.0)
_, candidates = tree.query(outside,k=8,workers=2)
delta = np.abs(outside[:,None,:]-centers[candidates])
delta = np.minimum(delta,1.0-delta)
distance = np.sqrt(np.sum(np.maximum(delta-0.5/nside,0.0)**2,axis=2)).min(axis=1)*2048.0
thresholds = [0.01,0.1,0.5,1.0,2.0]
report = {'scale_factor':a, 'mismatched_level11_particles':len(outside),
          'distance_units':'level-11 cell widths (25/2048 comoving Mpc/h)',
          'quantiles_upper_bound':{str(q):float(np.quantile(distance,q)) for q in (0,0.5,0.9,0.99,1)},
          'within_upper_bound_counts':{str(t):int(np.count_nonzero(distance<=t)) for t in thresholds},
          'method':'periodic distance to closest among eight nearest level-11 grid centres, subtracting box half-width componentwise',
          'limits':'Distances are upper bounds on distance to any level-11 box, not exact nearest-box distances. They cannot diagnose cause, rule out stale levels, or determine force error.'}
(out/'mismatch_distance.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
