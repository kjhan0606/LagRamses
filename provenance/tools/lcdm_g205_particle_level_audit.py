#!/usr/bin/env python3
"""Count particle integration levels by IC mass from a short RAMSES output."""
import json
import sys
from pathlib import Path

import h5py
import numpy as np

snapshot, destination = map(Path, sys.argv[1:])
expected_mass = (1.0 - 0.0493017 / 0.315) * 2.0**-33
counts = {}
with h5py.File(snapshot, 'r') as h:
    scale = float(np.asarray(h['header'].attrs['aexp']).item())
    if not 0.01009 <= scale <= 0.01011:
        raise ValueError(f'unexpected snapshot epoch {scale}')
    p = h['particles']
    if not all(x in p for x in ('mass', 'levelp')):
        raise ValueError('missing particle masses or levels')
    n = len(p['mass'])
    if n != len(p['levelp']):
        raise ValueError('inconsistent particle array lengths')
    for start in range(0, n, 500_000):
        mass = np.asarray(p['mass'][start:start+500_000])
        lev = np.asarray(p['levelp'][start:start+500_000])
        if not np.isfinite(mass).all() or np.any(mass <= 0) or np.any((lev < 1) | (lev > 11)):
            raise ValueError('invalid mass or particle level')
        fine = np.isclose(mass, expected_mass, rtol=1e-5, atol=0.0)
        for cls, mask in [('fine_mass', fine), ('other_mass', ~fine)]:
            hist = np.bincount(lev[mask].astype(np.int64), minlength=12)
            for level in range(1, 12):
                key = f'{cls}_level_{level}'
                counts[key] = counts.get(key, 0) + int(hist[level])
fine_count = sum(counts[f'fine_mass_level_{i}'] for i in range(1, 12))
if fine_count == 0 or sum(counts.values()) != n:
    raise ValueError('missing finest-mass particles or incomplete histogram')
summary = {'scale_factor': scale, 'particles_total': n,
           'finest_cdm_code_mass_expected': expected_mass,
           'finest_mass_particles': fine_count,
           'finest_mass_by_particle_level': {str(i): counts[f'fine_mass_level_{i}'] for i in range(1,12)},
           'finest_mass_fraction_level11': counts['fine_mass_level_11']/fine_count,
           'counts': counts,
           'limits': 'levelp is the saved particle integration/tree level; this is not an independent AMR leaf-location or force-accuracy measurement, not a final halo or contamination count.'}
destination.write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
