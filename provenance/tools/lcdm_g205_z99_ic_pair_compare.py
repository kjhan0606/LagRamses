#!/usr/bin/env python3
"""Compare same-coordinate GRAFIC fields; does not claim spectral convergence."""
import json
import os
import struct
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

base = Path('/gpfs/kjhan/RAR_QUMOND/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084')
other = Path(os.environ['RAR_HIGHK_IC'])
out = Path(os.environ['RAR_COMPARE_OUTPUT'])

def read(path):
    with FortranFile(path, 'r') as f:
        header = f.read_record(np.uint8).tobytes()
        if len(header) != 48:
            raise ValueError(f'bad GRAFIC header: {path}')
        nx, ny, nz, *_ = struct.unpack('<3i9f', header)
        if min(nx, ny, nz) < 1 or max(nx, ny, nz) > 1024:
            raise ValueError('unexpected dimensions')
        a = np.empty((nz, ny, nx), dtype=np.float32)
        for z in range(nz):
            a[z] = f.read_reals(np.float32).reshape(ny, nx)
    if not np.isfinite(a).all():
        raise ValueError('nonfinite input')
    return header, a

fields = ['ic_deltab', 'ic_velbx', 'ic_velby', 'ic_velbz',
          'ic_poscx', 'ic_poscy', 'ic_poscz', 'ic_velcx', 'ic_velcy', 'ic_velcz']
report = {'baseline': str(base), 'highk': str(other),
          'scope': 'same-coordinate raw field comparison; not AMR leaf, Fourier spectrum, or evolved galaxy', 'levels': {}}
for level in (8, 9, 10, 11):
    ad = base / 'ics' / f'level_{level:03d}'
    bd = other / 'ics' / f'level_{level:03d}'
    ha, maska = read(ad / 'ic_refmap')
    hb, maskb = read(bd / 'ic_refmap')
    if ha != hb or not np.array_equal(maska, maskb):
        raise ValueError(f'different headers or masks at level {level}')
    mask = maska > 0.5
    group = {'refmap_one_count': int(mask.sum()), 'fields': {}}
    del maska, maskb
    for field in fields:
        ha, a = read(ad / field)
        hb, b = read(bd / field)
        if ha != hb:
            raise ValueError(f'different header: {level} {field}')
        metrics = {}
        for name, selected in [('all_cells', None), ('refmap_one_cells', mask)]:
            x = a if selected is None else a[selected]
            y = b if selected is None else b[selected]
            x = x.astype(np.float64).ravel()
            y = y.astype(np.float64).ravel()
            denominator = np.linalg.norm(x)
            if not len(x) or denominator == 0:
                raise ValueError(f'empty or zero reference: {level} {field} {name}')
            metrics[name] = {'n': int(len(x)),
                             'relative_l2': float(np.linalg.norm(y-x)/denominator),
                             'baseline_rms': float(np.sqrt(np.mean(x*x))),
                             'highk_rms': float(np.sqrt(np.mean(y*y)))}
        group['fields'][field] = metrics
    report['levels'][str(level)] = group
(out / 'field_comparison.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
