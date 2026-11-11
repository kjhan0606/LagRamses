#!/usr/bin/env python3
"""Inventory conservative gas tracer and mass on stored AMR leaf cells."""
import json
import os
from pathlib import Path

import h5py
import numpy as np

snapshot = Path(os.environ['RAR_SNAPSHOT'])
out = Path(os.environ['RAR_ANALYSIS_OUTPUT'])
report = {'scope': 'single-epoch leaf-cell gas and passive-tracer inventory; not a conservation test or final galaxy', 'levels': {}}
with h5py.File(snapshot,'r') as h:
    a = float(np.asarray(h['header'].attrs['aexp']).item())
    if abs(a-0.02) > 1e-6:
        raise ValueError('unexpected snapshot epoch')
    report['scale_factor'] = a
    for lev in range(1,8):
        path = f'amr/level_{lev}/son_flag'
        if path in h and np.any(np.asarray(h[path]) != 1):
            raise ValueError(f'unexpected leaf below levelmin: {lev}')
    for lev in range(8,12):
        flags = h[f'amr/level_{lev}/son_flag']
        gas = h[f'hydro/level_{lev}/uold_1']
        tracer = h[f'hydro/level_{lev}/uold_18']
        if len(flags) != len(gas) or len(gas) != len(tracer):
            raise ValueError('AMR and hydro dataset size mismatch')
        nleaf = 0
        mass = tracer_mass = selected_mass = 0.0
        fraction_min, fraction_max = np.inf, -np.inf
        below, above = 0, 0
        volume = 2.0**(-3*lev)
        for lo in range(0,len(flags),500_000):
            hi = min(lo+500_000,len(flags))
            f = np.asarray(flags[lo:hi])
            if not np.isin(f,[0,1]).all():
                raise ValueError('unexpected son flag')
            leaf = f==0
            d = np.asarray(gas[lo:hi])[leaf]
            q = np.asarray(tracer[lo:hi])[leaf]
            if not np.isfinite(d).all() or not np.isfinite(q).all() or np.any(d<=0):
                raise ValueError('invalid leaf density or tracer')
            fraction = q/d
            if len(d):
                nleaf += len(d)
                mass += float(np.sum(d,dtype=np.float64))*volume
                tracer_mass += float(np.sum(q,dtype=np.float64))*volume
                selected_mass += float(np.sum(d[fraction>0.5],dtype=np.float64))*volume
                fraction_min = min(fraction_min,float(fraction.min()))
                fraction_max = max(fraction_max,float(fraction.max()))
                below += int(np.count_nonzero(fraction < -1e-10))
                above += int(np.count_nonzero(fraction > 1+1e-10))
        if nleaf==0:
            raise ValueError('missing leaf cells')
        report['levels'][str(lev)] = {'leaf_cells':nleaf,'gas_mass_code':mass,
            'tracer_mass_code':tracer_mass,'gas_mass_fraction_gt_0_5_code':selected_mass,
            'tracer_fraction_min':fraction_min,'tracer_fraction_max':fraction_max,
            'fraction_below_minus_1e_10_count':below,'fraction_above_1_plus_1e_10_count':above}
report['total_leaf_gas_mass_code'] = sum(x['gas_mass_code'] for x in report['levels'].values())
report['total_leaf_tracer_mass_code'] = sum(x['tracer_mass_code'] for x in report['levels'].values())
report['level11_fraction_of_total_tracer_mass'] = report['levels']['11']['tracer_mass_code']/report['total_leaf_tracer_mass_code']
report['limits'] = 'Tracer mass is an advected passive diagnostic, not an independently measured target-halo gas mass. One epoch cannot establish conservation, mixing accuracy, physical metallicity, or the later galaxy.'
baseline_path = Path(os.environ['RAR_BASELINE_INVENTORY'])
baseline = json.loads(baseline_path.read_text())
if abs(baseline['scale_factor']-0.0101) > 1e-6:
    raise ValueError('wrong baseline inventory epoch')
comparison = {'early_a':baseline['scale_factor'],'late_a':a,
              'early_gas_mass':baseline['total_leaf_gas_mass_code'],
              'late_gas_mass':report['total_leaf_gas_mass_code'],
              'early_tracer_mass':baseline['total_leaf_tracer_mass_code'],
              'late_tracer_mass':report['total_leaf_tracer_mass_code'],
              'gas_mass_relative_change':report['total_leaf_gas_mass_code']/baseline['total_leaf_gas_mass_code']-1.0,
              'tracer_mass_relative_change':report['total_leaf_tracer_mass_code']/baseline['total_leaf_tracer_mass_code']-1.0,
              'limits':'Separate runs from the same input IC have different output/timestep schedules; changes do not isolate conservation from time-stepping, AMR and floating-point differences.'}
(out/'gas_tracer_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
(out/'gas_tracer_inventory.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
