#!/usr/bin/env python3
"""Measure any cooling-operator energy change outside accepted thermal records."""
import argparse,collections,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('log',type=Path);p.add_argument('inventory',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
groups=collections.defaultdict(list);lines=iter(a.log.read_text().splitlines())
for line in lines:
    if 'RAR_COOLING_COMMIT_V1' not in line:continue
    t=line.split('RAR_COOLING_COMMIT_V1',1)[1].split()
    while len(t)<9:t.extend(next(lines).split())
    if len(t)!=9:raise ValueError('bad thermal record')
    v=[float(x.replace('D','E')) for x in t]
    if not all(map(math.isfinite,v)):raise ValueError('nonfinite thermal record')
    groups[round(v[3],12)].append(v)
result=[]
for r in json.loads(a.inventory.read_text())['epochs']:
    rows=groups.pop(round(r['before']['aexp'],12))
    if len(rows)!=32768 or r['before']['leaf_cells']!=32768:raise ValueError('requires uniform cell coverage')
    accepted=math.fsum(v[7] for v in rows);rhs=math.fsum(v[6] for v in rows)
    measured=r['gas_total_energy_delta_code'];residual=measured-accepted
    scale=max(abs(r[p]['gas_energy_code']) for p in ('before','after'))
    tolerance=256*math.ulp(max(scale,abs(accepted)))+1e-10*abs(accepted)
    result.append({'step':r['step'],'aexp':r['before']['aexp'],'records':len(rows),
       'effective_rhs_code':rhs,'accepted_thermal_code':accepted,'whole_operator_code':measured,
       'unaccounted_final_correction_code':residual,'tolerance':tolerance,'within_tolerance':abs(residual)<=tolerance})
if groups:raise ValueError('unmatched thermal epochs')
a.output.write_text(json.dumps({'schema':'rar_cooling_whole_operator_v1','epochs':result,
  'limits':'Effective thermal source and whole cooling call, including final floors; no cosmic gravitational or stellar energy closure.'},indent=2,allow_nan=False)+'\n')
if not result or not all(r['within_tolerance'] for r in result):raise ValueError('unexplained final cooling correction')
