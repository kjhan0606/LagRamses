#!/usr/bin/env python3
"""Check physical mass across selected operators, excluding sink-cloud masses."""
import argparse,json,math
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('log',type=Path);p.add_argument('output',type=Path);p.add_argument('--before-phase',default='before_agn');p.add_argument('--after-phase',default='after_agn');p.add_argument('--activity',choices=['auto','any_energy'],default='auto');a=p.parse_args()
s=a.log.read_text()
if 'Run completed' not in s:raise ValueError('incomplete evolution')
for bad in ('failed to converge','NaN/Inf detected','segmentation fault','FATAL:'):
    if bad.lower() in s.lower():raise ValueError('failed runtime physics gate')
records={};lines=iter(s.splitlines())
for line in lines:
    marker='RAR_GLOBAL_INVENTORY_V1'
    if marker not in line:continue
    t=line.split(marker,1)[1].split()
    while len(t)<4:t.extend(next(lines).split())
    phase=t[0];step=int(t[1]);epoch=float(t[2].replace('D','E'));nvar=int(t[3])
    if nvar<5:raise ValueError('bad header')
    while len(t)<nvar+16:t.extend(next(lines).split())
    if len(t)!=nvar+16:raise ValueError('bad record length')
    if phase not in (a.before_phase,a.after_phase):continue
    x=[float(v.replace('D','E')) for v in t[4:]]
    if not all(map(math.isfinite,x)) or epoch<=0:raise ValueError('nonfinite inventory')
    if x[nvar+8] or x[nvar+9]:raise ValueError('invalid or duplicated physical particle list')
    if x[nvar]<=0 or min(x[0],x[nvar+2],x[nvar+4],x[nvar+7])<0:raise ValueError('invalid mass/count')
    mass=math.fsum(x[i] for i in (0,nvar+2,nvar+4,nvar+7,nvar+10))
    key=(step,phase)
    if key in records:raise ValueError('duplicate global inventory')
    records[key]={'aexp':epoch,'gas_mass':x[0],'star_mass':x[nvar+4],'sink_mass':x[nvar+7],
                  'cloud_mass_excluded':x[nvar+6],'physical_mass':mass,'gas_energy_code':x[4],
                  'leaf_cells':x[nvar],'nvar':nvar,'raw_totals':x}
if not records:raise ValueError('no inventory')
result=[];active=False
for step in sorted({k[0] for k in records}):
    before=records[(step,a.before_phase)];after=records[(step,a.after_phase)]
    if before['aexp']!=after['aexp'] or before['nvar']!=after['nvar']:raise ValueError('phase mismatch')
    delta=after['physical_mass']-before['physical_mass']
    if not math.isclose(before['physical_mass'],after['physical_mass'],rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError('operator global physical mass mismatch')
    energy=after['gas_energy_code']-before['gas_energy_code']
    if a.activity=='any_energy':active|=energy!=0
    elif a.before_phase=='before_agn':active|=energy>0
    else:active|=before['star_mass']>after['star_mass'] and after['gas_mass']>before['gas_mass']
    result.append({'step':step,'before':before,'after':after,'physical_mass_delta':delta,
                   'gas_total_energy_delta_code':energy})
if not active:raise ValueError('no active operator exchange')
a.output.write_text(json.dumps({'schema':'rar_global_operator_inventory_v1','epochs':result,'before_phase':a.before_phase,'after_phase':a.after_phase,
    'status':'operator_mass_checked_energy_and_full_step_budgets_require_interpretation',
    'limits':'Gas passive fields measured raw; stellar nuclear production, sink accretion/radiation and global energy closure not certified.'},indent=2,allow_nan=False)+'\n')
