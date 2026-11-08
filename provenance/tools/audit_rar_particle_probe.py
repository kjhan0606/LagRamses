#!/usr/bin/env python3
"""Measure the initial native AMR particle self force and symmetric pair force."""
import json,subprocess,sys
from pathlib import Path
import numpy as np
from rar_native_particles import epoch

root=Path(sys.argv[1]);results={};fields={}
for gravity in ('newtonian','qumond'):
 for fixture in ('a','b','pair'):
  name=f'{gravity}_{fixture}';case=root/name
  leaves=np.loadtxt(case/'composite_epoch_00000_leaf.tsv')
  rows,rho,report=epoch(case/'particles_epoch_00000.tsv',leaves)
  if np.max(abs(leaves[:,6]-1-rho))>1e-13:raise ValueError('uniform gas plus cloud source mismatch')
  if np.any(rows[:,7:10]!=0):raise ValueError('initial-force probe acquired a velocity')
  subprocess.run([sys.executable,str(root/'verify.py'),str(root),str(root/f'verify_{name}.json'),
   '--case',name,'--leaf-file','composite_epoch_00000_leaf.tsv','--snapshots','0','--evolving','--particles'],check=True)
  acceleration=rows[:,10:13];mass=rows[:,3];forces=mass[:,None]*acceleration
  report.update(acceleration=acceleration.tolist(),mass_weighted_net_force=forces.sum(axis=0).tolist(),
   net_force_fraction=float(np.linalg.norm(forces.sum(axis=0))/max(np.linalg.norm(forces,axis=1).sum(),1e-300)))
  if fixture!='pair':
   # G=1 in these physical, noncosmological code units. This dimensionless
   # measure uses the declared cloud scale, not a fitted galaxy acceleration.
   report['self_acceleration_over_Gm_s2']=float(np.linalg.norm(acceleration[0])/(mass[0]/report['cloud_width']**2))
  results[name]=report;fields[name]=leaves
for name,field in fields.items():
 if not np.array_equal(field[:,:6],fields['newtonian_a'][:,:6]):raise ValueError('probe topologies differ')
pair=fields['newtonian_pair'][:,10:13]
superposed=fields['newtonian_a'][:,10:13]+fields['newtonian_b'][:,10:13]
residual=float(np.max(abs(pair-superposed))/max(np.max(abs(pair)),1e-300))
if residual>1e-8:raise ValueError('Newtonian particle source superposition failed')
(root/'result.json').write_text(json.dumps({'schema':'rar_native_particle_initial_force_v1',
 'cases':results,'newtonian_field_superposition_relative_max_error':residual,
 'scope':'Time-zero native fixed AMR forces for one/two synthetic clouds in uniform gas; explicit one-step executions with no snapshots. Evolved fields are retained but not used in this initial-force measurement.',
 'limits':'A symmetric pair and two interface positions are a bounded measurement, not a spatial convergence study. Matched maps do not imply zero self force. QUMOND is nonlinear and its single-particle fields are not superposed. No calibrated galaxy softening or long-time momentum-conservation claim.'},indent=2,allow_nan=False)+'\n')
