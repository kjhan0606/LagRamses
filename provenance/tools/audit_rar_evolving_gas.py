#!/usr/bin/env python3
"""Read actual moving-gas epochs, gravity transactions and saved endpoints."""
import argparse,hashlib,json,re,subprocess,sys
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('out',type=Path);p.add_argument('verifier',type=Path);p.add_argument('--cooling',action='store_true');p.add_argument('--thermal-accuracy',action='store_true');p.add_argument('--thermal-energy-tolerance',type=float);a=p.parse_args()
def records(text,tag,n):
 result=[];lines=iter(text.splitlines())
 for line in lines:
  if tag not in line:continue
  tokens=line.split(tag,1)[1].split()
  while len(tokens)<n:tokens.extend(next(lines).split())
  if len(tokens)!=n:raise ValueError('bad record '+tag)
  row=np.array([float(x.replace('D','E')) for x in tokens])
  if not np.isfinite(row).all():raise ValueError('nonfinite '+tag)
  result.append(row)
 return np.array(result)
def snap_state(path,reference):
 lookup={(int(row[1]),*np.rint(row[2:5]*2**(int(row[1])+1)).astype(int)):i for i,row in enumerate(reference)}
 result=np.empty((len(reference),5));seen=np.zeros(len(reference),bool)
 offsets=np.array([[((c>>d)&1)-.5 for d in range(3)] for c in range(8)])
 with h5py.File(path) as h:
  for lev in sorted(set(reference[:,1].astype(int))):
   prefix=f'amr/level_{lev}';grids=np.stack([h[f'{prefix}/xg_{d}'][:] for d in (1,2,3)],axis=1)
   centers=np.repeat(grids,8,axis=0)+np.tile(offsets,(len(grids),1))*2.**(-int(lev));mask=h[f'{prefix}/son_flag'][:]==0
   indices=np.array([lookup[(int(lev),*np.rint(c*2**(int(lev)+1)).astype(int))] for c in centers[mask]])
   result[indices]=np.stack([h[f'hydro/level_{lev}/uold_{d}'][:][mask] for d in range(1,6)],axis=1);seen[indices]=True
  time=float(np.asarray(h['header'].attrs['time']).item());dt=np.asarray(h['header'].attrs['dtnew'])
 if not seen.all() or not np.isfinite(result).all():raise ValueError('bad snapshot leaf join')
 return result,time,dt
summaries={};initials={};finals={};reference=None;all_inputs={}
cases=[('newtonian',6,.0002),('qumond',6,.0002),('qumond_half',12,.0001)]
if a.cooling:cases.extend([('qumond_off',6,.0002),('qumond_heating',6,.0002)])
if a.thermal_accuracy:cases.append(('qumond_accuracy',6,.0002))
for mode,steps,dt in cases:
 case=a.run/mode;text=(case/'run.log').read_text()
 starts=records(text,'RAR_QUMOND_EVOLVE_BEGIN_V1',5);ends=records(text,'RAR_QUMOND_EVOLVE_END_V1',5)
 if len(starts)!=steps+1 or len(ends)!=steps or not np.array_equal(ends[:,0],np.arange(steps)):raise ValueError('incomplete step coverage')
 if not np.array_equal(starts[:,0],np.arange(steps+1)):raise ValueError('duplicate/missing begin')
 if not np.allclose(starts[:,1],np.arange(steps+1)*dt,rtol=0,atol=1e-16) or not np.all(starts[:,2]==dt):raise ValueError('common timestep/clock mismatch')
 if not np.array_equal(starts[:,1],starts[:,3]):raise ValueError('start field source time mismatch')
 if not np.allclose(ends[:,1],np.arange(steps)*dt,rtol=0,atol=1e-16) or not np.allclose(ends[:,2],np.arange(1,steps+1)*dt,rtol=0,atol=1e-16) or not np.all(ends[:,3]==dt):raise ValueError('end clock mismatch')
 paths=sorted(case.glob('composite_epoch_*_leaf.tsv'))
 if len(paths)!=steps+1:raise ValueError('missing/duplicate field epoch')
 states=[np.loadtxt(f) for f in paths];first=states[0];volume=first[:,5]**3;n=len(first)
 if reference is None:reference=first
 if not np.array_equal(reference[:,:6],first[:,:6]):raise ValueError('control topology mismatch')
 for k,(path,state) in enumerate(zip(paths,states)):
  header=path.read_text().splitlines()[1].split();time=float(header[3])
  if abs(time-k*dt)>1e-16 or not np.array_equal(first[:,:6],state[:,:6]):raise ValueError('source time/topology changed')
  if not np.isfinite(state).all() or np.any(state[:,6]<=0):raise ValueError('bad evolving source')
  all_inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
 expected=1+.1*np.cos(2*np.pi*first[:,2])+.05*np.sin(2*np.pi*first[:,3])
 if np.max(abs(expected-first[:,6]))>1e-14:raise ValueError('smooth physical IC source mismatch')
 mass=np.array([np.sum(volume*v[:,6]) for v in states]);mass_error=float(np.max(abs(mass/mass[0]-1)))
 if mass_error>1e-12:raise ValueError('gas mass does not close across AMR fluxes')
 norm=lambda v:np.sqrt(np.sum(volume*v*v))
 evolution=float(norm(states[-1][:,6]-first[:,6])/norm(first[:,6]))
 if evolution<=1e-8:raise ValueError('gas density did not evolve measurably')
 thermal=None;clip_error=0.;clipping_total=0.
 nml=(case/'effective.nml').read_text()
 def param(name):return float(re.search(r'(?im)^\s*'+name+r'\s*=\s*([^\n]+)',nml)[1].replace('d','e').replace('D','E'))
 if a.cooling and mode!='qumond_off':
  from rar_evolving_thermal import thermal_accounting
  thermal=thermal_accounting(text,first,steps,dt,param('units_time'),records)
  if mode=='qumond_heating' and thermal['statistics']['heating_records']==0:raise ValueError('heating case has no heating')
  if mode!='qumond_heating' and thermal['statistics']['cooling_records']==0:raise ValueError('cooling case has no cooling')
 elif 'RAR_COOLING_COMMIT_V2' in text:raise ValueError('adiabatic case applied a thermal source')
 max_momentum_error=0.;max_energy_error=0.;max_phase_join=0.
 lookup={int(cell):i for i,cell in enumerate(first[:,0])}
 for k in range(steps):
  path=case/f'gravity_kicks_{k:05d}.tsv';rows=np.loadtxt(path)
  all_inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
  if rows.shape!=(2*n,17) or not np.isfinite(rows).all():raise ValueError('incomplete gravity transaction coverage')
  phases=[]
  for phase in (0,1):
   row=rows[rows[:,0]==phase]
   if len(row)!=n or len(np.unique(row[:,1]))!=n:raise ValueError('duplicate/missing leaf kick')
   ix=np.array([lookup[int(c)] for c in row[:,1]]);row=row[np.argsort(ix)];phases.append(row)
   old=states[k];new=states[k+1];source=old if phase==0 else new
   if not np.array_equal(row[:,2],first[:,1]) or not np.allclose(row[:,3],source[:,6],rtol=1e-14,atol=0) or not np.allclose(row[:,4],new[:,6],rtol=1e-14,atol=0):raise ValueError('kick density source/receiver epoch mismatch')
   if not np.all(row[:,5]==.5*dt) or not np.array_equal(row[:,6:9],source[:,10:13]):raise ValueError('kick timestep/force epoch mismatch')
   expected_dp=row[:,3,None]*row[:,5,None]*row[:,6:9];actual_dp=row[:,13:16]-row[:,9:12]
   expected_de=(np.sum(row[:,13:16]**2,axis=1)-np.sum(row[:,9:12]**2,axis=1))/(2*row[:,4]);actual_de=row[:,16]-row[:,12]
   mp_error=float(np.max(abs(actual_dp-expected_dp)));en_error=float(np.max(abs(actual_de-expected_de)))
   mp_scale=max(float(np.max(abs(row[:,9:12]))),float(np.max(abs(row[:,13:16]))),1e-300)
   en_scale=max(float(np.max(abs(row[:,12]))),float(np.max(abs(row[:,16]))),1e-300)
   if mp_error>64*np.finfo(float).eps*mp_scale or en_error>64*np.finfo(float).eps*en_scale:raise ValueError('gravity source accounting mismatch')
   max_momentum_error=max(max_momentum_error,mp_error);max_energy_error=max(max_energy_error,en_error)
  energy_join=phases[1][:,12]-phases[0][:,16]
  if thermal is not None:
   before=phases[0];rho=before[:,4];internal=before[:,16]-np.sum(before[:,13:16]**2,axis=1)/(2*rho)
   # Constants are those of selected patch/cuRamses/cooling_module.f90, used by units.f90.
   scale_t2=(1.6600000e-24/1.3806200e-16)*(param('units_length')/param('units_time'))**2
   temperature=(param('gamma')-1)*internal/rho*scale_t2
   clip=rho/(param('gamma')-1)/scale_t2*(np.clip(temperature,1e-2,param('T2max'))-temperature)
   reconstruction=thermal['accepted'][k]-thermal['effective'][k]-clip
   error=float(np.max(abs(reconstruction)));clip_error=max(clip_error,error);clipping_total+=float(np.sum(clip*volume))
   if error>128*np.finfo(float).eps*max(float(np.max(abs(before[:,16]))),1e-300):raise ValueError('unexplained thermal/clipping source change')
   energy_join-=thermal['accepted'][k]
  join=max(float(np.max(abs(phases[0][:,13:16]-phases[1][:,9:12]))),float(np.max(abs(energy_join))))
  if join>1e-13:raise ValueError('gas changed between old half-source and endpoint half-source')
  max_phase_join=max(max_phase_join,join)
 # The final update_time stops before another transport: a header-only kick file is expected.
 unused=case/f'gravity_kicks_{steps:05d}.tsv'
 if any(line.strip() and not line.startswith('#') for line in unused.read_text().splitlines()):raise ValueError('unused terminal step applied a kick')
 snapshots=[]
 for output,epoch in ((1,steps//2),(2,steps)):
  snapshot=case/f'output_{output:05d}'/f'data_{output:05d}.h5'
  gas,time,dts=snap_state(snapshot,first)
  if abs(time-epoch*dt)>1e-16 or not np.all(dts[3:5]==dt):raise ValueError('saved hierarchy time/dt mismatch')
  verification=a.out/f'verify_{mode}_{output}.json'
  subprocess.run([sys.executable,str(a.verifier),str(a.run),str(verification),'--case',mode,'--snapshots','1','--leaf-file',paths[epoch].name,'--snapshot-glob',f'output_{output:05d}/data_{output:05d}.h5','--evolving'],check=True)
  snapshots.append({'time':time,'verification':str(verification)})
  if output==2:finals[mode]=gas
 initials[mode]=first[:,6]
 summaries[mode]={'steps':steps,'common_dt':dt,'endpoint':float(ends[-1,2]),'leaves':n,'density_relative_change':evolution,'mass_relative_residual':mass_error,'maximum_momentum_source_residual':max_momentum_error,'maximum_energy_source_residual':max_energy_error,'maximum_half_phase_state_join':max_phase_join,'snapshots':snapshots}
 if thermal is not None:
  summaries[mode]['thermal']=thermal['statistics'];summaries[mode]['thermal']['clipping_energy_code']=clipping_total;summaries[mode]['thermal']['max_rhs_clipping_reconstruction_residual']=clip_error
if not np.array_equal(initials['newtonian'],initials['qumond']) or not np.array_equal(initials['qumond'],initials['qumond_half']):raise ValueError('unmatched initial gas')
volume=reference[:,5]**3
relative=lambda x,y:float(np.sqrt(np.sum(volume[:,None]*(x-y)**2)/np.sum(volume[:,None]*y*y)))
halved={};contrasts={}
for field,cols in (('density',slice(0,1)),('momentum',slice(1,4)),('energy',slice(4,5))):
 halved[field]=relative(finals['qumond_half'][:,cols],finals['qumond'][:,cols]);contrasts[field]=relative(finals['newtonian'][:,cols],finals['qumond'][:,cols])
if halved['density']>=summaries['qumond']['density_relative_change']:raise ValueError('timestep replay differs more than actual density evolution')
cooling_effect={}
if a.cooling:
 for field,cols in (('density',slice(0,1)),('momentum',slice(1,4)),('energy',slice(4,5))):cooling_effect[field]=relative(finals['qumond'][:,cols],finals['qumond_off'][:,cols])
thermal_precision={}
if a.thermal_accuracy:
 for field,cols in (('density',slice(0,1)),('momentum',slice(1,4)),('energy',slice(4,5))):thermal_precision[field]=relative(finals['qumond_accuracy'][:,cols],finals['qumond'][:,cols])
 if a.thermal_energy_tolerance is None:raise ValueError('precision replay needs a declared energy tolerance')
 if max(halved['energy'],thermal_precision['energy'])>a.thermal_energy_tolerance:raise ValueError(f'thermal precision unresolved: dt={halved["energy"]} substep={thermal_precision["energy"]}')
(a.out/'result.json').write_text(json.dumps({'schema':'rar_evolving_native_gas_v1','cases':summaries,'halved_timestep_relative_state_difference':halved,'newtonian_vs_qumond_relative_state_difference':contrasts,'cooling_on_vs_off_relative_state_difference':cooling_effect,'tightened_thermal_substep_relative_state_difference':thermal_precision,'declared_energy_tolerance':a.thermal_energy_tolerance,'input_sha256':all_inputs,
 'scope':'Actual gas transport with old-density/old-force and new-density/new-force half sources, common hierarchy timestep and fixed native topology. Independent source-field snapshots and per-leaf momentum/kinetic-work accounting. Thermal cases additionally reconcile original cooling/heating, initial temperature clipping and complete AMR call budgets.',
 'limits':'Noncosmological periodic smooth gas at one spatial resolution and short endpoint. Cooling uses a level-local Lie split, not demonstrated second-order full physics. No total QUMOND field-plus-gas energy conservation, exact self-force cancellation, particles, stellar/AGN feedback, enrichment, MPI/regrid/restart or galaxy-formation readiness.'},indent=2,allow_nan=False)+'\n')
