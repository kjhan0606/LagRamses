#!/usr/bin/env python3
"""Matched prescribed-disc rotation curves from native leaf force exports."""
import argparse,csv,hashlib,json,math,re
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('low',type=Path);p.add_argument('high',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
modes=('newtonian_baryon','newtonian_halo','qumond_baryon')
def scalar(text,key):
 values=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',text)
 if len(values)!=1:raise ValueError('ambiguous '+key)
 return float(values[0].replace('d','e').replace('D','E'))
def load(run):
 states={};inputs={};nmls={};constants={};commits=set()
 for mode in modes:
  case=run/mode;nml=(case/'effective.nml').read_text();nmls[mode]=nml
  constants[mode]=np.array((case/'composite_leaf_state.tsv').read_text().splitlines()[0].split()[2:],float)
  states[mode]=np.loadtxt(case/'composite_leaf_state.tsv')
  commits.add((case/'source_commit.txt').read_text().strip())
  for name in ('source_commit.txt','effective.nml','composite_leaf_state.tsv','composite_faces.tsv'):
   f=case/name;inputs[str(f)]=hashlib.sha256(f.read_bytes()).hexdigest()
 if len(commits)!=1:raise ValueError('unmatched source revisions within suite')
 ref=states[modes[0]];bar=ref[:,6]
 for mode in modes:
  if not np.array_equal(ref[:,:6],states[mode][:,:6]):raise ValueError('unmatched hierarchy')
  if inputs[str(run/mode/'composite_faces.tsv')]!=inputs[str(run/modes[0]/'composite_faces.tsv')]:raise ValueError('unmatched shared subfaces')
 if not np.array_equal(bar,states['qumond_baryon'][:,6]):raise ValueError('baryons differ between gravity models')
 if np.any(states['newtonian_halo'][:,6]<bar) or not np.any(states['newtonian_halo'][:,6]>bar):raise ValueError('missing/negative CDM halo')
 difference=states['qumond_baryon'][:,7]-ref[:,7]
 pn_error=np.linalg.norm(difference)/np.linalg.norm(ref[:,7])
 if pn_error>2e-8:raise ValueError('QUMOND first stage differs from baryon-only Newtonian solution')
 lev=int(scalar(nmls[modes[0]],'levelmax'));h=2.**(-lev)
 grids={}
 for mode in modes:
  v=states[mode];mask=(v[:,1]==lev)&np.isclose(abs(v[:,4]-.5),.5*h,rtol=0,atol=1e-14)
  xy=np.floor(v[mask,2:4]/h).astype(int);grid=np.zeros((2**lev,2**lev,3));count=np.zeros(grid.shape[:2],int)
  np.add.at(grid,(xy[:,0],xy[:,1]),v[mask,10:13]);np.add.at(count,(xy[:,0],xy[:,1]),1)
  grid[count>0]/=count[count>0,None];grids[mode]=(grid,count)
 meta={'leaves':len(ref),'levelmax':lev,'source_commit':next(iter(commits)),'baryon_source_bitwise_equal':True,'q_first_stage_relative_difference':pn_error,
       'sampled_source_integrals':{mode:float(np.sum(v[:,6]*v[:,5]**3)) for mode,v in states.items()},'inputs':inputs}
 return grids,constants['qumond_baryon'][1],nmls[modes[0]],meta
low,a0l,nmll,metal=load(a.low);high,a0h,nmlh,metah=load(a.high)
if metal['source_commit']!=metah['source_commit'] or metah['levelmax']!=metal['levelmax']+1:raise ValueError('unmatched source or grid refinement across resolutions')
if a0l!=a0h:raise ValueError('MOND parameter differs with resolution')
for key in ('units_length','units_time','units_density','gravity_params'):
 if re.findall(r'(?im)^\s*'+key+r'\s*=([^\n]+)',nmll)!=re.findall(r'(?im)^\s*'+key+r'\s*=([^\n]+)',nmlh):raise ValueError('physical source/units differ with resolution')
unit_l=scalar(nmlh,'units_length');unit_t=scalar(nmlh,'units_time');unit_d=scalar(nmlh,'units_density')
if not math.isclose(6.67430e-8*unit_t**2*unit_d,1,rel_tol=1e-12):raise ValueError('physical/code G convention mismatch')
radii=np.linspace(4*2.**(-metal['levelmax']),.25,32);angles=np.arange(128)*2*np.pi/128
rx=radii[:,None]*np.cos(angles);ry=radii[:,None]*np.sin(angles)
def sample(grids,lev):
 h=2.**(-lev);qx=(.5+rx)/h-.5;qy=(.5+ry)/h-.5;i=np.floor(qx).astype(int);j=np.floor(qy).astype(int);wx=qx-i;wy=qy-j
 result={}
 for mode,(grid,count) in grids.items():
  force=np.zeros((*rx.shape,3))
  for di,dj,weight in ((0,0,(1-wx)*(1-wy)),(1,0,wx*(1-wy)),(0,1,(1-wx)*wy),(1,1,wx*wy)):
   if np.any(count[i+di,j+dj]!=2):raise ValueError('polar sample lacks symmetric fine midplane cells')
   force+=weight[...,None]*grid[i+di,j+dj]
  inward=-(force[:,:,0]*rx+force[:,:,1]*ry)/radii[:,None]
  mean=inward.mean(axis=1)
  if not np.isfinite(inward).all() or np.any(mean<=0):raise ValueError('nonfinite/outward mean rotation force')
  result[mode]={'mean':mean,'azimuth_std':inward.std(axis=1),'vertical_rms':np.sqrt(np.mean(force[:,:,2]**2,axis=1))}
 return result
lc=sample(low,metal['levelmax']);hc=sample(high,metah['levelmax']);gbar=hc['newtonian_baryon']['mean'];gq=hc['qumond_baryon']['mean'];gcdm=hc['newtonian_halo']['mean']
galgebraic=.5*(1+np.sqrt(1+4*a0h/gbar))*gbar
delta=np.log10(gq/galgebraic);changes={mode:hc[mode]['mean']/lc[mode]['mean']-1 for mode in modes}
columns={'radius_code':radii,'radius_kpc':radii*unit_l/3.0856775814913673e21,'gbar_cgs':gbar*unit_l/unit_t**2,'g_lcdm_cgs':gcdm*unit_l/unit_t**2,'g_qumond_cgs':gq*unit_l/unit_t**2,'log10_qumond_over_algebraic_simple_nu':delta}
for mode in modes:
 columns['v_'+mode+'_kms']=np.sqrt(radii*hc[mode]['mean'])*unit_l/unit_t/1e5
 columns['fractional_force_resolution_change_'+mode]=changes[mode]
 columns['azimuth_relative_std_'+mode]=hc[mode]['azimuth_std']/hc[mode]['mean']
with (a.out/'rotation_rar.csv').open('w') as f:
 writer=csv.writer(f);writer.writerow(columns);writer.writerows(zip(*columns.values()))
result={'schema':'rar_native_disc_comparison_v1','low':metal,'high':metah,'radius_range_kpc':[float(columns['radius_kpc'][0]),float(columns['radius_kpc'][-1])],
 'a0_cgs':a0h*unit_l/unit_t**2,'force_resolution_changes':{mode:{'relative_rms':float(np.sqrt(np.mean(v*v))),'relative_max_abs':float(np.max(abs(v)))} for mode,v in changes.items()},
 'qumond_vs_algebraic_simple_nu_dex':{'mean':float(delta.mean()),'rms':float(np.sqrt(np.mean(delta*delta))),'max_abs':float(np.max(abs(delta)))},
 'lcdm_vs_qumond_force_ratio_range':[float((gcdm/gq).min()),float((gcdm/gq).max())],
 'sampling':'32 common radii, 128 azimuths each; symmetric fine z planes averaged to the midplane, followed by bilinear xy interpolation. Scatter is numerical azimuthal variation, not an intrinsic population RAR scatter.',
 'scope':'Same prescribed stellar/gas discs. Newtonian baryon+truncated NFW is the CDM-halo comparison; Newtonian baryon-only provides gbar. QUMOND uses the solved field, not an algebraic boost.',
 'limits':'Synthetic fixed periodic 10 kpc source; finite thickness, resolution and image/background effects. Halo is prescribed, not fitted or cosmologically formed. No observational model selection, self-consistent galaxy formation, population RAR scatter, MPI/evolving/restart claim.'}
(a.out/'result.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
import subprocess,sys
subprocess.run([sys.executable,str(Path(__file__).with_name('plot_rar_native_disc.py')),str(a.out)],check=True)
