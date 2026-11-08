#!/usr/bin/env python3
"""Independent NumPy reconstruction of native subface operators and saved fields."""
import argparse,json,math,re
from pathlib import Path
import h5py
import numpy as np
p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('output',type=Path);p.add_argument('--case',default='new_on');p.add_argument('--snapshots',type=int,default=2)
p.add_argument('--leaf-file',default='composite_leaf_state.tsv');p.add_argument('--snapshot-glob',default='output_*/data_*.h5');p.add_argument('--evolving',action='store_true');p.add_argument('--particles',action='store_true');a=p.parse_args()
case=a.run/a.case;leaf=np.loadtxt(case/a.leaf_file);face=np.loadtxt(case/'composite_faces.tsv')
constants=np.array((case/a.leaf_file).read_text().splitlines()[0].split()[2:],float)
fourpi,a0=constants[:2];external=constants[2:]
if leaf.ndim!=2 or leaf.shape[1]!=13 or face.ndim!=2 or face.shape[1]!=5:raise ValueError('invalid export shape')
if not np.isfinite(leaf).all() or not np.isfinite(face).all():raise ValueError('nonfinite export')
levels=leaf[:,1].astype(int);center=leaf[:,2:5];width=leaf[:,5];volume=width**3;n=len(leaf)
l=face[:,0].astype(int)-1;r=face[:,1].astype(int)-1;axis=face[:,2].astype(int)-1;area=face[:,3];distance=face[:,4]
if min(l.min(),r.min())<0 or max(l.max(),r.max())>=n or np.any(l==r) or np.any((axis<0)|(axis>2)):raise ValueError('invalid connectivity')
if len(np.unique(face[:,:3],axis=0))!=len(face):raise ValueError('duplicate subface')
if not np.allclose(width,2.**(-levels),rtol=0,atol=1e-15) or abs(volume.sum()-1)>1e-12:raise ValueError('invalid leaf volume')
delta=center[r]-center[l];delta-=np.rint(delta)
if not np.allclose(delta[np.arange(len(face)),axis],distance,rtol=0,atol=1e-14):raise ValueError('face orientation/distance')
if not np.allclose(area,np.minimum(width[l],width[r])**2,rtol=0,atol=1e-15):raise ValueError('subface area')
coverage=np.zeros((n,3,2));np.add.at(coverage,(l,axis,np.zeros(len(l),int)),area);np.add.at(coverage,(r,axis,np.ones(len(r),int)),area)
if not np.allclose(coverage,width[:,None,None]**2,rtol=1e-13,atol=0):raise ValueError('incomplete face coverage')
w=area/distance;matrix=np.zeros((n,3,3))
for i in range(3):
 for j in range(3):
  term=w*delta[:,i]*delta[:,j]
  matrix[:,i,j]=np.bincount(l,term,minlength=n)+np.bincount(r,term,minlength=n)
def gradient(phi):
 rhs=np.empty((n,3))
 for d in range(3):
  term=w*delta[:,d]*(phi[r]-phi[l])
  rhs[:,d]=np.bincount(l,term,minlength=n)+np.bincount(r,term,minlength=n)
 return np.linalg.solve(matrix,rhs[...,None])[...,0]
tangent=delta.copy();tangent[np.arange(len(face)),axis]=0
coarse_l=width[l]>width[r];coarse_r=width[r]>width[l]
def faces(phi,g):
 pl=phi[l]+np.einsum('ij,ij->i',g[l],tangent)*coarse_l
 pr=phi[r]-np.einsum('ij,ij->i',g[r],tangent)*coarse_r
 out=.5*(g[l]+g[r]);out[np.arange(len(face)),axis]=(pr-pl)/distance
 return out
def divergence(normal):
 term=area*normal
 return (np.bincount(l,term,minlength=n)-np.bincount(r,term,minlength=n))/volume
def norm(v):return np.sqrt(np.sum(volume*v*v))
source,pn,pq,ph=leaf[:,6],leaf[:,7],leaf[:,8],leaf[:,9]
particle_density=np.zeros(n);particle_report=None
if a.particles:
 from rar_native_particles import epoch
 index=re.fullmatch(r'composite_epoch_(\d+)_leaf.tsv',a.leaf_file)
 if not a.evolving or index is None:raise ValueError('particles require an evolving field epoch')
 _,particle_density,particle_report=epoch(case/f'particles_epoch_{int(index[1]):05d}.tsv',leaf)
gn=gradient(pn);gq=gradient(pq);fn=faces(pn,gn);fq=faces(pq,gq)
rhsn=fourpi*(source-np.sum(volume*source)/volume.sum())
resn=norm(divergence(fn[np.arange(len(face)),axis])-rhsn)/norm(rhsn)
gfield=-fn+external;mag=np.linalg.norm(gfield,axis=1);boost=np.zeros_like(mag)
active=mag>0;boost[active]=2*a0/(np.sqrt(mag[active]**2+4*a0*mag[active])+mag[active])
ph_ref=-divergence(boost*gfield[np.arange(len(face)),axis])/fourpi
ph_error=norm(ph-ph_ref)/max(norm(ph_ref),1e-300)
rhsq=rhsn+fourpi*ph_ref;rhsq-=np.sum(volume*rhsq)/volume.sum()
resq=norm(divergence(fq[np.arange(len(face)),axis])-rhsq)/norm(rhsq)
force_error=np.sqrt(np.sum(volume[:,None]*(gq+leaf[:,10:13])**2)/np.sum(volume[:,None]*gq*gq))
if max(resn,resq)>1e-9 or ph_error>1e-11 or force_error>1e-11:raise ValueError(f'independent operator mismatch {resn} {resq} {ph_error} {force_error}')
nml=(case/'effective.nml').read_text()
def scalar(name):
 v=re.findall(r'(?im)^\s*'+name+r'\s*=\s*([^\n]+)',nml)
 if len(v)!=1:raise ValueError('missing/ambiguous '+name)
 return float(v[0].replace('d','e').replace('D','E'))
params=re.findall(r'(?im)^\s*gravity_params\s*=\s*([^\n]+)',nml)
params=np.array([float(x.replace('d','e').replace('D','E')) for x in params[0].split(',')])
params=np.pad(params,(0,10-len(params)))
kind=int(scalar('gravity_type'))
if kind==0 and a.evolving:
 # The independently saved gas snapshot supplies the source oracle below.
 expected=source.copy()
elif kind==-1:
 rad=np.linalg.norm(center-params[1:4],axis=1)
 cap=1/(params[0]*(1+params[0])**2)
 expected=scalar('d_region')+np.minimum(1/(np.maximum(rad,1e-300)*(1+rad)**2),cap)
elif kind in (-2,-3):
 xyz=center-.5*scalar('boxlen');rad=np.linalg.norm(xyz[:,:2],axis=1)
 expected=np.full(n,scalar('d_region'))
 discs=np.zeros(n)
 for mass,scale,height in (params[:3],params[3:6]):
  discs+=mass/(4*np.pi*scale**2*height)*np.exp(-rad/scale)/np.cosh(xyz[:,2]/height)**2
 if params[9]>0:
  def taper(r):return .5*(1+np.cos(np.pi*np.clip((r-.8*params[9])/(.2*params[9]),0,1)))
  discs*=taper(rad)*taper(abs(xyz[:,2]))
 expected+=discs
 if kind==-3:
  rr=np.linalg.norm(xyz,axis=1);c=params[8]/params[7];u=rr/params[7]
  normalization=params[6]/(4*np.pi*params[7]**3*(np.log1p(c)-c/(1+c)))
  expected+=np.where(rr<=params[8],normalization/(u*(1+u)**2),0)
else:raise ValueError('unsupported independent source oracle')
source_error=np.max(abs(source-expected))/np.max(abs(expected))
if source_error>1e-13:raise ValueError('physical analytic source mismatch')
mond=re.findall(r'(?im)^\s*use_mond\s*=\s*(\.[a-z]+\.)',nml)
if len(mond)!=1:raise ValueError('missing gravity model declaration')
expected_a0=scalar('a0_mond')*scalar('units_time')**2/scalar('units_length') if mond[0].lower()=='.true.' else 0.
if not math.isclose(a0,expected_a0,rel_tol=1e-14):raise ValueError('MOND units/model mismatch')
if expected_a0==0 and (np.any(ph!=0) or np.any(pn!=pq) or np.any(external!=0)):raise ValueError('Newtonian control acquired a MOND contribution')
keys=[(int(lev),*np.rint(c*2**(int(lev)+1)).astype(int)) for lev,c in zip(levels,center)]
lookup={k:i for i,k in enumerate(keys)}
if len(lookup)!=n:raise ValueError('duplicate leaf positions')
snapshots=[]
offsets=np.array([[((child>>d)&1)-.5 for d in range(3)] for child in range(8)])
for snapshot in sorted(case.glob(a.snapshot_glob)):
 count=0;max_phi=0.;max_force=0.;covered_phi=0.;covered_force=0.;max_momentum=0.;energy_error=0.
 with h5py.File(snapshot) as h:
  saved={}
  for lev in sorted(set(levels)):
   lev=int(lev);prefix=f'amr/level_{lev}'
   grids=np.stack([h[f'{prefix}/xg_{d}'][:] for d in (1,2,3)],axis=1)
   centers=np.repeat(grids,8,axis=0)+np.tile(offsets,(len(grids),1))*2.**(-lev)
   leaf_mask=h[f'{prefix}/son_flag'][:]==0
   phi=h[f'gravity/level_{lev}/phi'][:];force=np.stack([h[f'gravity/level_{lev}/f_{d}'][:] for d in (1,2,3)],axis=1)
   indices=np.array([lookup[(lev,*np.rint(c*2**(lev+1)).astype(int))] for c in centers[leaf_mask]])
   count+=len(indices)
   max_phi=max(max_phi,float(np.max(abs(phi[leaf_mask]-pq[indices]))))
   max_force=max(max_force,float(np.max(abs(force[leaf_mask]-leaf[indices,10:13]))))
   density_saved=h[f'hydro/level_{lev}/uold_1'][:]
   if a.evolving:
    if kind!=0:raise ValueError('evolving source oracle only supports gas self-gravity')
    if not np.isfinite(density_saved).all() or np.any(density_saved<=0):raise ValueError('invalid evolving density')
    source_error=max(source_error,float(np.max(abs(density_saved[leaf_mask]+particle_density[indices]-source[indices]))/np.max(abs(source[indices]))))
    if source_error>1e-13:raise ValueError('evolving density/source epoch mismatch')
   elif not np.allclose(density_saved,scalar('d_region'),rtol=1e-13,atol=0):raise ValueError('gas density evolved')
   for d in (2,3,4):
    max_momentum=max(max_momentum,float(np.max(abs(h[f'hydro/level_{lev}/uold_{d}'][:]))))
   if not a.evolving:
    if max_momentum>1e-14:raise ValueError('gas received gravity kicks')
    initial_energy=scalar('p_region')/(scalar('gamma')-1)
    energy_error=max(energy_error,float(np.max(abs(h[f'hydro/level_{lev}/uold_5'][:]-initial_energy)))/initial_energy)
    if energy_error>1e-13:raise ValueError('frozen gas energy changed')
   else:
    momenta=np.stack([h[f'hydro/level_{lev}/uold_{d}'][:] for d in (2,3,4)],axis=1)
    internal=h[f'hydro/level_{lev}/uold_5'][:]-np.sum(momenta**2,axis=1)/(2*density_saved)
    if not np.isfinite(internal).all() or np.any(internal<=0):raise ValueError('nonpositive evolving internal energy')
   saved[lev]=(grids,centers,leaf_mask,phi,force)
  lo,hi=min(saved),max(saved)
  grids,_,_,phi_hi,force_hi=saved[hi]
  fine_map={tuple(np.rint(c*2**(lo+1)).astype(int)):i for i,c in enumerate(grids)}
  _,centers,mask,phi_lo,force_lo=saved[lo]
  child_phi=phi_hi.reshape(-1,8).mean(axis=1);child_force=force_hi.reshape(-1,8,3).mean(axis=1)
  for idx in np.flatnonzero(~mask):
   j=fine_map[tuple(np.rint(centers[idx]*2**(lo+1)).astype(int))]
   covered_phi=max(covered_phi,abs(phi_lo[idx]-child_phi[j]));covered_force=max(covered_force,float(np.max(abs(force_lo[idx]-child_force[j]))))
 if count!=n or max_phi>1e-11*max(np.max(abs(pq)),1) or max_force>1e-11*max(np.max(abs(leaf[:,10:13])),1):raise ValueError('published leaf fields differ')
 if max(covered_phi,covered_force)>1e-12:raise ValueError('covered-cell restriction mismatch')
 snapshots.append({'snapshot':str(snapshot),'leaves':count,'phi_max_difference':max_phi,'force_max_difference':max_force,'covered_phi_difference':covered_phi,'covered_force_difference':covered_force,'gas_max_momentum':max_momentum,'gas_energy_relative_difference':None if a.evolving else energy_error})
if len(snapshots)!=a.snapshots:raise ValueError('incorrect snapshot coverage')
a.output.write_text(json.dumps({'schema':'rar_composite_publication_independent_v1','leaves':n,'subfaces':len(face),'newtonian_relative_residual':resn,'qumond_relative_residual':resq,'phantom_relative_difference':ph_error,'force_relative_difference':force_error,'source_relative_difference':source_error,'snapshots':snapshots,
 'gas_state_mode':'evolving' if a.evolving else 'fixed','particle_source':particle_report,
 'limits':'Independent native operator/source/field snapshot reconstruction. Transport/source accounting, continuum accuracy, particles and evolving MPI/restart require their own evidence.'},indent=2,allow_nan=False)+'\n')
