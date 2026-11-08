#!/usr/bin/env python3
"""Audit the compiled FR02 NGC 2403 fixed-source periodic fields."""
import argparse,hashlib,json,math
from pathlib import Path
import h5py
import numpy as np
import torch
from analyze_qumond_spherical_reference import read_uniform_snapshot
from analyze_qumond_source_diagnostic import parse_log

parser=argparse.ArgumentParser();parser.add_argument('root',type=Path);parser.add_argument('n',type=int,nargs='?',default=512);parser.add_argument('--output',type=Path)
args=parser.parse_args();root=args.root;n=args.n
output=args.output if args.output is not None else root/'analysis.json'
if output.exists():raise FileExistsError(output)
if n not in (128,512): raise ValueError('FR02 audit requires 128^3 or 512^3')
level=int(math.log2(n))
# Check historical log compatibility before allocating the large field arrays.
# This run predates post-second-Poisson telemetry. Its saved force/source
# equations are independently reconstructed below; absent telemetry is not
# relabelled as a measured post-solve source trace.
diags,post_diags,fft_diags=parse_log((root/'qumond/run.log').read_text(),allow_source_only=True)
if len(diags)!=2 or [(d['step'],d['level']) for d in diags]!=[(0,level),(1,level)]:
    raise ValueError('unexpected fixed-source diagnostic call coverage')
for d in diags:
    if not (d['cells']==n**3 and d['updated_cells']==n**3 and d['skipped_cells']==0 and d['missing_face_links']==0):
        raise ValueError('incomplete phantom source cell coverage')
    if set(d['fields'])!={'newtonian_source','phantom_source','combined_source'}:
        raise ValueError('incomplete historical source statistics')
    if not all(math.isfinite(v) for field in d['fields'].values() for v in field.values()):
        raise ValueError('nonfinite historical source statistics')
    if d['fields']['phantom_source']['rms']<=1e-6:raise ValueError('missing nonzero phantom source')
level=int(math.log2(n)); T=lambda x:torch.as_tensor(x,device='cuda',dtype=torch.float64)
if not torch.cuda.is_available():raise SystemExit('GPU required')
inputs=root/'star_profile_code.dat',root/'hi_profile_code.dat',root/'source_contract.json'
meta=json.loads(inputs[2].read_text())
if (meta.get('stellar_rows'),meta.get('gas_rows')) != (73,44): raise ValueError('unexpected FR02 tables')
if not math.isclose(float(meta['stellar_mass_to_light']),0.5,rel_tol=0.0,abs_tol=1e-12): raise ValueError('unexpected stellar M/L')
if not math.isclose(float(meta['gas_helium_factor']),1.33,rel_tol=0.0,abs_tol=1e-12): raise ValueError('unexpected helium factor')
if not math.isclose(float(meta['box_kpc']),64.0,rel_tol=0.0,abs_tol=1e-12): raise ValueError('unexpected source box')
def table(p):
    a=np.loadtxt(p,skiprows=1);assert len(a)==int(p.read_text().splitlines()[0]);return T(a)
st=table(inputs[0]);hi=table(inputs[1]);dx=1/n
a=(torch.arange(n,device='cuda',dtype=torch.float64)+.5)*dx-.5
r=torch.sqrt(a[:,None]**2+a[None,:]**2);z=torch.abs(a)
def interpolate(x,p):
    i=torch.searchsorted(p[:,0].contiguous(),x).clamp(1,len(p)-1);u=((x-p[i-1,0])/(p[i,0]-p[i-1,0])).clamp(0,1)
    return p[i-1,1]*(1-u)+p[i,1]*u
def vert(h):return (torch.exp(-(z-dx/2)/h)-torch.exp(-(z+dx/2)/h))/(2*dx)
star=interpolate(r,st);star=torch.where(r<=st[-1,0],star,torch.zeros_like(star))
gas=interpolate(r,hi)*((meta['gas_taper_end_code']-r)/(meta['gas_taper_end_code']-meta['gas_taper_start_code'])).clamp(0,1)
rho=star[:,:,None]*vert(meta['star_exp_scaleheight_code'])[None,None,:]+gas[:,:,None]*vert(meta['gas_exp_scaleheight_code'])[None,None,:]
idx=torch.arange(n,device='cuda',dtype=torch.float64);idz=torch.arange(n//2+1,device='cuda',dtype=torch.float64)
ev=2*(torch.cos(2*math.pi*idx[:,None,None]/n)+torch.cos(2*math.pi*idx[None,:,None]/n)+torch.cos(2*math.pi*idz[None,None,:]/n)-3)
green=torch.where(ev!=0,dx**2/ev,torch.zeros_like(ev))
def poisson(d):return torch.fft.irfftn(torch.fft.rfftn(4*math.pi*(d-d.mean()))*green,s=(n,n,n))
def force(p):
    return torch.stack([(2/(3*dx))*(torch.roll(p,1,i)-torch.roll(p,-1,i))-(1/(12*dx))*(torch.roll(p,2,i)-torch.roll(p,-2,i)) for i in range(3)])
newt=force(poisson(rho));mag=torch.linalg.vector_norm(newt,dim=0)
a0=1.2e-8/(meta['unit_length_cm']/(3.8707679655184575e15**2))
nu=.5*(torch.sqrt(1+4/(mag/a0).clamp(min=1e-30))-1)
vec=newt*nu
div=sum((torch.roll(vec[i],-1,i)-torch.roll(vec[i],1,i))/(2*dx) for i in range(3))
phantom=-div/(4*math.pi);qref=force(poisson(rho+phantom))
del vec,nu,div,mag,phantom,ev
def sha256_file(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(32*1024*1024),b''):digest.update(block)
    return digest.hexdigest()
rec={}
for label,reference in [('newtonian',newt),('qumond',qref)]:
    path=root/label/'output_00001/data_00001.h5'
    measured,coords,box,mesh=read_uniform_snapshot(path);assert mesh==n and box==1
    v=T(measured); norm=torch.linalg.vector_norm(reference,dim=0)
    R=torch.sqrt(a[:,None,None]**2+a[None,:,None]**2+a[None,None,:]**2)
    # Include the selected SPARC radii, which reach 20.87 kpc in the 64 kpc box.
    mask=(R>=.03)&(R<=.35)&(norm>norm.max()*1e-4)
    errors=torch.linalg.vector_norm(v-reference,dim=0)[mask]/norm[mask]
    # torch.quantile rejects this >2**24-element mask. NumPy's exact linear
    # quantile has the same interpolation convention and no such size limit.
    # Preserve the original lower-median statistic from torch.median.
    rec[label]={'median_relative_force_error':float(torch.median(errors)),
                'p90_relative_force_error':float(np.quantile(errors.detach().cpu().numpy(),.9,method='linear')),
                'force_sample_count':int(errors.numel()),'quantile_method':'numpy_exact_linear',
                'snapshot':str(path),'snapshot_sha256':sha256_file(path)}
    if rec[label]['median_relative_force_error']>1e-6 or rec[label]['p90_relative_force_error']>1e-6:raise ValueError(f'{label} force gate failed')
    if label=='newtonian':
        with h5py.File(path,'r') as f:
            am=f[f'amr/level_{level}'];gr=f[f'gravity/level_{level}']
            xg=np.column_stack([np.asarray(am[f'xg_{i}']) for i in (1,2,3)])
            child=np.arange(8);xyz=np.column_stack((child%2,(child//2)%2,child//4))
            pos=(xg[:,None,:]+(xyz[None,:,:]-.5)*dx)%1
            ijk=np.rint(pos/dx-.5).astype(np.int64)%n
            p=np.empty((n,n,n));p[ijk[...,0],ijk[...,1],ijk[...,2]]=np.asarray(gr['phi']).reshape(len(xg),8)
        pt=T(p);lap=sum((torch.roll(pt,1,i)+torch.roll(pt,-1,i)-2*pt)/dx**2 for i in range(3))/(4*math.pi)
        diff=lap-(rho-rho.mean());metric=float(torch.sqrt(torch.mean(diff**2))/torch.sqrt(torch.mean((rho-rho.mean())**2)))
        rec[label]['source_rms_relative_error']=metric
        if metric>1e-6:raise ValueError(f'source gate failed: {metric}')
        del xg,child,xyz,pos,ijk,p,pt,lap,diff
    del measured,coords,v,norm,R,mask,errors
result={'status':'PERIODIC_SOURCE_AND_FORCE_GATE_PASS','cases':rec,'qumond_diagnostic_calls':len(diags),
        'diagnostic_scope':'source_only' if not post_diags else 'source_and_post_second_poisson',
        'post_second_poisson_diagnostic_calls':len(post_diags),'fft_diagnostic_calls':len(fft_diags),
        'scope':'NGC 2403 graph-digitized FR02 gas plus fixed stellar source; uniform periodic QUMOND, no CDM halo or cosmological evolution',
        'grid':{'nmesh':n,'box_kpc':64.0,'cell_width_kpc':64.0/n},
        'stellar_ml':float(meta['stellar_mass_to_light']),
        'source_files_sha256':{p.name:sha256_file(p) for p in inputs}}
with output.open('x') as stream:stream.write(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
