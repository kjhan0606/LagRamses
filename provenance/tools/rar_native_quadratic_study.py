#!/usr/bin/env python3
"""Native-leaf quadratic reconstruction study; not a Fortran validation claim."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from rar_native_particles import epoch

p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('out',type=Path);args=p.parse_args()
leaf=np.loadtxt(args.run/'newtonian_a/composite_epoch_00000_leaf.tsv')
face=np.loadtxt(args.run/'newtonian_a/composite_faces.tsv')
n=len(leaf);nf=len(face);center=leaf[:,2:5];width=leaf[:,5];volume=width**3
l=face[:,0].astype(int)-1;r=face[:,1].astype(int)-1;axis=face[:,2].astype(int)-1
area=face[:,3];distance=face[:,4];delta=center[r]-center[l];delta-=np.rint(delta)
constant=np.ones(n)
D=sparse.csr_matrix((np.r_[area/volume[l],-area/volume[r]],(np.r_[l,r],np.tile(np.arange(nf),2))),shape=(n,nf))
J=sparse.csr_matrix((np.r_[-np.ones(nf),np.ones(nf)],(np.tile(np.arange(nf),2),np.r_[l,r])),shape=(nf,n))
diag=lambda values:sparse.diags(values,dtype=float,format='csr')
norm=lambda x:float(np.sqrt(np.sum(volume*x*x)))


def linear_coefficients():
    matrices=np.zeros((n,3,3));w=area/distance
    for a in range(3):
        for b in range(3):
            term=w*delta[:,a]*delta[:,b]
            matrices[:,a,b]=np.bincount(l,term,minlength=n)+np.bincount(r,term,minlength=n)
    inverse=np.linalg.inv(matrices);result=[]
    for a in range(3):
        entries=[]
        for cells in (l,r):
            entries.append(w*np.einsum('ij,ij->i',inverse[cells,a,:],delta))
        # Each adjacent cell uses the same oriented potential difference.
        C=sparse.csr_matrix((np.r_[entries[0],entries[1]],(np.r_[l,r],np.tile(np.arange(nf),2))),shape=(n,nf))
        result.append((C@J).tocsr())
    return result,None,{}


def quadratic_coefficients():
    graph=[set() for _ in range(n)]
    for i,j in zip(l,r):graph[i].add(j);graph[j].add(i)
    row=[];column=[];coefficients=[];max_condition=0.;max_reproduction=0.;minimum_neighbors=n
    for i in range(n):
        neighbors=set(graph[i])
        for j in graph[i]:neighbors.update(graph[j])
        neighbors.discard(i);neighbors=np.array(sorted(neighbors),int)
        s=center[neighbors]-center[i];s-=np.rint(s);s/=width[i]
        basis=np.column_stack((s,.5*s*s,s[:,0]*s[:,1],s[:,0]*s[:,2],s[:,1]*s[:,2]))
        weights=1/np.sum(s*s,axis=1)
        normal=basis.T@(weights[:,None]*basis)
        condition=np.linalg.cond(normal);max_condition=max(max_condition,float(condition))
        mapping=np.linalg.solve(normal,basis.T*weights)
        max_reproduction=max(max_reproduction,float(np.max(abs(mapping@basis-np.eye(9)))))
        minimum_neighbors=min(minimum_neighbors,len(neighbors))
        # Convert dimensionless polynomial coefficients to physical derivatives.
        mapping[:3]/=width[i];mapping[3:]/=width[i]**2
        row.extend([i]*(len(neighbors)+1));column.extend([*neighbors,i])
        coefficients.append(np.column_stack((mapping,-mapping.sum(axis=1))))
    data=np.concatenate(coefficients,axis=1)
    C=[sparse.csr_matrix((data[a],(row,column)),shape=(n,n)) for a in range(9)]
    H=[[None]*3 for _ in range(3)]
    for a in range(3):H[a][a]=C[3+a]
    for a,b,k in ((0,1,6),(0,2,7),(1,2,8)):H[a][b]=H[b][a]=C[k]
    if max_reproduction>1e-11 or not np.isfinite(max_condition):raise ValueError('quadratic stencil rank/consistency')
    return C[:3],H,{'max_normal_condition':max_condition,'max_polynomial_reproduction_error':max_reproduction,'minimum_neighbors':minimum_neighbors}


def operator(G,H):
    coarse_l=width[l]>width[r];coarse_r=width[r]>width[l];interface=coarse_l|coarse_r
    tangent=delta.copy();tangent[np.arange(nf),axis]=0
    difference=J.copy()
    for d in range(3):
        difference-=diag(tangent[:,d]*coarse_l)@G[d][l]
        difference-=diag(tangent[:,d]*coarse_r)@G[d][r]
    if H is not None:
        for d in range(3):
            for e in range(3):
                difference-=diag(.5*tangent[:,d]*tangent[:,e]*coarse_l)@H[d][e][l]
                difference+=diag(.5*tangent[:,d]*tangent[:,e]*coarse_r)@H[d][e][r]
    normal=diag(1/distance)@difference
    h_l=width[l]/2;h_r=width[r]/2
    if H is not None:
        for d in range(3):
            mask=interface&(axis==d)
            normal+=diag(.5*h_l*h_l/distance*mask)@H[d][d][l]
            normal-=diag(.5*h_r*h_r/distance*mask)@H[d][d][r]
    F=[]
    for d in range(3):
        if H is None:
            tangential=.5*(G[d][l]+G[d][r])
        else:
            gl=G[d][l].copy();gr=G[d][r].copy()
            for e in range(3):
                gl+=diag(tangent[:,e]*coarse_l)@H[d][e][l]
                gr-=diag(tangent[:,e]*coarse_r)@H[d][e][r]
            tangential=diag(h_r/distance)@gl+diag(h_l/distance)@gr
        F.append((diag(axis==d)@normal+diag(axis!=d)@tangential).tocsr())
    L=(D@normal).tocsr();zero=float(np.max(abs(L@constant)))
    conserved=float(np.max(abs(volume@L)))
    if zero>1e-9 or conserved>1e-12:raise ValueError(f'constant/conservation failure {zero} {conserved}')
    return L,F,{'constant_residual':zero,'volume_conservation_residual':conserved,'nonzeros':L.nnz}


def evaluate(name,G,H,diagnostics):
    L,F,checks=operator(G,H);diagnostics.update(checks)
    # Fix a gauge for the sparse direct study. Native Fortran uses FGMRES;
    # a successful study says nothing about that compiled solver yet.
    factor=splu(L[1:,1:].tocsc())
    def solve(rhs):
        rhs=rhs-np.sum(rhs*volume)/volume.sum()
        phi=np.r_[0.,factor.solve(rhs[1:])];phi-=np.sum(phi*volume)/volume.sum()
        residual=norm(L@phi-rhs)/max(norm(rhs),1e-300)
        if residual>1e-9:raise ValueError(f'{name} true residual {residual}')
        return phi,residual
    def mond(pn,a0,external):
        gf=-np.stack([f@pn for f in F],axis=1)+external
        mag=np.linalg.norm(gf,axis=1);boost=np.zeros(nf);active=mag>0
        boost[active]=2*a0/(np.sqrt(mag[active]**2+4*a0*mag[active])+mag[active])
        ph=-(D@(boost*gf[np.arange(nf),axis]))/(4*np.pi)
        return ph
    k=2*np.pi;exact=np.cos(k*center[:,0])/k
    pn,rn=solve(-k*k*exact);expected_n=np.zeros((n,3));expected_n[:,0]=np.sin(k*center[:,0])
    ph=mond(pn,1.,np.array([0.,1.,0.]));pq,rq=solve(-k*k*exact+4*np.pi*ph)
    expected_q=expected_n.copy();expected_q[:,0]*=.5*(1+np.sqrt(1+4/np.sqrt(1+expected_n[:,0]**2)))
    fn=-np.stack([g@pn for g in G],axis=1);fq=-np.stack([g@pq for g in G],axis=1)
    relative=lambda a,b:float(np.sqrt(np.sum(volume[:,None]*(a-b)**2)/np.sum(volume[:,None]*b*b)))
    diagnostics['continuum']={'newtonian_force_relative_l2':relative(fn,expected_n),'qumond_force_relative_l2':relative(fq,expected_q),'N_residual':rn,'Q_residual':rq}
    cases={}
    for fixture in ('a','b','pair'):
        path=args.run/f'newtonian_{fixture}';old=np.loadtxt(path/'composite_epoch_00000_leaf.tsv')
        rows,rho,_=epoch(path/'particles_epoch_00000.tsv',old)
        if np.max(abs(old[:,6]-1-rho))>1e-13:raise ValueError('source mismatch')
        rhs=4*np.pi*(rho-np.sum(rho*volume));pn,rn=solve(rhs)
        qpath=args.run/f'qumond_{fixture}/composite_epoch_00000_leaf.tsv'
        constants=np.array(qpath.read_text().splitlines()[0].split()[2:],float);a0=constants[1]
        ph=mond(pn,a0,constants[2:]);pq,rq=solve(rhs+4*np.pi*ph)
        for model,potential,res in (('newtonian',pn,rn),('qumond',pq,rq)):
            reconstructed=old.copy();reconstructed[:,10:13]=-np.stack([g@potential for g in G],axis=1)
            # Independent separable integrations of each physical cubic cloud.
            s=.0625;lo=center-width[:,None]/2;hi=center+width[:,None]/2;weights=[]
            for x in rows[:,4:7]:
                length=np.zeros_like(center)
                for shift in (-1.,0.,1.):length+=np.maximum(0,np.minimum(hi,x+shift+s/2)-np.maximum(lo,x+shift-s/2))
                weights.append(np.prod(length,axis=1)/s**3)
            acceleration=np.asarray(weights)@reconstructed[:,10:13]
            observed=np.loadtxt(args.run/f'{model}_{fixture}/particles_epoch_00000.tsv',ndmin=2)
            observed=observed[np.argsort(observed[:,1])]
            replay=float(np.max(abs(acceleration-observed[:,10:13])))
            if name=='linear' and replay>1e-7:raise ValueError('baseline does not reproduce compiled force')
            mass=rows[:,3];forces=mass[:,None]*acceleration
            cases[f'{model}_{fixture}']={'acceleration':acceleration.tolist(),'true_residual':res,
                'relative_net_pair_force':float(np.linalg.norm(forces.sum(axis=0))/max(np.linalg.norm(forces,axis=1).sum(),1e-300)),
                'compiled_baseline_max_difference':replay}
            if fixture!='pair':cases[f'{model}_{fixture}']['self_over_Gm_s2']=float(np.linalg.norm(acceleration[0])/(mass[0]/s**2))
    return {'geometry':diagnostics,'cases':cases}

result={}
for name,builder in (('linear',linear_coefficients),('quadratic',quadratic_coefficients)):
    G,H,diagnostics=builder();result[name]=evaluate(name,G,H,diagnostics)
ratios={}
for model in ('newtonian','qumond'):
    ratios[model]=result['quadratic']['cases'][f'{model}_a']['self_over_Gm_s2']/result['linear']['cases'][f'{model}_a']['self_over_Gm_s2']
result['self_force_ratio_quadratic_over_linear']=ratios
result['scope']='Native geometry sparse prototype using the same source, cloud, face conservation and QUMOND equation. Direct solve rather than compiled FGMRES; no implementation readiness claim.'
args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'ratios':ratios,'continuum':{k:v['geometry']['continuum'] for k,v in result.items() if k in ('linear','quadratic')}},indent=2))
