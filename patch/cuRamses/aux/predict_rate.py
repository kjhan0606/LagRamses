#!/usr/bin/env python3
"""
Exact offline prediction of the SIDM MC scattering rate for a UNIGRID
(levelmin=levelmax=6) run, compared against the code's own per-step
'SIDM level 6: scattered=..., pairs=...' log lines.

  pairs_pred  = sum_cells floor(N/2)                      (pure geometry)
  events_pred = (sigma/m) * m_p * dt_phys / V_cell *
                sum_cells sum_{i<j} |v_i - v_j|           (full pair sum;
                equals the MC expectation floor(N/2)*(N-1)*<v_pair> with
                the odd-N correction N/(2*floor(N/2)) applied in code)

Units come straight from the info_XXXXX.txt of the snapshot.
"""
import numpy as np, re, sys, glob

outdir  = sys.argv[1]              # e.g. output_00002
logfile = sys.argv[2] if len(sys.argv) > 2 else 'run.log'
sigma_m = float(sys.argv[3]) if len(sys.argv) > 3 else 1000.0
NC = 64                            # 2^levelmin cells per side

# ---------- read info file ----------
info = {}
for line in open(glob.glob(f'{outdir}/info_?????.txt')[0]):
    m = re.match(r'\s*(\w+)\s*=\s*([\dEe\.\+\-]+)', line)
    if m: info[m.group(1)] = float(m.group(2))
aexp, unit_l, unit_d, unit_t = info['aexp'], info['unit_l'], info['unit_d'], info['unit_t']

# ---------- read particle files ----------
def read_records(fn):
    recs = []
    with open(fn,'rb') as f:
        while True:
            h = f.read(4)
            if len(h) < 4: break
            n = np.frombuffer(h,'<i4')[0]
            recs.append(f.read(n)); f.read(4)
    return recs

X=[];V=[];M=[]
for fn in sorted(glob.glob(f'{outdir}/part_?????.out?????')):
    recs = read_records(fn)
    npart = np.frombuffer(recs[2],'<i4')[0]
    big = [r for r in recs if len(r)==8*npart]
    # order: x,y,z,vx,vy,vz,m,idp,(potp)
    x = np.stack([np.frombuffer(big[i],'<f8') for i in range(3)],axis=1)
    v = np.stack([np.frombuffer(big[i],'<f8') for i in range(3,6)],axis=1)
    m = np.frombuffer(big[6],'<f8')
    X.append(x); V.append(v); M.append(m)
X=np.concatenate(X); V=np.concatenate(V); M=np.concatenate(M)
print(f'{outdir}: read {len(X)} particles, aexp={aexp:.4f}')

# ---------- bin into 64^3 cells; pair counts and pair |dv| sums ----------
ci = np.clip((X*NC).astype(np.int64),0,NC-1)
cell = (ci[:,0]*NC+ci[:,1])*NC+ci[:,2]
order = np.argsort(cell); cell=cell[order]; Vs=V[order]
bnd = np.searchsorted(cell, np.arange(NC**3+1)-0.5)  # cell boundaries
# faster: unique
uc, start, cnt = np.unique(cell, return_index=True, return_counts=True)
pairs_pred = int(np.sum(cnt//2))
sum_vrel = 0.0
for s,c in zip(start,cnt):
    if c < 2: continue
    vv = Vs[s:s+c]
    if c <= 512:
        d = vv[:,None,:]-vv[None,:,:]
        vm = np.sqrt((d**2).sum(-1))
        sum_vrel += vm[np.triu_indices(c,1)].sum()
    else:
        for i in range(c-1):
            sum_vrel += np.sqrt(((vv[i+1:]-vv[i])**2).sum(-1)).sum()

mp_code = np.median(M)
mp_g    = mp_code*unit_d*unit_l**3
# events per unit code-dt:  (sigma/m)*mp*<sum v_phys>*unit_t/V_cell_phys
#   v_phys = v_code*unit_l/unit_t ; V = (unit_l/NC)^3 ; dt_phys = dt_code*unit_t
K = sigma_m*mp_g*NC**3/unit_l**2   # so that E = K*sum_vrel*dt_code
print(f'pairs_pred = {pairs_pred}')
print(f'mp = {mp_g:.4e} g   K*sum_vrel = {K*sum_vrel:.6e} events per code-dt')

# ---------- parse log: per coarse step (a, dt, scattered, pairs) ----------
steps=[]
a=dt=None; import io
for line in open(logfile):
    m = re.search(r'Fine step=\s*(\d+).*?dt=\s*([\dEe\.\+\-]+)\s+a=\s*([\dEe\.\+\-]+)',line)
    if m:
        steps.append(dict(n=int(m.group(1)),dt=float(m.group(2)),a=float(m.group(3)),sc=None,pr=None))
    m = re.search(r'SIDM level\s*\d+: scattered=\s*(\d+)\s*pairs=\s*(\d+)',line)
    if m and steps and steps[-1]['sc'] is None:
        steps[-1]['sc']=int(m.group(1)); steps[-1]['pr']=int(m.group(2))
# SIDM line appears before its Fine step summary? attach to nearest step:
# fallback: pair sequentially if above left gaps
sc_list=[(s['n'],s['a'],s['dt'],s['sc'],s['pr']) for s in steps if s['sc'] is not None]

# steps closest to the snapshot aexp
near = [s for s in sc_list if abs(s[1]-aexp)/aexp < 0.03 and s[1]>=aexp*0.995]
if not near:
    near = sorted(sc_list,key=lambda s:abs(s[1]-aexp))[:6]
near = near[:6]
print(f'\n step      a        dt_code    measured_sc  pred_sc   meas_pairs  pred_pairs')
tot_m=tot_p=0.0
for n,aa,dtc,sc,pr in near:
    pred = K*sum_vrel*dtc*(aexp/aa)**2
    tot_m+=sc; tot_p+=pred
    print(f'{n:5d} {aa:9.5f} {dtc:11.4e} {sc:9d} {pred:10.1f} {pr:11d} {pairs_pred:10d}')
print(f'\nTOTAL measured/predicted = {tot_m}/{tot_p:.1f} = {tot_m/tot_p:.4f}'
      f'   (Poisson ~ +-{1/np.sqrt(max(tot_m,1)):.3f})')
