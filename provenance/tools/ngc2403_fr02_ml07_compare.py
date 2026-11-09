"""Compare a recomputed FR02 M/L=0.7 field to the saved matched 0.5 field."""
import csv,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ngc2403_lagr_rar import read_uniform,M_PER_KPC
from ngc2403_compiled512_cdm_nfw_holdout import sample_radial_acceleration
run,baseline,out=map(Path,sys.argv[1:])
assert (run/'COMPLETE').is_file()
out.mkdir(exist_ok=False)
with baseline.open() as f: rows=list(csv.DictReader(f))
r=np.array([float(x['radius_kpc']) for x in rows]); obs=np.array([float(x['vobs_km_s']) for x in rows]); err=np.array([float(x['vobs_error_km_s']) for x in rows]); old=np.array([float(x['vqumond_fr02_km_s']) for x in rows])
f,meta=read_uniform(run/'qumond')
assert meta['nmesh']==512 and np.isclose(meta['dx_kpc'],.125)
g=np.array([sample_radial_acceleration(f,x,meta)['median'] for x in r]); del f
v=np.sqrt(g*r*M_PER_KPC)/1000
assert np.all(np.isfinite(v)) and len(v)==64
with (out/'comparison.csv').open('w') as stream:
 w=csv.writer(stream); w.writerow(['radius_kpc','vobs_km_s','error_km_s','qumond_ml05_km_s','qumond_ml07_km_s']); w.writerows(zip(r,obs,err,old,v))
stats={k:{'velocity_rms_km_s':float(np.sqrt(np.mean((a-obs)**2))),'mean_velocity_residual_km_s':float(np.mean(a-obs))} for k,a in [('ml05',old),('ml07',v)]}
(out/'summary.json').write_text(json.dumps({'statistics':stats,'force_run':str(run),'baseline':str(baseline),'limits':'Fixed FR02 gas, periodic 64 kpc box, uniform 512 cubed; not a fitted M/L or cosmological evolution.'},indent=2))
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
ax.errorbar(r,obs,err,fmt='.',color='black',label='SPARC observed')
for a,label,color in [(old,'QUMOND M/L=0.5','#247f7f'),(v,'QUMOND M/L=0.7','#7951a0')]:
 ax.plot(r,a,label=label,color=color);bx.plot(r,a-obs,label=label,color=color)
ax.set(xlabel='Radius (kpc)',ylabel='Circular speed (km/s)');bx.set(xlabel='Radius (kpc)',ylabel='Model minus observed (km/s)');bx.axhline(0,color='.6',lw=.8)
for a,title in [(ax,'NGC 2403: fixed FR02 gas'),(bx,'Matched-radius residuals')]:
 a.text(.03,.97,title,transform=a.transAxes,va='top',fontsize=10);a.legend(loc='lower right',fontsize=10,frameon=False)
fig.savefig(out/'ml_comparison.pdf');fig.savefig(out/'ml_comparison.png',dpi=180)
