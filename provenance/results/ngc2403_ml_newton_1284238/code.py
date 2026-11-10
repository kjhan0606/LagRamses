from pathlib import Path
import csv,sys,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ngc2403_lagr_rar import read_uniform,M_PER_KPC
from ngc2403_compiled512_cdm_nfw_holdout import sample_radial_acceleration
out=Path(sys.argv[1]); out.mkdir(parents=True,exist_ok=False)
base=Path('/home/kjhan/BackUP/rar/reports/ngc2403_fr02_ml07_comparison_1284190/comparison.csv')
with base.open() as f: rows=list(csv.DictReader(f))
r=np.array([float(x['radius_kpc']) for x in rows]); obs=np.array([float(x['vobs_km_s']) for x in rows]); err=np.array([float(x['error_km_s']) for x in rows]); q05=np.array([float(x['qumond_ml05_km_s']) for x in rows]);q07=np.array([float(x['qumond_ml07_km_s']) for x in rows]); assert len(r)==64
paths=[Path('/gpfs/kjhan/RAR_QUMOND/runs/ngc2403_fr02_uniform512_456162'),Path('/gpfs/kjhan/RAR_QUMOND/runs/ngc2403_fr02_ml07_uniform512_457508')]
vs=[];metas=[]
for p in paths:
 assert (p/'COMPLETE').is_file()
 f,m=read_uniform(p/'newtonian'); assert m['nmesh']==512 and np.isclose(m['dx_kpc'],.125)
 g=np.array([sample_radial_acceleration(f,x,m)['median'] for x in r]);del f
 vs.append(np.sqrt(g*r*M_PER_KPC)/1000);metas.append(m)
assert all(metas[0][k]==metas[1][k] for k in ['unit_length_cm','unit_time_s','boxlen_code'])
with (out/'comparison.csv').open('w',newline='') as h:
 w=csv.writer(h);w.writerow(['radius_kpc','vobs','vobs_err','newton_ml05','newton_ml07','qumond_ml05','qumond_ml07']);w.writerows(zip(r,obs,err,*vs,q05,q07))
summary={'points':len(r),'radius_kpc':[float(r.min()),float(r.max())],'qumond_velocity_rms_km_s':{'0.5':float(np.sqrt(np.mean((q05-obs)**2))),'0.7':float(np.sqrt(np.mean((q07-obs)**2)))},'newtonian_fields':metas,'scope':'static prescribed stellar and FR02 gas source; periodic 64 kpc box; not a live halo or cosmological simulation'}
(out/'summary.json').write_text(json.dumps(summary,indent=2))
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
ax.errorbar(r,obs,yerr=err,fmt='.',color='black',ms=3,lw=.6,label='SPARC observed')
for v,ls,c,label in [(q05,'-','#247f7f','QUMOND M/L=0.5'),(q07,'-','#7951a0','QUMOND M/L=0.7'),(vs[0],'--','#247f7f','Newton M/L=0.5'),(vs[1],'--','#7951a0','Newton M/L=0.7')]:
 ax.plot(r,v,ls=ls,color=c,label=label);bx.plot(r,v-obs,ls=ls,color=c,label=label)
for a,title in [(ax,'NGC 2403, FR02 gas'),(bx,'Residuals at SPARC radii')]:
 a.text(.03,.97,title,transform=a.transAxes,va='top',fontsize=10);a.set_xlabel('Radius (kpc)');a.tick_params(labelsize=10)
ax.set_ylabel('Circular speed (km/s)');bx.set_ylabel('Model minus observed (km/s)');bx.axhline(0,color='0.65',lw=.8)
ax.legend(loc='lower right',frameon=False,fontsize=9);bx.legend(loc='lower left',frameon=False,fontsize=9)
fig.savefig(out/'comparison.pdf');fig.savefig(out/'comparison.png',dpi=180)
