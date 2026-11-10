import csv,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
source=Path(sys.argv[1]); out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
with source.open() as h: data=list(csv.DictReader(h))
r=np.array([float(x['radius_kpc']) for x in data]);obs=np.array([float(x['vobs']) for x in data]);err=np.array([float(x['vobs_err']) for x in data]);
fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
ax.errorbar(r,obs,yerr=err,fmt='.',color='black',ms=3,lw=.6,label='SPARC observed')
for key,ls,c,label in [('qumond_ml05','-','#247f7f','QUMOND M/L=0.5'),('qumond_ml07','-','#7951a0','QUMOND M/L=0.7'),('newton_ml05','--','#247f7f','Newton M/L=0.5'),('newton_ml07','--','#7951a0','Newton M/L=0.7')]:
 v=np.array([float(x[key]) for x in data]);ax.plot(r,v,ls=ls,color=c,label=label);bx.plot(r,v-obs,ls=ls,color=c,label=label)
for a,title in [(ax,'NGC 2403, FR02 gas'),(bx,'Residuals at SPARC radii')]:
 a.text(.03,.97,title,transform=a.transAxes,va='top',fontsize=10);a.set_xlabel('Radius (kpc)');a.tick_params(labelsize=10)
ax.set_ylabel('Circular speed (km/s)');ax.set_ylim(25,190);bx.set_ylabel('Model minus observed (km/s)');bx.set_ylim(-100,62);bx.axhline(0,color='0.65',lw=.8)
ax.legend(loc='upper right',frameon=False,fontsize=9);bx.legend(loc='upper right',frameon=False,fontsize=9)
fig.savefig(out/'comparison.pdf');fig.savefig(out/'comparison.png',dpi=180)
