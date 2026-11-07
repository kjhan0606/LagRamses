#!/usr/bin/env python3
"""Render the saved native-disc CSV without repeating field calculations."""
import argparse,csv
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('out',type=Path);a=p.parse_args()
with (a.out/'rotation_rar.csv').open() as f:
 rows=list(csv.DictReader(f))
columns={key:np.array([float(row[key]) for row in rows]) for key in rows[0]}
modes=('newtonian_baryon','newtonian_halo','qumond_baryon')
changes={mode:columns['fractional_force_resolution_change_'+mode] for mode in modes}
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11})
fig,axes=plt.subplots(1,3,figsize=(15,4.5));colors={'newtonian_baryon':'0.4','newtonian_halo':'tab:blue','qumond_baryon':'tab:red'};labels={'newtonian_baryon':'Baryons (diagnostic)','newtonian_halo':'Baryons + CDM halo','qumond_baryon':'QUMOND (baryons)'}
for mode in modes:
 axes[0].plot(columns['radius_kpc'],columns['v_'+mode+'_kms'],color=colors[mode],label=labels[mode])
 axes[2].plot(columns['radius_kpc'],100*changes[mode],color=colors[mode])
axes[0].set(xlabel='R [kpc]',ylabel='Circular speed [km/s]',ylim=(35,125));axes[0].legend(loc='center right',bbox_to_anchor=(.99,.57),fontsize=11)
axes[1].loglog(columns['gbar_cgs']/100,columns['g_lcdm_cgs']/100,color='tab:blue',label='CDM halo')
axes[1].loglog(columns['gbar_cgs']/100,columns['g_qumond_cgs']/100,color='tab:red',label='QUMOND')
axes[1].loglog(columns['gbar_cgs']/100,columns['gbar_cgs']/100,':',color='0.4',label='Baryons only')
axes[1].set(xlabel=r'$g_{\rm bar}$ [m/s$^2$]',ylabel=r'$g_{\rm tot}$ [m/s$^2$]');axes[1].legend(loc='lower right',fontsize=11)
axes[2].axhline(0,color='0.6',ls=':');axes[2].set(xlabel='R [kpc]',ylabel='Force change: fine/coarse - 1 [%]',ylim=(-.1,4.0))
for ax,title in zip(axes,('Prescribed disc; periodic field','Radial acceleration relation','Resolution comparison')):
 ax.text(.04,.96,title,transform=ax.transAxes,va='top',fontsize=11);ax.grid(alpha=.15)
fig.tight_layout();fig.savefig(a.out/'rotation_rar.pdf');fig.savefig(a.out/'rotation_rar.png',dpi=150);plt.close(fig)
