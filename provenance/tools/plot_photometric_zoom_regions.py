#!/home/kjhan/miniconda3/bin/python3
#SBATCH --job-name=rar-photo-regions
#SBATCH --partition=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=00:03:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/photometric_regions_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/photometric_regions_%j.err
"""Show the mass-selected pilot regions at the same comoving scale."""
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np


def plot(run, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for line in (run/'results.sha256').read_text().splitlines():
        sha, name = line.split()
        if hashlib.sha256((run/name).read_bytes()).hexdigest() != sha:
            raise ValueError(f'changed trace product {name}')
    with (run/'shortlist.csv').open(newline='') as f:
        groups = [int(row['group_id']) for row in csv.DictReader(f) if row['point_file']]
    candidates = {r['group_id']:r for r in json.loads((run/'summary.json').read_text())['candidates']}
    fig, axes = plt.subplots(len(groups),3,figsize=(7.8,2.65*len(groups)),sharex=True,sharey=True,
                             squeeze=False,gridspec_kw={'wspace':0,'hspace':0})
    with np.load(run/'candidate_traces.npz') as ar:
        width = max(max(candidates[g]['traces'][1]['periodic_minimum_bbox_width_xyz']) for g in groups)*25.
        for row,g in enumerate(groups):
            obj = candidates[g]['traces'][1]
            center = np.array(obj['periodic_minimum_bbox_start_xyz']) + .5*np.array(obj['periodic_minimum_bbox_width_xyz'])
            cells = ar[f'g{g}_r3_cell_xyz']
            xyz = ((cells+.5)/256 - center + .5) % 1 - .5
            np.testing.assert_allclose(np.ptp(xyz,axis=0),obj['periodic_minimum_bbox_width_xyz'])
            xyz *= 25.
            inner = np.isin(ar[f'g{g}_r3_ids'],ar[f'g{g}_r2_ids'])
            for ax,(i,j),tag in zip(axes[row],[(0,1),(0,2),(1,2)],['xy','xz','yz']):
                ax.scatter(xyz[~inner,i],xyz[~inner,j],s=3,color='#bc693d',alpha=.65,lw=0,
                           label='Additional 3R cells',rasterized=True)
                ax.scatter(xyz[inner,i],xyz[inner,j],s=3,color='#185f83',alpha=.55,lw=0,
                           label='2R cells',rasterized=True)
                ax.text(.04,.96,f'{g}   {tag}',transform=ax.transAxes,va='top')
                ax.set(xlim=(-.7*width,.7*width),ylim=(-.7*width,.7*width),aspect='equal')
    axes[0,-1].legend(loc='lower right',fontsize=9,frameon=False,markerscale=2)
    fig.supxlabel('Comoving offset [Mpc/h]',fontsize=10)
    fig.supylabel('Comoving offset [Mpc/h]',fontsize=10)
    fig.subplots_adjust(left=.092,right=.99,bottom=.11,top=.99)
    fig.savefig(out/'regions.pdf')
    fig.savefig(out/'regions.png',dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    base = Path('/gpfs/kjhan/RAR_QUMOND/analysis')
    run = base/'lcdm_photometric_zoom_1578501'
    out = base/('lcdm_photometric_regions_'+os.environ['SLURM_JOB_ID'])
    out.mkdir()
    os.environ['MPLCONFIGDIR'] = str(out/'mplconfig')
    plot(run,out)
    (out/'code_revision.txt').write_text(sys.argv[1]+'\n')
