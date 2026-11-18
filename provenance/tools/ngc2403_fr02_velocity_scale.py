#!/usr/bin/env python3
#SBATCH --job-name=rar-fr02-vscale
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem=512M
#SBATCH --time=00:03:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/fr02_velocity_scale_%j.csv
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/fr02_velocity_scale_%j.err
"""Compare the velocity normalizations required by the inner and outer rings."""
import argparse
import csv
import hashlib
import math
from pathlib import Path
import sys


def measure(population_run):
    curves = population_run / 'rings.csv'
    sha, = [row.split()[0] for row in (population_run/'results.sha256').read_text().splitlines()
            if row.split()[1] == 'rings.csv']
    if not (population_run/'COMPLETE').is_file() or hashlib.sha256(curves.read_bytes()).hexdigest() != sha:
        raise ValueError('incomplete or changed rotation curves')
    with curves.open(newline='') as f:
        rows = list(csv.DictReader(f))
    writer = csv.writer(sys.stdout, lineterminator='\n')
    writer.writerow(['stellar_ml', 'n_inner', 'n_outer', 'scale_inner', 'scale_outer',
                     'outer_over_inner_scale', 'inner_minus_outer_mean_log_g_ratio_dex'])
    for ml in sorted({float(row['stellar_ml']) for row in rows}):
        rings = [row for row in rows if float(row['stellar_ml']) == ml]
        inner = [row for row in rings if 1 <= float(row['radius_kpc']) < 5]
        outer = [row for row in rings if 10 <= float(row['radius_kpc']) <= 20.87]
        if (len(rings), len(inner), len(outer)) != (64, 32, 22):
            raise ValueError('SPARC radial sampling changed')
        scales, ratios = [], []
        for region in (inner, outer):
            vq = [float(row['vqumond_km_s']) for row in region]
            vo = [float(row['vobs_km_s']) for row in region]
            scales.append(sum(vo) / sum(vq))
            ratios.append(sum(2 * math.log10(q/o) for q, o in zip(vq, vo)) / len(vo))
        writer.writerow([ml, len(inner), len(outer), *scales, scales[1] / scales[0],
                         ratios[0] - ratios[1]])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('population_run', type=Path)
    measure(p.parse_args().population_run)
