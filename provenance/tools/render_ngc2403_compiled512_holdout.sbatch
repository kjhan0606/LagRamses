#!/bin/bash
# Rerender only existing tabulated results; no force read or NFW refit.
#SBATCH --job-name=rar-holdout-figure
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_holdout_figure_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_holdout_figure_%j.err
set -Eeuo pipefail
umask 077
project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
run=$base/analysis/ngc2403_compiled512_cdm_nfw_holdout_${RAR_HOLDOUT_RUN_ID:?}
out=$base/analysis/ngc2403_holdout_figure_${SLURM_JOB_ID:?}
report=$project/reports/ngc2403_holdout_figure_${SLURM_JOB_ID:?}
revision=${RAR_HOLDOUT_PLOT_REVISION:?}
repo=$project/build/rar_qumond_amr_history_20261107
test -f "$run/COMPLETE"
mkdir "$out"
for name in ngc2403_compiled512_cdm_nfw_holdout.py ngc2403_lagr_rar.py ngc2403_cpu512_radial_holdout.py;do
 git -C "$repo" show "$revision:provenance/tools/$name" > "$out/$name"
done
cp "$run/holdout.csv" "$out/holdout.csv"
printf '%s\n' "$revision" > "$out/source_revision.txt"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH PYTHONPATH="$out" /home/kjhan/miniconda3/bin/python3 - "$out" <<'PY'
import csv,sys
from pathlib import Path
from ngc2403_compiled512_cdm_nfw_holdout import save_plot
p=Path(sys.argv[1])
with (p/'holdout.csv').open() as f:rows=sorted(csv.DictReader(f),key=lambda r:float(r['radius_kpc']))
def values(key):return [float(r[key]) for r in rows]
observed=[{'gbar_m_s2':r['gbar_sparc_m_s2']} for r in rows]
save_plot(p,observed,values('radius_kpc'),values('gbar_compiled_newtonian_m_s2'),
          values('gqumond_compiled_m_s2'),values('gnfw_heldout_m_s2'),values('gobs_m_s2'))
PY
(cd "$out" && sha256sum *.py *.csv *.pdf *.png source_revision.txt job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'EXISTING_HOLDOUT_RERENDERED_NO_REFIT' > "$out/COMPLETE"
mkdir "$report"
cp "$out"/*.py "$out"/*.csv "$out"/*.pdf "$out"/*.png \
 "$out/source_revision.txt" "$out/job.sbatch" "$out/artifacts.sha256" "$out/COMPLETE" "$report/"
(cd "$report" && sha256sum -c artifacts.sha256)
