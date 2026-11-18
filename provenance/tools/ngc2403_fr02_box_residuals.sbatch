#!/bin/bash
#SBATCH --job-name=rar-fr02-box-radii
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/fr02_box_radii_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/fr02_box_radii_%j.err
set -euo pipefail
umask 077
project=/home/kjhan/BackUP/rar
repo=$project/build/rar_qumond_amr_history_20261107
analysis=/gpfs/kjhan/RAR_QUMOND/analysis
rev=$1
shift
test "$#" -gt 0
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
out=$analysis/ngc2403_fr02_box_radii_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
trap 'touch FAILED' ERR
git -C "$repo" show "$rev:provenance/tools/ngc2403_fr02_box_residuals.py" > report.py
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
sha256sum "$project/data/sparc/SPARC_Lelli2016c.mrt" > catalogue.sha256
awk '$1 == "NGC2403"' "$project/data/sparc/SPARC_Lelli2016c.mrt" > sparc_row.txt
roots=()
for id in "$@"; do
    [[ $id =~ ^[0-9]+$ ]]
    run=$analysis/ngc2403_fr02_periodic_boxes_$id
    (cd "$run" && test -f COMPLETE && sha256sum -c inputs.sha256 --status && sha256sum -c results.sha256 --status)
    roots+=("$run")
done
printf '%s\n' "$@" > upstream_job_ids.txt
sha256sum report.py job.sbatch sparc_row.txt catalogue.sha256 upstream_job_ids.txt code_revision.txt > inputs.sha256
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$out/mplconfig"
/home/kjhan/miniconda3/bin/python report.py . sparc_row.txt "${roots[@]}"
sha256sum -c inputs.sha256 --status
sha256sum residuals.pdf residuals.png zones.csv rings.csv summary.json > results.sha256
touch COMPLETE
