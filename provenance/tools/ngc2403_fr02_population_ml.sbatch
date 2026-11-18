#!/bin/bash
#SBATCH --job-name=rar-fr02-photo-ml
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/fr02_population_ml_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/fr02_population_ml_%j.err
set -euo pipefail
umask 077
project=/home/kjhan/BackUP/rar
repo=$project/build/rar_qumond_amr_history_20261107
analysis=/gpfs/kjhan/RAR_QUMOND/analysis
rev=$1
[[ $2 =~ ^[0-9]+$ ]]
up=$analysis/ngc2403_fr02_periodic_boxes_$2
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
(cd "$up" && test -f COMPLETE && sha256sum --status -c inputs.sha256 && sha256sum --status -c results.sha256)
out=$analysis/ngc2403_fr02_population_ml_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
trap 'touch FAILED' ERR
for file in ngc2403_fr02_periodic_boxes.py ngc2403_fr02_population_ml.py; do
    git -C "$repo" show "$rev:provenance/tools/$file" > "$file"
done
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
printf '%s\n' "$up" > upstream.txt
sha256sum "$up/inputs.sha256" "$up/results.sha256" > upstream.sha256
sha256sum "$project/data/sparc/SPARC_Lelli2016c.mrt" > catalogue.sha256
awk '$1 == "NGC2403"' "$project/data/sparc/SPARC_Lelli2016c.mrt" > sparc_row.txt
sha256sum ngc2403_fr02_*.py job.sbatch code_revision.txt upstream.txt upstream.sha256 catalogue.sha256 sparc_row.txt > inputs.sha256
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLCONFIGDIR="$out/mplconfig"
/home/kjhan/miniconda3/bin/python ngc2403_fr02_population_ml.py "$up" sparc_row.txt . --workers "$SLURM_CPUS_PER_TASK" > analysis.log
sha256sum --status -c inputs.sha256
sha256sum --status -c upstream.sha256
sha256sum masses.csv rings.csv zones.csv summary.json rotation_curves.pdf rotation_curves.png > results.sha256
touch COMPLETE
