#!/bin/bash
#SBATCH --job-name=rar-fr02-boxes
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/fr02_boxes_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/fr02_boxes_%j.err
set -euo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
rev=$1
shift
if [ "$#" -eq 0 ]; then set -- 64 96 128 192; fi
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
archive=provenance/results/rar_ngc2403_fr02_recovery_20261108
baseline=provenance/results/ngc2403_ml_newton_1284238
out=/gpfs/kjhan/RAR_QUMOND/analysis/ngc2403_fr02_periodic_boxes_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
mkdir inputs
trap 'touch FAILED' ERR
for file in star_profile_code.dat hi_profile_code.dat source_contract.json; do
    git -C "$repo" show "$rev:$archive/gate512/$file" > "inputs/$file"
done
git -C "$repo" show "$rev:$archive/comparison/results/comparison.csv" > inputs/heldout.csv
git -C "$repo" show "$rev:$baseline/comparison.csv" > inputs/compiled.csv
git -C "$repo" show "$rev:$baseline/summary.json" > inputs/units.json
git -C "$repo" show "$rev:provenance/tools/ngc2403_fr02_periodic_boxes.py" > compute.py
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
sha256sum compute.py job.sbatch code_revision.txt inputs/* > inputs.sha256
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python -c 'import scipy,numpy; print(scipy.__version__, numpy.__version__)' > versions.txt
/home/kjhan/miniconda3/bin/python compute.py inputs . --workers "$SLURM_CPUS_PER_TASK" --boxes "$@" > analysis.log
sha256sum -c inputs.sha256 --status
sha256sum grid.json rings_*.csv analysis.log versions.txt > results.sha256
touch COMPLETE
