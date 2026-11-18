#!/bin/bash
#SBATCH --job-name=rar-ngc2403-halo
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/fr02_external_halo_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/fr02_external_halo_%j.err
set -euo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
rev=$1
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
out=/gpfs/kjhan/RAR_QUMOND/analysis/ngc2403_external_halo_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
trap 'touch FAILED' ERR
mkdir inputs
data=provenance/results/ngc2403_fr02_boxes_20261118/population
for file in rings.csv masses.csv summary.json results.sha256; do
    git -C "$repo" show "$rev:$data/$file" > "inputs/$file"
done
(cd inputs && awk '$2 == "rings.csv" || $2 == "masses.csv" || $2 == "summary.json"' results.sha256 > used.sha256 &&
    test "$(wc -l < used.sha256)" -eq 3 && sha256sum -c used.sha256)
for file in ngc2403_external_halo.py ngc2403_cpu512_radial_holdout.py; do
    git -C "$repo" show "$rev:provenance/tools/$file" > "$file"
done
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
sha256sum job.sbatch *.py code_revision.txt inputs/* > inputs.sha256
export MPLCONFIGDIR="$out/mplconfig" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python ngc2403_external_halo.py inputs .
sha256sum -c inputs.sha256 --status
sha256sum comparison.pdf comparison.png rings.csv zones.csv halo.csv scenarios.csv > results.sha256
touch COMPLETE
