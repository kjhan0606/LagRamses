#!/bin/bash
#SBATCH --job-name=rar-g205-highk-cmp
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_highk_compare_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_highk_compare_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_z99_highk_compare.py
out="$base/analysis/lcdm_g205_z99_highk_compare_${SLURM_JOB_ID:?}"
test ! -e "$out"
for i in 0 1; do
    d="$base/transfer_probes/lcdm_g205_z99_highk_1290351_$i"
    test -f "$d/COMPLETE"
    sha256sum -c "$d/artifacts.sha256" --status
done
old="$base/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084/transfer.dat"
printf '%s  %s\n' adcac18fd4e3801c7eb553ac6e50be81e5912c23e4bf01ede8ad18812b03720f "$old" | sha256sum -c --status
mkdir "$out"
cp -p "$source" "$out/compare.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_COMPARE_OUTPUT="$out" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/compare.py" > "$out/compare.log" 2>&1
sha256sum "$out/compare.py" "$out/job.sbatch" "$out/comparison.json" "$out/compare.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'TRANSFER_TABLE_COMPARISON_REQUIRES_SCIENTIFIC_REVIEW' > "$out/COMPLETE"
