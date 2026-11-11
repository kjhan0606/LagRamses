#!/bin/bash
#SBATCH --job-name=rar-g205-ic-compare
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_ic_compare_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_ic_compare_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
old="$base/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084"
new="${RAR_HIGHK_IC:?}"
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_z99_ic_pair_compare.py
out="$base/analysis/lcdm_g205_z99_ic_compare_${SLURM_JOB_ID:?}"
test ! -e "$out"
for d in "$old" "$new"; do
    test -f "$d/COMPLETE"
    sha256sum -c "$d/artifacts.sha256" --status
    sha256sum -c "$d/inputs.sha256" --status
    sha256sum -c "$d/ic_fields.sha256" --status
done
mkdir "$out"
cp -p "$source" "$out/compare.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_HIGHK_IC="$new" RAR_COMPARE_OUTPUT="$out" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/compare.py" > "$out/compare.log" 2>&1
sha256sum "$out/compare.py" "$out/job.sbatch" "$out/field_comparison.json" "$out/compare.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'IC_FIELD_COMPARISON_REQUIRES_REVIEW' > "$out/COMPLETE"
