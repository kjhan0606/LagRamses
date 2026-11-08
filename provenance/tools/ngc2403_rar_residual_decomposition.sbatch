#!/bin/bash
#SBATCH --job-name=rar-residual-split
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_residual_split_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_residual_split_%j.err
set -Eeuo pipefail
umask 077
project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
repo=$project/build/rar_qumond_amr_history_20261107
revision=${RAR_RESIDUAL_SOURCE_REVISION:?}
input=$base/analysis/ngc2403_compiled512_cdm_nfw_holdout_1279415
out=$base/analysis/ngc2403_rar_residual_split_${SLURM_JOB_ID:?}
report=$project/reports/ngc2403_rar_residual_split_${SLURM_JOB_ID:?}
test -f "$input/COMPLETE"
test ! -e "$out" && test ! -e "$report"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 33554432
mkdir "$out"
git -C "$repo" show "$revision:provenance/tools/ngc2403_rar_residual_decomposition.py" > "$out/compute.py"
cp "$input/holdout.csv" "$out/holdout.csv"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
printf 'purpose=postprocessing 64 saved radii; full_snapshots=0; reserved_bytes=33554432; free_bytes=%s\n' "$free_bytes" > "$out/output_policy.txt"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/compute.py" "$out/holdout.csv" "$out" > "$out/compute.log" 2>&1
(cd "$out" && sha256sum compute.py holdout.csv decomposition.csv summary.json job.sbatch source_revision.txt output_policy.txt compute.log > artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'SAVED_RESIDUAL_IDENTITY_VERIFIED_NOT_CAUSAL_ATTRIBUTION' > "$out/COMPLETE"
mkdir "$report"
cp "$out"/* "$report/"
(cd "$report" && sha256sum -c artifacts.sha256)
