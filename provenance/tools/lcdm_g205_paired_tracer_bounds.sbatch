#!/bin/bash
#SBATCH --job-name=rar-g205-paired-bounds
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_paired_bounds_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_paired_bounds_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
src=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_paired_tracer_bounds.py
run="$base/runs/lcdm_g205_z99_z49_paired_tracer_1290972"
early="$run/output_00001/data_00001.h5"
late="$run/output_00002/data_00002.h5"
out="$base/analysis/lcdm_g205_paired_tracer_bounds_${SLURM_JOB_ID:?}"
test ! -e "$out"
grep -Fxq 'PAIRED_TRACER_EVOLUTION_REVIEW_REQUIRED' "$run/outcome.txt"
grep -Fxq '0' "$run/srun_exit_code.txt"
test -f "$run/output_00001/COMPLETE" && test -f "$run/output_00002/COMPLETE"
sha256sum -c "$run/artifacts.sha256" --status
mkdir "$out"
sha256sum "$early" "$late" > "$out/snapshots.sha256"
cp -p "$src" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_EARLY_SNAPSHOT="$early" RAR_LATE_SNAPSHOT="$late" RAR_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1
sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/paired_tracer_bounds.json" "$out/snapshots.sha256" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'PAIRED_TRACER_BOUNDS_REVIEW_REQUIRED' > "$out/COMPLETE"
