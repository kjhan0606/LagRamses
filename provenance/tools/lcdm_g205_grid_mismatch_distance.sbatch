#!/bin/bash
#SBATCH --job-name=rar-g205-boundary
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_boundary_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_boundary_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
run="$base/runs/lcdm_g205_z99_cover_snap_1290834"
previous="$base/analysis/lcdm_g205_early_proximity_1290874"
leaf="$base/analysis/lcdm_g205_leaf_coverage_1290893"
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_grid_mismatch_distance.py
out="$base/analysis/lcdm_g205_boundary_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$run/output_00001/COMPLETE" && test -f "$previous/COMPLETE" && test -f "$leaf/COMPLETE"
sha256sum -c "$previous/artifacts.sha256" --status
sha256sum -c "$previous/snapshot.sha256" --status
sha256sum -c "$leaf/artifacts.sha256" --status
mkdir "$out"
cp -p "$source" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_SNAPSHOT="$run/output_00001/data_00001.h5" RAR_ANALYSIS_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1
sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/mismatch_distance.json" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'BOUNDARY_DISTANCE_UPPER_BOUND_ONLY' > "$out/COMPLETE"
