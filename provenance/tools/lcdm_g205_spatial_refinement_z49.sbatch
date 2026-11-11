#!/bin/bash
#SBATCH --job-name=rar-g205-z49-spatial
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_z49_spatial_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_z49_spatial_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
adaptive_run="$base/runs/lcdm_g205_z99_z49_paired_tracer_1290972"
fixed_run="$base/runs/lcdm_g205_z99_z49_fixedmap_tracer_1292078"
adaptive_analysis="$base/analysis/lcdm_g205_paired_tracer_bounds_1291945"
fixed_analysis="$base/analysis/lcdm_g205_fixedmap_tracer_bounds_1292186"
adaptive_snapshot="$adaptive_run/output_00002/data_00002.h5"
fixed_snapshot="$fixed_run/output_00002/data_00002.h5"
source="$project/provenance/tools/lcdm_g205_spatial_refinement_z49.py"
root="$base/analysis/lcdm_g205_spatial_refinement_z49_${SLURM_JOB_ID:?}"

module purge
test -f "$adaptive_run/output_00002/COMPLETE"
test -f "$fixed_run/output_00002/COMPLETE"
test -f "$adaptive_analysis/COMPLETE"
test -f "$fixed_analysis/COMPLETE"
test -s "$adaptive_snapshot"
test -s "$fixed_snapshot"
test ! -e "$root"
mkdir "$root"
cp -p "$source" "$root/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
git -C "$project" rev-parse HEAD > "$root/source_revision.txt"

grep -F "$adaptive_snapshot" "$adaptive_analysis/snapshots.sha256" > "$root/adaptive_snapshot.sha256"
grep -F "$fixed_snapshot" "$fixed_analysis/snapshots.sha256" > "$root/fixed_map_snapshot.sha256"
test "$(wc -l < "$root/adaptive_snapshot.sha256")" -eq 1
test "$(wc -l < "$root/fixed_map_snapshot.sha256")" -eq 1
(cd / && sha256sum -c "$root/adaptive_snapshot.sha256")
(cd / && sha256sum -c "$root/fixed_map_snapshot.sha256")
adaptive_sha256=$(awk 'NR==1 {print $1}' "$root/adaptive_snapshot.sha256")
fixed_sha256=$(awk 'NR==1 {print $1}' "$root/fixed_map_snapshot.sha256")

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$root/analysis.py" \
    --adaptive "$adaptive_snapshot" \
    --fixed-map "$fixed_snapshot" \
    --adaptive-sha256 "$adaptive_sha256" \
    --fixed-map-sha256 "$fixed_sha256" \
    --output "$root/spatial_refinement_z49.json" > "$root/analysis.log" 2>&1

sha256sum "$root/analysis.py" "$root/job.sbatch" "$root/source_revision.txt" \
    "$root/adaptive_snapshot.sha256" "$root/fixed_map_snapshot.sha256" \
    "$root/spatial_refinement_z49.json" "$root/analysis.log" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'Z49_FINEST_MASS_SPATIAL_LEVEL_REVIEW_REQUIRED' > "$root/COMPLETE"
