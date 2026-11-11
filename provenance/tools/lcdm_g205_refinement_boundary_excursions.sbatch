#!/bin/bash
#SBATCH --job-name=rar-g205-iface
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_iface_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_iface_%j.err
set -Eeuo pipefail
umask 077
module purge

base=/gpfs/kjhan/RAR_QUMOND
src=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_refinement_boundary_excursions.py
adaptive_run="$base/runs/lcdm_g205_z99_z49_paired_tracer_1290972"
fixed_run="$base/runs/lcdm_g205_z99_z49_fixedmap_tracer_1292078"
adaptive_audit="$base/analysis/lcdm_g205_paired_tracer_bounds_1291945"
fixed_audit="$base/analysis/lcdm_g205_fixedmap_tracer_bounds_1292186"
adaptive_snapshot="$adaptive_run/output_00002/data_00002.h5"
fixed_snapshot="$fixed_run/output_00002/data_00002.h5"
adaptive_bounds="$adaptive_audit/paired_tracer_bounds.json"
fixed_bounds="$fixed_audit/paired_tracer_bounds.json"
out="$base/analysis/lcdm_g205_refinement_boundary_excursions_${SLURM_JOB_ID:?}"

# Reconfirm the two completed-run and audit records, then rehash the large HDF5
# inputs on this Grammar compute job before the spatial analysis reads them.
grep -Fxq 'PAIRED_TRACER_EVOLUTION_REVIEW_REQUIRED' "$adaptive_run/outcome.txt"
grep -Fxq 'FIXEDMAP_TRACER_CONTROL_REVIEW_REQUIRED' "$fixed_run/outcome.txt"
grep -Fxq '0' "$adaptive_run/srun_exit_code.txt"
grep -Fxq '0' "$fixed_run/srun_exit_code.txt"
test -f "$adaptive_run/output_00002/COMPLETE"
test -f "$fixed_run/output_00002/COMPLETE"
sha256sum -c "$adaptive_run/artifacts.sha256" --status
sha256sum -c "$fixed_run/artifacts.sha256" --status
grep -Fxq 'PAIRED_TRACER_BOUNDS_REVIEW_REQUIRED' "$adaptive_audit/COMPLETE"
grep -Fxq 'FIXEDMAP_TRACER_BOUNDS_REVIEW_REQUIRED' "$fixed_audit/COMPLETE"
sha256sum -c "$adaptive_audit/artifacts.sha256" --status
sha256sum -c "$fixed_audit/artifacts.sha256" --status
sha256sum -c "$adaptive_audit/snapshots.sha256" --status
sha256sum -c "$fixed_audit/snapshots.sha256" --status
test ! -e "$out"
mkdir "$out"
cp -p "$src" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
cat > "$out/input_provenance.txt" <<EOF
adaptive_run=1290972
adaptive_bounds_audit=1291945
fixed_map_run=1292078
fixed_map_bounds_audit=1292186
adaptive_snapshot=$adaptive_snapshot
fixed_map_snapshot=$fixed_snapshot
analysis_scope=late a=0.02 level-11 gas-tracer cell-face association with level-10 AMR refinement boundaries
calculation_type=post-processing only; no RAMSES evolution and no input mutation
EOF
sha256sum "$adaptive_snapshot" "$fixed_snapshot" > "$out/snapshots.sha256"

export RAR_ADAPTIVE_SNAPSHOT="$adaptive_snapshot"
export RAR_ADAPTIVE_BOUNDS_JSON="$adaptive_bounds"
export RAR_FIXED_SNAPSHOT="$fixed_snapshot"
export RAR_FIXED_BOUNDS_JSON="$fixed_bounds"
export RAR_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1

sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/input_provenance.txt" \
    "$out/snapshots.sha256" "$out/refinement_boundary_excursions.json" \
    "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'REFINEMENT_BOUNDARY_EXCURSIONS_REVIEW_REQUIRED' > "$out/COMPLETE"
