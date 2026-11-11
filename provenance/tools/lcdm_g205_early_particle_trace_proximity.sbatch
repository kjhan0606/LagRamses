#!/bin/bash
#SBATCH --job-name=rar-g205-proximity
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_proximity_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_proximity_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
run="$base/runs/lcdm_g205_z99_cover_snap_1290834"
trace="$base/initial_conditions/lcdm_z0_group205_r3_ellipsoid_l11_1273554"
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_early_particle_trace_proximity.py
out="$base/analysis/lcdm_g205_early_proximity_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$run/output_00001/COMPLETE"
test -f "$trace/GEOMETRY_AUDIT_COMPLETE"
sha256sum -c "$run/artifacts.sha256" --status
mkdir "$out"
awk '$2=="./candidate_traces.npz" {print}' "$trace/artifacts.sha256" > "$out/trace.sha256"
test "$(awk 'END {print NR}' "$out/trace.sha256")" -eq 1
(cd "$trace" && sha256sum -c "$out/trace.sha256" --status)
cp -p "$source" "$out/proximity.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_SNAPSHOT="$run/output_00001/data_00001.h5" RAR_PARENT_TRACE="$trace/candidate_traces.npz" RAR_ANALYSIS_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sha256sum "$RAR_SNAPSHOT" > "$out/snapshot.sha256"
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/proximity.py" > "$out/analysis.log" 2>&1
sha256sum "$out/proximity.py" "$out/job.sbatch" "$out/trace.sha256" "$out/snapshot.sha256" "$out/proximity.json" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'EARLY_EULERIAN_PROXIMITY_ONLY_REVIEW_REQUIRED' > "$out/COMPLETE"
