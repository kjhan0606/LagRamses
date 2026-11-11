#!/bin/bash
#SBATCH --job-name=rar-g205-refine-edge
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_refine_edge_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_refine_edge_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_tracer_refine_edge.py
early="$base/runs/lcdm_g205_z99_cover_snap_1290834/output_00001/data_00001.h5"
late="$base/runs/lcdm_g205_z49_adapt_1290886/output_00001/data_00001.h5"
run="$base/runs/lcdm_g205_z49_adapt_1290886"
out="$base/analysis/lcdm_g205_tracer_refine_edge_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$run/output_00001/COMPLETE"
sha256sum -c "$run/artifacts.sha256" --status
sha256sum -c "$base/analysis/lcdm_g205_early_proximity_1290874/snapshot.sha256" --status
test "$(grep -Eic 'ivar_refine[[:space:]]*=[[:space:]]*18' "$run/effective.nml")" -eq 1
grep -Eiq 'var_cut_refine[[:space:]]*=[[:space:]]*0\.5' "$run/effective.nml"
mkdir "$out"
sha256sum "$early" "$late" > "$out/snapshots.sha256"
cp -p "$source" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_EARLY_SNAPSHOT="$early" RAR_LATE_SNAPSHOT="$late" RAR_OUTPUT="$out"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1
sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/tracer_refine_edge.json" "$out/snapshots.sha256" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'TRACER_REFINE_EDGE_DIAGNOSTIC_REVIEW_REQUIRED' > "$out/COMPLETE"
