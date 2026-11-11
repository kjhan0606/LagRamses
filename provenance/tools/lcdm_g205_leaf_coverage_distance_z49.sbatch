#!/bin/bash
#SBATCH --job-name=rar-g205-z49-dist
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_z49_dist_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_z49_dist_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
snapshot="$base/runs/lcdm_g205_z99_z49_paired_tracer_1290972/output_00002/data_00002.h5"
coverage="$base/analysis/lcdm_g205_leaf_coverage_z49_1292344"
source="$project/provenance/tools/lcdm_g205_leaf_coverage_distance_z49.py"
root="$base/analysis/lcdm_g205_leaf_coverage_distance_z49_${SLURM_JOB_ID:?}"

module purge
test -f "$coverage/COMPLETE"
test "$(cat "$coverage/COMPLETE")" = "Z49_LEVEL11_PARTICLE_BOX_COVERAGE_REVIEW_REQUIRED"
test -s "$coverage/leaf_coverage_z49.json"
(cd "$coverage" && sha256sum -c artifacts.sha256 --status)
grep -F '"outside_level11_box": 8326' "$coverage/leaf_coverage_z49.json" >/dev/null
test -s "$snapshot"
test ! -e "$root"
mkdir "$root"
cp -p "$source" "$root/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
git -C "$project" rev-parse HEAD > "$root/source_revision.txt"
cp -p "$coverage/adaptive_snapshot.sha256" "$root/snapshot.sha256"
(cd / && sha256sum -c "$root/snapshot.sha256")
snapshot_sha256=$(awk 'NR==1 {print $1}' "$root/snapshot.sha256")

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$root/analysis.py" \
    --snapshot "$snapshot" \
    --snapshot-sha256 "$snapshot_sha256" \
    --expected-outside 8326 \
    --output "$root/mismatch_distance_z49.json" > "$root/analysis.log" 2>&1

sha256sum "$root/analysis.py" "$root/job.sbatch" "$root/source_revision.txt" \
    "$root/snapshot.sha256" "$root/mismatch_distance_z49.json" \
    "$root/analysis.log" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'Z49_LEVEL11_PARTICLE_BOX_DISTANCE_REVIEW_REQUIRED' > "$root/COMPLETE"
