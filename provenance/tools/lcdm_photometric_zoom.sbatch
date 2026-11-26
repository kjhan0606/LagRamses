#!/bin/bash
#SBATCH --job-name=rar-photo-trace
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/photometric_zoom_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/photometric_zoom_%j.err
set -euo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
rev=$1
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
old=$base/analysis/lcdm_z0_candidate_trace_1272106
parent=$base/runs/lcdm_z0_parent_direct_1271100
hop=$base/analysis/hop_lcdm_z0_1271118
so=$base/analysis/lcdm_z0_so_m200_1271672
map=$base/analysis/dmo_lagrangian_id_1269366
snap=$parent/output_00003/data_00003.h5
test -f "$old/COMPLETE" && test -f "$parent/COMPLETE"
test "$(readlink -f "$parent/ic")" = "$(readlink -f "$base/initial_conditions/dmo_parent_z49_n256_1269226/ics/level_008")"
out=$base/analysis/lcdm_photometric_zoom_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
trap 'touch FAILED' ERR
(cd "$old" && sha256sum -c artifacts.sha256 --status && sha256sum -c upstream_inputs.sha256 --status)
(cd "$map" && sha256sum -c artifacts.sha256 --status)
cp "$old/upstream_inputs.sha256" parent_inputs.sha256
sha256sum "$map/audit.json" >> parent_inputs.sha256
for name in lcdm_z0_candidate_lagrangian_trace lcdm_z0_center_sensitivity lcdm_photometric_zoom_points; do
    git -C "$repo" show "$rev:provenance/tools/$name.py" > "$name.py"
done
src=provenance/results/ngc2403_external_halo_1324385
git -C "$repo" show "$rev:$src/halo.csv" > halo.csv
git -C "$repo" show "$rev:$src/results.sha256" | awk '$2=="halo.csv"' > halo.sha256
sha256sum -c halo.sha256
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
sha256sum *.py job.sbatch code_revision.txt halo.csv halo.sha256 parent_inputs.sha256 > inputs.sha256
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
py=/home/kjhan/miniconda3/bin/python3
"$py" lcdm_z0_candidate_lagrangian_trace.py "$snap" "$hop/catalog.csv" "$hop/groups.tag" \
    "$hop/converted/id_by_row.npy" "$so/candidate_m200c.csv" "$so/summary.json" \
    "$map/cell_xyz_by_id.npy" "$map/audit.json" "$out" halo.csv 3.16 > trace.log
"$py" lcdm_z0_center_sensitivity.py "$snap" "$hop/catalog.csv" "$hop/groups.tag" \
    "$so/candidate_m200c.csv" summary.json "$out" > centers.log
"$py" lcdm_photometric_zoom_points.py "$out"
sha256sum -c inputs.sha256 --status
shopt -s nullglob
sha256sum candidate_traces.npz summary.json center_sensitivity.json shortlist.csv points_g*.txt > results.sha256
touch COMPLETE
