#!/bin/bash
#SBATCH --job-name=rar-g205-z99-geom
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_geom_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_z99_geom_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l11_z99_diag_1290084"
parent="$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246"
trace="$base/initial_conditions/lcdm_z0_group205_r3_ellipsoid_l11_1273554"
audit=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_z99_geometry_audit.py
out="$base/analysis/g205_z99_geometry_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$ic/COMPLETE"
test -f "$parent/COMPLETE"
test -f "$trace/GEOMETRY_AUDIT_COMPLETE"
mkdir "$out"
cp -p "$audit" "$out/audit.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' "$ic" "$parent" "$trace" > "$out/input_paths.txt"
sha256sum -c "$ic/artifacts.sha256" --status
sha256sum -c "$ic/inputs.sha256" --status
sha256sum -c "$ic/ic_fields.sha256" --status
awk -v path="$parent/ics_lcdm/level_008/ic_deltab" '$2==path {print}' "$parent/artifacts.sha256" > "$out/parent_input.sha256"
test "$(awk 'END {print NR}' "$out/parent_input.sha256")" -eq 1
(cd "$parent" && sha256sum -c "$out/parent_input.sha256" --status)
awk '$2=="./summary.json" || $2=="./candidate_traces.npz" {print}' "$trace/artifacts.sha256" > "$out/trace_inputs.sha256"
test "$(awk 'END {print NR}' "$out/trace_inputs.sha256")" -eq 2
(cd "$trace" && sha256sum -c "$out/trace_inputs.sha256" --status)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$ic" "$parent" "$trace" "$out" > "$out/audit.log"
sha256sum "$out/audit.py" "$out/job.sbatch" "$out/input_paths.txt" "$out/parent_input.sha256" "$out/trace_inputs.sha256" "$out/geometry_audit.json" "$out/audit.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'Z99_GEOMETRY_AUDIT_COMPLETE_NO_EVOLUTION' > "$out/COMPLETE"
