#!/bin/bash
#SBATCH --job-name=rar-g205-ic-scalars
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_ic_scalars_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_ic_scalars_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l11_1286049"
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/g205_ic_scalar_audit.py
out="$base/analysis/g205_ic_scalars_${SLURM_JOB_ID:?}"
test -f "$ic/GENERATION_ONLY"
test -f "$ic/GEOMETRY_AUDIT_COMPLETE"
test ! -e "$out"
mkdir "$out"
cp -p "$source" "$out/audit.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' "$ic" > "$out/input_ic_path.txt"
awk '$2 ~ /^\.\/ics\/level_00(8|9)\/ic_(deltab|pvar_00001|refmap)$/ || $2 ~ /^\.\/ics\/level_01(0|1)\/ic_(deltab|pvar_00001|refmap)$/ {print}' "$ic/artifacts.sha256" > "$out/selected_inputs.sha256"
test "$(awk 'END {print NR}' "$out/selected_inputs.sha256")" -eq 12
(cd "$ic" && sha256sum -c "$out/selected_inputs.sha256" --status)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$ic" "$out/summary.json"
sha256sum "$out/audit.py" "$out/job.sbatch" "$out/input_ic_path.txt" "$out/selected_inputs.sha256" "$out/summary.json" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'IC_SCALAR_AUDIT_COMPLETE' > "$out/COMPLETE"
