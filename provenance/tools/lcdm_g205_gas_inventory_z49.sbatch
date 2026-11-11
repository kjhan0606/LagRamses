#!/bin/bash
#SBATCH --job-name=rar-g205-gas-z49
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_gas_z49_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_gas_z49_%j.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
run="$base/runs/lcdm_g205_z49_adapt_1290886"
previous="$base/analysis/lcdm_g205_gas_inventory_1290910"
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_gas_inventory_z49.py
out="$base/analysis/lcdm_g205_gas_z49_${SLURM_JOB_ID:?}"
test ! -e "$out"
test -f "$run/output_00001/COMPLETE" && test -f "$previous/COMPLETE"
sha256sum -c "$run/artifacts.sha256" --status
sha256sum -c "$previous/artifacts.sha256" --status

mkdir "$out"
sha256sum "$run/output_00001/data_00001.h5" > "$out/snapshot.sha256"
cp -p "$source" "$out/analysis.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_SNAPSHOT="$run/output_00001/data_00001.h5" RAR_ANALYSIS_OUTPUT="$out" RAR_BASELINE_INVENTORY="$previous/gas_tracer_inventory.json"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/analysis.py" > "$out/analysis.log" 2>&1
sha256sum "$out/analysis.py" "$out/job.sbatch" "$out/gas_tracer_inventory.json" "$out/gas_tracer_comparison.json" "$out/snapshot.sha256" "$out/analysis.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'TWO_RUN_GAS_INVENTORY_REVIEW_REQUIRED' > "$out/COMPLETE"
