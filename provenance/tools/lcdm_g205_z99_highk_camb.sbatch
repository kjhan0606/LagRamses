#!/bin/bash
#SBATCH --job-name=rar-g205-highk
#SBATCH --partition=normal
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_highk_%A_%a.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_highk_%A_%a.err
set -Eeuo pipefail
umask 077
module purge
base=/gpfs/kjhan/RAR_QUMOND
source=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107/provenance/tools/lcdm_g205_z99_highk_camb.py
out="$base/transfer_probes/lcdm_g205_z99_highk_${SLURM_ARRAY_JOB_ID:?}_${SLURM_ARRAY_TASK_ID:?}"
case "$SLURM_ARRAY_TASK_ID" in 0) sampling=100;; 1) sampling=200;; *) exit 1;; esac
test ! -e "$out"
test "$(df -B1 "$base" | awk 'END {print $4}')" -ge 10737418240
mkdir "$out"
cp -p "$source" "$out/generate.py"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
export RAR_CAMB_OUTPUT="$out" RAR_CAMB_SAMPLING="$sampling"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/generate.py" > "$out/run.log" 2>&1
sha256sum "$out/generate.py" "$out/job.sbatch" "$out/transfer_z99.dat" "$out/settings.json" "$out/run.log" > "$out/artifacts.sha256"
sha256sum -c "$out/artifacts.sha256" --status
printf '%s\n' 'HIGHK_TRANSFER_GENERATED_NOT_CONVERGENCE_VALIDATED' > "$out/COMPLETE"
