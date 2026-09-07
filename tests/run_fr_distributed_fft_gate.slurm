#!/bin/bash
#SBATCH --job-name=p1_fr_fft512
#SBATCH --partition=normal
#SBATCH --exclude=grammar005,grammar006,grammar007,grammar012,grammar023,grammar077,grammar078,grammar081,grammar112
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=320G
#SBATCH --time=03:00:00
set -euo pipefail
release=${FR_DISTRIBUTED_RELEASE:?set audited test root}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset OMP_PLACES I_MPI_PIN_PROCESSOR_LIST I_MPI_FABRICS
cd "$release"
python3 "$release/source/run_fr_distributed_fft_gate.py" verify --root "$release"
printf 'job=%s nodes=%s start=%s\n' "$SLURM_JOB_ID" "$SLURM_JOB_NODELIST" "$(date --iso-8601=seconds)"
srun --ntasks=2 --ntasks-per-node=1 --cpus-per-task=1 hostname
python3 "$release/source/run_fr_distributed_fft_gate.py" run --root "$release"
