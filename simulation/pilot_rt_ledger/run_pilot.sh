#!/usr/bin/env bash
#SBATCH --job-name=p3_rt_ledger
#SBATCH --partition=normal
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=pilot_%j.log
#SBATCH --error=pilot_%j.err
#SBATCH --exclusive

set -euo pipefail
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH="/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}"
export UCX_RC_TIMEOUT=30s
export UCX_RC_RETRY_COUNT=7
export UCX_LOG_LEVEL=warn
export UCX_IB_REG_METHODS=direct
export UCX_RNDV_THRESH=65536

[[ -L output_00016 ]] || { echo 'run prepare_pilot.sh first' >&2; exit 1; }
[[ -x ramses_final3d ]] || { echo 'missing ramses_final3d' >&2; exit 1; }
srun --mpi=pmi2 ./ramses_final3d cosmo.nml
