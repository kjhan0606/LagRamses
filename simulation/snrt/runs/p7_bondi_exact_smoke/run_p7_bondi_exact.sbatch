#!/bin/bash
#SBATCH -J p7_bondi_exact
#SBATCH -p normal
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=10
#SBATCH -t 01:00:00
#SBATCH -o p7_bondi_exact.log
#SBATCH -e p7_bondi_exact.err

set -euo pipefail
export OMP_NUM_THREADS=10
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:${LD_LIBRARY_PATH:-}
cd "$SLURM_SUBMIT_DIR"
mpirun -np 12 ./ramses_final3d cosmo.nml
