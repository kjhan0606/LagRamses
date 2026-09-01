#!/bin/bash
#SBATCH -J p7_sinkprops
#SBATCH -p normal
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=10
#SBATCH -t 01:00:00
#SBATCH -o p7_sinkprops.log
#SBATCH -e p7_sinkprops.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
export OMP_NUM_THREADS=10
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:${LD_LIBRARY_PATH:-}
cd "$SLURM_SUBMIT_DIR"
mpirun -np 12 ./ramses_final3d cosmo.nml
