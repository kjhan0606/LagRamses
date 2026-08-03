#!/bin/bash
#SBATCH --job-name=dynpart_cdm
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/dynamic_particle_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/dynamic_particle_%j.err

set -euo pipefail
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_dir=/home/kjhan/BACKUP/FDM/dynamic_particle_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp /home/kjhan/BACKUP/lagRamses/tests/dynamic_particle_cdm.nml "$run_dir/cosmo.nml"
: > "$run_dir/jobcontrol.txt"
cd "$run_dir"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  /home/kjhan/BACKUP/lagRamses/bin/ramses_final3d cosmo.nml | tee run.log
grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
