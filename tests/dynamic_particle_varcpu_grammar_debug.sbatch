#!/bin/bash
#SBATCH --job-name=part_varcpu
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/dynamic_particle_varcpu_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/dynamic_particle_varcpu_%j.err

set -euo pipefail
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
run_dir=/home/kjhan/BACKUP/FDM/dynamic_particle_varcpu_${SLURM_JOB_ID}
mkdir -p "$run_dir"
cp /home/kjhan/BACKUP/lagRamses/tests/dynamic_particle_cdm.nml "$run_dir/cosmo.nml"
sed -i "s/'hdf5'/'original'/g" "$run_dir/cosmo.nml"
: > "$run_dir/jobcontrol.txt"
cd "$run_dir"

# Produce an 8-rank checkpoint, then restore it with 16 ranks.  The second
# launch exercises redistribute_particles_by_position instead of only the
# ordinary same-ncpu IC path.
sed -i 's/^nstepmax=.*/nstepmax=1/' cosmo.nml
export OMP_NUM_THREADS=4
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 8 -c 4 --cpu-bind=cores \
  /home/kjhan/BACKUP/lagRamses/bin/ramses_final3d cosmo.nml | tee stage1.log
grep -q 'Run completed' stage1.log
test -f output_00001/COMPLETE
sed -i 's/^nrestart=.*/nrestart=1/; s/^nstepmax=.*/nstepmax=2/' cosmo.nml
export OMP_NUM_THREADS=2
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  /home/kjhan/BACKUP/lagRamses/bin/ramses_final3d cosmo.nml | tee stage2.log
grep -q 'Particle redistribution done' stage2.log
grep -q 'Run completed' stage2.log
