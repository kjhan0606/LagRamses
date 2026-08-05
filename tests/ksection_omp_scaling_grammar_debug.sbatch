#!/bin/bash
#SBATCH --job-name=ksec_omp
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/ksection_omp_%A_%a.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/ksection_omp_%A_%a.err

set -euo pipefail

threads=(1 2 4)
export OMP_NUM_THREADS=${threads[$SLURM_ARRAY_TASK_ID]}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_dir=/home/kjhan/BACKUP/FDM/ksection_omp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "${run_dir}"
cp "${repo}/tests/dynamic_particle_cdm.nml" "${run_dir}/cosmo.nml"
sed -i "s/^ordering=.*/ordering='ksection'/" "${run_dir}/cosmo.nml"
sed -i 's/^nstepmax=.*/nstepmax=1/' "${run_dir}/cosmo.nml"
sed -i 's/^pic=.*/pic=.false./' "${run_dir}/cosmo.nml"
sed -i 's/^poisson=.*/poisson=.false./' "${run_dir}/cosmo.nml"
: > "${run_dir}/jobcontrol.txt"

echo "Testing ksection OMP_NUM_THREADS=${OMP_NUM_THREADS} with 16 fixed MPI ranks"
cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 4 --cpu-bind=cores \
  "${repo}/bin/ramses_final3d" cosmo.nml | tee run.log

grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
grep -q 'LB compact stages:' run.log
grep -q 'cmp_new_cpu_map stages \[ksection\]:' run.log
