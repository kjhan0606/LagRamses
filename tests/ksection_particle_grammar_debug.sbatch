#!/bin/bash
#SBATCH --job-name=ksec_part
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/ksection_particle_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/ksection_particle_%j.err

set -euo pipefail

export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_dir=/home/kjhan/BACKUP/FDM/ksection_particle_${SLURM_JOB_ID}
mkdir -p "${run_dir}"
cp "${repo}/tests/dynamic_particle_cdm.nml" "${run_dir}/cosmo.nml"
sed -i "s/^ordering=.*/ordering='ksection'/" "${run_dir}/cosmo.nml"
sed -i 's/^nstepmax=.*/nstepmax=1/' "${run_dir}/cosmo.nml"
sed -i 's/^poisson=.*/poisson=.false./' "${run_dir}/cosmo.nml"
# Ksection ownership of this zoom IC is intentionally less uniform than
# Hilbert.  Give each rank ample particle capacity so the test measures the
# common particle-aware cost/remap path rather than npartmax exhaustion.
sed -i 's/^nparttot=.*/nparttot=16000000/' "${run_dir}/cosmo.nml"
sed -i '/^nparttot=/a npartmax=4000000' "${run_dir}/cosmo.nml"
: > "${run_dir}/jobcontrol.txt"

cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 4 --cpu-bind=cores \
  "${repo}/bin/ramses_final3d" cosmo.nml | tee run.log

grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
grep -q 'LB compact stages:' run.log
grep -q 'cmp_new_cpu_map stages \[ksection\]:' run.log
if grep -Eq 'No more free|Increase npartmax|npartmax should|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'Ksection particle/remap diagnostic failure' >&2
  exit 2
fi
