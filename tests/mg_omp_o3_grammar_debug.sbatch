#!/bin/bash
#SBATCH --job-name=mg_omp_o3
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/mg_omp_o3_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/mg_omp_o3_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_dir=/home/kjhan/BACKUP/FDM/mg_omp_o3_${SLURM_JOB_ID}
mkdir -p "${run_dir}"

cd "${repo}/bin"
make clean
make -j1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_dir}/ramses_final3d.o3"
sha256sum "${run_dir}/ramses_final3d.o3"

cp "${repo}/tests/fdm_refine_mass.nml" "${run_dir}/cosmo.nml"
: > "${run_dir}/jobcontrol.txt"
cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 4 --cpu-bind=cores \
  ./ramses_final3d.o3 cosmo.nml | tee run.log

grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'MG O3 smoke diagnostic failure' >&2
  exit 2
fi
awk '/FDM: M_tot=/{mass=$3} END {
  if(mass=="") exit 1
  delta=mass-1.0
  if(delta<0) delta=-delta
  if(delta>1.0e-6) exit 1
}' run.log
echo MG_OMP_O3_PASS
