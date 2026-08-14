#!/bin/bash
#SBATCH --job-name=runtime_lb
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/runtime_lb_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/runtime_lb_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_root=/home/kjhan/BACKUP/FDM/runtime_lb_${SLURM_JOB_ID}
mkdir -p "${run_root}/cdm" "${run_root}/fdm_alias"

cd "${repo}/bin"
make clean
make HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_root}/ramses_final3d"

cp "${repo}/tests/runtime_lb_cdm.nml" "${run_root}/cdm/cosmo.nml"
cp "${repo}/tests/runtime_lb.jobcontrol" "${run_root}/cdm/jobcontrol.txt"
cd "${run_root}/cdm"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  ../ramses_final3d cosmo.nml | tee run.log

grep -Fq 'Job control: load-balance mode memory -> timed at step 1' run.log
grep -Fq 'Job control: load-balance mode timed -> work at step 3' run.log
grep -Fq 'Job control: load-balance mode work -> memory at step 4' run.log
grep -Fq 'Job control: forced load balance requested at step 5' run.log
test "$(grep -Fc 'Runtime load balance: forcing remap at coarse step' run.log)" -eq 4
grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'Runtime load-balance CDM regression failure' >&2
  exit 2
fi

# Verify that legacy FDM inputs actually select the generic timed work model.
cp "${repo}/tests/dynamic_fdm.nml" "${run_root}/fdm_alias/cosmo.nml"
: > "${run_root}/fdm_alias/jobcontrol.txt"
cd "${run_root}/fdm_alias"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  ../ramses_final3d cosmo.nml | tee run.log
grep -Fq 'fdm_cost_mode=1 (deprecated alias): selecting timed work balance' run.log
grep -Fq 'effective LB=timed work' run.log
grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'FDM timed-work alias regression failure' >&2
  exit 3
fi

echo "RUNTIME_LB_PASS binary=${run_root}/ramses_final3d"
