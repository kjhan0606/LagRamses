#!/bin/bash
#SBATCH --job-name=fdm_refmass
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_refine_mass_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_refine_mass_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_dir=/home/kjhan/BACKUP/FDM/fdm_refine_mass_${SLURM_JOB_ID}
mkdir -p "${run_dir}"

cd "${repo}/bin"
make clean
make FDMDEBUG=1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_dir}/ramses_final3d.debug"
cp "${repo}/tests/fdm_refine_mass.nml" "${run_dir}/cosmo.nml"
: > "${run_dir}/jobcontrol.txt"

cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
  ./ramses_final3d.debug cosmo.nml | tee run.log

grep -q 'Run completed' run.log
test -f output_00001/COMPLETE
awk '/Level[[:space:]]+7 has/{if($4+0>0) found=1} END{exit(found?0:1)}' run.log
if rg -q 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'FDM refinement diagnostic failure' >&2
  exit 2
fi

last_mass=$(sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log | tail -n 1)
nmass=$(rg -c 'FDM: M_tot=' run.log)
test "${nmass}" -ge 2
# This includes one complete CN wave step in addition to the refinement.  The
# observed corrected error is 1.27e-8; the production defect was 1.3e-2.
awk -v m="${last_mass}" 'BEGIN{d=m-1.0;if(d<0)d=-d; if(d>2e-8) exit 1}'
echo "FDM_REFINE_MASS_PASS last_mass=${last_mass}"
