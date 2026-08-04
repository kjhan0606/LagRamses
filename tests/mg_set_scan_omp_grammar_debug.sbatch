#!/bin/bash
#SBATCH --job-name=mg_scan_omp
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/mg_set_scan_omp_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/mg_set_scan_omp_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
suite=/home/kjhan/BACKUP/FDM/mg_set_scan_omp_${SLURM_JOB_ID}
ntasks=${MG_NTASKS:-16}
if [[ ! "${ntasks}" =~ ^[1-9][0-9]*$ ]] || (( ntasks > SLURM_NTASKS )); then
  echo "invalid MG_NTASKS=${ntasks} (allocation has ${SLURM_NTASKS})" >&2
  exit 2
fi
mkdir -p "${suite}"

cd "${repo}/bin"
if [[ "${MG_INCREMENTAL_BUILD:-0}" != "1" ]]; then
  make clean
fi
make -j1 FDMDEBUG=1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${suite}/ramses_final3d.debug"

for threads in 1 4; do
  run_dir="${suite}/omp${threads}"
  mkdir -p "${run_dir}"
  cp "${repo}/tests/fdm_refine_mass.nml" "${run_dir}/cosmo.nml"
  : > "${run_dir}/jobcontrol.txt"

  cd "${run_dir}"
  export OMP_NUM_THREADS=${threads}
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n "${ntasks}" -c 4 --cpu-bind=cores \
    "${suite}/ramses_final3d.debug" cosmo.nml | tee run.log

  grep -q 'Run completed' run.log
  test -f output_00001/COMPLETE
  grep -q 'MG set_scan profile:' run.log
  if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
    echo "MG OMP${threads} diagnostic failure" >&2
    exit 2
  fi
  grep -E 'MG (set_scan|smooth|residual) profile:|Time elapsed since last coarse step:|FDM: M_tot=' run.log \
    > "${suite}/omp${threads}.summary"
  grep -E '==> Level=|FDM: M_tot=|rho_max=|Level  7 has|Morton \[step\]|MCHK:.*post-refine' run.log \
    > "${suite}/omp${threads}.state"
done

diff -u "${suite}/omp1.state" "${suite}/omp4.state"

t1=$(sed -n 's/.*MG set_scan profile:.*wall=[[:space:]]*\([0-9.]*\)[[:space:]]*s.*/\1/p' \
  "${suite}/omp1/run.log" | awk '{sum+=$1; n++} END{if(n==0) exit 1; printf "%.9f",sum/n}')
t4=$(sed -n 's/.*MG set_scan profile:.*wall=[[:space:]]*\([0-9.]*\)[[:space:]]*s.*/\1/p' \
  "${suite}/omp4/run.log" | awk '{sum+=$1; n++} END{if(n==0) exit 1; printf "%.9f",sum/n}')
speedup=$(awk -v t1="${t1}" -v t4="${t4}" 'BEGIN{if(t4<=0) exit 1; printf "%.3f",t1/t4}')
echo "MG_SET_SCAN_OMP_PASS t1=${t1}s t4=${t4}s speedup=${speedup}"
awk -v speedup="${speedup}" 'BEGIN{if(speedup<1.5) exit 1}'
