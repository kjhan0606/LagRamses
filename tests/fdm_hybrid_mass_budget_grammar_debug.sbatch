#!/bin/bash
#SBATCH --job-name=fdm_mass_budget
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_mass_budget_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_mass_budget_%j.err

# Production-topology miniature: fluid base L6, wave L7--8, with AMR
# subcycling and periodic remaps.  FDMDEBUG attributes any total-mass change
# to HJM stages, CN interface exchange, refinement, or restriction.

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u

export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/FDM/lagramses_science_clean
suite=/home/kjhan/BACKUP/FDM/fdm_hybrid_mass_budget
run_root=${suite}/runs/${SLURM_JOB_ID}
mkdir -p "${run_root}"

git -C "${repo}" rev-parse HEAD > "${run_root}/source.commit"
git -C "${repo}" status --short > "${run_root}/git-status.txt"
git -C "${repo}" diff > "${run_root}/source.patch"
sha256sum "${repo}/patch/lagRamses/fdm_step.f90" \
  "${repo}/patch/lagRamses/fdm_hjm.f90" \
  "${repo}/patch/lagRamses/amr_step.jaehyun.f90" > "${run_root}/source.sha256"
sha256sum "${repo}/tests/analyze_fdm_mass_budget.py" \
  "${repo}/tests/fdm_hybrid_mass_budget_grammar_debug.sbatch" > "${run_root}/test.sha256"

cd "${repo}/bin"
make clean
make FDMDEBUG=1 HDF5=1 USE_FFTW=1
test -x ./ramses_final3d
cp ./ramses_final3d "${run_root}/ramses_final3d.debug"
sha256sum "${run_root}/ramses_final3d.debug" > "${run_root}/binary.sha256"

run_case () {
  arm=$1
  use_hjm=$2
  first_wave=$3
  courant=$4
  method=$5
  case_dir=${run_root}/${arm}

  mkdir -p "${case_dir}"
  cp "${repo}/tests/dynamic_fdm.nml" "${case_dir}/cosmo.nml"
  : > "${case_dir}/jobcontrol.txt"
  sed -i 's/nstepmax=2/nstepmax=40/' "${case_dir}/cosmo.nml"
  sed -i 's/nremap=1/nremap=5/' "${case_dir}/cosmo.nml"
  sed -i 's/levelmax=9/levelmax=8/' "${case_dir}/cosmo.nml"
  sed -i 's/nsubcycle=1,1,1,1,1,1,1,1,1/nsubcycle=1,1,1,1,1,1,2,2,2/' \
    "${case_dir}/cosmo.nml"
  sed -i "s/fdm_courant=0.3/fdm_courant=${courant}/" "${case_dir}/cosmo.nml"
  sed -i "s/fdm_use_hjm=.true./fdm_use_hjm=${use_hjm}/" "${case_dir}/cosmo.nml"
  sed -i "s/fdm_first_wave_level=99/fdm_first_wave_level=${first_wave}/" \
    "${case_dir}/cosmo.nml"
  sed -i "/&RUN_PARAMS/a exchange_method='${method}'" "${case_dir}/cosmo.nml"
  sed -i '/&FDM_PARAMS/a fdm_refine_matched=.true.\nfdm_ghost2=.true.\nfdm_ghost2_rev=.true.' \
    "${case_dir}/cosmo.nml"
  sed -i '/&POISSON_PARAMS/a maxiter_fine=40\nrestart_phi_warm_start=.false.\nabort_on_mg_nonconvergence=.true.' \
    "${case_dir}/cosmo.nml"

  cd "${case_dir}"
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
    "${run_root}/ramses_final3d.debug" cosmo.nml | tee run.log
}

failures=0
run_case hybrid_l7_target_budget .true. 7 0.30 p2p || failures=$((failures+1))

python3 "${repo}/tests/analyze_fdm_mass_budget.py" "${run_root}" \
  --output "${run_root}/report.json" | tee "${run_root}/report.txt"
test "${failures}" -eq 0
echo "FDM_HYBRID_MASS_BUDGET_DONE $(date --iso-8601=seconds)"
