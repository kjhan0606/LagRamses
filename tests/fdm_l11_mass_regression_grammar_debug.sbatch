#!/bin/bash
#SBATCH --job-name=fdm_l11_mass
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --time=01:30:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_l11_mass_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_l11_mass_%j.err

# Exercise the same base-to-first-refined-level transition as production with
# a 64^3 -> level-7 mesh that fits on grammar-debug.  The three sequential
# cases differ only in their wave/HJM representation.

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u

export OMP_NUM_THREADS=2
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/FDM/lagramses_science_clean
suite=/home/kjhan/BACKUP/FDM/fdm_l11_regression
run_root=${suite}/runs/${SLURM_JOB_ID}
mkdir -p "${run_root}"

git -C "${repo}" status --short > "${run_root}/git-status.txt"
git -C "${repo}" diff > "${run_root}/source.patch"
sha256sum "${repo}/patch/lagRamses/fdm_step.f90" \
  "${repo}/patch/lagRamses/fdm_hjm.f90" \
  "${repo}/patch/lagRamses/amr_step.jaehyun.f90" > "${run_root}/source.sha256"

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
  max_level=${4:-7}
  refine_threshold=${5:-1.1}
  case_dir=${run_root}/${arm}

  mkdir -p "${case_dir}"
  cp "${repo}/tests/fdm_refine_mass.nml" "${case_dir}/cosmo.nml"
  : > "${case_dir}/jobcontrol.txt"
  sed -i 's/nstepmax=2/nstepmax=3/' "${case_dir}/cosmo.nml"
  sed -i "s/levelmax=7/levelmax=${max_level}/" "${case_dir}/cosmo.nml"
  sed -i "s/m_refine=9\*1.1/m_refine=9*${refine_threshold}/" "${case_dir}/cosmo.nml"
  sed -i "s/fdm_use_hjm=.false./fdm_use_hjm=${use_hjm}/" "${case_dir}/cosmo.nml"
  sed -i "s/fdm_first_wave_level=6/fdm_first_wave_level=${first_wave}/" "${case_dir}/cosmo.nml"
  sed -i '/&FDM_PARAMS/a fdm_ghost2=.true.\nfdm_ghost2_rev=.true.' "${case_dir}/cosmo.nml"
  sed -i '/&POISSON_PARAMS/a maxiter_fine=20\nrestart_phi_warm_start=.false.\nabort_on_mg_nonconvergence=.true.' \
    "${case_dir}/cosmo.nml"

  cd "${case_dir}"
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n 16 -c 2 --cpu-bind=cores \
    "${run_root}/ramses_final3d.debug" cosmo.nml | tee run.log

  grep -q 'Run completed' run.log
  awk -v target="${max_level}" \
    '$1=="Level" && $2+0==target && $3=="has" {if($4+0>0) found=1} END{exit(found?0:1)}' \
    run.log
  if grep -Eiq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
    echo "${arm}: fatal diagnostic" >&2
    return 2
  fi

  # Every individual refine_fine call must conserve composite leaf mass.
  awk '
    /MCHK:.*pre-refine/  {x=$0; sub(/^.*total=[[:space:]]*/,"",x); sub(/[[:space:]].*$/,"",x); pre=x}
    /MCHK:.*post-refine/ {x=$0; sub(/^.*total=[[:space:]]*/,"",x); sub(/[[:space:]].*$/,"",x);
                           if(pre!="") {d=x-pre; if(d<0)d=-d; if(d>2e-10) bad=1}}
    END{exit(bad?1:0)}' run.log

  # The complete operator sequence must conserve composite leaf mass, not only
  # mesh creation.  This catches base FFT/coarse-fine transfer, HJM floors,
  # fine-level CN drift, kick, and restriction independently.
  awk '
    /MCHK:/ {x=$0; sub(/^.*total=[[:space:]]*/,"",x); sub(/[[:space:]].*$/,"",x);
              if(!seen) {ref=x; seen=1}
              d=x-ref; if(d<0)d=-d; if(d>2e-10) bad=1}
    END{exit(!seen || bad)}' run.log

  grep 'MCHK:' run.log > mass-phases.log
  sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log > mass-total.log
  echo "${arm}: PASS full-phase mass conservation and diagnostics"
}

run_case wave .false. 6
run_case hybrid .true. 7
run_case fluid .true. 99
run_case wave_l8 .false. 6 8 0.1

echo "FDM_L11_MASS_REGRESSION_PASS $(date --iso-8601=seconds)"
