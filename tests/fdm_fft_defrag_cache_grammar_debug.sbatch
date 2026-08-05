#!/bin/bash
#SBATCH --job-name=fdm_fftcache
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_fft_defrag_cache_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_fft_defrag_cache_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
run_dir=/home/kjhan/BACKUP/FDM/fdm_fft_defrag_cache_${SLURM_JOB_ID}
mkdir -p "${run_dir}"

cd "${repo}/bin"
make clean
make FDMDEBUG=1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_dir}/ramses_final3d.debug"
cp "${repo}/tests/fdm_refine_mass.nml" "${run_dir}/cosmo.nml"
: > "${run_dir}/jobcontrol.txt"

# A 256^3 uniform base level selects the distributed FDM FFT.  Keep one AMR
# level above it: amr_step intentionally disables periodic load balancing when
# levelmin==nlevelmax, so levelmax=9 is required to exercise a second defrag
# between cached FFT uses.  Remapping every step then forces the epoch change.
sed -i 's/levelmin=6/levelmin=8/' "${run_dir}/cosmo.nml"
sed -i 's/levelmax=7/levelmax=9/' "${run_dir}/cosmo.nml"
sed -i 's/nremap=1000000/nremap=1/' "${run_dir}/cosmo.nml"
sed -i 's/nstepmax=2/nstepmax=3/' "${run_dir}/cosmo.nml"
# Refinement can flag most of the 256^3 base cells at the deliberately low
# test threshold.  Reserve enough for the worst-case 16.8M level-9 grids plus
# the base hierarchy and load-balance headroom.
sed -i 's/ngridtot=4000000/ngridtot=32000000/' "${run_dir}/cosmo.nml"
sed -i 's@/home/kjhan/BACKUP/FDM/ksection_zoom_test/ic_zoom/level_006@/gpfs/kjhan/Hydro/FDM/ic_lv8_z30_15Mpc/level_008@' "${run_dir}/cosmo.nml"

cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 --cpu-bind=cores \
  ./ramses_final3d.debug cosmo.nml | tee run.log

grep -q 'Run completed' run.log
test "$(grep -c 'FDM FFT map rebuild: epoch' run.log)" -ge 2
if grep -Eq 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'FDM FFT/defrag cache diagnostic failure' >&2
  exit 2
fi

first_mass=$(sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log | head -n 1)
last_mass=$(sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log | tail -n 1)
test -n "${first_mass}"
test -n "${last_mass}"
# Refinement itself must be conservative to roundoff.  The intentionally
# abrupt creation of millions of L9 grids can leave the 10-cycle fine MG
# solve short of convergence, so allow a small solver-scale evolution drift
# while still rejecting the historical percent-level corruption.
awk '
  /MCHK:.*pre-refine/  {pre=$0; sub(/^.*total=[[:space:]]*/,"",pre); sub(/[[:space:]].*$/,"",pre)}
  /MCHK:.*post-refine/ {post=$0; sub(/^.*total=[[:space:]]*/,"",post); sub(/[[:space:]].*$/,"",post);
                         if(pre!="" && post!="") {d=post-pre; if(d<0)d=-d; if(d>2e-10) bad=1}}
  END{exit(bad?1:0)}' run.log
awk -v a="${first_mass}" -v b="${last_mass}" 'BEGIN{d=b-a;if(d<0)d=-d; if(d>5e-5) exit 1}'
echo "FDM_FFT_DEFRAG_CACHE_PASS first=${first_mass} last=${last_mass}"
