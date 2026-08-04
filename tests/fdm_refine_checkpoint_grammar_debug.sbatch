#!/bin/bash
#SBATCH --job-name=fdm_refckpt
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_refine_checkpoint_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_refine_checkpoint_%j.err

set -eo pipefail
source /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh >/dev/null 2>&1
set -u
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

repo=/home/kjhan/BACKUP/lagRamses
source_run=/gpfs/kjhan/Hydro/FDM/run_3way_wave_matched_v2
run_dir=/home/kjhan/BACKUP/FDM/fdm_refine_checkpoint_${SLURM_JOB_ID}
mkdir -p "${run_dir}"

cd "${repo}/bin"
make clean
make FDMTRACE=1 HDF5=1 USE_FFTW=1
cp ramses_final3d "${run_dir}/ramses_final3d.trace"
cp "${source_run}/cosmo.nml" "${run_dir}/cosmo.nml"
ln -s "${source_run}/output_00006" "${run_dir}/output_00006"
ln -s "${source_run}/yield_table.asc" "${run_dir}/yield_table.asc"
: > "${run_dir}/jobcontrol.txt"

# Keep the same per-rank grid capacity while fitting the 96-rank checkpoint on
# one 64-core debug node.  Stop immediately after the historical first-jump
# window and suppress all scientific output.
sed -i 's/ngridtot=1300000000/ngridtot=184000000/' "${run_dir}/cosmo.nml"
# The checkpoint has exactly zero particles (all 96 part files are header-only).
# Reserve enough for diagnostics/possible star creation, not the production
# 40-million-particle capacity, which only wastes memory in this short replay.
sed -i 's/nparttot=40000000/nparttot=1000000/' "${run_dir}/cosmo.nml"
sed -i -E 's/^nrestart=.*/nrestart=6/' "${run_dir}/cosmo.nml"
sed -i '/ncontrol=1/a nstepmax=273' "${run_dir}/cosmo.nml"
sed -i 's/foutput=100/foutput=1000000/' "${run_dir}/cosmo.nml"

cd "${run_dir}"
srun --kill-on-bad-exit=1 --mpi=pmi2 -n 64 -c 1 --cpu-bind=cores \
  ./ramses_final3d.trace cosmo.nml | tee run.log

grep -q 'Run completed' run.log
awk '/Fine step=/{step=$3+0} END{exit(step>=272?0:1)}' run.log
awk '/Level[[:space:]]+11 has/{if($4+0>0) found=1} END{exit(found?0:1)}' run.log
grep -q 'FDM FFT map rebuild: epoch' run.log
grep -q 'FDM FFT map check:.*mismatch=0' run.log
if rg -q 'No more free|FATAL|mismatch=[[:space:]]*[1-9]|NaN_CHK.*=[[:space:]]*[1-9]' run.log; then
  echo 'FDM checkpoint refinement diagnostic failure' >&2
  exit 2
fi

# Refinement itself must not change leaf mass.  Production evolution may have
# solver-scale roundoff, so compare paired pre/post-refine MCHK lines directly.
awk '
  /MCHK:.*pre-refine/  {pre=$0; sub(/^.*total=[[:space:]]*/,"",pre); sub(/[[:space:]].*$/,"",pre)}
  /MCHK:.*post-refine/ {post=$0; sub(/^.*total=[[:space:]]*/,"",post); sub(/[[:space:]].*$/,"",post);
                         if(pre!="" && post!="") {d=post-pre; if(d<0)d=-d; if(d>2e-10) bad=1}}
  END{exit(bad?1:0)}' run.log

# The complete evolution through the historical corruption window must also
# conserve the composite leaf mass, not merely each refine_fine call.
first_mass=$(sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log | head -n 1)
last_mass=$(sed -n 's/.*FDM: M_tot=[[:space:]]*\([0-9.E+-]*\).*/\1/p' run.log | tail -n 1)
awk -v a="${first_mass}" -v b="${last_mass}" 'BEGIN{d=b-a;if(d<0)d=-d; if(d>1e-6) exit 1}'
echo 'FDM_REFINE_CHECKPOINT_PASS'
