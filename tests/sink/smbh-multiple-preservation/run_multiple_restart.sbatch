#!/bin/bash
#SBATCH --job-name=smbh-multiple-restart
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-restart-%j.out
#SBATCH --error=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-restart-%j.err

# One MPI step from a copied, complete checkpoint. Source outputs are immutable.
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 FFTW_NUM_THREADS=1
if [[ -r /usr/lib64/libpmi2.so ]]; then
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
fi
repo=${1:?Provide the isolated lagRamses checkout}
continuous=${2:?Provide the complete two-output smoke directory}
restart_phi_policy=${3:-predictor}
fixture=${repo}/tests/sink/smbh-multiple-preservation
binary=${repo}/bin/ramses_final3d
test -z "$(git -C "${repo}" status --porcelain)"
test -x "${binary}"
test -f "${continuous}/output_00001/COMPLETE"
test -f "${continuous}/output_00002/COMPLETE"
test -f "${continuous}/binary.sha256"
test "$(sha256sum "${binary}" | cut -d' ' -f1)" = \
     "$(cut -d' ' -f1 "${continuous}/binary.sha256")"
# Test-only commits after the baseline build are permitted; the selected
# lagRamses patch sources must still be byte-identical to that build.
baseline_revision=$(<"${continuous}/source-revision.txt")
test -z "$(git -C "${repo}" diff --name-only "${baseline_revision}" HEAD -- patch/lagRamses)"
root=/gpfs/kjhan/FDM_WRITER_BUILD/MULTIPLE_RUNTIME
run_dir=${root}/${SLURM_JOB_ID}
mkdir "${run_dir}"
cp -a "${continuous}/output_00001" "${run_dir}/output_00001"
cp "${continuous}/cosmo.nml" "${run_dir}/cosmo.nml"
if [[ "${restart_phi_policy}" == warm ]]; then
  sed -i '/^&POISSON_PARAMS$/a restart_phi_warm_start=.true.' "${run_dir}/cosmo.nml"
elif [[ "${restart_phi_policy}" != predictor ]]; then
  printf 'unsupported restart phi policy: %s\n' "${restart_phi_policy}" >&2
  exit 2
fi
printf '%s\n' "$(git -C "${repo}" rev-parse HEAD)" > "${run_dir}/test-revision.txt"
sha256sum "${run_dir}/cosmo.nml" > "${run_dir}/inputs.sha256"
printf 'effective_namelist=%s\n' "${run_dir}/cosmo.nml"
printf 'command_line_nrestart=1\n'
printf 'restart_phi_policy=%s\n' "${restart_phi_policy}"
printf 'run_purpose=restart_short_evolution_test\n'
printf 'noutput=0; aout=none; foutput=1; fbackup=1; nstepmax=2\n'
printf 'expected_copied_dumps=1; expected_new_dumps=1; estimated_max_total_storage=200MiB\n'
df -h /gpfs
cd "${run_dir}"
: > jobcontrol.txt
srun --kill-on-bad-exit=1 --mpi=pmi2 -n2 -c1 --cpu-bind=cores \
  "${binary}" cosmo.nml 1 > run.log 2>&1
grep -Fq 'Run completed' run.log
test -f output_00002/COMPLETE
test ! -e output_00003
python "${repo}/patch/lagRamses/aux/validate_smbh_capture_ledger.py" \
  smbh_capture_ledger_v1.jsonl > capture-validation.json
python "${fixture}/verify_multiple_restart.py" "${continuous}" "${run_dir}" \
  > restart-verification.json
test -z "$(git -C "${repo}" status --porcelain)"
printf 'multiple_restart_status=passed\n'
printf 'run_directory=%s\n' "${run_dir}"
