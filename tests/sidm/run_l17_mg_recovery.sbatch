#!/bin/bash
# Short operator/restart gate. These inherited z~4 checkpoints are diagnostic
# inputs, not certified science seeds. No production continuation is released.
# Submit with a frozen snapshot of this runner, analyze_mg_convergence_log.py,
# and audit_poisson_binary.py. Never modify a previously submitted snapshot:
# sbatch --export=ALL,SOURCE_REV=<commit>,VALIDATION_SNAPSHOT=<absolute-dir> ...
#SBATCH --job-name=sidm_l17_mg_recovery
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --array=0-1%1
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/mg_recovery_%A_%a.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/mg_recovery_%A_%a.err

set -euo pipefail
: "${SLURM_JOB_ID:?submit this calculation through Slurm}"
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17

code=/home/kjhan/BACKUP/lagRamses-SIDM
revision=${SOURCE_REV:?set an explicit immutable source commit}
validation_snapshot=${VALIDATION_SNAPSHOT:?set a frozen runner/analyzer directory}
base=/gpfs/kjhan/Hydro/Sidm/validation
case "${SLURM_ARRAY_TASK_ID:?}" in
  0) label=cdm; seed=${base}/z4_force_l17/cdm_323683/output_00002 ;;
  1) label=sidm3; seed=${base}/z4_force_l17/sidm3_323771/output_00002 ;;
  *) exit 2 ;;
esac
root=${base}/l17_mg_recovery/${SLURM_ARRAY_JOB_ID}/${label}
if [[ -e ${root} ]]; then
  printf 'refusing to reuse %s\n' "${root}" >&2
  exit 3
fi
test -f "${seed}/info_00002.txt"
test -f "${validation_snapshot}/analyze_mg_convergence_log.py"
test -f "${validation_snapshot}/audit_poisson_binary.py"
git -C "${code}" cat-file -e "${revision}^{commit}"
mkdir -p "${root}/source" "${root}/snapshot"

# Source-only archive and build happen on the allocated compute node. Include
# no paper edits or executable from another run. Override git build metadata
# because an extracted archive is deliberately not a Git checkout.
git -C "${code}" archive "${revision}" | tar -xf - -C "${root}/source"
cp -p "${validation_snapshot}/analyze_mg_convergence_log.py" "${root}/snapshot/"
cp -p "${validation_snapshot}/audit_poisson_binary.py" "${root}/snapshot/"
cp -p "${validation_snapshot}/run_l17_mg_recovery.sbatch" "${root}/snapshot/"
module -t list >"${root}/snapshot/modules.txt" 2>&1
printf '%s\n' "${revision}" >"${root}/snapshot/source_revision.txt"
# The legacy Makefile lacks complete Fortran module dependencies: even -j2
# races amr_commons against amr_parameters on an empty source archive.
make -C "${root}/source/bin" -j1 HDF5=1 USE_FFTW=1 \
  GITHASH="${revision}" GITBRANCH=frozen-sidm-recovery \
  GITREPO=git@github.com:kjhan0606/LagRamses.git >"${root}/build.log" 2>&1
binary=${root}/source/bin/ramses_final3d
test -x "${binary}"
sha256sum "${binary}" >"${root}/snapshot/binary.sha256"
ldd "${binary}" >"${root}/snapshot/ldd.txt"
if grep -q 'not found' "${root}/snapshot/ldd.txt"; then exit 4; fi
python3 "${root}/snapshot/audit_poisson_binary.py" "${binary}" \
  --json "${root}/snapshot/poisson_binary_gate.json" --require-fftw
coarse=$(awk -F= '/^nstep_coarse[[:space:]]*=/{gsub(/[[:space:]]/,"",$2);print $2}' \
  "${seed}/info_00002.txt")
[[ ${coarse} =~ ^[0-9]+$ ]]
stop_step=$((coarse + 2))

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_FABRICS=shm
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

gate_failed=0
for cap in 40 80; do
  case_root=${root}/cap${cap}
  mkdir "${case_root}"
  ln -s "${seed}" "${case_root}/output_00002"
  cp -p "${root}/source/tests/sidm/z0_l17_continue_${label}.nml.in" "${case_root}/template.nml"
  # All replacements target the explicit source template. A future dummy
  # clock plus the absolute coarse-step cap suppress every new full dump.
  sed -e 's/__NRESTART__/2/' \
      -e "/^nrestart=2/a nstepmax=${stop_step}" \
      -e 's/^noutput=5$/noutput=1/' \
      -e 's/^aout=.*/aout=1.1\ntout=1.0d100/' \
      -e 's/^foutput=.*/foutput=1000000000\nfbackup=1000000000/' \
      -e "/^&POISSON_PARAMS/a epsilon=1.0d-4\nmaxiter_fine=${cap}\nabort_on_mg_nonconvergence=.true.\nrestart_phi_warm_start=.false." \
      "${case_root}/template.nml" >"${case_root}/run.nml"
  : >"${case_root}/jobcontrol.txt"
  printf 'effective_namelist=%s noutput=1 dummy_aout=1.1 dummy_z=-0.090909 nstepmax=%s expected_new_dumps=0\n' \
    "${case_root}/run.nml" "${stop_step}"
  if (
    cd "${case_root}"
    mpirun -np "${SLURM_NTASKS}" "${binary}" run.nml >run.log 2>&1
  ); then runtime_rc=0; else runtime_rc=$?; fi
  if python3 "${root}/snapshot/analyze_mg_convergence_log.py" "${case_root}/run.log" \
    --epsilon 1e-4 --maxiter "${cap}" --require-converged --require-fftw \
    --json "${case_root}/mg_summary.json"; then analyzer_rc=0; else analyzer_rc=$?; fi
  if [[ ${runtime_rc} -ne 0 || ${analyzer_rc} -ne 0 ]] || \
     ! grep -Fq 'Run completed' "${case_root}/run.log"; then gate_failed=1; fi
  printf 'cap=%s runtime_exit=%s analyzer_exit=%s\n' "${cap}" "${runtime_rc}" "${analyzer_rc}" \
    >"${case_root}/exit_status.txt"
  # Only the inherited seed link is allowed. There are no tree walks here.
  for output in "${case_root}"/output_*; do
    [[ ${output} == "${case_root}/output_00002" ]] || { printf 'unexpected dump: %s\n' "${output}" >&2; exit 5; }
  done
done
if [[ ${gate_failed} -ne 0 ]]; then
  printf 'L17_MG_RECOVERY_RESULT=failed label=%s (inspect both cap summaries)\n' "${label}" >&2
  exit 6
fi
printf 'L17_MG_RECOVERY_RESULT=passed label=%s caps=40,80\n' "${label}"
