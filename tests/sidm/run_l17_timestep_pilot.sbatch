#!/bin/bash
# Diagnostic evolution only: a legacy seed is not a certified science seed.
# Submit with --dependency=afterok:<paired recovery array>, a frozen PILOT_BUNDLE
# and RECOVERY_JOB. This job cannot release a production continuation.
#SBATCH --job-name=sidm_l17_dt_pilot
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/timestep_pilot_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/timestep_pilot_%j.err

set -euo pipefail
: "${SLURM_JOB_ID:?submit through Slurm}"
bundle=${PILOT_BUNDLE:?set a frozen diagnostic bundle}
recovery_job=${RECOVERY_JOB:?set the successful paired operator gate}
[[ ${recovery_job} =~ ^[0-9]+$ ]]
base=/gpfs/kjhan/Hydro/Sidm/validation
root=${base}/l17_timestep_pilot/${SLURM_JOB_ID}
seed=${base}/z0_production_l17/sidm3/output_00005
[[ ! -e ${root} ]] || { printf 'refusing to reuse %s\n' "${root}" >&2; exit 2; }
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
# Slurm afterok is the scheduler gate. The manifest independently checks all
# four numerical outcomes and the frozen inputs before any binary is executed.
python3 "${bundle}/validate_pilot_bundle.py" "${bundle}" \
  --recovery-root "${base}/l17_mg_recovery/${recovery_job}"
mkdir -p "${root}/snapshot"
cp -p "${bundle}/run.nml" "${root}/run.nml"
cp -p "${bundle}/manifest.json" "${root}/snapshot/manifest.json"
cp -p "${bundle}/analyze_timestep_diagnostics.py" \
  "${bundle}/analyze_mg_convergence_log.py" "${root}/snapshot/"
module -t list >"${root}/snapshot/modules.txt" 2>&1
ln -s "${seed}" "${root}/output_00005"
: >"${root}/jobcontrol.txt"
printf 'purpose=short_diagnostic_evolution effective_namelist=%s noutput=1 dummy_aout=1.1 dummy_z=-0.090909 nstepmax=858 expected_new_dumps=0 seed_gib=6.680\n' \
  "${root}/run.nml"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export I_MPI_FABRICS=shm
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
if (cd "${root}"; mpirun -np "${SLURM_NTASKS}" "${bundle}/ramses_dt_pilot" run.nml >run.log 2>&1); then
  runtime_rc=0
else runtime_rc=$?; fi
printf 'runtime_exit=%s\n' "${runtime_rc}" >"${root}/snapshot/runtime_exit.txt"
for output in "${root}"/output_*; do
  [[ ${output} == "${root}/output_00005" ]] || { printf 'unexpected dump: %s\n' "${output}" >&2; exit 3; }
done
python3 "${root}/snapshot/analyze_mg_convergence_log.py" "${root}/run.log" \
  --epsilon 1e-4 --maxiter 80 --require-converged --require-fftw \
  --json "${root}/mg_summary.json"
[[ ${runtime_rc} == 0 ]] || exit 4
python3 "${root}/snapshot/analyze_timestep_diagnostics.py" "${root}/run.log" \
  --pilot-seed-info "${seed}/info_00005.txt" --pilot-steps 20 --warmup-steps 2 \
  --json "${root}/pilot_summary.json"
printf 'L17_TIMESTEP_PILOT_RESULT=passed scope=diagnostic_only\n'
