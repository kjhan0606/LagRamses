#!/bin/bash
#SBATCH --job-name=smbh-multiple-h5restart
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-h5restart-%j.out
#SBATCH --error=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-h5restart-%j.err

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 FFTW_NUM_THREADS=1
if [[ -r /usr/lib64/libpmi2.so ]]; then
  export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
fi
repo=${1:?Provide the isolated HDF5 lagRamses checkout}
continuous=${2:?Provide the complete HDF5 smoke directory}
expectation=${3:-pass}
restart_index=${4:-2}
if [[ "${restart_index}" != 1 && "${restart_index}" != 2 ]]; then
  printf 'unsupported restart output index: %s\n' "${restart_index}" >&2
  exit 2
fi
printf -v source_name 'output_%05d' "${restart_index}"
printf -v source_number '%05d' "${restart_index}"
target_index=$((restart_index+1))
printf -v target_name 'output_%05d' "${target_index}"
test -z "$(git -C "${repo}" status --porcelain)"
binary=${repo}/bin/ramses_final3d
test -x "${binary}"
test -f "${continuous}/${source_name}/COMPLETE"
test -f "${continuous}/${target_name}/COMPLETE"
test "$(sha256sum "${binary}" | cut -d' ' -f1)" = \
     "$(cut -d' ' -f1 "${continuous}/binary.sha256")"
root=/gpfs/kjhan/FDM_WRITER_BUILD/MULTIPLE_HDF5_RUNTIME
run_dir=${root}/${SLURM_JOB_ID}
mkdir "${run_dir}"
cp -a "${continuous}/${source_name}" "${run_dir}/${source_name}"
cp "${continuous}/cosmo.nml" "${run_dir}/cosmo.nml"
# Stop after exactly one evolved step, even when the uninterrupted baseline
# contains later outputs used by another restart gate.
sed -i "s/^nstepmax=[0-9][0-9]*$/nstepmax=${target_index}/" "${run_dir}/cosmo.nml"
grep -Fqx "nstepmax=${target_index}" "${run_dir}/cosmo.nml"
sha256sum "${run_dir}/cosmo.nml" > "${run_dir}/inputs.sha256"
printf 'effective_namelist=%s\n' "${run_dir}/cosmo.nml"
printf 'command_line_nrestart=%s\n' "${restart_index}"
printf 'run_purpose=restart_short_hdf5_evolution_test; expectation=%s\n' "${expectation}"
printf 'noutput=0; aout=none; foutput=1; fbackup=1; %s\n' \
  "$(grep '^nstepmax=' "${run_dir}/cosmo.nml")"
if [[ "${expectation}" == reject-legacy-bounds ]]; then
  printf 'expected_copied_dumps=1; expected_new_dumps=0; estimated_max_total_storage=100MiB\n'
else
  printf 'expected_copied_dumps=1; expected_new_dumps=1; estimated_max_total_storage=200MiB\n'
fi
df -h /gpfs
cd "${run_dir}"
: > jobcontrol.txt
if [[ "${expectation}" == reject-legacy-bounds ]]; then
  python "${repo}/tests/sink/smbh-multiple-preservation/remove_hdf5_split_bounds.py" \
    "${source_name}/data_${source_number}.h5"
fi
if [[ "${expectation}" == reject-legacy-bounds ]]; then
  if srun --kill-on-bad-exit=1 --mpi=pmi2 -n2 -c1 --cpu-bind=cores \
      "${binary}" cosmo.nml "${restart_index}" > run.log 2>&1; then
    printf '%s\n' 'ERROR: unsafe MPI HDF5 sink restart was accepted' >&2
    exit 1
  fi
  grep -Fq 'ERROR: MPI sink HDF5 restart requires lossless Hilbert bound_key' run.log
  test ! -e "${target_name}"
  printf 'multiple_hdf5_restart_status=correctly_rejected\n'
  printf 'run_directory=%s\n' "${run_dir}"
  exit 0
fi
test "${expectation}" = pass
srun --kill-on-bad-exit=1 --mpi=pmi2 -n2 -c1 --cpu-bind=cores \
  "${binary}" cosmo.nml "${restart_index}" > run.log 2>&1
grep -Fq 'Run completed' run.log
test -f "${target_name}/COMPLETE"
printf -v unexpected_name 'output_%05d' "$((target_index+1))"
test ! -e "${unexpected_name}"
python "${repo}/patch/lagRamses/aux/validate_smbh_capture_ledger.py" \
  smbh_capture_ledger_v1.jsonl > capture-validation.json
python "${repo}/tests/sink/smbh-multiple-preservation/verify_multiple_hdf5.py" \
  "${continuous}" "${run_dir}" --target-output "${target_index}" \
  > restart-verification.json
test -z "$(git -C "${repo}" status --porcelain)"
printf 'multiple_hdf5_restart_status=passed\n'
printf 'run_directory=%s\n' "${run_dir}"
