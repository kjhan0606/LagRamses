#!/bin/bash
#SBATCH --job-name=smbh-multiple-hdf5
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-hdf5-%j.out
#SBATCH --error=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-hdf5-%j.err

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
fixture_mode=${2:-two-step}
steps=2
if [[ "${fixture_mode}" == three-step ]]; then
  steps=3
elif [[ "${fixture_mode}" == three-step-sink-refined ]]; then
  steps=3
elif [[ "${fixture_mode}" != two-step ]]; then
  printf 'unsupported fixture mode: %s\n' "${fixture_mode}" >&2
  exit 2
fi
test -z "$(git -C "${repo}" status --porcelain)"
revision=$(git -C "${repo}" rev-parse HEAD)
fixture=${repo}/tests/sink/smbh-multiple-preservation
binary=${repo}/bin/ramses_final3d
test -x "${binary}"
root=/gpfs/kjhan/FDM_WRITER_BUILD/MULTIPLE_HDF5_RUNTIME
mkdir -p "${root}"
run_dir=${root}/${SLURM_JOB_ID}
mkdir "${run_dir}"
cp "${fixture}/smbh-multiple.nml" "${run_dir}/cosmo.nml"
if [[ "${fixture_mode}" == three-step-sink-refined ]]; then
  sed -i 's/^levelmax=4$/levelmax=5/' "${run_dir}/cosmo.nml"
  sed -i '/^&HYDRO_PARAMS$/i\&REFINE_PARAMS\nm_refine=1.0d6\nmass_sph=1.0d0\nsink_refine=.true.\n/\n' \
    "${run_dir}/cosmo.nml"
  sed -i 's/^nlevelmax_sink=4$/nlevelmax_sink=5/' "${run_dir}/cosmo.nml"
fi
if [[ "${steps}" == 3 ]]; then
  sed -i 's/^nstepmax=2$/nstepmax=3/' "${run_dir}/cosmo.nml"
fi
sed -i '/^&OUTPUT_PARAMS$/a outformat="hdf5"\ninformat="hdf5"' "${run_dir}/cosmo.nml"
cp "${fixture}/ic_sink" "${run_dir}/ic_sink"
printf '%s\n' "${revision}" > "${run_dir}/source-revision.txt"
sha256sum "${binary}" > "${run_dir}/binary.sha256"
sha256sum "${run_dir}/cosmo.nml" "${run_dir}/ic_sink" > "${run_dir}/inputs.sha256"
printf 'effective_namelist=%s\n' "${run_dir}/cosmo.nml"
printf 'run_purpose=short_hdf5_evolution_test\n'
printf 'fixture_mode=%s; noutput=0; aout=none; foutput=1; fbackup=1; nstepmax=%s\n' \
  "${fixture_mode}" "${steps}"
printf 'expected_full_dumps=%s; estimated_max_total_storage=%sMiB\n' \
  "${steps}" "$((steps*100))"
df -h /gpfs
cd "${run_dir}"
: > jobcontrol.txt
srun --kill-on-bad-exit=1 --mpi=pmi2 -n2 -c1 --cpu-bind=cores \
  "${binary}" cosmo.nml > run.log 2>&1
grep -Fq 'Run completed' run.log
if [[ "${fixture_mode}" == three-step-sink-refined ]]; then
  awk '/Level +5 has/ {if ($4 > 0 && $4 < 4096) populated=1; if ($4 >= 4096) full=1} END {exit !populated || full}' run.log
fi
test -f output_00001/COMPLETE
test -f output_00002/COMPLETE
if [[ "${steps}" == 3 ]]; then
  test -f output_00003/COMPLETE
  test ! -e output_00004
else
  test ! -e output_00003
fi
python "${repo}/patch/lagRamses/aux/validate_smbh_capture_ledger.py" \
  smbh_capture_ledger_v1.jsonl > capture-validation.json
python "${fixture}/verify_multiple_hdf5.py" "${run_dir}" \
  --target-output "${steps}" > multiple-verification.json
test "$(git -C "${repo}" rev-parse HEAD)" = "${revision}"
test -z "$(git -C "${repo}" status --porcelain)"
printf 'multiple_hdf5_smoke_status=passed\n'
printf 'run_directory=%s\n' "${run_dir}"
