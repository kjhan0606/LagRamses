#!/bin/bash
#SBATCH --job-name=smbh-multiple-build
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-build-%j.out
#SBATCH --error=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-build-%j.err

# Compilation only: no simulation, snapshots, or physical acceptance.
set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 gnu13/13.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
repo=${1:?Provide the isolated lagRamses checkout}
test -f "${repo}/patch/lagRamses/sink_particle.kjhan.f90"
test -z "$(git -C "${repo}" status --porcelain)"
build_revision=$(git -C "${repo}" rev-parse HEAD)
run_dir=/gpfs/kjhan/FDM_WRITER_BUILD/multiple-build-${SLURM_JOB_ID}
mkdir "${run_dir}"
module list 2>&1 | tee "${run_dir}/modules.log"
printf '%s\n' "${build_revision}" | tee "${run_dir}/source-revision.txt"
make --directory="${repo}/bin" -j4 ramses 2>&1 | tee "${run_dir}/build.log"
test "$(git -C "${repo}" rev-parse HEAD)" = "${build_revision}"
test -z "$(git -C "${repo}" status --porcelain)"
test -x "${repo}/bin/ramses_final3d"
sha256sum "${repo}/bin/ramses_final3d" | tee "${run_dir}/binary.sha256"
printf '%s\n' 'compile_status=passed'
