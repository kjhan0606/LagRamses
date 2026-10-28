#!/bin/bash
#SBATCH --job-name=fof-equivalence
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/fof-equivalence-%j.out
#SBATCH --error=/gpfs/kjhan/fof-equivalence-%j.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?submit this script from the lagRamses repository root}
test_dir="$repo_root/tests/sink/smbh-capture-restart"
source_file="$repo_root/patch/lagRamses/kjhan.f90"
work_base=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
work_dir=$(mktemp -d "$work_base/fof-equivalence-${SLURM_JOB_ID:?}.XXXXXX")

[[ -f "$source_file" && -f "$test_dir/fof_backend_equivalence_driver.f90" ]]
[[ ${SLURM_NTASKS:?} -eq 2 ]]
printf 'node=%s\nwork_dir=%s\n' "$(hostname -s)" "$work_dir"

module purge
module load intel/tbb/2022.3 intel/umf/1.0.2 \
  intel/compiler-rt/2025.3.0 intel/compiler/2025.3.0 intel/mpi/2021.17

sha256sum "$source_file" "$test_dir/fof_backend_equivalence_driver.f90"
mpiifx -fpp -qopenmp -D_KJHAN_TEST -DNDIM=3 -O0 -g -check bounds \
  "$source_file" "$test_dir/fof_backend_equivalence_driver.f90" \
  -o "$work_dir/fof_backend_equivalence"

export OMP_NUM_THREADS=2
export OMP_DYNAMIC=FALSE
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
srun --mpi=pmi2 --ntasks=2 --cpus-per-task=2 --cpu-bind=cores \
  "$work_dir/fof_backend_equivalence"
