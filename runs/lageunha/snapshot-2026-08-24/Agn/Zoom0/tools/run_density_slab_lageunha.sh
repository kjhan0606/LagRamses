#!/usr/bin/env bash
set -euo pipefail

run_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm
figure_dir=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/figures/zstart_comparison
executable=/home/kjhan/BACKUP/lagRamses-de-nonstd/bin/ramses_density_slab3d
prefix=${figure_dir}/amr_cdm_zstart100_z3p97
pid_file=${figure_dir}/amr_cdm_zstart100_z3p97.pid

mkdir -p "${figure_dir}"
echo "$$" > "${pid_file}"
trap 'rm -f "${pid_file}"' EXIT

cd "${run_dir}"
export OMP_NUM_THREADS=1
export OMP_STACKSIZE=128M
export I_MPI_FABRICS=shm
export RAMSES_LBOX_MPC_H=128
export RAMSES_SLAB_PREFIX="${prefix}"

mpirun -np 32 "${executable}" zoom_cdm_restart3.nml 4
