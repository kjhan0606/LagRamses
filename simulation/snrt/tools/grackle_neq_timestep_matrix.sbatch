#!/usr/bin/env bash
# Focused chemistry accuracy check using the exact bridge object from a fresh build.
# Submit with GRACKLE_BRIDGE_OBJECT=/path/to/fresh/grackle_neq_bridge.o.
#SBATCH --job-name=grackle-dt-matrix
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-dt-matrix-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
bridge=${GRACKLE_BRIDGE_OBJECT:?set GRACKLE_BRIDGE_OBJECT to a fresh-build bridge object}
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
run=$root/.grackle-neq-timestep-matrix-${SLURM_JOB_ID}
test -r "$bridge" && test -r "$uvdata" && test ! -e "$run"
mkdir -p "$run"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 OMP_STACKSIZE=64M OMP_PROC_BIND=close

mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_timestep_convergence_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/timestep_convergence_test"
sha256sum "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_timestep_convergence_test.c" \
  "$bridge" "$uvdata" > "$run/assets.sha256"
srun --ntasks=1 --cpus-per-task=1 "$run/timestep_convergence_test" "$uvdata"
