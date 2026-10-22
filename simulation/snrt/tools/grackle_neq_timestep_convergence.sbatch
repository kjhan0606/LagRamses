#!/usr/bin/env bash
# Focused Grackle stiff-regime timestep-convergence measurement.
#SBATCH --job-name=grackle-neq-dtconv
#SBATCH --partition=a40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:A40:1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-neq-dtconv-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
bridge=$root/.grackle-neq-build-lageunha-20261006T1943/bin/grackle_neq_bridge.o
run=$root/.grackle-neq-timestep-convergence-${SLURM_JOB_ID}
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
test -r "$bridge" && test -r "$uvdata" && test ! -e "$run"
mkdir "$run"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_timestep_convergence_test.c" \
  "$bridge" -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/grackle_neq_timestep_convergence_test"
sha256sum "$run/grackle_neq_timestep_convergence_test" "$bridge" \
  /home/kjhan/local/lib64/libgrackle.a > "$run/assets.sha256"
export OMP_NUM_THREADS=1 OMP_STACKSIZE=64M
TIMEFORMAT='wall=%3R user=%3U sys=%3S'
{ time srun --ntasks=1 --cpus-per-task=1 "$run/grackle_neq_timestep_convergence_test" "$uvdata"; } \
  > "$run/convergence.log" 2> "$run/wall.txt"
printf 'GRACKLE_NEQ_TIMESTEP_CONVERGENCE_COMPLETED %s\n' "$run"
