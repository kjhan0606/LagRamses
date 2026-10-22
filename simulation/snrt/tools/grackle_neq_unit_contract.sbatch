#!/usr/bin/env bash
# Compare the adapter with Grackle's direct Local API at two cosmological
# epochs; cross-epoch physical changes are diagnostic, not an equality gate.
# CPU-only regression; the GPU request is the cluster allocation token.
#SBATCH --job-name=grackle-unit
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-unit-contract-%j.log
set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
run=$root/.grackle-unit-contract-${SLURM_JOB_ID}
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
test -r "$uvdata" && test ! -e "$run"
mkdir -p "$run"

module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=1 OMP_STACKSIZE=64M OMP_PROC_BIND=close

mpicc -O2 -fopenmp -I/home/kjhan/local/include \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_unit_contract_test.c" \
  "$root/patch/lagRamses/grackle_neq_bridge.c" \
  -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib \
  -lhdf5 -lz -lgfortran -lquadmath -lm \
  -Wl,-rpath,/home/kjhan/local/lib64 -Wl,-rpath,/home/kjhan/local/lib \
  -o "$run/unit_contract_test"
sha256sum "$run/unit_contract_test" \
  "$root/simulation/snrt/tests/fixtures/phase0/grackle_neq_unit_contract_test.c" \
  "$root/patch/lagRamses/grackle_neq_bridge.c" "$uvdata" > "$run/assets.sha256"
srun --ntasks=1 --cpus-per-task=1 "$run/unit_contract_test" "$uvdata" \
  | tee "$run/results.txt"
grep -q 'Grackle cosmological unit contract: PASS' "$run/results.txt"
printf 'GRACKLE_NEQ_UNIT_CONTRACT_PASS %s\n' "$run"
