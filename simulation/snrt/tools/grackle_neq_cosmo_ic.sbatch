#!/usr/bin/env bash
# One-step actual cosmological IC test; chemistry remains CPU based.
#SBATCH --job-name=grackle-ic128
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-ic128-%j.log
set -euo pipefail
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
ic=/gpfs/kjhan/LRD_JWST/.galaxy-preenriched.IJ8aeX/ics/level_007
run=$root/.grackle-neq-cosmo-ic-${SLURM_JOB_ID}
test ! -e "$run"
for name in ic_deltab ic_velcx ic_velcy ic_velcz ic_poscx ic_poscy ic_poscz \
            ic_velbx ic_velby ic_velbz ic_pvar_00001 ic_pvar_00012; do
  test -r "$ic/$name"
done
mkdir "$run"
cp "$root/simulation/snrt/config/grackle_neq_cosmo_ic_smoke.nml" "$run/effective.nml"
stat -c '%n %s bytes' "$ic/ic_deltab" "$ic/ic_pvar_00001" "$ic/ic_pvar_00012"
sha256sum "$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d" > "$run/binary.sha256"
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
  intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close
cd "$run"
cp effective.nml effective.nml.launched
TIMEFORMAT='wall=%3R user=%3U sys=%3S'
time timeout 2400s srun --ntasks=1 --cpus-per-task=8 \
  "$root/.grackle-neq-build-lageunha-20261006T1943/bin/ramses_grackle3d" effective.nml \
  > run.log 2>&1
grep -E 'Grackle NEQ:|Reading file|Initializing|Run completed|Total running time|NaN_CHK|Fine step=' run.log
if grep -E 'ERROR:|Aborting|failed|NaN_CHK.*=[[:space:]]*[1-9][0-9]*|WARN: Fine multigrid Poisson failed' run.log; then exit 1; fi
test -d output_00001 && { echo 'Unexpected raw snapshot'; exit 1; }
grep -q 'Run completed' run.log
