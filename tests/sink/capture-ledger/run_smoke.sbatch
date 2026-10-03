#!/bin/bash
# One-rank, two-step structural test of native MULTIPLE capture output.
#SBATCH --job-name=smbh_capture_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_smoke_%j.log

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export I_MPI_FABRICS=shm

code=/home/kjhan/BACKUP/lagRamses-capture-ledger
input="$code/tests/sink/capture-ledger"
run=/gpfs/kjhan/lagramses_capture_smoke_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp "$input/capture-ledger.nml" "$run/effective.nml"
cp "$input/ic_sink" "$run/ic_sink"
cp "$input/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$code/bin/ramses_final3d" "$run/ramses_final3d"
cp "$code/patch/lagRamses/aux/validate_smbh_capture_ledger.py" "$run/validate_ledger.py"
echo "purpose=two-step hydrodynamic structural capture smoke"
echo "effective_namelist=$run/effective.nml"
echo "run_directory=$run"
echo "git_commit=$(git -C "$code" rev-parse HEAD)"
sha256sum "$run/effective.nml" "$run/ic_sink" "$run/legacy_yields.dat" \
    "$run/ramses_final3d" "$run/validate_ledger.py"
cd "$run"
srun --kill-on-bad-exit=1 -n 1 ./ramses_final3d effective.nml >run.log 2>&1
python3 validate_ledger.py capture.jsonl >validation.json
python3 -c 'import json; d=json.load(open("validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["unique_events"]==1 and d["committed_batches"]==1, d'
echo 'CAPTURE_SMOKE_PASS'
