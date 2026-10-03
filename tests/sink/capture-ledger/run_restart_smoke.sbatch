#!/bin/bash
# One initial checkpoint, then the same MULTIPLE capture before and after restart.
#SBATCH --job-name=smbh_capture_restart
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_restart_%j.log

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
run=/gpfs/kjhan/lagramses_capture_restart_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp "$input/ic_sink" "$run/ic_sink"
cp "$input/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$code/bin/ramses_final3d" "$run/ramses_final3d"
cp "$code/patch/lagRamses/aux/validate_smbh_capture_ledger.py" "$run/validate_ledger.py"
cp "$input/verify_restart_smoke.py" "$run/verify_restart_smoke.py"
# Mechanical materialization of the committed source namelist.  In read_params,
# positive tend overwrites tout, so both values must change together.  tout(1)=0
# requests one full initial output; tout(2) is an unreachable dummy.  Explicitly
# setting both entries prevents the default tout(2)=0 from creating another dump.
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=0.0d0,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    "$input/capture-ledger.nml" >"$run/fresh.nml"
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=0.0d0,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    -e 's/nrestart=0/nrestart=1/' \
    "$input/capture-ledger.nml" >"$run/restart.nml"
echo "purpose=one-rank hydrodynamic checkpoint/restart structural smoke"
echo "fresh_effective_namelist=$run/fresh.nml"
echo "restart_effective_namelist=$run/restart.nml"
echo "run_directory=$run"
echo "git_commit=$(git -C "$code" rev-parse HEAD)"
sha256sum "$run/fresh.nml" "$run/restart.nml" "$run/ic_sink" \
    "$run/legacy_yields.dat" "$run/ramses_final3d" \
    "$run/validate_ledger.py" "$run/verify_restart_smoke.py"
cd "$run"
srun --kill-on-bad-exit=1 -n 1 ./ramses_final3d fresh.nml >fresh.log 2>&1
test -f output_00001/COMPLETE
test ! -e output_00002
python3 validate_ledger.py capture.jsonl >fresh_validation.json
python3 -c 'import json; d=json.load(open("fresh_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["committed_batches"]==1, d'
srun --kill-on-bad-exit=1 -n 1 ./ramses_final3d restart.nml >restart.log 2>&1
test ! -e output_00002
python3 validate_ledger.py capture.jsonl >restart_validation.json
python3 -c 'import json; d=json.load(open("restart_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["superseded_batches"]==1 and d["run_attempts"]==2, d'
python3 verify_restart_smoke.py "$run" >restart_evaluation.json
echo 'CAPTURE_RESTART_SMOKE_PASS'
