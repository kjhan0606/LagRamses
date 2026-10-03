#!/bin/bash
# Two-rank step-2 legacy checkpoint ownership discriminator, not a physics pass.
#SBATCH --job-name=smbh_capture_legacy_probe
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_legacy_probe_%j.log

set -euo pipefail
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export I_MPI_FABRICS=shm
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export I_MPI_PIN=0

code=/home/kjhan/BACKUP/lagRamses-capture-ledger
input="$code/tests/sink/capture-ledger"
binary=${CAPTURE_TEST_BINARY:?set CAPTURE_TEST_BINARY to the legacy probe executable}
run=/gpfs/kjhan/lagramses_capture_evolved_legacy_probe_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
test -x "$binary"
mkdir "$run"
cp "$input/ic_sink" "$input/legacy_yields.dat" "$run/"
cp "$binary" "$run/ramses_final3d"
cp "$code/patch/lagRamses/aux/validate_smbh_capture_ledger.py" "$run/validate_ledger.py"
cp "$input/evaluate_restart_ownership.py" "$run/evaluate_restart_ownership.py"
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=3.0d-2,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    "$input/capture-ledger-perturbed.nml" >"$run/fresh.nml"
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=3.0d-2,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    -e 's/nrestart=0/nrestart=1/' \
    "$input/capture-ledger-perturbed.nml" >"$run/restart.nml"
if [ "${CAPTURE_GROUPED_IO:-0}" = 1 ]; then
    sed -i '/^fbackup=100000$/a IOGROUPSIZEREP=1' "$run/fresh.nml" "$run/restart.nml"
    checkpoint_info=output_00001/group_00001/info_00001.txt
    particle_files=(output_00001/group_00001/part_00001.out00001 \
        output_00001/group_00002/part_00001.out00002)
else
    checkpoint_info=output_00001/info_00001.txt
    particle_files=(output_00001/part_00001.out00001 output_00001/part_00001.out00002)
fi
echo "purpose=two-rank evolving-gas legacy restart ownership discriminator; grouped_io=${CAPTURE_GROUPED_IO:-0}"
echo 'full_dump_policy=noutput=2 tout=0.03,1e100 aout=1.1,1.1 foutput=fbackup=100000; only one reachable dump'
echo 'expected_storage=one small legacy dump, below 15 MiB based on the HDF5 fixture; free space checked before submission'
echo "fresh_effective_namelist=$run/fresh.nml"
echo "restart_effective_namelist=$run/restart.nml"
echo "binary_source=$binary"
echo "source_commit=$(git -C "$code" rev-parse HEAD)"
sha256sum "$run/fresh.nml" "$run/restart.nml" "$run/ramses_final3d" "$run/ic_sink" "$run/legacy_yields.dat"
cd "$run"
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 ./ramses_final3d fresh.nml >fresh.log 2>&1
test -f output_00001/COMPLETE
test ! -e output_00002
grep -Eq '^nstep_coarse=[[:space:]]*2[[:space:]]*$' "$checkpoint_info"
if [ "${CAPTURE_GROUPED_IO:-0}" = 1 ]; then
    test -f output_00001/group_00001/sink_00001.out
    test ! -e output_00001/group_00002/sink_00001.out
fi
sha256sum "${particle_files[@]}" >checkpoint_particles.sha256
python3 validate_ledger.py capture.jsonl >fresh_validation.json
python3 -c 'import json; d=json.load(open("fresh_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["run_attempts"]==1, d'
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 ./ramses_final3d restart.nml >restart.log 2>&1
sha256sum -c checkpoint_particles.sha256
test ! -e output_00002
python3 validate_ledger.py capture.jsonl >restart_validation.json
python3 -c 'import json; d=json.load(open("restart_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["committed_batches"]==1 and d["run_attempts"]==2, d'
python3 evaluate_restart_ownership.py fresh.log restart.log >ownership_comparison.json
echo 'CAPTURE_EVOLVED_LEGACY_LEDGER_STRUCTURAL_PASS'
if ! python3 -c 'import json,sys; d=json.load(open("ownership_comparison.json")); sys.exit(0 if d["owned_cloud_count_equal"] and d["sink_stat_equal"] else 1)'; then
    echo 'CAPTURE_EVOLVED_LEGACY_OWNERSHIP_REPLAY_FAIL'
    exit 1
fi
echo 'CAPTURE_EVOLVED_LEGACY_OWNERSHIP_REPLAY_EXACT'
