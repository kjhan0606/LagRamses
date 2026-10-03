#!/bin/bash
# A pre-schema multi-rank sink checkpoint must not silently replay its partial statistic.
#SBATCH --job-name=smbh_capture_old_sink_stat
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_old_sink_stat_%j.log

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

source_run=${CAPTURE_SOURCE_RUN:?set CAPTURE_SOURCE_RUN to a pre-schema legacy run}
binary=${CAPTURE_TEST_BINARY:?set CAPTURE_TEST_BINARY to the schema-aware legacy executable}
test -f "$source_run/output_00001/COMPLETE"
test -f "$source_run/output_00001/sink_00001.out"
test -f "$source_run/restart.nml"
test -x "$binary"
run=/gpfs/kjhan/lagramses_capture_old_sink_stat_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp -a "$source_run/output_00001" "$run/output_00001"
cp "$source_run/restart.nml" "$run/restart.nml"
cp "$source_run/capture.jsonl" "$run/capture.jsonl"
cp "$source_run/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$source_run/ic_sink" "$run/ic_sink"
cp "$binary" "$run/ramses_final3d"
echo 'purpose=negative two-rank pre-schema sink_stat restart gate'
echo "source_checkpoint=$source_run/output_00001"
echo "effective_namelist=$run/restart.nml"
echo 'full_dump_policy=restart should fail before any new dump; source noutput=2 tout=0.03,1e100 foutput=fbackup=100000'
echo 'expected_storage=one copied small legacy checkpoint, below 15 MiB; no new dump'
sha256sum "$run/restart.nml" "$run/ramses_final3d" "$run/output_00001/sink_00001.out"
cd "$run"
ledger_before=$(sha256sum capture.jsonl | awk '{print $1}')
if srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 \
    ./ramses_final3d restart.nml >restart.log 2>&1; then
    echo 'FATAL: pre-schema multi-rank sink_stat checkpoint was accepted'
    exit 1
fi
grep -Fq 'FATAL: legacy sink_stat lacks global schema for multi-rank restart' restart.log
test ! -e output_00002
ledger_after=$(sha256sum capture.jsonl | awk '{print $1}')
test "$ledger_before" = "$ledger_after"
echo "unchanged_capture_ledger_sha256=$ledger_after"
echo 'LEGACY_OLD_SINK_STAT_REJECTED_PASS'
