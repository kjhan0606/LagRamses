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
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export I_MPI_PIN=0

code=/home/kjhan/BACKUP/lagRamses-capture-ledger
input="$code/tests/sink/capture-ledger"
format=${CAPTURE_FORMAT:-original}
case_id=${CAPTURE_CASE:-standard}
case "$case_id" in
    standard) source_nml="$input/capture-ledger.nml" ;;
    uniform) source_nml="$input/capture-ledger-uniform.nml" ;;
    *) echo "FATAL: unsupported CAPTURE_CASE=$case_id"; exit 2 ;;
esac
case "$format" in
    original)
        if [ "$case_id" = standard ]; then
            run=/gpfs/kjhan/lagramses_capture_restart_${SLURM_JOB_ID:?}
        else
            run=/gpfs/kjhan/lagramses_capture_uniform_restart_${SLURM_JOB_ID:?}
        fi ;;
    hdf5)
        if [ "$case_id" = standard ]; then
            run=/gpfs/kjhan/lagramses_capture_hdf5_restart_${SLURM_JOB_ID:?}
        else
            run=/gpfs/kjhan/lagramses_capture_uniform_hdf5_restart_${SLURM_JOB_ID:?}
        fi ;;
    *) echo "FATAL: unsupported CAPTURE_FORMAT=$format"; exit 2 ;;
esac
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp "$input/ic_sink" "$run/ic_sink"
cp "$input/legacy_yields.dat" "$run/legacy_yields.dat"
cp "$code/bin/ramses_final3d" "$run/ramses_final3d"
if [ "$format" = hdf5 ]; then
    ldd "$run/ramses_final3d" | grep -F libhdf5_fortran || {
        echo 'FATAL: requested HDF5 restart with a non-HDF5 binary'
        exit 2
    }
fi
cp "$code/patch/lagRamses/aux/validate_smbh_capture_ledger.py" "$run/validate_ledger.py"
cp "$input/verify_restart_smoke.py" "$run/verify_restart_smoke.py"
# Mechanical materialization of the committed source namelist.  In read_params,
# positive tend overwrites tout, so both values must change together.  tout(1)=0
# requests one full initial output; tout(2) is an unreachable dummy.  Explicitly
# setting both entries prevents the default tout(2)=0 from creating another dump.
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=0.0d0,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    "$source_nml" >"$run/fresh.nml"
sed -e 's/noutput=1/noutput=2/' -e 's/aout=1.1/aout=1.1,1.1/' \
    -e 's/tout=1.0d100/tout=0.0d0,1.0d100/' -e 's/tend=1.0d100/tend=0.0d0/' \
    -e 's/nrestart=0/nrestart=1/' \
    "$source_nml" >"$run/restart.nml"
if [ "$format" = hdf5 ]; then
    for nml in "$run/fresh.nml" "$run/restart.nml"; do
        awk '/^fbackup=100000$/ {print; print "informat=\"hdf5\""; print "outformat=\"hdf5\""; next} {print}' \
            "$nml" >"$nml.tmp"
        mv -- "$nml.tmp" "$nml"
    done
fi
echo "purpose=${SLURM_NTASKS}-rank hydrodynamic checkpoint/restart structural smoke"
echo "checkpoint_format=$format"
echo "capture_case=$case_id"
echo "source_namelist=$source_nml"
echo "fresh_effective_namelist=$run/fresh.nml"
echo "restart_effective_namelist=$run/restart.nml"
echo "run_directory=$run"
echo "git_commit=$(git -C "$code" rev-parse HEAD)"
sha256sum "$run/fresh.nml" "$run/restart.nml" "$run/ic_sink" \
    "$run/legacy_yields.dat" "$run/ramses_final3d" \
    "$run/validate_ledger.py" "$run/verify_restart_smoke.py"
cd "$run"
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none \
    -n "${SLURM_NTASKS}" ./ramses_final3d fresh.nml >fresh.log 2>&1
test -f output_00001/COMPLETE
if [ "$format" = hdf5 ]; then
    test -f output_00001/data_00001.h5
else
    test ! -e output_00001/data_00001.h5
fi
test ! -e output_00002
python3 validate_ledger.py capture.jsonl >fresh_validation.json
python3 -c 'import json; d=json.load(open("fresh_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["committed_batches"]==1, d'
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none \
    -n "${SLURM_NTASKS}" ./ramses_final3d restart.nml >restart.log 2>&1
test ! -e output_00002
python3 validate_ledger.py capture.jsonl >restart_validation.json
python3 -c 'import json; d=json.load(open("restart_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["superseded_batches"]==1 and d["run_attempts"]==2, d'
python3 verify_restart_smoke.py "$run" >restart_evaluation.json
if [ "$case_id" = uniform ]; then
    python3 -c 'import json; d=json.load(open("restart_evaluation.json")); assert d["uniform_static_roundoff_gate_passed"] and d["hydro_replay_within_tolerance"] and d["active_multiple_events"]==1, d'
    echo 'CAPTURE_UNIFORM_STATIC_ROUNDOFF_PASS_NOT_GALAXY_MODEL'
elif [ "$format" = hdf5 ]; then
    python3 -c 'import json; d=json.load(open("restart_evaluation.json")); assert d["checkpoint_format"]=="hdf5" and d["replayed_event_records_identical"] and d["active_multiple_events"]==1, d'
    echo 'CAPTURE_HDF5_LEDGER_STRUCTURAL_PASS_HYDRO_NOT_ADMITTED'
else
    echo 'CAPTURE_RESTART_SMOKE_PASS'
fi
