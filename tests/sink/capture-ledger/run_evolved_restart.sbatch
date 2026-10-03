#!/bin/bash
# Step-2 HDF5 checkpoint and step-3/4 state replay of mildly perturbed gas.
#SBATCH --job-name=smbh_capture_evolved_h5
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_evolved_hdf5_%j.log

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
binary=${CAPTURE_TEST_BINARY:-$code/bin/ramses_final3d}
run=/gpfs/kjhan/lagramses_capture_evolved_hdf5_restart_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp "$input/ic_sink" "$run/ic_sink"
cp "$input/legacy_yields.dat" "$run/legacy_yields.dat"
test -x "$binary"
cp "$binary" "$run/ramses_final3d"
ldd "$run/ramses_final3d" | grep -F libhdf5_fortran
cp "$code/patch/lagRamses/aux/validate_smbh_capture_ledger.py" "$run/validate_ledger.py"
cp "$input/compare_evolved_hdf5.py" "$run/compare_evolved_hdf5.py"

# The first reachable time is after step 1 sets the conservation baseline.
# Two comparison slots produce dumps at steps 3 and 4; periodic clocks are inactive.
sed -e 's/noutput=1/noutput=3/' -e 's/aout=1.1/aout=1.1,1.1,1.1/' \
    -e 's/tout=1.0d100/tout=3.0d-2,5.5d-2,7.5d-2/' -e 's/tend=1.0d100/tend=0.0d0/' \
    -e 's/nstepmax=4/nstepmax=5/' \
    -e '/^fbackup=100000$/a\informat="hdf5"\noutformat="hdf5"' \
    "$input/capture-ledger-perturbed.nml" >"$run/fresh.nml"
sed -e 's/noutput=1/noutput=3/' -e 's/aout=1.1/aout=1.1,1.1,1.1/' \
    -e 's/tout=1.0d100/tout=3.0d-2,5.5d-2,7.5d-2/' -e 's/tend=1.0d100/tend=0.0d0/' \
    -e 's/nstepmax=4/nstepmax=5/' \
    -e 's/nrestart=0/nrestart=1/' \
    -e '/^fbackup=100000$/a\informat="hdf5"\noutformat="hdf5"' \
    "$input/capture-ledger-perturbed.nml" >"$run/restart.nml"
echo 'purpose=two-rank moving-gas step-2 checkpoint capture restart diagnostic'
echo "fresh_effective_namelist=$run/fresh.nml"
echo "restart_effective_namelist=$run/restart.nml"
echo "run_directory=$run"
echo "binary_source=$binary"
echo "source_commit=$(git -C "$code" rev-parse HEAD)"
echo 'full_dump_policy=noutput=3 tout=0.03,0.055,0.075 aout=1.1,1.1,1.1 foutput=fbackup=100000'
echo 'expected_storage=three HDF5 dumps, about 21 MiB based on the prior fixture; free space checked before submission'
sha256sum "$run/fresh.nml" "$run/restart.nml" "$run/ic_sink" \
    "$run/legacy_yields.dat" "$run/ramses_final3d" "$run/validate_ledger.py" "$run/compare_evolved_hdf5.py"
cd "$run"
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 \
    ./ramses_final3d fresh.nml >fresh.log 2>&1
test -f output_00001/COMPLETE
test -f output_00001/data_00001.h5
test -f output_00002/COMPLETE
test -f output_00002/data_00002.h5
test -f output_00003/COMPLETE
test -f output_00003/data_00003.h5
grep -Eq '^nstep_coarse=[[:space:]]*3[[:space:]]*$' output_00002/info_00002.txt
grep -Eq '^nstep_coarse=[[:space:]]*4[[:space:]]*$' output_00003/info_00003.txt
grep -Eq '^nstep_coarse=[[:space:]]*2[[:space:]]*$' output_00001/info_00001.txt
python3 -c 'import h5py; h=h5py.File("output_00001/data_00001.h5"); a=h["/header"].attrs; assert a["nstep_coarse"].item()==2 and a["mass_tot_0"].item()>0 and a["const"].item()!=0; print("checkpoint_step=2 mass_tot_0=",a["mass_tot_0"]," const=",a["const"])'
python3 -c 'import h5py,numpy as np; h=h5py.File("output_00001/data_00001.h5"); names=("bhspin_1","bhspin_2","bhspin_3","spinmag"); assert all(np.isfinite(h["/sinks/"+n][...]).all() for n in names); print("checkpoint_sink_spin_finite=True")'
sha256sum output_00001/data_00001.h5 >checkpoint.sha256
python3 validate_ledger.py capture.jsonl >fresh_validation.json
python3 -c 'import json; d=json.load(open("fresh_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["committed_batches"]==1 and d["run_attempts"]==1, d'
mv output_00002 fresh_output_00002
mv output_00003 fresh_output_00003
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 \
    ./ramses_final3d restart.nml >restart.log 2>&1
sha256sum -c checkpoint.sha256
test -f output_00002/COMPLETE
test -f output_00002/data_00002.h5
test -f output_00003/COMPLETE
test -f output_00003/data_00003.h5
grep -Eq '^nstep_coarse=[[:space:]]*3[[:space:]]*$' output_00002/info_00002.txt
grep -Eq '^nstep_coarse=[[:space:]]*4[[:space:]]*$' output_00003/info_00003.txt
python3 validate_ledger.py capture.jsonl >restart_validation.json
python3 -c 'import json; d=json.load(open("restart_validation.json")); assert d["status"]=="valid" and d["multiple_events"]==1 and d["committed_batches"]==1 and d["superseded_batches"]==0 and d["run_attempts"]==2, d'
python3 compare_evolved_hdf5.py fresh_output_00002/data_00002.h5 output_00002/data_00002.h5 >evolved_state_comparison.json
python3 compare_evolved_hdf5.py fresh_output_00003/data_00003.h5 output_00003/data_00003.h5 >evolved_state_comparison_step4.json
echo 'CAPTURE_EVOLVED_HDF5_LEDGER_STRUCTURAL_PASS'
if ! python3 -c 'import json,sys; files=("evolved_state_comparison.json","evolved_state_comparison_step4.json"); sys.exit(0 if all(json.load(open(p))["active_state_exact"] for p in files) else 1)'; then
    echo 'CAPTURE_EVOLVED_HDF5_ACTIVE_STATE_REPLAY_FAIL'
    exit 1
fi
echo 'CAPTURE_EVOLVED_HDF5_ACTIVE_STATE_REPLAY_EXACT'
