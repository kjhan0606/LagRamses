#!/bin/bash
# Synthetic cosmological DMO gate for the pre-update_time restart energy history.
#SBATCH --job-name=smbh_capture_cosmo_energy
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_cosmo_energy_%j.log

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
binary=${CAPTURE_TEST_BINARY:?set CAPTURE_TEST_BINARY to the HDF5 probe executable}
test -x "$binary"
run=/gpfs/kjhan/lagramses_capture_cosmo_energy_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing run directory $run"
    exit 2
fi
mkdir "$run"
cp "$binary" "$run/ramses_final3d"
cp "$input/compare_cosmo_energy_restart.py" "$run/compare_cosmo_energy_restart.py"
python3 "$code/tests/make_cuda_ndgp_ic.py" "$run/ic" \
    --profile matched-coarse-v3 --omega-b-header 0.05 >"$run/ic_generation.log"
for mode in fresh restart; do
    restart_index=0
    if [ "$mode" = restart ]; then restart_index=2; fi
    sed -e "s|__IC_ROOT__|$run/ic|g" \
        -e "s|__NRESTART__|$restart_index|g" \
        "$input/cosmo-energy-restart.nml.in" >"$run/$mode.nml"
done
echo 'purpose=synthetic cosmological DMO restart energy-history regression'
echo "fresh_effective_namelist=$run/fresh.nml"
echo "restart_effective_namelist=$run/restart.nml"
echo "binary_source=$binary"
echo 'full_dump_policy=noutput=1 unreachable aout=1.1/tout=1e100; foutput=1 outputs coarse steps 1-5; fbackup=100000'
echo 'expected_storage=at most five fresh and three restarted small two-level HDF5 dumps, below 1 GiB; free space checked before submission'
sha256sum "$run/fresh.nml" "$run/restart.nml" "$run/ramses_final3d" \
    "$run/compare_cosmo_energy_restart.py" "$run/ic/MANIFEST.sha256"
cd "$run"
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 \
    ./ramses_final3d fresh.nml >fresh.log 2>&1
for number in 00001 00002 00003 00004 00005; do
    test -f "output_${number}/COMPLETE"
    test -f "output_${number}/data_${number}.h5"
done
grep -Eq '^nstep_coarse=[[:space:]]*2[[:space:]]*$' output_00002/info_00002.txt
grep -Fq ' IC/namelist omega_b consistent:' fresh.log
if grep -Fq 'WARN: Fine multigrid Poisson failed to converge' fresh.log; then
    echo 'FATAL: fine Poisson nonconvergence in fresh leg'
    exit 1
fi
sha256sum output_00002/data_00002.h5 >checkpoint.sha256
for number in 00003 00004 00005; do
    mv "output_${number}" "fresh_output_${number}"
done
srun --kill-on-bad-exit=1 --mpi=pmi2 --cpu-bind=none -n 2 \
    ./ramses_final3d restart.nml >restart.log 2>&1
sha256sum -c checkpoint.sha256
if grep -Fq 'WARN: Fine multigrid Poisson failed to converge' restart.log; then
    echo 'FATAL: fine Poisson nonconvergence in restart leg'
    exit 1
fi
for number in 00003 00004 00005; do
    test -f "output_${number}/COMPLETE"
    python3 compare_cosmo_energy_restart.py \
        "fresh_output_${number}/data_${number}.h5" \
        "output_${number}/data_${number}.h5" \
        >"comparison_${number}.json"
done
python3 - <<'PY'
import json
import os
from pathlib import Path

passes = []
for number in ('00003', '00004', '00005'):
    report = json.loads(Path(f'comparison_{number}.json').read_text())
    passes.append(report['replay_within_roundoff'])
    print(number, 'replay_within_roundoff=', report['replay_within_roundoff'],
          'epot_tot_int=', report['headers']['epot_tot_int']['reference'][0])
if os.environ.get('CAPTURE_EXPECT_REPLAY_FAIL') == '1':
    assert not all(passes), 'old binary unexpectedly passed the energy replay gate'
else:
    assert all(passes), 'cosmological restart differs beyond roundoff'
PY
if [ "${CAPTURE_EXPECT_REPLAY_FAIL:-0}" = 1 ]; then
    echo 'CAPTURE_OLD_COSMO_ENERGY_REPLAY_REJECTED'
else
    echo 'CAPTURE_COSMO_ENERGY_HISTORY_WITHIN_ROUNDOFF'
fi
