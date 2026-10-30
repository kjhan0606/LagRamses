#!/usr/bin/env bash
# Bounded CPU/HDF5 integration smoke for delayed-cooling restart identity.
# Tests same-model restart acceptance and mismatched-t_diss rejection on a tiny grid.
# It does not test feedback yields.
#SBATCH --job-name=dc-hdf5-smoke
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.dc-restart-smoke-%j.log

set -euo pipefail

root=/gpfs/kjhan/chimes-cvode-gpu-20261001
binary=$root/bin/ramses_final3d
source_nml=$root/simulation/snrt/config/grackle_neq_delayed_validation.nml
run=$root/.delayed-cooling-restart-smoke-${SLURM_JOB_ID}

test -x "$binary"
test -r "$source_nml"
test ! -e "$run"
mkdir "$run"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=2 OMP_STACKSIZE=64M OMP_PROC_BIND=close

sed \
  -e 's/nstepmax=8, ncontrol=1/nstepmax=4, ncontrol=1/' \
  -e "s/cooling_method='grackle_neq'/cooling_method='original'/" \
  -e '/grackle_neq_data_file=/d' \
  -e "s/delayed_cooling=.true., neq_chem=.false., dust_mass_enabled=.false./delayed_cooling=.true., delayed_cooling_model='first_event_fixed', t_diss=20d0, neq_chem=.false., dust_mass_enabled=.false./" \
  "$source_nml" > "$run/effective-first.nml"
sed \
  -e 's/nstepmax=4, ncontrol=1/nstepmax=5, ncontrol=1/' \
  -e 's/foutput=4, fbackup=1000000/foutput=1, fbackup=1000000/' \
  -e '/hydro=.true.,/a\ nrestart=1' \
  "$run/effective-first.nml" > "$run/effective-restart.nml"
sed 's/t_diss=20d0/t_diss=21d0/' \
  "$run/effective-restart.nml" > "$run/effective-mismatch.nml"

grep -q "cooling_method='original'" "$run/effective-first.nml"
grep -q "delayed_cooling_model='first_event_fixed', t_diss=20d0" "$run/effective-first.nml"
grep -q 'nrestart=1' "$run/effective-restart.nml"

cd "$run"
mpirun -np 1 "$binary" effective-first.nml > first.log 2>&1
grep -q 'Run completed' first.log
grep -q 'Delayed cooling model: first_event_fixed' first.log
test -s output_00001/data_00001.h5

mpirun -np 1 "$binary" effective-restart.nml > restart.log 2>&1
grep -q 'HDF5 restore AMR:.*output_00001/data_00001.h5' restart.log
grep -q 'Restarting at t=' restart.log
grep -q 'Run completed' restart.log
test -s output_00002/data_00002.h5

set +e
mpirun -np 1 "$binary" effective-mismatch.nml > mismatch.log 2>&1
mismatch_rc=$?
set -e
if [[ $mismatch_rc -eq 0 ]]; then
  echo 'ERROR: mismatched delayed-cooling restart unexpectedly succeeded' >&2
  exit 1
fi
grep -q 'delayed-cooling model/t_diss restart identity mismatch' mismatch.log

if grep -E 'ERROR:|Aborting|Too many errors|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' first.log restart.log; then
  exit 1
fi
sha256sum "$binary" "$source_nml" "$run/effective-first.nml" \
  "$run/effective-restart.nml" "$run/effective-mismatch.nml" "$0" > assets.sha256
cat > summary.txt <<EOF
scope=delayed-cooling HDF5 restart identity smoke; tiny fixed-grid CPU run
build_sha256=$(sha256sum "$binary" | awk '{print $1}')
valid_restart=PASS
mismatched_t_diss_rejected=PASS
mismatch_exit_code=$mismatch_rc
outputs_evaluated=output_00001 output_00002
EOF
cat > cleanup_manifest.txt <<'EOF'
Evaluated raw simulation output directories removed after same-model restart
and mismatched-t_diss rejection checks passed:
  output_00001/
  output_00002/
Retained effective namelists, Slurm stdout, run logs, checksums and summary.
EOF
rm -r -- output_00001 output_00002
printf 'DELAYED_COOLING_HDF5_RESTART_SMOKE_PASS %s\n' "$run"
