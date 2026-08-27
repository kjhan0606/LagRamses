#!/bin/bash
#SBATCH --job-name=dmo_wait4
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_wait4_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_wait4_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
exe=${RAMSES_TEST_EXE:-${code}/bin/ramses_final3d}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_wait4_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

for arm in default_a active_a active_b default_b; do
  unset OMP_WAIT_POLICY KMP_BLOCKTIME
  if [[ "${arm}" == active_* ]]; then
    export OMP_WAIT_POLICY=active
    export KMP_BLOCKTIME=infinite
  fi

  run=${root}/${arm}
  mkdir -p "${run}"
  ln -s "${restart}" "${run}/output_00002"
  cp "${paper}/j02_species_boundary256_one.nml" "${run}/cosmo.nml"
  sed -i 's/nrestart=0/nrestart=2/; s/nstepmax=100000/nstepmax=66/; s/nremap=10/nremap=1000000/' \
    "${run}/cosmo.nml"
  sed -i '/use_fftw=.true./a timer_report_interval=1' "${run}/cosmo.nml"
  sed -i 's/aout=0.07692307692307692/aout=2.0/; s/match_aout=.true./match_aout=.false./' \
    "${run}/cosmo.nml"
  sed -i 's/foutput=1000000/foutput=100000000/; s/fbackup=1000000/fbackup=100000000/' \
    "${run}/cosmo.nml"
  : > "${run}/jobcontrol.txt"

  cd "${run}"
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n 4 -c 4 --cpu-bind=cores \
    "${exe}" cosmo.nml > run.log
  grep -q 'Run completed' run.log
done

python3 - "${root}" <<'PY'
import pathlib
import re
import statistics
import sys

root = pathlib.Path(sys.argv[1])
pattern = re.compile(r"Time elapsed since last coarse step:\s+([0-9.]+) s")
values = {}
for arm in ("default_a", "active_a", "active_b", "default_b"):
    matches = pattern.findall((root / arm / "run.log").read_text(errors="replace"))
    if len(matches) != 1:
        raise SystemExit(f"{arm}: expected one timed coarse step, found {len(matches)}")
    values[arm] = float(matches[0])
default = statistics.mean(values[name] for name in ("default_a", "default_b"))
active = statistics.mean(values[name] for name in ("active_a", "active_b"))
gain = (default / active - 1.0) * 100.0
with (root / "summary.txt").open("w") as stream:
    for arm, value in values.items():
        stream.write(f"{arm}_step_seconds {value:.6f}\n")
    stream.write(f"default_mean_seconds {default:.6f}\n")
    stream.write(f"active_mean_seconds {active:.6f}\n")
    stream.write(f"active_wall_gain_percent {gain:.3f}\n")
print("DMO_OMP4_WAIT_POLICY_PASS")
PY
