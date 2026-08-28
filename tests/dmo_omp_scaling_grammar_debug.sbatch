#!/bin/bash
#SBATCH --job-name=dmo_omp8
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --mem=80G
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp8_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp8_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
exe=${RAMSES_TEST_EXE:-${code}/bin/ramses_final3d}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp8_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_arm () {
  arm=$1
  threads=$2
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

  export OMP_NUM_THREADS=${threads}
  cd "${run}"
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n 4 -c "${threads}" --cpu-bind=cores \
    "${exe}" cosmo.nml > run.log
  grep -q 'Run completed' run.log
  grep -q 'FFTW3 CPU direct Poisson solver enabled' run.log
  if grep -Eiq 'FATAL|SIG(SEGV|ABRT)|problem in cic' run.log; then
    echo "${arm}: runtime failure" >&2
    exit 2
  fi
}

# ABBA order brackets the OMP=8 measurements with OMP=4 controls.
run_arm t4_a 4
run_arm t8_a 8
run_arm t8_b 8
run_arm t4_b 4
run_arm t1 1

python3 - "${root}" <<'PY'
import pathlib
import re
import statistics
import sys

root = pathlib.Path(sys.argv[1])
elapsed_re = re.compile(r"Time elapsed since last coarse step:\s+([0-9.]+) s")
rho_re = re.compile(r"^\s*([0-9.]+)\s+[0-9.]+\s+[0-9.]+.*\brho\s*$", re.M)
total_re = re.compile(r"^\s*([0-9.]+)\s+100\.0\s+TOTAL\s*$", re.M)
state_re = re.compile(
    r"Mtot=\s*([^ ]+)\s+Ekin_tot=\s*([^ ]+)\s+Epop_tot=\s*([^ ]+)\s+Eint_tot=\s*([^ ]+)"
)
arms = ("t1", "t4_a", "t8_a", "t8_b", "t4_b")
data = {}
for arm in arms:
    text = (root / arm / "run.log").read_text(errors="replace")
    elapsed = [float(value) for value in elapsed_re.findall(text)]
    rho = [float(value) for value in rho_re.findall(text)]
    totals = [float(value) for value in total_re.findall(text)]
    states = state_re.findall(text)
    if len(elapsed) != 1 or not rho or not totals or not states:
        raise SystemExit(f"{arm}: incomplete timing/state evidence")
    data[arm] = {
        "timer_total": totals[-1], "restart_elapsed": elapsed[-1],
        "rho": rho[-1], "state": states[-1]
    }

states = {data[name]["state"] for name in arms}
if len(states) != 1:
    raise SystemExit(f"printed final conservation state differs: {states}")

t1 = data["t1"]["timer_total"]
t4 = statistics.mean(data[name]["timer_total"] for name in ("t4_a", "t4_b"))
t8 = statistics.mean(data[name]["timer_total"] for name in ("t8_a", "t8_b"))
rho4 = statistics.mean(data[name]["rho"] for name in ("t4_a", "t4_b"))
rho8 = statistics.mean(data[name]["rho"] for name in ("t8_a", "t8_b"))

with (root / "summary.txt").open("w") as stream:
    for arm in arms:
        stream.write(f"{arm}_timer_total_seconds {data[arm]['timer_total']:.6f}\n")
        stream.write(f"{arm}_restart_elapsed_seconds {data[arm]['restart_elapsed']:.6f}\n")
        stream.write(f"{arm}_rho_seconds {data[arm]['rho']:.6f}\n")
    stream.write(f"omp1_timer_total_seconds {t1:.6f}\n")
    stream.write(f"omp4_mean_timer_total_seconds {t4:.6f}\n")
    stream.write(f"omp8_mean_timer_total_seconds {t8:.6f}\n")
    stream.write(f"omp4_speedup {t1 / t4:.6f}\n")
    stream.write(f"omp8_speedup {t1 / t8:.6f}\n")
    stream.write(f"omp4_efficiency_percent {100.0 * t1 / (4.0 * t4):.3f}\n")
    stream.write(f"omp8_efficiency_percent {100.0 * t1 / (8.0 * t8):.3f}\n")
    stream.write(f"omp4_to_omp8_efficiency_percent {100.0 * t4 / (2.0 * t8):.3f}\n")
    stream.write(f"omp4_to_omp8_speedup {t4 / t8:.6f}\n")
    stream.write(f"omp4_to_omp8_wall_reduction_percent {(1.0 - t8 / t4) * 100.0:.3f}\n")
    stream.write(f"rho4_mean_seconds {rho4:.6f}\n")
    stream.write(f"rho8_mean_seconds {rho8:.6f}\n")
    stream.write("final_conservation_state " + " ".join(next(iter(states))) + "\n")
print("DMO_OMP8_SCALING_PASS")
PY
