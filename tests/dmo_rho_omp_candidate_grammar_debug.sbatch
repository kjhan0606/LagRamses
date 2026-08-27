#!/bin/bash
#SBATCH --job-name=dmo_rho4
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_rho4_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_rho4_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
: "${RHO_OMP_BASELINE_EXE:?set RHO_OMP_BASELINE_EXE to the pre-change executable}"
baseline=${RHO_OMP_BASELINE_EXE}
candidate=${RHO_OMP_CANDIDATE_EXE:-${code}/bin/ramses_final3d}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_rho4_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_arm () {
  arm=$1
  executable=$2
  threads=$3
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
    "${executable}" cosmo.nml > run.log
  grep -q 'Run completed' run.log
  grep -q 'FFTW3 CPU direct Poisson solver enabled' run.log
  if grep -Eiq 'FATAL|SIG(SEGV|ABRT)|problem in cic' run.log; then
    echo "${arm}: runtime failure" >&2
    exit 2
  fi
}

run_arm candidate_t4_a "${candidate}" 4
run_arm baseline_t4_a "${baseline}" 4
run_arm candidate_t4_b "${candidate}" 4
run_arm baseline_t4_b "${baseline}" 4
run_arm candidate_t1 "${candidate}" 1

python3 - "${root}" <<'PY'
import pathlib
import re
import statistics
import sys

root = pathlib.Path(sys.argv[1])
step_re = re.compile(r"Time elapsed since last coarse step:\s+([0-9.]+) s")
rho_re = re.compile(r"^\s*([0-9.]+)\s+[0-9.]+\s+[0-9.]+.*\brho\s*$", re.M)
state_re = re.compile(
    r"Mtot=\s*([^ ]+)\s+Ekin_tot=\s*([^ ]+)\s+Epop_tot=\s*([^ ]+)\s+Eint_tot=\s*([^ ]+)"
)
arms = ("candidate_t1", "baseline_t4_a", "candidate_t4_a", "candidate_t4_b", "baseline_t4_b")
data = {}
for arm in arms:
    text = (root / arm / "run.log").read_text(errors="replace")
    steps = [float(value) for value in step_re.findall(text)]
    rho = [float(value) for value in rho_re.findall(text)]
    states = state_re.findall(text)
    if len(steps) != 1 or len(rho) < 1 or not states:
        raise SystemExit(f"{arm}: incomplete timing/state evidence")
    data[arm] = {
        "step": steps[-1],
        "rho": rho[-1],
        "state": states[-1],
    }

base4 = statistics.mean(data[name]["step"] for name in ("baseline_t4_a", "baseline_t4_b"))
cand4 = statistics.mean(data[name]["step"] for name in ("candidate_t4_a", "candidate_t4_b"))
base_rho = statistics.mean(data[name]["rho"] for name in ("baseline_t4_a", "baseline_t4_b"))
cand_rho = statistics.mean(data[name]["rho"] for name in ("candidate_t4_a", "candidate_t4_b"))
cand1 = data["candidate_t1"]["step"]

states = {data[name]["state"] for name in arms}
if len(states) != 1:
    raise SystemExit(f"printed final conservation state differs: {states}")

with (root / "summary.txt").open("w") as stream:
    for arm in arms:
        stream.write(f"{arm}_mean_step_seconds {data[arm]['step']:.6f}\n")
        stream.write(f"{arm}_mean_rho_seconds {data[arm]['rho']:.6f}\n")
    stream.write(f"baseline_t4_mean_step_seconds {base4:.6f}\n")
    stream.write(f"candidate_t4_mean_step_seconds {cand4:.6f}\n")
    stream.write(f"candidate_wall_gain_percent {(base4 / cand4 - 1.0) * 100.0:.3f}\n")
    stream.write(f"candidate_rho_gain_percent {(base_rho / cand_rho - 1.0) * 100.0:.3f}\n")
    stream.write(f"candidate_omp4_efficiency_percent {100.0 * cand1 / (4.0 * cand4):.3f}\n")
    stream.write("final_conservation_state " + " ".join(next(iter(states))) + "\n")
print("DMO_RHO_OMP_CANDIDATE_PASS")
PY
