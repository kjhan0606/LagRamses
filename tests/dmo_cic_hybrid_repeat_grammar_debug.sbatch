#!/bin/bash
#SBATCH --job-name=cic_repeat
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --mem=80G
#SBATCH --time=00:12:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_repeat_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_repeat_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
atomic_exe=${ATOMIC_EXE:?ATOMIC_EXE is required}
colored_exe=${COLORED_EXE:?COLORED_EXE is required}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_repeat_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_arm () {
  variant=$1
  nmpi=$2
  threads=$3
  repeat=$4
  arm=${variant}_r${nmpi}t${threads}_${repeat}
  run=${root}/${arm}
  if [[ ${variant} == atomic ]]; then
    exe=${atomic_exe}
  else
    exe=${colored_exe}
  fi
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
  srun --exclusive --kill-on-bad-exit=1 --mpi=pmi2 -n "${nmpi}" -c "${threads}" \
    --cpu-bind=cores "${exe}" cosmo.nml > run.log
  grep -q 'Run completed' run.log
  grep -q 'FFTW3 CPU direct Poisson solver enabled' run.log
  if grep -Eiq 'FATAL|SIG(SEGV|ABRT)|problem in cic' run.log; then
    echo "${arm}: runtime failure" >&2
    exit 2
  fi
}

for config in '8 8' '16 4' '32 2'; do
  read -r nmpi threads <<< "${config}"
  # ABBA ordering brackets both candidate measurements with the reference.
  run_arm atomic "${nmpi}" "${threads}" a
  run_arm colored "${nmpi}" "${threads}" a
  run_arm colored "${nmpi}" "${threads}" b
  run_arm atomic "${nmpi}" "${threads}" b
done

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
configs = ((8, 8), (16, 4), (32, 2))
variants = ("atomic", "colored")
repeats = ("a", "b")
data = {}
for nmpi, threads in configs:
    for variant in variants:
        for repeat in repeats:
            arm = f"{variant}_r{nmpi}t{threads}_{repeat}"
            text = (root / arm / "run.log").read_text(errors="replace")
            elapsed = [float(value) for value in elapsed_re.findall(text)]
            rho = [float(value) for value in rho_re.findall(text)]
            total = [float(value) for value in total_re.findall(text)]
            states = state_re.findall(text)
            if len(elapsed) != 1 or not rho or not total or not states:
                raise SystemExit(f"{arm}: incomplete timing/state evidence")
            data[(variant, nmpi, threads, repeat)] = {
                "elapsed": elapsed[-1], "rho": rho[-1],
                "total": total[-1], "state": states[-1]
            }

states = {entry["state"] for entry in data.values()}
if len(states) != 1:
    raise SystemExit(f"printed final conservation state differs: {states}")

with (root / "summary.txt").open("w") as stream:
    for nmpi, threads in configs:
        means = {}
        for variant in variants:
            totals = [data[(variant, nmpi, threads, repeat)]["total"] for repeat in repeats]
            rhos = [data[(variant, nmpi, threads, repeat)]["rho"] for repeat in repeats]
            means[variant] = (statistics.mean(totals), statistics.mean(rhos))
            stream.write(
                f"{variant}_r{nmpi}t{threads}_timer_total_mean_s {means[variant][0]:.6f}\n"
                f"{variant}_r{nmpi}t{threads}_timer_total_range_s "
                f"{min(totals):.6f} {max(totals):.6f}\n"
                f"{variant}_r{nmpi}t{threads}_rho_mean_s {means[variant][1]:.6f}\n"
                f"{variant}_r{nmpi}t{threads}_rho_range_s "
                f"{min(rhos):.6f} {max(rhos):.6f}\n"
            )
        stream.write(
            f"r{nmpi}t{threads}_timer_total_speedup "
            f"{means['atomic'][0] / means['colored'][0]:.6f}\n"
            f"r{nmpi}t{threads}_rho_speedup "
            f"{means['atomic'][1] / means['colored'][1]:.6f}\n"
        )
    stream.write("final_conservation_state " + " ".join(next(iter(states))) + "\n")
print("DMO_CIC_HYBRID_REPEAT_PASS")
PY
