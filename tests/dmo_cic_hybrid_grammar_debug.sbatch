#!/bin/bash
#SBATCH --job-name=cic_hybrid
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --exclusive
#SBATCH --mem=80G
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_hybrid_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_hybrid_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
atomic_exe=${ATOMIC_EXE:?ATOMIC_EXE is required}
colored_exe=${COLORED_EXE:?COLORED_EXE is required}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_hybrid_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

run_arm () {
  variant=$1
  nmpi=$2
  threads=$3
  label=r${nmpi}t${threads}
  run=${root}/${variant}_${label}
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
    echo "${variant}_${label}: runtime failure" >&2
    exit 2
  fi
}

configs=('8 1' '8 4' '8 8' '16 1' '16 2' '16 4' '32 1' '32 2')
index=0
for config in "${configs[@]}"; do
  read -r nmpi threads <<< "${config}"
  # Alternate execution order so monotonic node drift does not favor one path.
  if (( index % 2 == 0 )); then
    run_arm atomic "${nmpi}" "${threads}"
    run_arm colored "${nmpi}" "${threads}"
  else
    run_arm colored "${nmpi}" "${threads}"
    run_arm atomic "${nmpi}" "${threads}"
  fi
  index=$((index+1))
done

python3 - "${root}" <<'PY'
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
step_re = re.compile(r"Time elapsed since last coarse step:\s+([0-9.]+) s")
rho_re = re.compile(r"^\s*([0-9.]+)\s+[0-9.]+\s+[0-9.]+.*\brho\s*$", re.M)
total_re = re.compile(r"^\s*([0-9.]+)\s+100\.0\s+TOTAL\s*$", re.M)
state_re = re.compile(
    r"Mtot=\s*([^ ]+)\s+Ekin_tot=\s*([^ ]+)\s+Epop_tot=\s*([^ ]+)\s+Eint_tot=\s*([^ ]+)"
)
configs = ((8, 1), (8, 4), (8, 8), (16, 1), (16, 2), (16, 4), (32, 1), (32, 2))
variants = ("atomic", "colored")
data = {}
for variant in variants:
    for nmpi, threads in configs:
        arm = f"{variant}_r{nmpi}t{threads}"
        text = (root / arm / "run.log").read_text(errors="replace")
        steps = [float(value) for value in step_re.findall(text)]
        rho = [float(value) for value in rho_re.findall(text)]
        totals = [float(value) for value in total_re.findall(text)]
        states = state_re.findall(text)
        if len(steps) != 1 or not rho or not totals or not states:
            raise SystemExit(f"{arm}: incomplete timing/state evidence")
        data[(variant, nmpi, threads)] = {
            "restart_elapsed": steps[-1], "timer_total": totals[-1],
            "rho": rho[-1], "state": states[-1]
        }

states = {entry["state"] for entry in data.values()}
if len(states) != 1:
    raise SystemExit(f"printed final conservation state differs: {states}")

with (root / "summary.tsv").open("w") as stream:
    stream.write(
        "variant\tnmpi\tthreads\tcores\trestart_elapsed_s\ttimer_total_s\t"
        "rho_s\ttimer_omp_efficiency_percent\trho_omp_efficiency_percent\n"
    )
    for variant in variants:
        for nmpi, threads in configs:
            entry = data[(variant, nmpi, threads)]
            base = data[(variant, nmpi, 1)]
            timer_efficiency = 100.0 * base["timer_total"] / (
                threads * entry["timer_total"]
            )
            rho_efficiency = 100.0 * base["rho"] / (threads * entry["rho"])
            stream.write(
                f"{variant}\t{nmpi}\t{threads}\t{nmpi * threads}\t"
                f"{entry['restart_elapsed']:.6f}\t{entry['timer_total']:.6f}\t"
                f"{entry['rho']:.6f}\t{timer_efficiency:.3f}\t{rho_efficiency:.3f}\n"
            )

with (root / "comparison.txt").open("w") as stream:
    for nmpi, threads in configs:
        old = data[("atomic", nmpi, threads)]
        new = data[("colored", nmpi, threads)]
        stream.write(
            f"r{nmpi}t{threads}_timer_total_speedup "
            f"{old['timer_total'] / new['timer_total']:.6f}\n"
            f"r{nmpi}t{threads}_timer_total_reduction_percent "
            f"{100.0 * (1.0 - new['timer_total'] / old['timer_total']):.3f}\n"
            f"r{nmpi}t{threads}_rho_speedup {old['rho'] / new['rho']:.6f}\n"
        )
    for variant in variants:
        for cores, pairs in (
            (32, ((8, 4), (16, 2), (32, 1))),
            (64, ((8, 8), (16, 4), (32, 2))),
        ):
            values = " ".join(
                f"r{nmpi}t{threads}="
                f"{data[(variant, nmpi, threads)]['timer_total']:.6f}"
                for nmpi, threads in pairs
            )
            stream.write(f"{variant}_fixed_{cores}_cores {values}\n")
    stream.write("final_conservation_state " + " ".join(next(iter(states))) + "\n")
print("DMO_CIC_HYBRID_PASS")
PY
