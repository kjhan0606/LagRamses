#!/bin/bash
#SBATCH --job-name=dmo_omp4
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=00:20:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp4_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp4_%j.err
set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-DE/code
paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
exe=${RAMSES_TEST_EXE:-${code}/bin/ramses_final3d}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_omp4_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

for arm in mg_t1 mg_t4 fftw_t1 fftw_t4; do
  case "${arm}" in
    *_t1) threads=1 ;;
    *_t4) threads=4 ;;
  esac
  case "${arm}" in
    mg_*) fftw=.false. ;;
    fftw_*) fftw=.true. ;;
  esac

  run=${root}/${arm}
  mkdir -p "${run}"
  ln -s "${restart}" "${run}/output_00002"
  cp "${paper}/j02_species_boundary256_one.nml" "${run}/cosmo.nml"
  sed -i 's/nrestart=0/nrestart=2/; s/nstepmax=100000/nstepmax=68/; s/nremap=10/nremap=1000000/' \
    "${run}/cosmo.nml"
  sed -i "s/use_fftw=.true./use_fftw=${fftw}/" "${run}/cosmo.nml"
  sed -i '/use_fftw=/a timer_report_interval=1' "${run}/cosmo.nml"
  sed -i 's/aout=0.07692307692307692/aout=2.0/; s/match_aout=.true./match_aout=.false./' \
    "${run}/cosmo.nml"
  sed -i 's/foutput=1000000/foutput=100000000/; s/fbackup=1000000/fbackup=100000000/' \
    "${run}/cosmo.nml"
  : > "${run}/jobcontrol.txt"

  export OMP_NUM_THREADS=${threads}
  cd "${run}"
  start=$(date +%s.%N)
  srun --kill-on-bad-exit=1 --mpi=pmi2 -n 4 -c "${threads}" --cpu-bind=cores \
    "${exe}" cosmo.nml > run.log
  stop=$(date +%s.%N)
  awk -v start="${start}" -v stop="${stop}" 'BEGIN{printf "ELAPSED %.6f\n",stop-start}' > elapsed.txt
  grep -q 'Run completed' run.log
  if grep -Eiq 'FATAL|NaN_CHK.*=[[:space:]]*[1-9]|No more free|SIG(SEGV|ABRT)' run.log; then
    echo "${arm}: numerical or runtime failure" >&2
    exit 2
  fi
done

python3 - "${root}" <<'PY'
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
step_re = re.compile(r"Time elapsed since last coarse step:\s+([0-9.]+) s")
main_re = re.compile(r"Main step=\s*(\d+).*mcons=\s*([^ ]+).*econs=\s*([^ ]+)")
results = {}
for arm in ("mg_t1", "mg_t4", "fftw_t1", "fftw_t4"):
    log = (root / arm / "run.log").read_text(errors="replace")
    steps = [float(value) for value in step_re.findall(log)]
    mains = main_re.findall(log)
    if len(steps) < 3 or not mains:
        raise SystemExit(f"{arm}: incomplete step timing/conservation output")
    results[arm] = (sum(steps[-3:]) / 3.0, mains[-1])

mg_eff = results["mg_t1"][0] / (4.0 * results["mg_t4"][0])
fftw_eff = results["fftw_t1"][0] / (4.0 * results["fftw_t4"][0])
base_speedup = results["mg_t4"][0] / results["fftw_t4"][0]
with (root / "summary.txt").open("w") as stream:
    for arm, (step, state) in results.items():
        stream.write(f"{arm}_mean_step_seconds {step:.6f}\n")
        stream.write(f"{arm}_final {' '.join(state)}\n")
    stream.write(f"mg_omp4_parallel_efficiency_percent {100.0 * mg_eff:.3f}\n")
    stream.write(f"fftw_omp4_parallel_efficiency_percent {100.0 * fftw_eff:.3f}\n")
    stream.write(f"fftw_vs_mg_omp4_speedup {base_speedup:.3f}\n")
print("DMO_OMP4_BASE_SOLVER_PASS")
PY
