#!/bin/bash
#SBATCH --job-name=cic_omp1
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --exclusive
#SBATCH --mem=80G
#SBATCH --time=00:10:00
#SBATCH --output=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_omp1_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_omp1_%j.err
set -euo pipefail

paper=/home/kjhan/BACKUP/lagRamses-DE/papers/paper-Ib_jwst_early_structure
restart=${paper}/output_00002
exe=${COLORED_EXE:?COLORED_EXE is required}
root=/home/kjhan/BACKUP/lagRamses-DE/dmo_cic_omp1_${SLURM_JOB_ID}

mkdir -p "${root}"
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close
unset OMP_WAIT_POLICY KMP_BLOCKTIME
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

for nmpi in 4 8 16 32; do
  run=${root}/r${nmpi}t1
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
  srun --exclusive --kill-on-bad-exit=1 --mpi=pmi2 -n "${nmpi}" -c 1 \
    --cpu-bind=cores "${exe}" cosmo.nml > run.log
  grep -q 'Run completed' run.log
  grep -q 'FFTW3 CPU direct Poisson solver enabled' run.log
done

python3 - "${root}" <<'PY'
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
rho_re = re.compile(r"^\s*([0-9.]+)\s+[0-9.]+\s+[0-9.]+.*\brho\s*$", re.M)
total_re = re.compile(r"^\s*([0-9.]+)\s+100\.0\s+TOTAL\s*$", re.M)
state_re = re.compile(
    r"Mtot=\s*([^ ]+)\s+Ekin_tot=\s*([^ ]+)\s+Epop_tot=\s*([^ ]+)\s+Eint_tot=\s*([^ ]+)"
)
states = set()
with (root / "summary.txt").open("w") as stream:
    for nmpi in (4, 8, 16, 32):
        text = (root / f"r{nmpi}t1" / "run.log").read_text(errors="replace")
        rho = [float(value) for value in rho_re.findall(text)]
        total = [float(value) for value in total_re.findall(text)]
        state = state_re.findall(text)
        if not rho or not total or not state:
            raise SystemExit(f"r{nmpi}t1: incomplete timing/state evidence")
        states.add(state[-1])
        stream.write(f"r{nmpi}t1_timer_total_s {total[-1]:.6f}\n")
        stream.write(f"r{nmpi}t1_rho_s {rho[-1]:.6f}\n")
    if len(states) != 1:
        raise SystemExit(f"printed final conservation state differs: {states}")
    stream.write("final_conservation_state " + " ".join(next(iter(states))) + "\n")
print("DMO_CIC_OMP1_MPI_PASS")
PY
