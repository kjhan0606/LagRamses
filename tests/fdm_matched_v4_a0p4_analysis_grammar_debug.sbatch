#!/bin/bash
#SBATCH --job-name=fdm_v4_a04
#SBATCH --partition=debug
#SBATCH --nodelist=grammar-debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/home/kjhan/BACKUP/FDM/fdm_matched_v4_a0p4_%j.slurm
#SBATCH --error=/home/kjhan/BACKUP/FDM/fdm_matched_v4_a0p4_%j.err

# First preregistered matched-v4 analysis epoch.  This is post-processing only.
# RAMSES itself is never launched by this job.

set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

repo=/home/kjhan/BACKUP/FDM/lagramses_science_clean
analysis=/home/kjhan/BACKUP/FDM/matched_v4_analysis/a0p4_ngrid256
script=${repo}/patch/cuRamses/aux/analyze_fdm_matched.py

wave=/gpfs/kjhan/Hydro/FDM/run_3way_wave_matched_v4_7ba393f
hybrid=/gpfs/kjhan/Hydro/FDM/run_3way_hybrid_matched_v4_7ba393f
fluid=/gpfs/kjhan/Hydro/FDM/run_3way_fluid_matched_v4_7ba393f

mkdir -p "${analysis}"
sha256sum "${script}" "$0" > "${analysis}/analysis.sha256"
git -C "${repo}" rev-parse HEAD > "${analysis}/source.commit"
git -C "${repo}" status --short > "${analysis}/git-status.txt"

scratch_parent=${SLURM_TMPDIR:-/tmp}
python3 "${script}" \
  --arm wave "${wave}/output_00039" "${wave}/run_405304.log" \
  --arm hybrid "${hybrid}/output_00024" "${hybrid}/run_405305.log" \
  --arm fluid "${fluid}/output_00023" "${fluid}/run_405306.log" \
  --output "${analysis}" \
  --ngrid 256 \
  --nkbin 48 \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --scratch "${scratch_parent}"

echo "FDM_MATCHED_V4_A0P4_ANALYSIS_COMPLETE $(date --iso-8601=seconds)"
