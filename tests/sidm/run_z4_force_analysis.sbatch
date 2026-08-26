#!/bin/bash
#SBATCH --job-name=sidm_z4_l17_analysis
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17/analysis_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17/analysis_%j.err

set -euo pipefail

source_root=/home/kjhan/BACKUP/lagRamses-SIDM
validation_root=/gpfs/kjhan/Hydro/Sidm/validation/z4_force_l17
force_job=${FORCE_JOB_ID:-}
cdm_force_job=${CDM_FORCE_JOB_ID:-${force_job}}
sidm3_force_job=${SIDM3_FORCE_JOB_ID:-${force_job}}
if [[ -z "${cdm_force_job}" || -z "${sidm3_force_job}" ]]; then
  printf '%s\n' 'set FORCE_JOB_ID or both CDM_FORCE_JOB_ID and SIDM3_FORCE_JOB_ID' >&2
  exit 2
fi
analysis_root=${validation_root}/analysis_${SLURM_JOB_ID}
mkdir -p "${analysis_root}/snapshot"

git -C "${source_root}" status --short >"${analysis_root}/snapshot/git_status.txt"
git -C "${source_root}" diff --binary >"${analysis_root}/snapshot/source.diff"
cp -p "${source_root}/paper/analyse_zoom_sidm_profiles.py" \
  "${analysis_root}/snapshot/"
cp -p "${source_root}/paper/analyse_zoom_sidm_force_convergence.py" \
  "${analysis_root}/snapshot/"
cp -p "${source_root}/tests/sidm/analyze_mg_convergence_log.py" \
  "${analysis_root}/snapshot/"
cp -p "${source_root}/tests/sidm/check_z4_force_acceptance.py" \
  "${analysis_root}/snapshot/"
cp -p "${source_root}/tests/sidm/run_z4_force_analysis.sbatch" \
  "${analysis_root}/snapshot/"

select_target_output() {
  python3 - "$1" <<'PY'
import re
import sys
from pathlib import Path

run_root = Path(sys.argv[1])
candidates = []
for snapshot in sorted(run_root.glob("output_*")):
    number = snapshot.name.removeprefix("output_")
    info = snapshot / f"info_{number}.txt"
    if not info.is_file():
        continue
    match = re.search(r"^aexp\s*=\s*([^\s]+)", info.read_text(), re.MULTILINE)
    if match is None:
        continue
    aexp = float(match.group(1).replace("D", "E"))
    candidates.append((abs(aexp - 0.2), aexp, snapshot))
if not candidates:
    raise SystemExit(f"no readable outputs below {run_root}")
distance, aexp, snapshot = min(candidates)
if distance > 0.005:
    raise SystemExit(f"closest output {snapshot} has aexp={aexp}, not the a=0.2 target")
print(snapshot)
PY
}

cdm_snapshot=$(select_target_output "${validation_root}/cdm_${cdm_force_job}")
sidm3_snapshot=$(select_target_output "${validation_root}/sidm3_${sidm3_force_job}")

python3 "${analysis_root}/snapshot/analyze_mg_convergence_log.py" \
  "${validation_root}/cdm_${cdm_force_job}/run.log" \
  --json "${analysis_root}/mg_convergence_cdm.json"
python3 "${analysis_root}/snapshot/analyze_mg_convergence_log.py" \
  "${validation_root}/sidm3_${sidm3_force_job}/run.log" \
  --json "${analysis_root}/mg_convergence_sidm3.json"

python3 "${analysis_root}/snapshot/analyse_zoom_sidm_force_convergence.py" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --l17-cdm "${cdm_snapshot}" \
  --l17-sidm3 "${sidm3_snapshot}" \
  --output "${analysis_root}/force_convergence.json"

python3 "${analysis_root}/snapshot/check_z4_force_acceptance.py" \
  "${analysis_root}/force_convergence.json"
