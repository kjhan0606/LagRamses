#!/bin/bash
# Run after both snapshots in PROFILE_CASE have reached a=1.  Example:
#   sbatch --export=ALL,PROFILE_CASE=l17_sidm3 run_lowz_profile_analysis.sbatch
# The script refuses to use an output farther than 0.005 in a from z=0.
#SBATCH --job-name=sidm_lowz_profile
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/slurm_lowz_profile_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/slurm_lowz_profile_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
validation_root=/gpfs/kjhan/Hydro/Sidm/validation
profile_case=${PROFILE_CASE:?PROFILE_CASE must be l18_sidm1 or l17_sidm3}
analysis_root=${validation_root}/z0_profiles/${profile_case}_analysis_${SLURM_JOB_ID}

case "${profile_case}" in
  l18_sidm1)
    cdm_root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_cdm
    sidm_root=${validation_root}/z0_production_l18/sidm1
    sidm_option=--sidm1
    ;;
  l17_sidm3)
    cdm_root=${validation_root}/z0_production_l17/cdm
    sidm_root=${validation_root}/z0_production_l17/sidm3
    sidm_option=--sidm3
    ;;
  *)
    printf 'unsupported PROFILE_CASE: %s\n' "${profile_case}" >&2
    exit 2
    ;;
esac

select_z0_output() {
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
    candidates.append((abs(aexp - 1.0), aexp, snapshot))
if not candidates:
    raise SystemExit(f"no readable outputs below {run_root}")
distance, aexp, snapshot = min(candidates)
if distance > 0.005:
    raise SystemExit(f"closest output {snapshot} has aexp={aexp}, not the z=0 target")
print(snapshot)
PY
}

cdm_snapshot=$(select_z0_output "${cdm_root}")
sidm_snapshot=$(select_z0_output "${sidm_root}")
mkdir -p "${analysis_root}/snapshot"
git -C "${code}" status --short >"${analysis_root}/snapshot/git_status.txt"
git -C "${code}" diff --binary >"${analysis_root}/snapshot/source.diff"
cp -p "${code}/paper/analyse_lowz_sidm_profiles.py" "${analysis_root}/snapshot/"
cp -p "${code}/paper/analyse_zoom_sidm_profiles.py" "${analysis_root}/snapshot/"
cp -p "${code}/tests/sidm/run_lowz_profile_analysis.sbatch" "${analysis_root}/snapshot/"

python3 "${analysis_root}/snapshot/analyse_lowz_sidm_profiles.py" \
  --cdm "${cdm_snapshot}" \
  "${sidm_option}" "${sidm_snapshot}" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --data "${analysis_root}/lowz_profiles.json" \
  --figure "${analysis_root}/lowz_profiles.pdf"
