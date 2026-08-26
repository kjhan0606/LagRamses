#!/bin/bash
#SBATCH --job-name=adm_temp_sens
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:25:00
#SBATCH --output=/gpfs/kjhan/Hydro/Sidm/validation/adm_temperature_sensitivity/slurm_%j.out
#SBATCH --error=/gpfs/kjhan/Hydro/Sidm/validation/adm_temperature_sensitivity/slurm_%j.err

set -euo pipefail

code=/home/kjhan/BACKUP/lagRamses-SIDM
root=/gpfs/kjhan/Hydro/Sidm/validation/adm_temperature_sensitivity/${SLURM_JOB_ID:?}
snapshot=${root}/snapshot
mkdir -p "${snapshot}"

git -C "${code}" status --short >"${snapshot}/git_status.txt"
git -C "${code}" diff --binary >"${snapshot}/source.diff"
for file in \
  tests/sidm/generate_amr_estimator_ic.py \
  tests/adm/adm_thermal_history.py \
  tests/adm/adm_temperature_sensitivity.nml.in \
  tests/adm/analyze_adm_temperature_sensitivity.py \
  tests/adm/run_adm_temperature_sensitivity.sbatch; do
  cp -p "${code}/${file}" "${snapshot}/$(basename "${file}")"
done

"${code}/bin/build_lageunha.sh"
cp -p "${code}/bin/build_lageunha.log" "${snapshot}/build_lageunha.log"
grep -Fq '=== make exit=0' "${snapshot}/build_lageunha.log"
sha256sum "${code}/bin/ramses_final3d" >"${snapshot}/binary.sha256"

python3 "${snapshot}/adm_thermal_history.py" \
  --xi 0.5 --z-init 100 --z-kd 100 300 1000 3000 10000 \
  --output "${root}/thermal_history.json" \
  | tee "${root}/thermal_history_summary.txt"

# initfile is a fixed-width Fortran string. Keep this IC path short enough
# for both the directory and its /ic_part suffix to survive namelist parsing.
ic_directory=/gpfs/kjhan/Hydro/Sidm/validation/admts/${SLURM_JOB_ID}/ic
if (( ${#ic_directory} + 8 > 80 )); then
  echo "IC directory exceeds the initfile width: ${ic_directory}" >&2
  exit 2
fi
mkdir -p "${ic_directory}"
python3 "${snapshot}/generate_amr_estimator_ic.py" "${ic_directory}" \
  --multiplicity 4 --coincident-copies >"${root}/ic_generation.log"

python3 - "${root}/thermal_history.json" >"${root}/candidate_cases.txt" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as stream:
    candidates = json.load(stream)["candidates"]
for candidate in candidates:
    if candidate["above_code_floor"]:
        print(candidate["label"], format(candidate["temperature_k"], ".16e"))
PY

# Keep the candidate list in a regular file. MPI launchers can inherit and
# consume a loop's process-substitution file descriptor.
while read -r label temperature; do
  case_directory=${root}/${label}
  mkdir -p "${case_directory}"
  sed -e "s|__IC_DIRECTORY__|${ic_directory}|" \
      -e "s|__T_INIT__|${temperature}|" \
    "${snapshot}/adm_temperature_sensitivity.nml.in" >"${case_directory}/adm_temperature_sensitivity.nml"
  : >"${case_directory}/jobcontrol.txt"
  cp -p "${code}/bin/ramses_final3d" "${case_directory}/ramses_final3d"
  (
    cd "${case_directory}"
    I_MPI_FABRICS=shm mpirun -np 1 \
      ./ramses_final3d adm_temperature_sensitivity.nml </dev/null >run.log 2>&1
  )
done <"${root}/candidate_cases.txt"

python3 "${snapshot}/analyze_adm_temperature_sensitivity.py" "${root}" \
  | tee "${root}/temperature_sensitivity_summary.txt"
