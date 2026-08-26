#!/bin/bash
#SBATCH --job-name=zoom0_validate_z0
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=validate_cdm_z0_%j.log
#SBATCH --error=validate_cdm_z0_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
run=${root}/zoom_run_cdm
validation=${root}/validation
tools=${root}/tools
work=${validation}/cdm_z0_catalog
sidm_job=386148
cdm_job=386315
cd "${validation}"

snapshot=$(find "${run}" -maxdepth 1 -type d -name 'output_*' -printf '%f\n' |
           sort -V | tail -1)
if [ -z "${snapshot}" ]; then
    printf 'No CDM output is available\n' >&2
    exit 2
fi
number=${snapshot#output_}
aexp=$(awk '$1=="aexp" {print $3; exit}' \
       "${run}/${snapshot}/info_${number}.txt")
awk -v value="${aexp}" 'BEGIN {exit !(value > 0.99)}'

if [ -e "${work}" ]; then
    printf 'Validation work directory already exists: %s\n' "${work}" >&2
    exit 3
fi
mkdir "${work}"

"${tools}/hop" -in "${run}/${snapshot}/part_${number}.out" -p 1.0 \
    -o "${work}/hop"
"${tools}/regroup" -root "${work}/hop" -douter 80.0 -dsaddle 200.0 \
    -dpeak 240.0 -f77 -o "${work}/grp"
"${tools}/poshalo" -inp "${run}/${snapshot}" -pre "${work}/grp" \
    -cut 4.656612873077393e-10 > "${work}/poshalo.stdout"

selection=$("${tools}/select_zoom_halo.py" "${work}/grp.pos" \
    --expected-center 0.4992375 0.4714 0.49128125 \
    --target-mass 7.9310408e13 \
    --json "${work}/matched_halo.json")
read -r group xc yc zc r200_box two_r200_box mass_hinv group_contam \
    <<< "${selection}"

"${tools}/zoom_particle_qc.py" "${run}/${snapshot}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --center "${xc}" "${yc}" "${zc}" \
    --radii-box "${r200_box}" "${two_r200_box}" \
    --json "${work}/particle_contamination.json" \
    > "${work}/particle_contamination.txt"

python3 - "${work}/matched_halo.json" \
          "${work}/particle_contamination.json" <<'PY'
import json
import sys

halo = json.load(open(sys.argv[1]))
particles = json.load(open(sys.argv[2]))
if halo["periodic_distance_from_expected"] > 0.02:
    raise SystemExit("Matched halo moved implausibly far from parent centre")
if not 0.5 <= halo["mass_ratio_to_target"] <= 2.0:
    raise SystemExit("Matched halo mass differs too much from parent target")
if halo["contamination_group"] > 1.0e-8:
    raise SystemExit("HOP group contains low-resolution mass")
if len(particles["spheres"]) != 2:
    raise SystemExit("Missing R200 contamination measurements")
for sphere in particles["spheres"]:
    fraction = sphere["coarse_mass_fraction"]
    if fraction is None or fraction > 1.0e-8:
        raise SystemExit(
            f"Low-resolution contamination at radius {sphere['radius_box']}"
        )
PY

if [ -s "${run}/zoom_cdm_restart3_${cdm_job}.err" ]; then
    printf 'CDM stderr is not empty\n' >&2
    exit 4
fi
if grep -Eiq 'segmentation fault|sigsegv|floating invalid|fatal error|nan detected|mpi_abort|out of memory' \
       "${run}/zoom_cdm_restart3_${cdm_job}.log"; then
    printf 'Fatal pattern found in CDM log\n' >&2
    exit 5
fi

touch CDM_VALIDATION_PASS
scontrol release "${sidm_job}"
printf 'CDM validation passed; released SIDM job %s at %s\n' \
       "${sidm_job}" "$(date)"
