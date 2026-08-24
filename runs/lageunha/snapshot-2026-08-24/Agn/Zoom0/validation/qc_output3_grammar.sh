#!/bin/bash
#SBATCH --job-name=zoom0_qc_z9
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=qc_output3_%j.log
#SBATCH --error=qc_output3_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
snapshot=${root}/zoom_run_cdm/output_00003
validation=${root}/validation
scanner=${root}/tools/zoom_particle_qc.py
cd "${validation}"

aexp=$(awk '$1=="aexp" {print $3; exit}' "${snapshot}/info_00003.txt")
awk -v value="${aexp}" 'BEGIN {exit !(value > 0.09 && value < 0.11)}'

"${scanner}" "${snapshot}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --center 0.4992375 0.4714 0.49128125 \
    --radii-box 0.0054586735 0.0109173470 \
    --json output_00003_particle_qc.json \
    > output_00003_particle_qc.txt

python3 - output_00003_particle_qc.json <<'PY'
import json
import sys

report = json.load(open(sys.argv[1]))
if report["files_scanned"] != 32:
    raise SystemExit("Expected 32 RAMSES particle files")
if report["highres_count"] <= 0:
    raise SystemExit("No level-11 particles found")
masses = [row["mass_code"] for row in report["mass_tiers"]]
expected = 1.0 / 2048.0**3
if not any(abs(mass / expected - 1.0) < 1.0e-7 for mass in masses):
    raise SystemExit("Expected level-11 particle mass is absent")
PY

touch OUTPUT3_QC_PASS
printf 'output_00003 QC passed at %s\n' "$(date)"
