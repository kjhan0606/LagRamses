#!/bin/bash
#SBATCH --job-name=z100_hr_cic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=make_z100_hr_cic_%j.log
#SBATCH --error=make_z100_hr_cic_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
snapshot=${root}/zoom_run_cdm/output_00004
figures=${root}/figures/zstart_comparison
mkdir -p "${figures}"

"${root}/tools/make_zoom_density_figure.py" "${snapshot}" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --assignment cic \
    --max-particle-mass 2.0e-10 \
    --zoom-bins 512 \
    --full-sigma 2.5 \
    --zoom-sigma 1.5 \
    --output "${figures}/cdm_zstart100_z3p97_hr_cic.png" \
    --title "CDM zoom, 2LPT start at z=100" \
    --label "Finest particle tier only; CIC projection"
