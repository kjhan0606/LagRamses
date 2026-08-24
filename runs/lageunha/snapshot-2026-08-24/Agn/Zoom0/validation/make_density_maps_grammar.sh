#!/bin/bash
#SBATCH --job-name=zoom0_density_maps
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=make_density_maps_%j.log
#SBATCH --error=make_density_maps_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
maker=${root}/tools/make_zoom_density_figure.py
figures=${root}/figures
mkdir -p "${figures}"

"${maker}" "${root}/zoom_run_cdm/output_00002" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --zoom-bins 384 \
    --full-sigma 2.5 \
    --zoom-sigma 3.0 \
    --output "${figures}/cdm_zoom_density_current_z90.png" \
    --title "Current level-11 CDM zoom" \
    --label "Latest complete zoom snapshot"

"${maker}" "${root}/parent_run/output_00008" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --shift 0.0859375 -0.125 -0.26171875 \
    --output "${figures}/cdm_zoom_target_parent_z0.png" \
    --title "Selected group environment in the parent volume" \
    --label "z=0 parent preview in the shifted zoom coordinate system"

printf 'Density figures completed at %s\n' "$(date)"
