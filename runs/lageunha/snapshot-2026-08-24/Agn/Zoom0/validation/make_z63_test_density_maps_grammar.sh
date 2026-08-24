#!/bin/bash
#SBATCH --job-name=z63_test_maps
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=make_z63_test_maps_%j.log
#SBATCH --error=make_z63_test_maps_%j.err

set -euo pipefail

root=/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0
run=${root}/zoom_run_cdm_z63_test
maker=${root}/tools/make_zoom_density_figure.py
figures=${root}/figures/zstart_comparison
mkdir -p "${figures}"

for number in 00001 00002 00003; do
    snapshot=${run}/output_${number}
    if [ ! -s "${snapshot}/info_${number}.txt" ]; then
        printf 'Missing complete snapshot %s\n' "${snapshot}" >&2
        exit 2
    fi
    "${maker}" "${snapshot}" \
        --workers "${SLURM_CPUS_PER_TASK}" \
        --assignment cic \
        --max-particle-mass 2.0e-10 \
        --zoom-bins 512 \
        --full-sigma 2.5 \
        --zoom-sigma 1.5 \
        --output "${figures}/cdm_zstart63_${number}_hr_cic.png" \
        --title "CDM zoom, 2LPT start at z=63" \
        --label "Finest particle tier only; CIC projection"
done

printf 'zstart=63 CIC figures completed at %s\n' "$(date --iso-8601=seconds)"
