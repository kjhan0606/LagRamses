#!/bin/bash
#SBATCH --job-name=validate_hhe_thermal
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/validate_hhe_thermal_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/validate_hhe_thermal_%j.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
cd /gpfs/kjhan/LRD_JWST/simulation/snrt

.venv/bin/python tools/validate_conservative_primordial_thermal.py \
  --thermal-atlas data/production_metal_thermal_atlas_v2.h5 \
  --output data/validate_conservative_primordial_thermal.json \
  --fixed-point-iterations 16 \
  --thermal-implicit-iterations 24
