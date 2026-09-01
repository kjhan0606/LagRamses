#!/bin/bash
#SBATCH --job-name=validate_h_stromgren
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/validate_h_stromgren_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/validate_h_stromgren_%j.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
cd /gpfs/kjhan/LRD_JWST/simulation/snrt

.venv/bin/python tools/validate_conservative_hydrogen_stromgren.py \
  --output data/validate_conservative_hydrogen_stromgren.json \
  --fixed-point-iterations 12
