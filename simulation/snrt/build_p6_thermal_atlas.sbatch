#!/bin/bash
#SBATCH --job-name=snrt_metal_atlas
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/build_p6_thermal_atlas_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/build_p6_thermal_atlas_%j.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
cd /gpfs/kjhan/LRD_JWST/simulation/snrt

.venv/bin/python tools/build_metal_thermal_atlas.py \
  --source-data /gpfs/kjhan/LRD_JWST/external/grackle/CloudyData_noUVB.h5 \
  --expected-source-sha256 0abe25cceeb5c0825381c5f17059982a9a2cdd27ce369a475c559fba6a8fa106 \
  --scale-factors config/p6_thermal_atlas_scale_factors.txt \
  --output data/production_metal_thermal_atlas_v2.h5
