#!/bin/bash
#SBATCH --job-name=p4_hdf5_stage
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/p4_hdf5_stage_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/p4_hdf5_stage_%j.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
readonly snapshot_path="${SNRT_P4_SNAPSHOT:-}"
if [[ -z "$snapshot_path" || ! -r "$snapshot_path" ]]; then
  echo "set SNRT_P4_SNAPSHOT to the external output_00017/data_00017.h5 before submission" >&2
  exit 2
fi
cd /gpfs/kjhan/LRD_JWST/simulation/snrt

.venv/bin/python tools/p4_stage_hdf5_level15.py \
  --snapshot "$snapshot_path" \
  --info ../pilot_rt_ledger/output_00017/info_00017.txt \
  --zoom-manifest data/p4_high_density_manifest.json \
  --field-map config/p4_hdf5_field_map_pilot.json \
  --thermal-atlas data/production_metal_thermal_atlas_v2.h5 \
  --source-ledger data/p4_pilot_agn_photon_ledger.csv \
  --output data/p4_coeval_static_rt_input.h5 \
  --metadata-output data/p4_coeval_static_rt_input.json \
  --analysis-level 15 \
  --scan-chunk-grids 1000000
