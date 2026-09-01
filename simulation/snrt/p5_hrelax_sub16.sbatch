#!/bin/bash
#SBATCH --job-name=p5_hrelax_sub16
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/p5_hrelax_sub16_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/p5_hrelax_sub16_%j.err

set -euo pipefail

readonly input_manifest="/gpfs/kjhan/LRD_JWST/manifests/lrd_jwst_external_assets.json"
[[ -s "$input_manifest" ]] || { echo "missing input manifest: $input_manifest" >&2; exit 1; }
cd /gpfs/kjhan/LRD_JWST/simulation/snrt

.venv/bin/python tools/p5_run_thermochemical_pilot.py \
  --input data/p4_coeval_static_rt_input.h5 \
  --photon-metadata data/p4_pilot_agn_photon_ledger.json \
  --thermal-atlas data/production_metal_thermal_atlas_v2.h5 \
  --scale-factor 0.208497764676753 \
  --output data/p5_coeval_hrelax_s8_sub16_6p37myr.h5 \
  --duration-myr 6.3697581 \
  --sn-order 8 \
  --reduced-light-fraction 0.01 \
  --courant 0.4 \
  --thermal-subcycles 16 \
  --time-averaged-absorption-iterations 4
