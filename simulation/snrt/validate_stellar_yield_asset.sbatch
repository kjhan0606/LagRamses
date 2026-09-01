#!/bin/bash
# Preflight only: never submit a production run with an unapproved yield asset.
# Set PHASE0_YIELD_TABLE to the exact candidate table before invoking sbatch.
# The audit itself is read-only; the JSON report is the only generated output.
#SBATCH --job-name=yield_asset_gate
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/yield_asset_gate_%j.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/yield_asset_gate_%j.err

set -euo pipefail

readonly project_root="/gpfs/kjhan/LRD_JWST"
readonly table_path="${PHASE0_YIELD_TABLE:-}"
readonly report_path="${YIELD_AUDIT_JSON:-${project_root}/simulation/snrt/data/stellar_yield_asset_audit.json}"
if [[ -z "${table_path}" || ! -r "${table_path}" ]]; then
  echo "set PHASE0_YIELD_TABLE to a readable canonical Phase-0 candidate" >&2
  exit 2
fi

cd "${project_root}/simulation/snrt"
.venv/bin/python tools/audit_stellar_yield_asset.py \
  "${table_path}" \
  --json-out "${report_path}"
