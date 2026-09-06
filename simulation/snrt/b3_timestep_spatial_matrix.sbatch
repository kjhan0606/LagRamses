#!/bin/bash
#SBATCH --job-name=snrt_b3_matrix
#SBATCH --partition=a10
#SBATCH --array=0-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/kjhan/LRD_JWST/simulation/snrt/data/b3_validation/b3_matrix_%A_%a.log
#SBATCH --error=/gpfs/kjhan/LRD_JWST/simulation/snrt/data/b3_validation/b3_matrix_%A_%a.err

set -euo pipefail

readonly project_root="/gpfs/kjhan/LRD_JWST"
readonly snrt_root="$project_root/simulation/snrt"
readonly expected_head="7e6dab63d87707dc4ee1749f242d3a809191cc00"
readonly coarse_input="$snrt_root/data/p4_coeval_static_rt_input_agn9.h5"
readonly refined_input="$snrt_root/data/b3_validation/p4_coeval_static_rt_input_agn9_refined2.h5"
readonly photon_metadata="$snrt_root/data/p4_pilot_agn_photon_ledger.json"
readonly thermal_atlas="$snrt_root/data/production_metal_thermal_atlas_v2.h5"
readonly output_dir="$snrt_root/data/b3_validation"

[[ "$(git -C "$project_root" rev-parse HEAD)" == "$expected_head" ]] || {
  echo "repository HEAD changed after B3 contract declaration" >&2
  exit 1
}
git -C "$project_root" diff --quiet HEAD -- \
  simulation/snrt/snrt_core \
  simulation/snrt/tools/p5_run_thermochemical_pilot.py \
  simulation/snrt/tools/refine_static_rt_input.py \
  simulation/snrt/data/p4_pilot_agn_photon_ledger.json || {
  echo "B3 runtime dependency has an uncommitted change" >&2
  exit 1
}

echo "c25720164947df6a7f64b01d898401cf00e1420e74281f233323d92d2316a8b2  $coarse_input" | sha256sum --check --status
echo "b786221c4ed9cbb4f0316542a18e5e69ff0b956d5a4b6e91c90097b0ca3a5f16  $refined_input" | sha256sum --check --status
echo "506bb47abd2481d7caab1f5ac3fa0451af10e0ec5907ec61ce03738639e9a53a  $photon_metadata" | sha256sum --check --status
echo "b1290d930b22ed049d6d3c5ed47ce56ecf3e0d2e693b39792740768e80fdf6ac  $thermal_atlas" | sha256sum --check --status

readonly case_names=(
  "n32_c0p1"
  "n32_c0p05"
  "n64_c0p1"
  "n64_c0p05"
)
readonly inputs=(
  "$coarse_input"
  "$coarse_input"
  "$refined_input"
  "$refined_input"
)
readonly courants=("0.1" "0.05" "0.1" "0.05")
readonly task_index="${SLURM_ARRAY_TASK_ID:?missing Slurm array index}"
readonly case_name="${case_names[$task_index]}"
readonly input_path="${inputs[$task_index]}"
readonly courant="${courants[$task_index]}"
readonly output_path="$output_dir/b3_matrix_${case_name}_0p1myr.h5"

[[ ! -e "$output_path" ]] || {
  echo "refusing to overwrite B3 output: $output_path" >&2
  exit 1
}

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PYTHONHASHSEED=0
cd "$snrt_root"

echo "B3_MATRIX_START case=$case_name input=$input_path courant=$courant cpus=$SLURM_CPUS_PER_TASK head=$expected_head"
.venv/bin/python tools/p5_run_thermochemical_pilot.py \
  --input "$input_path" \
  --photon-metadata "$photon_metadata" \
  --thermal-atlas "$thermal_atlas" \
  --scale-factor 0.208497764676753 \
  --output "$output_path" \
  --duration-myr 0.1 \
  --sn-order 4 \
  --reduced-light-fraction 0.01 \
  --courant "$courant" \
  --thermal-subcycles 1 \
  --source-cell-photons-per-neutral 0.25 \
  --source-deposition-mode point \
  --thermal-implicit-iterations 24 \
  --time-averaged-absorption-iterations 32 \
  --secondary-ionization fs2010 \
  --precision float64
echo "B3_MATRIX_DONE case=$case_name sha256=$(sha256sum "$output_path" | cut -d' ' -f1)"
