#!/bin/bash
#SBATCH --job-name=rar-fr02-source-compare
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_compare_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_compare_%j.err
set -Eeuo pipefail
umask 077
project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
repo=$project/build/rar_qumond_amr_history_20261107
revision=${RAR_COMPARISON_REVISION:?}
run_id=${RAR_FR02_FORCE_JOB:?}
[[ "$run_id" =~ ^[0-9]+$ ]]
run=$base/runs/ngc2403_fr02_uniform512_$run_id
baseline=$base/analysis/ngc2403_compiled512_cdm_nfw_holdout_1279415
observed=$base/analysis/ngc2403_input_448195/ngc2403_observed.csv
out=$base/analysis/ngc2403_fr02_comparison_${SLURM_JOB_ID:?}
report=$project/reports/ngc2403_fr02_comparison_${SLURM_JOB_ID:?}
test -f "$run/COMPLETE" && test -f "$baseline/COMPLETE"
test ! -e "$out" && test ! -e "$report"
(cd "$run" && sha256sum -c artifacts.sha256 --status)
(cd "$baseline" && sha256sum -c result.sha256 --status)
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 1073741824
mkdir "$out"
mkdir "$out/code" "$out/results" "$out/inputs"
for name in ngc2403_fr02_source_comparison.py ngc2403_compiled512_cdm_nfw_holdout.py ngc2403_lagr_rar.py ngc2403_cpu512_radial_holdout.py; do
 git -C "$repo" show "$revision:provenance/tools/$name" > "$out/code/$name"
done
printf '%s\n' "$revision" > "$out/analysis_revision.txt"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$baseline/holdout.csv" "$out/inputs/holdout.csv"
cp "$baseline/summary.json" "$out/inputs/baseline_summary.json"
cp "$observed" "$out/inputs/observed.csv"
cp "$run/analysis.json" "$out/inputs/force_analysis.json"
cp "$run/output_policy.json" "$run/run_identity.txt" "$run/source_contract.json" "$out/inputs/"
cp "$run/artifacts.sha256" "$out/inputs/force_artifacts.sha256"
(cd "$out" && sha256sum code/*.py inputs/* > inputs.sha256)
printf '%s\n' 'purpose=CPU comparison of compiled FR02 and THINGS-like force fields at the same resolution' "force_run=$run" "free_bytes=$free_bytes" 'new_full_dumps=0; output_budget_bytes=16777216; inputs_are_immutable' > "$out/output_policy.txt"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
srun --mpi=pmi2 env -u LD_LIBRARY_PATH PYTHONPATH="$out/code" \
 /home/kjhan/miniconda3/bin/python3 "$out/code/ngc2403_fr02_source_comparison.py" \
 "$run" "$out/inputs/holdout.csv" "$out/inputs/baseline_summary.json" \
 "$out/inputs/observed.csv" "$out/results" > "$out/compute.log" 2>&1
(cd "$out" && sha256sum -c inputs.sha256 --status)
(cd "$out" && find . -type f ! -name artifacts.sha256 ! -name COMPLETE -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256 && sha256sum -c artifacts.sha256 --status)
printf '%s\n' 'MATCHED_COMPILED_GAS_SOURCE_COMPARISON_COMPLETE_NOT_COSMOLOGY' > "$out/COMPLETE"
mkdir "$report"
cp -pr "$out"/. "$report/"
(cd "$report" && sha256sum -c artifacts.sha256 --status)
