#!/bin/bash
# Run on grammar only after the independent Syntax force analysis succeeds.
# Slurm cannot express a cross-cluster afterok dependency.
#SBATCH --job-name=rar-2403-c512-hold
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_compiled512_holdout_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_compiled512_holdout_%j.err

set -Eeuo pipefail
umask 077
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
project=/home/kjhan/BackUP/rar
analysis_repo=$project/build/rar_qumond_amr_history_20261107
analysis_revision=${RAR_HOLDOUT_REVISION:?}
base=/gpfs/kjhan/RAR_QUMOND
run_root=${RAR_QUMOND_512_RUN_ROOT:-$base/runs/ngc2403_uniform512_things_451793}
gate_root=${RAR_QUMOND_512_GATE_ROOT:-$run_root}
observed=$base/analysis/ngc2403_input_448195/ngc2403_observed.csv
observed_sha256_expected=852eb3bcbf10e0e088f70598a25075b991867eb64a0c7ca0cd46d6ff7eb6fb95
out=$base/analysis/ngc2403_compiled512_cdm_nfw_holdout_${SLURM_JOB_ID:?}
report=$project/reports/ngc2403_compiled512_cdm_nfw_holdout_${SLURM_JOB_ID:?}
runlog=$base/logs/ngc2403_compiled512_holdout_compute_${SLURM_JOB_ID:?}.log
code_dir=$base/analysis/ngc2403_compiled512_holdout_code_${SLURM_JOB_ID:?}

test -f "$gate_root/COMPLETE"
(cd "$gate_root" && sha256sum -c artifacts.sha256)
if [ "$gate_root" = "$run_root" ];then (cd "$run_root" && sha256sum -c analysis.sha256);fi
test -f "$observed"
test "$(sha256sum "$observed" | awk '{print $1}')" = "$observed_sha256_expected"
test ! -e "$out" && test ! -e "$report" && test ! -e "$code_dir"
free_bytes=$(df -B1 /gpfs | awk 'END {print $4}')
test "$free_bytes" -ge 1073741824
mkdir "$out"
mkdir "$code_dir"
for name in ngc2403_compiled512_cdm_nfw_holdout.py ngc2403_lagr_rar.py ngc2403_cpu512_radial_holdout.py;do
 git -C "$analysis_repo" show "$analysis_revision:provenance/tools/$name" > "$code_dir/$name"
done
printf '%s\n' "$analysis_revision" > "$code_dir/analysis_source_commit.txt"
(cd "$code_dir" && sha256sum *.py >source.sha256)
cp -p "$observed" "$out/ngc2403_observed.csv"
env -u LD_LIBRARY_PATH PYTHONPATH="$code_dir" \
  /home/kjhan/miniconda3/bin/python3 \
  "$code_dir/ngc2403_compiled512_cdm_nfw_holdout.py" "$out" "$run_root" \
  "$out/ngc2403_observed.csv" "$gate_root" \
  >"$runlog" 2>&1
printf '%s\n' "$free_bytes" > "$out/free_bytes_at_start.txt"
cp -p "$runlog" "$out/compute.log"
cp -p "$code_dir/ngc2403_compiled512_cdm_nfw_holdout.py" "$out/compute.py"
cp -p "$code_dir/ngc2403_lagr_rar.py" "$out/lagr_rar_provenance.py"
cp -p "$code_dir/ngc2403_cpu512_radial_holdout.py" "$out/nfw_holdout_provenance.py"
cp -p "$code_dir/analysis_source_commit.txt" "$out/analysis_source_commit.txt"
cp -p "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp -p "$gate_root/analysis.json" "$run_root/output_policy.json" "$run_root/source_contract.json" \
  "$run_root/run_identity.txt" "$out/"
cp -p "$gate_root/artifacts.sha256" "$out/force_analysis_artifacts.sha256"
(cd "$out" && sha256sum analysis.json >analysis.sha256)
printf '%s\n' "$gate_root" > "$out/force_analysis_root.txt"
printf '%s\n' "$SLURM_JOB_ID" > "$out/slurm_job_id.txt"
printf '%s\n' "$free_bytes" > "$out/free_bytes_at_start.txt"
(cd "$out" && sha256sum holdout.csv summary.json rar_comparison.pdf \
  rar_comparison.png compute.log compute.py lagr_rar_provenance.py \
  nfw_holdout_provenance.py job.sbatch analysis.json analysis.sha256 \
  output_policy.json source_contract.json run_identity.txt force_analysis_artifacts.sha256 force_analysis_root.txt \
  ngc2403_observed.csv slurm_job_id.txt free_bytes_at_start.txt analysis_source_commit.txt > result.sha256)
(cd "$code_dir" && sha256sum -c source.sha256)
(cd "$out" && sha256sum -c result.sha256)
mkdir "$report"
cp -p "$out"/{holdout.csv,summary.json,rar_comparison.pdf,rar_comparison.png,compute.log,compute.py,lagr_rar_provenance.py,nfw_holdout_provenance.py,job.sbatch,analysis.json,analysis.sha256,output_policy.json,source_contract.json,run_identity.txt,force_analysis_artifacts.sha256,force_analysis_root.txt,ngc2403_observed.csv,slurm_job_id.txt,free_bytes_at_start.txt,analysis_source_commit.txt,result.sha256} "$report/"
(cd "$report" && sha256sum -c result.sha256)
printf '%s\n' 'COMPILED_FIXED_SOURCE_NFW_HOLDOUT_COMPLETE_NOT_COSMOLOGICAL_LCDM' > "$out/COMPLETE"
cp "$out/COMPLETE" "$report/COMPLETE"
