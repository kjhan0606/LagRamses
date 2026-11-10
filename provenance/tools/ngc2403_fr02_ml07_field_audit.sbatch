#!/bin/bash
#SBATCH --job-name=rar-fr02-ml07-audit
#SBATCH --partition=a100_pcie
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_ml07_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_ml07_audit_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

project=/home/kjhan/BackUP/rar
analysis_repo="$project/build/rar_qumond_amr_history_20261107"
analysis_revision=${RAR_ANALYSIS_REVISION:?submit with the reviewed analysis commit}
data_root=/gpfs/kjhan/RAR_QUMOND
baseline_run="$data_root/runs/ngc2403_fr02_uniform512_456162"
current_run="$data_root/runs/ngc2403_fr02_ml07_uniform512_457508"
output_root="$data_root/analysis/ngc2403_fr02_ml07_field_audit_${SLURM_JOB_ID:?}"
minimum_free_bytes=1073741824

test -d "$analysis_repo"
git -C "$analysis_repo" cat-file -e "$analysis_revision^{commit}"
test -f "$baseline_run/COMPLETE"
test -f "$current_run/COMPLETE"
test ! -e "$output_root"
test -s "$baseline_run/analysis.json"
test -s "$baseline_run/star_profile_code.dat"
test -s "$baseline_run/hi_profile_code.dat"
test -s "$current_run/ml_override.json"
test -s "$current_run/star_profile_code.dat"
test -s "$current_run/hi_profile_code.dat"
for model in newtonian qumond; do
    test -s "$baseline_run/$model/output_00001/data_00001.h5"
    test -s "$current_run/$model/output_00001/data_00001.h5"
done
free_bytes=$(df -B1 "$data_root" | awk 'END {print $4}')
test "$free_bytes" -ge "$minimum_free_bytes"

gpu_record=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | sed -n '1p')
case "$gpu_record" in
    *A100*) ;;
    *) printf 'Expected an A100 GPU, found: %s\n' "$gpu_record" >&2; exit 1;;
esac

mkdir "$output_root"
trap 'rc=$?; if [ "$rc" -ne 0 ] && [ -d "$output_root" ]; then printf "ML07_FIELD_AUDIT_FAILED rc=%s\n" "$rc" > "$output_root/outcome.txt"; fi' ERR
for name in analyze_ngc2403_fr02_uniform.py analyze_qumond_spherical_reference.py analyze_qumond_source_diagnostic.py; do
    git -C "$analysis_repo" show "$analysis_revision:provenance/tools/$name" > "$output_root/$name"
done
cp -p "${BASH_SOURCE[0]}" "$output_root/job.sbatch"
printf '%s\n' "$analysis_revision" > "$output_root/analysis_source_commit.txt"
printf '%s\n' "$baseline_run" > "$output_root/baseline_run_path.txt"
printf '%s\n' "$current_run" > "$output_root/current_run_path.txt"
printf '%s\n' "$gpu_record" > "$output_root/runtime_gpu.txt"
printf 'purpose=read-only matched M/L=0.5 and M/L=0.7 periodic force audit\nsource_run=%s\nbaseline_run=%s\nanalysis_revision=%s\nexpected_outputs=analysis.json and provenance only; no RAMSES dumps\nfree_bytes_before=%s\n' \
    "$current_run" "$baseline_run" "$analysis_revision" "$free_bytes" > "$output_root/output_policy.txt"

export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE MKL_NUM_THREADS=1
env -u LD_LIBRARY_PATH PYTHONPATH="$output_root" \
    /home/kjhan/miniconda3/bin/python3 "$output_root/analyze_ngc2403_fr02_uniform.py" \
    "$current_run" 512 --baseline-root "$baseline_run" --output "$output_root/analysis.json" \
    > "$output_root/analysis.log" 2>&1
grep -Fq '"stellar_ml": 0.7' "$output_root/analysis.json"
grep -Fq '"stellar_ml_baseline": 0.5' "$output_root/analysis.json"
grep -Fq '"gas_profile_exactly_unchanged": true' "$output_root/analysis.json"
printf '%s\n' 'ML07_FIELD_AUDIT_COMPLETE' > "$output_root/outcome.txt"
sha256sum "$output_root/analysis.json" "$output_root/analysis.log" \
    "$output_root/analyze_ngc2403_fr02_uniform.py" \
    "$output_root/analyze_qumond_spherical_reference.py" \
    "$output_root/analyze_qumond_source_diagnostic.py" \
    "$output_root/analysis_source_commit.txt" "$output_root/baseline_run_path.txt" \
    "$output_root/current_run_path.txt" "$output_root/runtime_gpu.txt" \
    "$output_root/output_policy.txt" "$output_root/job.sbatch" \
    "$output_root/outcome.txt" > "$output_root/artifacts.sha256"
(cd "$output_root" && sha256sum -c artifacts.sha256 --status)
printf 'M/L=0.7 matched field audit complete. Analysis directory: %s\n' "$output_root"
