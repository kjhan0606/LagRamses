#!/bin/bash
#SBATCH --job-name=rar-ngc512-analysis
#SBATCH --partition=a100_pcie
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_512_analysis_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_512_analysis_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
analysis_repo=$project/build/rar_qumond_amr_history_20261107
analysis_revision=${RAR_ANALYSIS_REVISION:?}
run=$base/runs/ngc2403_uniform512_things_451793
out=$base/analysis/ngc2403_512_recovery_${SLURM_JOB_ID:?}
test -f "$run/newtonian/output_00001/data_00001.h5"
test -f "$run/qumond/output_00001/data_00001.h5"
test ! -e "$out"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
git -C "$analysis_repo" show "$analysis_revision:provenance/tools/analyze_ngc2403_lagr_512.py" > "$out/analyze.py"
for name in analyze_qumond_spherical_reference.py analyze_qumond_source_diagnostic.py;do
 git -C "$analysis_repo" show "$analysis_revision:provenance/tools/$name" > "$out/$name"
done
printf '%s\n' "$analysis_revision" > "$out/analysis_source_commit.txt"
cp "$run/run_identity.txt" "$out/simulation_identity.txt"
printf '%s\n' 'purpose=GPU postprocessing of existing two 512^3 force snapshots; no simulation launch, no new snapshots' "read_only_simulation_root=$run" 'gpu_cards=1; host_memory=128GiB; new_output_budget=16MiB; original failed wrapper and outputs preserved' > "$out/execution_policy.txt"
for model in newtonian qumond;do
 grep -Fq 'Run completed' "$run/$model/run.log"
 if grep -Eiq 'ERROR:|failed to converge|NaN/Inf detected|Too many errors|Aborting' "$run/$model/run.log";then exit 1;fi
 (cd "$run/$model" && sha256sum -c inputs.sha256)
done
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
env -u LD_LIBRARY_PATH PYTHONPATH="$out" /home/kjhan/miniconda3/bin/python3 "$out/analyze.py" "$run" 512 --output "$out/analysis.json" > "$out/analysis.log" 2>&1
(cd "$out" && sha256sum analysis.json analysis.log analyze.py analyze_qumond_spherical_reference.py analyze_qumond_source_diagnostic.py job.sbatch simulation_identity.txt analysis_source_commit.txt >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'EXISTING_512_SNAPSHOTS_REANALYZED_NO_SIMULATION_RERUN' > "$out/COMPLETE"
