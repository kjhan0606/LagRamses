#!/bin/bash
# Read-only independent code audit using Claude Code's current Opus alias.
#SBATCH --job-name=smbh_capture_opus_audit
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=00:12:00
#SBATCH --output=/gpfs/kjhan/lagramses_capture_opus_audit_%j.log

set -euo pipefail
module purge
module load nodejs
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
code=/home/kjhan/BACKUP/lagRamses-capture-ledger
run=/gpfs/kjhan/lagramses_capture_opus_audit_${SLURM_JOB_ID:?}
if [ -e "$run" ]; then
    echo "FATAL: refusing existing audit directory $run"
    exit 2
fi
mkdir "$run"
cp "$code/tests/sink/capture-ledger/OPUS_AUDIT_REQUEST_2026-10-03.md" "$run/request.md"
echo 'purpose=read-only independent Opus capture restart code audit'
echo "request=$run/request.md"
echo "source_commit=$(git -C "$code" rev-parse HEAD)"
sha256sum "$run/request.md" "$code/patch/lagRamses/sink_particle.kjhan.f90" \
    "$code/patch/lagRamses/init_sink.f90" "$code/patch/lagRamses/output_sink.f90" \
    "$code/tests/sink/capture-ledger/compare_evolved_hdf5.py"
cd "$code"
audit_prompt=$(<"$run/request.md")
timeout 8m /home/kjhan/.npm-global/bin/claude -p --model opus \
    --effort high --permission-mode plan --tools Read,Grep,Glob \
    --no-session-persistence --max-budget-usd 5 --output-format json \
    "$audit_prompt" >"$run/audit.json"
echo "audit_output=$run/audit.json"
echo 'CAPTURE_OPUS_AUDIT_FINISHED'
