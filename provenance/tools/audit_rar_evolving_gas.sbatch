#!/bin/bash
#SBATCH --job-name=rar-evolving-audit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_evolving_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_evolving_audit_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
run=$base/runs/rar_evolving_gas_${RAR_EVOLVE_RUN_ID:?}
out=$base/audits/rar_evolving_gas_${SLURM_JOB_ID:?}
test -f "$run/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp /home/kjhan/BackUP/rar/scripts/audit_rar_evolving_gas.py "$out/audit.py"
cp /home/kjhan/BackUP/rar/scripts/verify_rar_composite_publication.py "$out/verify.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$run" "$out" "$out/verify.py"
(cd "$out" && sha256sum result.json verify_*.json audit.py verify.py job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'EVOLVING_GAS_SOURCE_CLOCK_MASS_AND_KICK_BUDGET_CHECKED_NOT_PRODUCTION' > "$out/COMPLETE"
