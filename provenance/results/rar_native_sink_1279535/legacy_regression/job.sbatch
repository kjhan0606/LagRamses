#!/bin/bash
#SBATCH --job-name=rar-source-regression
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_source_regression_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_source_regression_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
run=$base/runs/rar_stellar_commit_${RAR_SOURCE_RUN_ID:?}
out=$base/audits/rar_source_regression_${SLURM_JOB_ID:?}
test -f "$run/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
for name in audit_rar_global_inventory.py reconcile_rar_cooling_global.py audit_rar_stellar_commit.py;do cp /home/kjhan/BackUP/rar/scripts/$name "$out/";done
sha256sum "$run/run.log" "$run/effective.nml" "$run/ramses_final3d" > "$out/inputs.sha256"
python3 "$out/audit_rar_global_inventory.py" "$run/run.log" "$out/agn_inventory.json"
python3 "$out/audit_rar_global_inventory.py" "$run/run.log" "$out/cooling_inventory.json" --before-phase before_cooling --after-phase after_cooling --activity any_energy
python3 "$out/reconcile_rar_cooling_global.py" "$run/run.log" "$out/cooling_inventory.json" "$out/cooling_result.json"
python3 "$out/audit_rar_stellar_commit.py" "$run/run.log" "$out/stellar_result.json" --require-elements
(cd "$out" && sha256sum *.json *.py job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'LEGACY_PIC_INVENTORIES_AND_SOURCE_TRANSACTIONS_REGRESSION_CHECKED' > "$out/COMPLETE"
