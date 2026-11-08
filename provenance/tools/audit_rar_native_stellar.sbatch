#!/bin/bash
#SBATCH --job-name=rar-native-source-audit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_stellar_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_stellar_audit_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
run=$base/runs/rar_native_stellar_${RAR_NATIVE_STELLAR_RUN_ID:?}
out=$base/audits/rar_native_stellar_${SLURM_JOB_ID:?}
project=/home/kjhan/BackUP/rar
test -f "$run/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$project/scripts/audit_rar_native_stellar.py" "$out/audit.py"
cp "$project/scripts/verify_rar_composite_publication.py" "$out/verify.py"
cp "$project/scripts/audit_rar_stellar_commit.py" "$out/stellar.py"
for name in rar_native_particles.py rar_evolving_thermal.py;do cp "$project/scripts/$name" "$out/$name";done
(cd "$out" && sha256sum *.py job.sbatch >audit_inputs.sha256)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$run" "$out"
(cd "$out" && sha256sum *.json >results.sha256 && sha256sum -c audit_inputs.sha256 && sha256sum -c results.sha256)
printf '%s\n' 'NATIVE_STELLAR_MASS_MOMENTUM_ELEMENTS_COOLING_KDK_AND_FIELDS_CHECKED_NOT_PRODUCTION' > "$out/COMPLETE"
