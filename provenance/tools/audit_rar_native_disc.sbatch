#!/bin/bash
#SBATCH --job-name=rar-disc-audit
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_disc_audit_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_disc_audit_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
low=$base/runs/rar_native_disc_${RAR_DISC_LOW_ID:?}
high=$base/runs/rar_native_disc_${RAR_DISC_HIGH_ID:?}
out=$base/audits/rar_native_disc_${SLURM_JOB_ID:?}
test -f "$low/COMPLETE";test -f "$high/COMPLETE"
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$project/scripts/verify_rar_composite_publication.py" "$out/verify.py"
cp "$project/scripts/analyze_rar_native_disc.py" "$out/analyze.py"
cp "$project/scripts/plot_rar_native_disc.py" "$out/plot_rar_native_disc.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python=/home/kjhan/miniconda3/bin/python3
for resolution in low high; do
 if [ "$resolution" = low ];then run=$low;else run=$high;fi
 for mode in newtonian_baryon newtonian_halo qumond_baryon;do
  "$python" "$out/verify.py" "$run" "$out/verify_${resolution}_${mode}.json" --case "$mode" --snapshots 1
 done
done
"$python" "$out/analyze.py" "$low" "$high" "$out"
(cd "$out" && sha256sum result.json rotation_rar.csv rotation_rar.pdf rotation_rar.png verify_*.json verify.py analyze.py plot_rar_native_disc.py job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'MATCHED_NATIVE_DISC_FORCE_AND_RAR_DIAGNOSTICS_CHECKED_NOT_GALAXY_FORMATION' > "$out/COMPLETE"
