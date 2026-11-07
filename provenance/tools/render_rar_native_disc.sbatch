#!/bin/bash
#SBATCH --job-name=rar-disc-render
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_disc_render_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_disc_render_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
parent=$base/audits/rar_native_disc_1275827
out=$base/audits/rar_native_disc_render_${SLURM_JOB_ID:?}
test -f "$parent/COMPLETE"
(cd "$parent" && sha256sum -c artifacts.sha256)
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$parent/rotation_rar.csv" "$parent/result.json" "$out/"
cp /home/kjhan/BackUP/rar/scripts/plot_rar_native_disc.py "$out/plot.py"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/kjhan/miniconda3/bin/python3 "$out/plot.py" "$out"
(cd "$out" && sha256sum result.json rotation_rar.csv rotation_rar.pdf rotation_rar.png plot.py job.sbatch >artifacts.sha256 && sha256sum -c artifacts.sha256)
printf '%s\n' 'FIGURE_LAYOUT_ONLY_PARENT_AUDIT_1275827_UNCHANGED' > "$out/COMPLETE"
