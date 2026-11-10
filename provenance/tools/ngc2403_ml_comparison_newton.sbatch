#!/bin/bash
#SBATCH --job-name=rar-ml-newton
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ml_newton_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ml_newton_%j.err
set -euo pipefail
project=/home/kjhan/BackUP/rar
out=$project/reports/ngc2403_ml_newton_${SLURM_JOB_ID:?}
test ! -e "$out"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
srun env -u LD_LIBRARY_PATH PYTHONPATH="$project/reports/ngc2403_fr02_comparison_1279603/code" /home/kjhan/miniconda3/bin/python3 "$project/scripts/ngc2403_ml_comparison_newton.py" "$out"
cp "$project/scripts/ngc2403_ml_comparison_newton.py" "$out/code.py"
cp "$project/scripts/ngc2403_ml_comparison_newton.sbatch" "$out/job.sbatch"
(cd "$out" && sha256sum code.py job.sbatch comparison.csv comparison.pdf comparison.png summary.json > SHA256SUMS && sha256sum -c SHA256SUMS)
touch "$out/COMPLETE"
