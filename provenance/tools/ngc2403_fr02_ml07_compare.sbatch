#!/bin/bash
#SBATCH --job-name=rar-ml07-compare
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ml07_compare_%j.out
set -euo pipefail
project=/home/kjhan/BackUP/rar
run=/gpfs/kjhan/RAR_QUMOND/runs/ngc2403_fr02_ml07_uniform512_457508
test -f "$run/COMPLETE"
cd "$run"
sha256sum -c artifacts.sha256 --status
out=$project/reports/ngc2403_fr02_ml07_comparison_${SLURM_JOB_ID:?}
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
srun env -u LD_LIBRARY_PATH PYTHONPATH="$project/reports/ngc2403_fr02_comparison_1279603/code" /home/kjhan/miniconda3/bin/python3 "$project/scripts/ngc2403_fr02_ml07_compare.py" "$run" "$project/reports/ngc2403_fr02_comparison_1279603/results/comparison.csv" "$out"
