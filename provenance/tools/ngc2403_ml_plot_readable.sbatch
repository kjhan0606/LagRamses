#!/bin/bash
#SBATCH --job-name=rar-ml-plot
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ml_plot_%j.out
set -euo pipefail
project=/home/kjhan/BackUP/rar
source=$project/reports/ngc2403_ml_newton_1284238
test -f "$source/COMPLETE"
(cd "$source" && sha256sum -c SHA256SUMS --status)
out=$project/reports/ngc2403_ml_plot_${SLURM_JOB_ID:?}
srun env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$project/scripts/ngc2403_ml_plot_readable.py" "$source/comparison.csv" "$out"
cp "$project/scripts/ngc2403_ml_plot_readable.py" "$out/code.py"
(cd "$out" && sha256sum code.py comparison.pdf comparison.png > SHA256SUMS && sha256sum -c SHA256SUMS)
touch "$out/COMPLETE"
