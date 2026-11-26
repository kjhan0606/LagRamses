#!/bin/bash
#SBATCH --job-name=rar-base-spectra-plot
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/base_spectra_plot_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/base_spectra_plot_%j.err
set -Eeuo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
source=/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_base_camb_spectra_1578580
out=/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_base_spectra_plot_${SLURM_JOB_ID:?}
revision=${1:?Pass a full committed revision}
test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test -f "$source/COMPLETE"
(cd "$source" && sha256sum -c artifacts.sha256 --status)
test ! -e "$out"
mkdir "$out"
trap 'touch "$out/FAILED"' ERR
git -C "$repo" show "$revision:provenance/tools/plot_lcdm_base_spectra.py" > "$out/plot.py"
cp -p "$0" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
export RAR_PLOT_OUTPUT="$out" MPLCONFIGDIR="$out/matplotlib"
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 "$out/plot.py"
(cd "$out" && sha256sum plot.py job.sbatch source_revision.txt base_spectra.pdf base_spectra.png > artifacts.sha256)
(cd "$out" && sha256sum -c artifacts.sha256 --status)
touch "$out/COMPLETE"
