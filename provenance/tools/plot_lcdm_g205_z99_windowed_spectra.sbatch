#!/bin/bash
#SBATCH --job-name=rar-g205-spectra-plot
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_spectra_plot_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_spectra_plot_%j.err
set -Eeuo pipefail
umask 077
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
source=/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_g205_z99_windowed_spectra_1578622
out=/gpfs/kjhan/RAR_QUMOND/analysis/lcdm_g205_z99_spectra_plot_${SLURM_JOB_ID:?}
revision=${1:?Pass a full committed revision}
test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test -f "$source/COMPLETE"
(cd "$source" && sha256sum -c artifacts.sha256 --status)
test ! -e "$out"
mkdir "$out"
trap 'touch "$out/FAILED"' ERR
git -C "$repo" show "$revision:provenance/tools/plot_lcdm_g205_z99_windowed_spectra.py" > "$out/plot.py"
cp -p "$0" "$out/job.sbatch"
printf '%s\n' "$revision" > "$out/source_revision.txt"
export RAR_PLOT_OUTPUT="$out" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for format in pdf png; do
    env -u LD_LIBRARY_PATH RAR_PLOT_FORMAT="$format" MPLCONFIGDIR="$out/matplotlib-$format" \
        /home/kjhan/miniconda3/bin/python3 "$out/plot.py" > "$out/$format.log" 2>&1 &
    pids+=("$!")
done
status=0
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
test "$status" -eq 0
(cd "$out" && sha256sum plot.py job.sbatch source_revision.txt fine_spectra.pdf fine_spectra.png > artifacts.sha256)
(cd "$out" && sha256sum -c artifacts.sha256 --status)
touch "$out/COMPLETE"
