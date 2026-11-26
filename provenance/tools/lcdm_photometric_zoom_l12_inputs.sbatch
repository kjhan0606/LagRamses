#!/bin/bash
#SBATCH --job-name=rar-photo-L12-input
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/photometric_l12_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/photometric_l12_%j.err
set -euo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
repo=/home/kjhan/BackUP/rar/build/rar_qumond_amr_history_20261107
rev=$1
test "$(git -C "$repo" rev-parse "$rev^{commit}")" = "$rev"
trace=$base/analysis/lcdm_photometric_zoom_1578501
parent=$base/initial_conditions/lcdm_g205_multifluid_l11_1286049
(cd "$trace" && test -f COMPLETE && sha256sum -c results.sha256 --status)
(cd "$parent" && sha256sum -c inputs.sha256 --status)
out=$base/initial_conditions/lcdm_photometric_l12_inputs_$SLURM_JOB_ID
mkdir "$out"
cd "$out"
trap 'touch FAILED' ERR
cp "$trace"/{points_g*.txt,shortlist.csv,summary.json} .
cp "$parent/transfer.dat" parent_transfer.dat
cp "$parent/effective.conf" base_l11.conf
for script in lcdm_g205_z99_highk_camb lcdm_photometric_zoom_l12_inputs; do
    git -C "$repo" show "$rev:provenance/tools/$script.py" > "$script.py"
done
cp "$0" job.sbatch
printf '%s\n' "$rev" > code_revision.txt
sha256sum *.py *.txt *.csv *.json *.sbatch *.conf *.dat > inputs.sha256
export RAR_CAMB_REDSHIFT=49 RAR_CAMB_KMAX_MPC=1024 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for n in 100 200; do
    mkdir "$n"
    RAR_CAMB_SAMPLING=$n RAR_CAMB_OUTPUT=$n /home/kjhan/miniconda3/bin/python3 lcdm_g205_z99_highk_camb.py > "$n/log.txt" 2>&1
done
/home/kjhan/miniconda3/bin/python3 lcdm_photometric_zoom_l12_inputs.py
sha256sum -c inputs.sha256 --status
sha256sum {100,200}/{settings.json,transfer_z49.dat} g*/{points.txt,music.conf} preparation.json transfer_checks.csv > results.sha256
touch COMPLETE
