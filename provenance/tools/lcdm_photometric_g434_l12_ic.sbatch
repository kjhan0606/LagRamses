#!/bin/bash
#SBATCH --job-name=rar-g434-mf-l12-ic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g434_multifluid_l12_ic_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g434_multifluid_l12_ic_%j.err
set -Eeuo pipefail
umask 077

project=/home/kjhan/BackUP/rar
repo=$project/build/rar_qumond_amr_history_20261107
base=/gpfs/kjhan/RAR_QUMOND
inputs=$base/initial_conditions/lcdm_photometric_l12_inputs_1578517
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
music=$project/build/lagmusic_qumond_cpu_20261103/MUSIC
music_sha=687220cdb68d6ee9fafd4da3a23b1ef751dcb10d3ce621a183ee1fd66f9fe415
revision=${1:?Pass the committed full Git revision}
root=$base/initial_conditions/lcdm_g434_multifluid_l12_${SLURM_JOB_ID:?}

test "$(git -C "$repo" rev-parse "$revision^{commit}")" = "$revision"
test -f "$inputs/COMPLETE"
(cd "$inputs" && sha256sum -c results.sha256 --status)
test -f "$parent/COMPLETE"
printf '%s  %s\n' "$music_sha" "$music" | sha256sum -c --status
test ! -e "$root"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 137438953472
mkdir "$root" "$root/g434" "$root/200" "$root/analysis"
trap 'touch "$root/FAILED"' ERR

cp -p "$inputs/g434/points.txt" "$root/g434/points.txt"
cp -p "$inputs/200/transfer_z49.dat" "$root/200/transfer_z49.dat"
cp -p "$inputs/200/settings.json" "$root/200/settings.json"
cp -p "$inputs/g434/music.conf" "$root/g434/original_music.conf"
sed 's|^filename = ics$|filename = ../ics|' "$root/g434/original_music.conf" > "$root/g434/music.conf"
cp -p "$music" "$root/MUSIC"
cp -p "$0" "$root/job.sbatch"
awk -v p="$parent" '$2 == p "/lcdm_n256.conf" ||
    $2 == p "/lcdm_control_transfer_out_z49.dat" ||
    $2 ~ /\/ics_lcdm\/level_008\/ic_(deltab|posc[xyz]|vel[bc][xyz])$/ {print}' \
    "$parent/artifacts.sha256" > "$root/parent_inputs.sha256"
test "$(wc -l < "$root/parent_inputs.sha256")" -eq 12
sha256sum -c "$root/parent_inputs.sha256" --status
for tool in lcdm_photometric_l12_geometry lcdm_photometric_base_coherence lcdm_base_spectra_camb; do
    git -C "$repo" show "$revision:provenance/tools/$tool.py" > "$root/$tool.py"
done
printf '%s\n' "$revision" > "$root/source_revision.txt"

for setting in 'boxlength = 25' 'zstart = 49' 'levelmin = 8' 'levelmax = 12' \
    'region = ellipsoid' 'region_point_level = 8' 'region_point_file = points.txt' \
    'baryons = yes' 'use_2LPT = no' 'Omega_m = 0.315' 'Omega_b = 0.0493017' \
    'seed[8] = 20261103' 'seed[9] = 20261104' 'seed[10] = 20261105' \
    'seed[11] = 20261106' 'seed[12] = 20261107' 'format = grafic2' \
    'transfer_file = ../200/transfer_z49.dat' 'filename = ../ics'; do
    grep -Fxq "$setting" "$root/g434/music.conf"
done
printf '%s\n' "$root/g434/music.conf" > "$root/effective_config_path.txt"
printf 'purpose=group434 two-fluid Lambda-CDM MUSIC IC generation and validation\nparent_box=25 Mpc/h; zstart=49; levels=8-12\nRAMSES_evolution=none; RAMSES_dumps=0\nprior_g205_L11_IC_and_work_bytes=4924539990\nL12_eightfold_scaling_planning_bytes=39396319920\nIC_and_work_budget_bytes=103079215104\nminimum_free_before_bytes=137438953472\nfree_before_bytes=%s\n' "$free_bytes" > "$root/storage_policy.txt"
(cd "$root" && sha256sum MUSIC g434/points.txt g434/music.conf g434/original_music.conf \
    200/transfer_z49.dat 200/settings.json job.sbatch source_revision.txt \
    parent_inputs.sha256 lcdm_photometric_l12_geometry.py lcdm_photometric_base_coherence.py \
    lcdm_base_spectra_camb.py > inputs.sha256)

export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
(cd "$root/g434" && srun --mpi=pmi2 --cpu-bind=cores "$root/MUSIC" music.conf > "$root/music.log" 2>&1)
if grep -Eiq 'ERROR:|FATAL|invalid conversion|conversion of .nan. failed|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$root/music.log"; then
    exit 1
fi
for level in 008 009 010 011 012; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz \
        ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
        test -s "$root/ics/level_$level/$field"
    done
done
output_bytes=$(du -sb "$root/ics" "$root/g434" | awk '{total += $1} END {print total}')
test "$output_bytes" -le 103079215104
printf 'actual_ic_and_work_bytes=%s\n' "$output_bytes" >> "$root/storage_policy.txt"

env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 \
    "$root/lcdm_photometric_l12_geometry.py" "$root" 434 2196 > "$root/geometry_audit.log"
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 \
    "$root/lcdm_photometric_base_coherence.py" "$root" "$parent" "$root/analysis" 434 > "$root/analysis.log"
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 \
    "$root/lcdm_base_spectra_camb.py" "$root/ics/level_008" "$root/200/transfer_z49.dat" \
    "$root/g434/input_powerspec.txt" "$root/g434/music.conf_log.txt" \
    "$root/analysis/camb" > "$root/camb_spectra.log"
(cd "$root" && sha256sum -c inputs.sha256 --status)
sha256sum -c "$root/parent_inputs.sha256" --status
(cd "$root" && find . -type f ! -name artifacts.sha256 ! -name COMPLETE \
    -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256)
(cd "$root" && sha256sum -c artifacts.sha256 --status)
printf 'generation=complete\ngeometry_parent_coherence_and_camb_spectra=measured\nRAMSES_evolution=none\n' > "$root/COMPLETE"
