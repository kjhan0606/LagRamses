#!/bin/bash
#SBATCH --job-name=rar-g205-mf-l11-ic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_ic_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_ic_%j.err
set -Eeuo pipefail
umask 077

project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
trace=$base/initial_conditions/lcdm_z0_group205_r3_ellipsoid_l11_1273554
music=$project/build/lagmusic_qumond_cpu_20261103/MUSIC
music_sha=687220cdb68d6ee9fafd4da3a23b1ef751dcb10d3ce621a183ee1fd66f9fe415
root=$base/initial_conditions/lcdm_g205_multifluid_l11_$SLURM_JOB_ID

for directory in "$parent" "$trace"; do
    test -f "$directory/artifacts.sha256"
    (cd "$directory" && sha256sum -c artifacts.sha256 --status)
done
test -f "$parent/COMPLETE"
test -f "$trace/GEOMETRY_AUDIT_COMPLETE"
printf '%s  %s\n' "$music_sha" "$music" | sha256sum -c -
test ! -e "$root"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 42949672960
mkdir "$root" "$root/work"

cp -p "$parent/lcdm_control_transfer_out_z49.dat" "$root/transfer.dat"
cp -p "$parent/lcdm_n256.conf" "$root/parent.conf"
cp -p "$trace/points.txt" "$root/points.txt"
cp -p "$music" "$root/MUSIC"
cp -p "$project/scripts/lcdm_g205_multifluid_l11_audit.py" "$root/audit.py"
cp -p "$0" "$root/job.sbatch"

sed \
    -e 's/levelmax = 8/levelmax = 11\nregion = ellipsoid\nregion_point_level = 8/' \
    -e "/region_point_level = 8/a region_point_file = $root/points.txt" \
    -e 's/seed\[8\] = 20261103/seed[8] = 20261103\nseed[9] = 20261104\nseed[10] = 20261105\nseed[11] = 20261106/' \
    -e "s|^transfer_file = .*|transfer_file = $root/transfer.dat|" \
    -e "s|^filename = .*|filename = $root/ics|" \
    "$root/parent.conf" > "$root/effective.conf"

for key in 'boxlength = 25' 'zstart = 49' 'levelmin = 8' 'levelmax = 11' \
    'region = ellipsoid' 'region_point_level = 8' \
    "region_point_file = $root/points.txt" 'baryons = yes' 'use_2LPT = no' \
    'Omega_m = 0.315' 'Omega_b = 0.0493017' \
    'seed[8] = 20261103' 'seed[9] = 20261104' \
    'seed[10] = 20261105' 'seed[11] = 20261106' 'format = grafic2' \
    "transfer_file = $root/transfer.dat" "filename = $root/ics"; do
    grep -Fxq "$key" "$root/effective.conf"
done

printf '%s\n' "$root/effective.conf" > "$root/effective_config_path.txt"
printf 'purpose=Group 205 LCDM two-fluid level-11 zoom IC generation and geometry audit; no RAMSES evolution\nmodel=Omega_m=0.315, Omega_b=0.0493017, baryons=yes\nbox=25 Mpc/h, zstart=49\nlevels=8,9,10,11\ntrace=4,536 initial cells within provisional 3R200c\nparent_seed=20261103\npreserved_level9_seed=20261104\nnew_level10_seed=20261105\nnew_level11_seed=20261106\nparticle_resolution=level-11 CDM mass estimate 1.991e5 Msun; provisional Group205 M200c predicts about 1.30e6 CDM particles inside R200c; evolved count is not established\nRAMSES_dumps=0\nIC_and_work_budget_bytes=34359738368\nminimum_free_headroom_bytes=8589934592\nfree_bytes_before=%s\ntransfer_normalization=not certified by geometry audit\n' "$free_bytes" > "$root/storage_policy.txt"

sha256sum "$root/MUSIC" "$root/effective.conf" "$root/points.txt" \
    "$root/transfer.dat" "$root/parent.conf" "$root/audit.py" \
    "$root/job.sbatch" > "$root/inputs.sha256"

export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
(cd "$root/work" && srun --mpi=pmi2 --cpu-bind=cores "$root/MUSIC" \
    "$root/effective.conf" > "$root/music.log" 2>&1)
if grep -Eiq 'ERROR:|FATAL|invalid conversion|conversion of .nan. failed|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' \
    "$root/music.log"; then
    exit 1
fi

for level in 008 009 010 011; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz \
        ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
        test -s "$root/ics/level_$level/$field"
    done
done
output_bytes=$(du -sb "$root/ics" "$root/work" | awk '{total += $1} END {print total}')
test "$output_bytes" -le 34359738368
env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 \
    "$root/audit.py" "$root" "$parent" "$trace" "$root" \
    > "$root/geometry_audit.log"
printf 'actual_ic_and_work_bytes=%s\n' "$output_bytes" >> "$root/storage_policy.txt"
sha256sum -c "$root/inputs.sha256" --status
(cd "$root" && find . -maxdepth 4 -type f \
    ! -name artifacts.sha256 ! -name GENERATION_ONLY \
    ! -name GEOMETRY_AUDIT_COMPLETE -print0 | sort -z | xargs -0 sha256sum \
    > artifacts.sha256)
sha256sum -c "$root/artifacts.sha256" --status
printf 'generation_state=GENERATION_ONLY\ngeometry_audit=passed\nRAMSES_evolution=not_run\n' \
    > "$root/GENERATION_ONLY"
touch "$root/GEOMETRY_AUDIT_COMPLETE"
printf 'Multi-fluid L11 IC geometry audit completed. Directory: %s\n' "$root"
