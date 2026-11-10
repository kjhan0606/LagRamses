#!/bin/bash
#SBATCH --job-name=rar-g205-twofluid-ic
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_twofluid_l9_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_twofluid_l9_%j.err
set -Eeuo pipefail
umask 077
project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
parent=$base/initial_conditions/qumond_lcdm_grafic_z49_n256_1260246
trace=$base/initial_conditions/lcdm_z0_group205_r3_ellipsoid_l11_1273554
music=$project/build/lagmusic_qumond_cpu_20261103/MUSIC
root=$base/initial_conditions/lcdm_g205_multifluid_l9_${SLURM_JOB_ID:?}
test -f "$parent/COMPLETE" && test -f "$trace/GEOMETRY_AUDIT_COMPLETE" && test ! -e "$root"
(cd "$parent" && sha256sum -c artifacts.sha256 --status)
(cd "$trace" && sha256sum -c artifacts.sha256 --status)
printf '%s  %s\n' 687220cdb68d6ee9fafd4da3a23b1ef751dcb10d3ce621a183ee1fd66f9fe415 "$music" | sha256sum -c -
free=$(df -B1 "$base" | awk 'END {print $4}'); test "$free" -ge 21474836480
mkdir "$root" "$root/work"
cp -p "$parent/lcdm_control_transfer_out_z49.dat" "$root/transfer.dat"
cp -p "$trace/points.txt" "$root/points.txt"
cp -p "$parent/lcdm_n256.conf" "$root/parent.conf"
cp -p "$music" "$root/MUSIC"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
sed -e 's/levelmax = 8/levelmax = 9\nregion = ellipsoid\nregion_point_level = 8/' -e "/region_point_level = 8/a region_point_file = $root/points.txt" -e 's/seed\[8\] = 20261103/seed[8] = 20261103\nseed[9] = 20261104/' -e "s|^transfer_file = .*|transfer_file = $root/transfer.dat|" -e "s|^filename = .*|filename = $root/ics|" "$root/parent.conf" > "$root/effective.conf"
for expected in 'boxlength = 25' 'zstart = 49' 'levelmin = 8' 'levelmax = 9' 'region = ellipsoid' 'baryons = yes' 'use_2LPT = no' 'Omega_m = 0.315' 'Omega_b = 0.0493017' 'seed[8] = 20261103' 'seed[9] = 20261104' 'format = grafic2'; do grep -Fxq "$expected" "$root/effective.conf"; done
printf 'purpose=LCDM Group 205 two-fluid nested IC generation feasibility, not production\nparent_group=205\ntrace=4536 initial cells, 3R200c of one-fluid DMO target\nlevels=8,9\nfull_simulation_dumps=0\nreserved_ic_and_work_bytes=12884901888\nfree_headroom_bytes=8589934592\nfree_bytes_before=%s\nreader_and_phase_validation=required\n' "$free" > "$root/storage_policy.txt"
(cd "$root" && sha256sum MUSIC effective.conf points.txt transfer.dat parent.conf job.sbatch > inputs.sha256)
export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
(cd "$root/work" && srun --mpi=pmi2 --cpu-bind=cores "$root/MUSIC" "$root/effective.conf" > "$root/music.log" 2>&1)
if grep -Eiq 'ERROR:|FATAL|invalid conversion|conversion of .nan. failed|(^|[^[:alpha:]])nan([^[:alpha:]]|$)' "$root/music.log"; then exit 1; fi
for level in 008 009; do for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz; do test -s "$root/ics/level_$level/$field"; done; done
test "$(du -sb "$root/ics" "$root/work" | awk '{s+=$1}END{print s}')" -le 12884901888
(cd "$root" && sha256sum -c inputs.sha256 --status && find . -type f ! -name artifacts.sha256 ! -name GENERATION_ONLY -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256 && sha256sum -c artifacts.sha256 --status)
printf 'Two-fluid IC generation only. No phase, geometry, reader or evolution validation.\n' > "$root/GENERATION_ONLY"
