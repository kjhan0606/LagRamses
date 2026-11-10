#!/bin/bash
#SBATCH --job-name=rar-g205-twofluid-reader
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/lcdm_g205_twofluid_reader_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/lcdm_g205_twofluid_reader_%j.err
set -Eeuo pipefail

project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l9_1284259"
template="$project/configs/lcdm_grafic_n256_reader_preflight_20261103.nml.in"
source_exe="$project/build/qumond_fix_cpu_20261103/bin/ramses_final3d"
exe_sha=e8a3eec8fc0b9823e6f7462a627866438623ee47b664268382c793a9f38e5b8b
root="$base/runs/lcdm_g205_multifluid_l9_reader_${SLURM_JOB_ID:?}"

test -f "$ic/GENERATION_ONLY"
test -f "$base/analysis/lcdm_g205_multifluid_l9_audit_1284269/GEOMETRY_ONLY_COMPLETE"
(cd "$base/analysis/lcdm_g205_multifluid_l9_audit_1284269" && sha256sum -c artifacts.sha256 --status)
(cd "$ic" && sha256sum -c artifacts.sha256)
printf '%s  %s\n' "$exe_sha" "$source_exe" | sha256sum -c -
test ! -e "$root"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 17179869184
mkdir "$root"
ln -s "$ic/ics/level_008" "$root/ic"
ln -s "$ic/ics/level_009" "$root/ic9"
test "${#root}" -lt 76
cp -p "$source_exe" "$root/ramses_final3d"
cp -p "$template" "$root/template.nml.in"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"

sed \
    -e "s|__IC_DIRECTORY__|$root/ic|g" \
    -e "/initfile(1)=/a\\  initfile(2)='$root/ic9'" \
    -e 's/levelmax=8/levelmax=9/' \
    -e 's/boxlen=25.0/boxlen=1.0/' \
    -e 's/nparttot=16777216/nparttot=18000000/' \
    -e '/&REFINE_PARAMS/a\  ivar_refine=0' \
    "$template" > "$root/effective.nml"

for key in 'cosmo=.true.' 'hydro=.true.' 'pic=.true.' 'poisson=.true.' \
    'use_fftw=.true.' 'use_mond=.false.' 'levelmin=8' 'levelmax=9' \
    'boxlen=1.0' 'ivar_refine=0' "initfile(2)='$root/ic9'" \
    'npartmax=20000000' 'nparttot=18000000' 'ngridmax=3000000' \
    'ngridtot=3000000' 'nstepmax=0' 'noutput=1' 'aout=1.1' \
    'tout=1.0d100' 'foutput=1000000000' 'fbackup=1000000000' \
    'abort_on_mg_nonconvergence=.true.' 'omega_b=0.04930169880390167'; do
    grep -Fq "$key" "$root/effective.nml"
done

for level in ic ic9; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz; do
        test -s "$root/$level/$field"
    done
done

printf '%s\n' "$root/effective.nml" > "$root/effective_namelist_path.txt"
printf 'purpose=standard Lambda-CDM gas and CDM multilevel GRAFIC reader and base FFTW/refined MG initialization preflight\ncosmology=Omega_m 0.315, Omega_b 0.0493017, Omega_Lambda 0.685, h 0.6736\ntransfer=CAMB LCDM input; independent field normalization validation pending\nlevelmin=8\nlevelmax=9\nivar_refine=0\nnstepmax=0\nnoutput=1\naout=1.1\ntout=1e100\nfoutput=fbackup=1000000000\nexpected_RAMSES_full_dumps=0\naccidental_dump_headroom_bytes=17179869184\nngridtot=3000000\nnparttot=18000000\nfree_bytes_before=%s\nhydro=true\ncooling=false\nmetal=false\nsink=false\nAGN=false\nrt=false\ndust=false\n' "$free_bytes" > "$root/output_policy.txt"
sha256sum "$root/effective.nml" "$root/ramses_final3d" "$root/template.nml.in" > "$root/inputs.sha256"
printf '%s\n' "$ic" > "$root/input_ic_path.txt"

export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
(cd "$root" && srun --mpi=pmi2 --ntasks=1 --cpu-bind=cores ./ramses_final3d effective.nml > run.log 2>&1)
for marker in 'Run completed' 'IC/namelist omega_b consistent' \
    'FFTW3 CPU direct Poisson solver enabled' 'FFT direct solve DONE'; do
    grep -Fq "$marker" "$root/run.log"
done
grep -Eq 'Level[[:space:]]+8 has[[:space:]]+2097152 grids' "$root/run.log"
grep -Eq 'Level[[:space:]]+9 has[[:space:]]+[1-9][0-9]* grids' "$root/run.log"
grep -Eq 'n1=[[:space:]]*256 n2=[[:space:]]*256 n3=[[:space:]]*256' "$root/run.log"
if grep -Eiq '(^|[[:space:]])(ERROR:|FATAL:|Aborting)|failed to converge|NaN/Inf detected|segmentation fault' "$root/run.log"; then
    exit 1
fi
if compgen -G "$root/output_[0-9]*" >/dev/null || compgen -G "$root/backup_[0-9]*" >/dev/null; then echo "Unexpected output or backup" >&2; exit 1; fi
sha256sum -c "$root/inputs.sha256"
(cd "$root" && sha256sum job.sbatch template.nml.in effective.nml \
    effective_namelist_path.txt input_ic_path.txt output_policy.txt inputs.sha256 \
    ramses_final3d run.log > artifacts.sha256 && sha256sum -c artifacts.sha256)
touch "$root/COMPLETE"
