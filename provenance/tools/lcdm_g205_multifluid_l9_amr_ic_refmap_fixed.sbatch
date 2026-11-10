#!/bin/bash
#SBATCH --job-name=rar-g205-mapfix
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_amr_ic_refmap_fixed_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_amr_ic_refmap_fixed_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l9_1284259"
audit="$base/analysis/lcdm_g205_multifluid_l9_audit_1284269"
template="$project/configs/lcdm_grafic_n256_reader_preflight_20261103.nml.in"
build_record="$base/builds/rar_amr_ic_refmap_fixed_20261111"
root="$base/runs/lcdm_g205_multifluid_l9_amr_ic_refmap_fixed_${SLURM_JOB_ID:?}"

module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
command -v mpiifx
test -f "$build_record/COMPLETE"
(cd "$build_record" && sha256sum -c PAYLOAD.sha256 --status)
(cd "$ic" && sha256sum -c artifacts.sha256 --status)
(cd "$audit" && sha256sum -c artifacts.sha256 --status)
test -f "$ic/GENERATION_ONLY"
test -f "$audit/GEOMETRY_ONLY_COMPLETE"
test ! -e "$root"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge 34359738368
mkdir "$root"
trap 'rc=$?; if [ "$rc" -ne 0 ] && [ -d "$root" ]; then printf "FIXED_IC_MAP_GATE_FAILED rc=%s\n" "$rc" > "$root/outcome.txt"; fi' ERR
ln -s "$ic/ics/level_008" "$root/ic"
ln -s "$ic/ics/level_009" "$root/ic9"
cp -p "$build_record/ramses_final3d" "$root/ramses_final3d"
cp -p "$template" "$root/template.nml.in"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
test "$(sha256sum "$root/ramses_final3d" | awk '{print $1}')" = \
     "$(awk 'NR==1 {print $1}' "$build_record/executable.sha256")"

sed \
    -e "s|__IC_DIRECTORY__|ic|g" \
    -e "/initfile(1)=/a\\  initfile(2)='ic9'" \
    -e 's/levelmax=8/levelmax=9/' \
    -e 's/boxlen=25.0/boxlen=1.0/' \
    -e 's/nparttot=16777216/nparttot=18000000/' \
    -e 's/m_refine=8\*8\./m_refine=100*-1./' \
    -e '/&AMR_PARAMS/a\\  ngridmax_auto=.false.' \
    "$template" | awk '
{
    if ($0 == "&REFINE_PARAMS") {
        print
        print "  ivar_refine=0"
        print "  ic_refmap_fixed=.true."
        next
    }
    print
}' > "$root/effective.nml"

for key in 'cosmo=.true.' 'hydro=.true.' 'pic=.true.' 'poisson=.true.' \
    'use_fftw=.true.' 'use_mond=.false.' 'levelmin=8' 'levelmax=9' \
    'ivar_refine=0' 'ic_refmap_fixed=.true.' \
    'ngridmax=3000000' 'ngridtot=3000000' 'ngridmax_auto=.false.' \
    'npartmax=20000000' 'nparttot=18000000' \
    'nstepmax=0' 'noutput=1' 'aout=1.1' 'tout=1.0d100' \
    'foutput=1000000000' 'fbackup=1000000000' \
    'abort_on_mg_nonconvergence=.true.' 'omega_b=0.04930169880390167'; do
    grep -Fq "$key" "$root/effective.nml"
done
test "$(grep -Fc 'ivar_refine=0' "$root/effective.nml")" -eq 1
test "$(grep -Fc 'ic_refmap_fixed=.true.' "$root/effective.nml")" -eq 1
grep -Fq 'm_refine=100*-1.' "$root/effective.nml"

for level in ic ic9; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz \
        ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
        test -s "$root/$level/$field"
    done
done

printf '%s\n' "$root/effective.nml" > "$root/effective_namelist_path.txt"
printf '%s\n' "$ic" > "$root/input_ic_path.txt"

printf 'purpose=one-coarse-step fixed-GRAFIC-IC-map AMR diagnostic; not production evolution\ncosmology=Omega_m 0.315, Omega_b 0.0493017, Omega_Lambda 0.685, h 0.6736\ntransfer=CAMB Lambda-CDM two-fluid; transfer normalization remains a separate gate\nlevels=8,9\nrefinement=ic_refmap_fixed=true, ivar_refine=0; initial GRAFIC mask applied and existing child topology preserved; density-derived map disabled\nsolver=FFTW base-level direct solve, refined MG, abort on MG nonconvergence\nphysics=hydro and CDM particles on; cooling, metal enrichment, star formation, sink/SMBH, AGN, radiation and dust off\nnstepmax=0; one coarse step may execute\nnoutput=1; requested aout=1.1 and tout=1e100; foutput=fbackup=1000000000\nexpected_RAMSES_dumps=0\nconservative_accidental_dump_estimate_bytes=17179869184\nexpected_total_output_bytes=0\nminimum_free_bytes_before=34359738368\nfree_bytes_before=%s\n' "$free_bytes" > "$root/output_policy.txt"

sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/template.nml.in" "$root/job.sbatch" \
    "$root/effective_namelist_path.txt" "$root/input_ic_path.txt" > "$root/inputs.sha256"

export OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_AMR_MAP_DIAGNOSTICS=1
set +e
(cd "$root" && srun --mpi=pmi2 --ntasks=1 --cpu-bind=cores \
    ./ramses_final3d effective.nml > run.log 2>&1)
run_rc=$?
set -e
printf '%s\n' "$run_rc" > "$root/srun_exit_code.txt"

rg -n '\[RAR_AMR_MAP\]|GRID_CAPACITY_FIXED|GRID_GROW|Level[[:space:]]+8 has|Level[[:space:]]+9 has|FFT direct solve DONE|WARN:' \
    "$root/run.log" > "$root/key_diagnostics.txt" || true
if [ "$run_rc" -ne 0 ]; then
    printf '%s\n' 'FIXED_IC_MAP_RUN_FAILED' > "$root/outcome.txt"
    exit 1
fi
for marker in 'Run completed' 'IC/namelist omega_b consistent' \
    'FFTW3 CPU direct Poisson solver enabled' 'FFT direct solve DONE' \
    'IC refinement map fixed' \
    '[RAR_AMR_MAP] source=ic_refmap level=8 active=16777216 flagged=12058'; do
    grep -Fq "$marker" "$root/run.log"
done
grep -Eq 'Reading file ic/ic_pvar_00001' "$root/run.log"
grep -Eq 'Reading file ic9/ic_pvar_00001' "$root/run.log"
grep -Eq 'Level[[:space:]]+8 has[[:space:]]+2097152 grids' "$root/run.log"

if grep -Fq '[RAR_AMR_MAP] source=density' "$root/run.log"; then
    printf '%s\n' 'unexpected density-derived refinement map' >&2
    printf '%s\n' 'FIXED_IC_MAP_DENSITY_MAP_REAPPEARED' > "$root/outcome.txt"
    exit 1
fi
if grep -Eiq 'ERROR:|FATAL:|Aborting|GRID_CAPACITY_FIXED|GRID_GROW|WARN: Fine multigrid Poisson failed to converge|NaN/Inf detected|NaN detected|Floating point exception|segmentation fault' "$root/run.log"; then
    printf '%s\n' 'FIXED_IC_MAP_PHYSICS_OR_RUNTIME_WARNING' > "$root/outcome.txt"
    exit 1
fi
if compgen -G "$root/output_[0-9]*" >/dev/null || \
   compgen -G "$root/backup_[0-9]*" >/dev/null; then
    printf '%s\n' 'FIXED_IC_MAP_UNEXPECTED_OUTPUT' > "$root/outcome.txt"
    exit 1
fi

awk '
/source=ic_refmap_preserved level=8 / {
    seen=1
    for (i=1;i<=NF;i++) if ($i ~ /^flagged=/) {
        split($i,a,"="); n=a[2]+0
    }
}
END { if (!seen || n<1 || n>100000) exit 1 }
' "$root/run.log"
last_level9_grids=$(awk '
/Level[[:space:]]+9 has/ {
    for (i=1;i<=NF;i++) if ($i=="has") {n=$(i+1)+0; seen=1}
}
END { if (!seen || n<1 || n>100000) exit 1; print n }
' "$root/run.log")

printf 'FIXED_IC_MASK_ONE_STEP_GATE_PASSED\nlast_level9_grids=%s\n' \
    "$last_level9_grids" > "$root/outcome.txt"
sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/output_policy.txt" "$root/run.log" "$root/job.sbatch" \
    "$root/effective_namelist_path.txt" "$root/input_ic_path.txt" \
    "$root/key_diagnostics.txt" "$root/outcome.txt" \
    "$root/srun_exit_code.txt" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'FIXED_IC_MAP_DIAGNOSTIC_CAPTURED' > "$root/FIXED_IC_MAP_DIAGNOSTIC_CAPTURED"
printf 'Fixed IC-map diagnostic captured. Run directory: %s\n' "$root"
