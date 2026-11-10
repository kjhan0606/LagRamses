#!/bin/bash
#SBATCH --job-name=rar-g205-mapdiag
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_amr_map_diag_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_amr_map_diag_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l9_1284259"
audit="$base/analysis/lcdm_g205_multifluid_l9_audit_1284269"
template="$project/configs/lcdm_grafic_n256_reader_preflight_20261103.nml.in"
build_record="$base/builds/rar_amr_map_diag_20261111_v2"
root="$base/runs/lcdm_g205_multifluid_l9_amrmapdiag_${SLURM_JOB_ID:?}"

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
test "$free_bytes" -ge 17179869184
mkdir "$root"
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
    -e '/&REFINE_PARAMS/a\\  ivar_refine=0' \
    -e '/&AMR_PARAMS/a\\  ngridmax_auto=.false.' \
    "$template" > "$root/effective.nml"

for key in 'cosmo=.true.' 'hydro=.true.' 'pic=.true.' 'poisson=.true.' \
    'use_fftw=.true.' 'use_mond=.false.' 'levelmin=8' 'levelmax=9' \
    'ivar_refine=0' 'ngridmax=3000000' 'ngridtot=3000000' \
    'ngridmax_auto=.false.' 'npartmax=20000000' 'nparttot=18000000' \
    'nstepmax=0' 'noutput=1' 'aout=1.1' 'tout=1.0d100' \
    'foutput=1000000000' 'fbackup=1000000000' \
    'abort_on_mg_nonconvergence=.true.' 'omega_b=0.04930169880390167'; do
    grep -Fq "$key" "$root/effective.nml"
done
grep -Fq 'm_refine=8*8.' "$root/effective.nml"

for level in ic ic9; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz \
        ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
        test -s "$root/$level/$field"
    done
done

printf 'purpose=one-coarse-step AMR refinement-map diagnosis; not production evolution\ncosmology=Omega_m 0.315, Omega_b 0.0493017, Omega_Lambda 0.685, h 0.6736\ntransfer=CAMB Lambda-CDM two-fluid; transfer normalization remains a separate gate\nlevels=8,9\nrefinement=original m_refine and ivar_refine settings; automatic grid growth disabled; capacity=3000000\nsolver=FFTW base-level direct solve, refined MG, abort on MG nonconvergence\nphysics=hydro and CDM particles on; cooling, metal enrichment, star formation, sink/SMBH, AGN, radiation and dust off\nnstepmax=0; one coarse step may execute\nnoutput=1; dummy aout=1.1 and tout=1e100; foutput=fbackup=1000000000\nexpected_RAMSES_dumps=0\nestimated_dump_bytes=0 for this bounded diagnostic; the future dummy schedule is not reached\nfree_bytes_before=%s\n' "$free_bytes" > "$root/output_policy.txt"

sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/template.nml.in" "$root/job.sbatch" > "$root/inputs.sha256"

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

grep -Fq '[RAR_AMR_MAP] source=ic_refmap level=9' "$root/run.log"
grep -Fq '[RAR_AMR_MAP] source=density level=8' "$root/run.log"
if grep -Eq 'WARN: Fine multigrid Poisson failed to converge|NaN detected|Floating point exception' "$root/run.log"; then
    exit 1
fi
if compgen -G "$root/output_[0-9]*" >/dev/null || \
   compgen -G "$root/backup_[0-9]*" >/dev/null; then
    exit 1
fi

if grep -Fq '[RESIZABLE] GRID_CAPACITY_FIXED' "$root/run.log"; then
    printf '%s\n' 'EXPECTED_CAPACITY_GUARD' > "$root/outcome.txt"
elif [ "$run_rc" -eq 0 ]; then
    printf '%s\n' 'ONE_STEP_COMPLETED_WITHIN_FIXED_CAPACITY' > "$root/outcome.txt"
else
    printf '%s\n' 'UNEXPECTED_RUN_FAILURE' > "$root/outcome.txt"
    exit 1
fi

rg -n '\[RAR_AMR_MAP\]|GRID_CAPACITY_FIXED|GRID_GROW|Level  8 has|Level  9 has|FFT direct solve DONE|WARN:' \
    "$root/run.log" > "$root/key_diagnostics.txt" || true
sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/output_policy.txt" "$root/run.log" "$root/job.sbatch" \
    "$root/outcome.txt" "$root/srun_exit_code.txt" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'MAP_DIAGNOSTIC_CAPTURED' > "$root/MAP_DIAGNOSTIC_CAPTURED"
printf 'Map diagnosis captured. Run directory: %s\n' "$root"
