#!/bin/bash
#SBATCH --job-name=rar-g205-l11-unitbox
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_mpi_unitbox_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/g205_multifluid_l11_mpi_unitbox_%j.err
set -Eeuo pipefail
umask 077
ulimit -c 0

project=/home/kjhan/BackUP/rar
base=/gpfs/kjhan/RAR_QUMOND
ic="$base/initial_conditions/lcdm_g205_multifluid_l11_1286049"
template="$project/configs/lcdm_grafic_n256_reader_preflight_20261103.nml.in"
build_record="$base/builds/rar_amr_ic_refmap_fixed_20261111"
root="$base/runs/lcdm_g205_multifluid_l11_mpi_unitbox_${SLURM_JOB_ID:?}"
expected_ic_status=NESTED_L11_LCDM_TWOFLUID_IC_GEOMETRY_PASS_NO_EVOLUTION
expected_template_sha=b8782fd4e0e06eae478c97d4206c22a10370c0b9747d2187eca9a15ac56658fd
expected_executable_sha=0537739ba9522bb923e53a3025f817968abff8b0d4f3794d4949e20831e9339d
minimum_free_bytes=206158430208

module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
command -v mpiifx
test -f "$build_record/COMPLETE"
(cd "$build_record" && sha256sum -c PAYLOAD.sha256 --status)
test "$(awk 'NR==1 {print $1}' "$build_record/executable.sha256")" = "$expected_executable_sha"
test "$(sha256sum "$template" | awk '{print $1}')" = "$expected_template_sha"
(cd "$ic" && sha256sum -c artifacts.sha256 --status)
test -f "$ic/GENERATION_ONLY"
test -f "$ic/GEOMETRY_AUDIT_COMPLETE"
grep -Fq "\"status\": \"$expected_ic_status\"" "$ic/geometry_audit.json"
test ! -e "$root"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -ge "$minimum_free_bytes"
mkdir "$root"
trap 'rc=$?; if [ "$rc" -ne 0 ] && [ -d "$root" ]; then printf "L11_READER_PREFLIGHT_FAILED rc=%s\n" "$rc" > "$root/outcome.txt"; fi' ERR

for level in 008 009 010 011; do
    level_num=$((10#$level))
    ln -s "$ic/ics/level_$level" "$root/ic$level_num"
done
cp -p "$build_record/ramses_final3d" "$root/ramses_final3d"
cp -p "$template" "$root/template.nml.in"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
test "$(sha256sum "$root/ramses_final3d" | awk '{print $1}')" = \
     "$(awk 'NR==1 {print $1}' "$build_record/executable.sha256")"

sed \
    -e 's|__IC_DIRECTORY__|ic8|g' \
    -e "/initfile(1)=/a\\\\  initfile(2)='ic9'" \
    -e "/initfile(1)=/a\\\\  initfile(3)='ic10'" \
    -e "/initfile(1)=/a\\\\  initfile(4)='ic11'" \
    -e 's/levelmax=8/levelmax=11/' \
    -e 's/boxlen=25.0/boxlen=1.0/' \
    -e 's/nsubcycle=1,2,2,2,2,2,2,2,2/nsubcycle=1,2,2,2,2,2,2,2,2,2,2/' \
    -e 's/m_refine=8\*8\./m_refine=100*-1./' \
    -e 's/ngridmax=3000000/ngridmax=4500000/' \
    -e 's/ngridtot=3000000/ngridtot=4500000/' \
    -e 's/npartmax=20000000/npartmax=32000000/' \
    -e 's/nparttot=16777216/nparttot=28000000/' \
    -e '/&AMR_PARAMS/a\\  ngridmax_auto=.false.' \
    "$root/template.nml.in" | awk '
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
    'use_fftw=.true.' 'use_mond=.false.' 'levelmin=8' 'levelmax=11' 'boxlen=1.0' \
    'ivar_refine=0' 'ic_refmap_fixed=.true.' \
    'm_refine=100*-1.' \
    'ngridmax=4500000' 'ngridtot=4500000' 'ngridmax_auto=.false.' \
    'npartmax=32000000' 'nparttot=28000000' \
    'nstepmax=0' 'noutput=1' 'aout=1.1' 'tout=1.0d100' \
    'foutput=1000000000' 'fbackup=1000000000' \
    'abort_on_mg_nonconvergence=.true.' 'omega_b=0.04930169880390167'; do
    grep -Fq "$key" "$root/effective.nml"
done
for index in 1 2 3 4; do
    grep -Fq "initfile($index)=" "$root/effective.nml"
done
test "$(grep -Fc 'ic_refmap_fixed=.true.' "$root/effective.nml")" -eq 1

for level in 8 9 10 11; do
    for field in ic_refmap ic_deltab ic_velbx ic_velby ic_velbz \
        ic_poscx ic_poscy ic_poscz ic_velcx ic_velcy ic_velcz ic_pvar_00001; do
        test -s "$root/ic$level/$field"
    done
done

printf '%s\n' "$root/effective.nml" > "$root/effective_namelist_path.txt"
printf '%s\n' "$ic" > "$root/input_ic_path.txt"
printf 'purpose=four-rank unit-box level-8-through-11 RAMSES IC reader and one-coarse-step preflight\nsource_revision=ca53ce39f455000868f0fbc25d810e45f6766a00\nlevels=8,9,10,11\ncode_boxlen=1.0; physical_box_length_from_grafic_header\nrefinement=ic_refmap_fixed=true; ivar_refine=0; no density-derived refinement\ninitial_grid_capacity_per_rank=4500000; prior_one_rank_final_global_grids=3249981\nparticle_capacity_per_rank=32000000; prior_one_rank_particles=23509074\nmemory_request_bytes=137438953472\nprior_level11_one_rank_maxrss=13336244K\nnstepmax=0; initialization may perform one coarse step\nnoutput=1; aout=1.1; tout=1.0d100; foutput=1000000000; fbackup=1000000000\nexpected_RAMSES_dumps=0\naccidental_single_dump_reserve_bytes=103079215104\nminimum_free_bytes_before=%s\nfree_bytes_before=%s\n' "$minimum_free_bytes" "$free_bytes" > "$root/output_policy.txt"
sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/template.nml.in" "$root/job.sbatch" \
    "$root/effective_namelist_path.txt" "$root/input_ic_path.txt" \
    "$root/output_policy.txt" > "$root/inputs.sha256"

export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_AMR_MAP_DIAGNOSTICS=1
set +e
(cd "$root" && srun --mpi=pmi2 --ntasks=4 --cpu-bind=cores \
    ./ramses_final3d effective.nml > run.log 2>&1)
run_rc=$?
set -e
printf '%s\n' "$run_rc" > "$root/srun_exit_code.txt"
grep -En '\[RAR_AMR_MAP\]|Level[[:space:]]+(8|9|10|11) has|npart=|FFT direct solve DONE|WARN:' \
    "$root/run.log" > "$root/key_diagnostics.txt" || true
if [ "$run_rc" -ne 0 ]; then
    printf '%s\n' 'L11_READER_PREFLIGHT_RUN_FAILED' > "$root/outcome.txt"
    exit 1
fi
for marker in 'Run completed' 'IC/namelist omega_b consistent' \
    'FFTW3 CPU direct Poisson solver enabled' 'FFT direct solve DONE' \
    'IC refinement map fixed'; do
    grep -Fq "$marker" "$root/run.log"
done
for level in 8 9 10 11; do
    grep -Eq "Reading file ic$level/ic_pvar_00001" "$root/run.log"
    grep -Eq "Level[[:space:]]+$level has[[:space:]]+[1-9][0-9]* grids" "$root/run.log"
done
if grep -Fq '[RAR_AMR_MAP] source=density' "$root/run.log"; then
    printf '%s\n' 'unexpected density-derived refinement map' >&2
    printf '%s\n' 'L11_READER_PREFLIGHT_DENSITY_MAP_REAPPEARED' > "$root/outcome.txt"
    exit 1
fi
if grep -Eiq 'ERROR:|FATAL:|Aborting|GRID_CAPACITY_FIXED|GRID_GROW|WARN: Fine multigrid Poisson failed to converge|NaN/Inf detected|NaN detected|Floating point exception|segmentation fault' "$root/run.log"; then
    printf '%s\n' 'L11_READER_PREFLIGHT_PHYSICS_OR_RUNTIME_WARNING' > "$root/outcome.txt"
    exit 1
fi
if compgen -G "$root/output_[0-9]*" >/dev/null || \
   compgen -G "$root/backup_[0-9]*" >/dev/null; then
    printf '%s\n' 'L11_READER_PREFLIGHT_UNEXPECTED_OUTPUT' > "$root/outcome.txt"
    exit 1
fi
python3 - "$root/run.log" "$root/particle_mesh_summary.txt" <<'PYVERIFY'
import re
import sys
from pathlib import Path

log = Path(sys.argv[1]).read_text(errors="replace")
particle_re = re.compile(r"npart=\s*(\d+)\s*/\s*(\d+)\s*for PE=\s*(\d+)")
particles = {}
for count, capacity, rank in particle_re.findall(log):
    rank = int(rank)
    if rank in particles:
        raise SystemExit(f"duplicate particle count for rank {rank}")
    particles[rank] = (int(count), int(capacity))
if set(particles) != {1, 2, 3, 4}:
    raise SystemExit(f"expected particle counts from ranks 1-4; found {sorted(particles)}")
if any(n <= 0 for n, cap in particles.values()):
    raise SystemExit(f'one or more MPI ranks received no particles: {particles}')
if any(n >= cap for n, cap in particles.values()):
    raise SystemExit("particle capacity exhausted on a rank")
total = sum(n for n, cap in particles.values())
if total != 23509074:
    raise SystemExit(f"global particle count {total} differs from prior reader 23509074")
meshes = {}
for level, count in re.findall(r"Level\s+(8|9|10|11) has\s+(\d+) grids", log):
    meshes[int(level)] = int(count)
if set(meshes) != {8, 9, 10, 11} or any(v <= 0 for v in meshes.values()):
    raise SystemExit(f"missing or empty mesh level: {meshes}")
if meshes[8] != 2097152 or sum(meshes.values()) > 4500000:
    raise SystemExit(f"unexpected global mesh extent: {meshes}")
nan_lines = [line for line in log.splitlines() if "NaN_CHK" in line]
if not nan_lines:
    raise SystemExit("missing NaN_CHK diagnostics")
for line in nan_lines:
    values = re.findall(r"=\s*([+-]?\d+)", line)
    if len(values) != 3 or any(int(v) != 0 for v in values):
        raise SystemExit(f"nonzero or malformed NaN diagnostic: {line}")
lines = [f"rank_{rank}_particles={particles[rank][0]}" for rank in sorted(particles)]
lines += [f"global_particles={total}", f"nan_check_lines={len(nan_lines)}"]
for level in sorted(meshes):
    lines.append(f"final_level_{level}_grids={meshes[level]}")
lines.append("single_rank_same_unitbox_reference=not_yet_run")
Path(sys.argv[2]).write_text("\n".join(lines) + "\n")
PYVERIFY
printf '%s\n' 'L11_UNITBOX_MPI_NUMERICAL_CHECKS_PASSED_REVIEW_REQUIRED' > "$root/outcome.txt"
sha256sum "$root/effective.nml" "$root/ramses_final3d" \
    "$root/output_policy.txt" "$root/run.log" "$root/job.sbatch" \
    "$root/effective_namelist_path.txt" "$root/input_ic_path.txt" \
    "$root/key_diagnostics.txt" "$root/particle_mesh_summary.txt" "$root/outcome.txt" \
    "$root/srun_exit_code.txt" > "$root/artifacts.sha256"
sha256sum -c "$root/artifacts.sha256" --status
printf '%s\n' 'L11_UNITBOX_MPI_READER_AWAITING_REVIEW' > "$root/PREFLIGHT_AWAITING_REVIEW"
printf 'Level-11 four-rank unit-box reader preflight captured; physical review required. Run directory: %s\n' "$root"
