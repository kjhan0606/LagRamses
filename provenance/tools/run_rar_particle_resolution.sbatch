#!/bin/bash
#SBATCH --job-name=rar-particle-resolution
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_particle_resolution_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_particle_resolution_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
build=$base/builds/rar_effective_ssp_${RAR_EVOLVE_BUILD_ID:?}
seed=$base/runs/rar_particle_probe_1279276
out=$base/runs/rar_particle_resolution_${SLURM_JOB_ID:?}
test -f "$build/COMPLETE"
test -f "$seed/COMPLETE"
(cd "$build/source/bin" && sha256sum -c "$build/executable.sha256")
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 8589934592
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$project/scripts/audit_rar_particle_resolution.py" "$out/audit.py"
cp "$project/scripts/rar_native_particles.py" "$out/rar_native_particles.py"
cp "$project/scripts/verify_rar_composite_publication.py" "$out/verify.py"
cp "$seed/result.json" "$out/original_self_pair.json"
IFS=',' read -r -a levels <<< "${RAR_PARTICLE_LEVELS:-4,5}"
for lev in "${levels[@]}";do [[ "$lev" =~ ^[456]$ ]];done
refinement=${RAR_PARTICLE_REFINEMENT:-sphere}
[[ "$refinement" = sphere || "$refinement" = fixed_box ]]
printf '%s\n' 'purpose=single-cloud initial force measurement at fixed physical width and changing native mesh resolution' 'one explicit evolution step per case; no snapshots; source mass=0.002; position=0.6873,0.51,0.52; cloud_width=0.0625; uniform gas density=1; periodic noncosmological CPU; cooling/feedback/dust/RT off' 'nstepmax=1;dt_cap=1e-6;noutput=1;aout=1.1;tout=1e100;foutput=1e9;fbackup=1e9;expected_full_dumps=0;reserved_logs_and_field_exports_bytes=4294967296' "base_levels=${levels[*]}" "refinement=$refinement" 'sphere: diameter=0.25, exponent=2, measure rasterized boundary; fixed_box: exponent=10, requested diameter=0.375-2*base_dx, require actual fine region [0.3125,0.6875]^3 after flag smoothing' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
set +e
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EVOLVE=1
export RAR_QUMOND_COMPOSITE_EXPORT=1 RAR_QUMOND_EVOLVE_DT=0.000001 RAR_NATIVE_PARTICLE_ROLE=diagnostic
unset RAR_QUMOND_GRAVITY_KICKS RAR_COOLING_COMMIT_LEDGER RAR_GLOBAL_INVENTORY
unset RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM
for lev in "${levels[@]}";do
 for gravity in newtonian qumond;do
  d=$out/L${lev}_$gravity
  mkdir "$d"
  cp "$build/source/bin/ramses_final3d" "$d/ramses_final3d"
  cp "$build/source_commit.txt" "$d/source_commit.txt"
  python3 - "$seed/newtonian_a/effective.nml" "$d/effective.nml" "$gravity" "$lev" "$refinement" <<'PY'
import re,sys
from pathlib import Path
s=Path(sys.argv[1]).read_text();lev=int(sys.argv[4]);assert lev in (4,5,6)
changes={'use_mond':'.true.' if sys.argv[3]=='qumond' else '.false.','levelmin':str(lev),'levelmax':str(lev+1),
 'void_refine_min_level':str(lev+1),'ngridtot':'200000' if lev==6 else '50000'}
if sys.argv[5]=='fixed_box':
 changes.update(r_refine=f'{.375-2*2.**(-lev):.17g}',exp_refine='10.0')
else:
 assert sys.argv[5]=='sphere'
for key,value in changes.items():
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1
for key,value in {'nstepmax':'1','noutput':'1','aout':'1.1','tout':'1.0d100','foutput':'1000000000','fbackup':'1000000000','pic':'.true.','cosmo':'.false.','hydro':'.true.','poisson':'.true.','sink':'.false.','gravity_type':'0','d_region':'1.0','boxlen':'1.0'}.items():
 vals=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s);assert len(vals)==1 and vals[0].strip().lower()==value
assert not re.search(r'(?im)^\s*(cooling|star)\s*=\s*\.true\.',s)
target=Path(sys.argv[2]);target.write_text(s);ic=target.parent/'ic';ic.mkdir()
(ic/'ic_part').write_text(' '.join(f'{v:.17g}' for v in (.6873-.5,.51-.5,.52-.5,0.,0.,0.,.002))+'\n')
PY
  if [ "$lev" = 4 ];then unset RAR_NATIVE_PARTICLE_CLOUD_WIDTH;else export RAR_NATIVE_PARTICLE_CLOUD_WIDTH=0.0625;fi
  printf '%s\n' "effective_namelist=$d/effective.nml" "cloud_width_setting=${RAR_NATIVE_PARTICLE_CLOUD_WIDTH:-base_cell_default}" >> "$out/effective_input_paths.txt"
  (cd "$d" && sha256sum effective.nml ic/ic_part ramses_final3d source_commit.txt >inputs.sha256 && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1)
  grep -Fq 'Run completed' "$d/run.log"
  grep -Fq 'RAR_NATIVE_PARTICLE_END_V1' "$d/run.log"
  test -f "$d/composite_epoch_00001_leaf.tsv"
  test ! -e "$d/composite_epoch_00002_leaf.tsv"
  test "$(find "$d" -maxdepth 1 -type d -name 'output_*' | wc -l)" -eq 0
  if grep -Eiq 'failed to converge|FATAL:|ERROR:|Native composite particles failed|QUMOND composite publication failed' "$d/run.log";then exit 1;fi
 done
done
/home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$out"
(cd "$out" && sha256sum result.json verify_*.json audit.py verify.py rar_native_particles.py job.sbatch */composite*.tsv */particles*.tsv */run.log >artifacts.sha256)
printf '%s\n' 'FIXED_CLOUD_SPATIAL_MEASUREMENTS_COMPLETE_NOT_PRODUCTION' > "$out/COMPLETE"
