#!/bin/bash
#SBATCH --job-name=rar-particle-self-pair
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_particle_probe_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_particle_probe_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
build=$base/builds/rar_effective_ssp_${RAR_EVOLVE_BUILD_ID:?}
seed=$base/runs/rar_evolving_gas_${RAR_PARTICLE_RUN_ID:?}
out=$base/runs/rar_particle_probe_${SLURM_JOB_ID:?}
test -f "$build/COMPLETE"
test -f "$seed/COMPLETE"
(cd "$build/source/bin" && sha256sum -c "$build/executable.sha256")
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 1073741824
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$project/scripts/audit_rar_particle_probe.py" "$out/audit.py"
cp "$project/scripts/rar_native_particles.py" "$out/rar_native_particles.py"
cp "$project/scripts/verify_rar_composite_publication.py" "$out/verify.py"
printf '%s\n' 'purpose=initial-field self/pair-force measurement with one explicit evolution step; only time-zero force is analyzed' 'six cases; uniform gas plus one/two synthetic particles; native fixed two-level CPU gravity; no cooling, feedback, enrichment, dust or RT' 'nstepmax=1;dt_cap=0.0002;noutput=1;aout=1.1;tout=1e100;foutput=1e9;fbackup=1e9;expected_full_dumps=0;reserve_for_logs_and_field_exports_bytes=536870912' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
set +e
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EVOLVE=1
export RAR_QUMOND_COMPOSITE_EXPORT=1 RAR_QUMOND_EVOLVE_DT=0.0002 RAR_NATIVE_PARTICLE_ROLE=diagnostic
unset RAR_QUMOND_GRAVITY_KICKS RAR_COOLING_COMMIT_LEDGER RAR_GLOBAL_INVENTORY
unset RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM
for gravity in newtonian qumond;do
 for particles in a b pair;do
  d=$out/${gravity}_${particles}
  mkdir "$d"
  cp "$build/source/bin/ramses_final3d" "$d/ramses_final3d"
  cp "$build/source_commit.txt" "$d/source_commit.txt"
  python3 - "$seed/newtonian/effective.nml" "$d/effective.nml" "$gravity" "$particles" <<'PY'
import re,sys
from pathlib import Path
s=Path(sys.argv[1]).read_text()
changes={'hydro_initial_condition':"'region'",'use_mond':'.true.' if sys.argv[3]=='qumond' else '.false.',
 'nstepmax':'1','foutput':'1000000000','u_region':'0.0','v_region':'0.0','w_region':'0.0'}
for key,value in changes.items():
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1
for key,value in {'noutput':'1','aout':'1.1','tout':'1.0d100','fbackup':'1000000000','pic':'.true.','cosmo':'.false.','hydro':'.true.','poisson':'.true.','sink':'.false.','gravity_type':'0','d_region':'1.0'}.items():
 vals=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s);assert len(vals)==1 and vals[0].strip().lower()==value
assert not re.search(r'(?im)^\s*(cooling|star)\s*=\s*\.true\.',s)
target=Path(sys.argv[2]);target.write_text(s);ic=target.parent/'ic';ic.mkdir()
positions={'a':[(.6873,.51,.52)],'b':[(.3127,.49,.48)],'pair':[(.6873,.51,.52),(.3127,.49,.48)]}[sys.argv[4]]
(ic/'ic_part').write_text(''.join(' '.join(f'{v:.17g}' for v in (*[x-.5 for x in p],0.,0.,0.,.002))+'\n' for p in positions))
PY
  printf '%s\n' "effective_namelist=$d/effective.nml" >> "$out/effective_input_paths.txt"
  (cd "$d" && sha256sum effective.nml ic/ic_part ramses_final3d source_commit.txt >inputs.sha256 && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1)
  grep -Fq 'Run completed' "$d/run.log"
  test -f "$d/particles_epoch_00000.tsv"
  test -f "$d/composite_epoch_00001_leaf.tsv"
  test ! -e "$d/composite_epoch_00002_leaf.tsv"
  test "$(find "$d" -maxdepth 1 -type d -name 'output_*' | wc -l)" -eq 0
  if grep -Eiq 'failed to converge|FATAL:|ERROR:|Native composite particles failed|QUMOND composite publication failed' "$d/run.log";then exit 1;fi
 done
done
/home/kjhan/miniconda3/bin/python3 "$out/audit.py" "$out"
(cd "$out" && sha256sum result.json verify_*.json audit.py verify.py rar_native_particles.py job.sbatch */composite*.tsv */particles*.tsv */run.log >artifacts.sha256)
printf '%s\n' 'INITIAL_NATIVE_SELF_PAIR_FORCES_MEASURED_NOT_A_CALIBRATION' > "$out/COMPLETE"
