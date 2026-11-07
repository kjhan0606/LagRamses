#!/bin/bash
#SBATCH --job-name=rar-stage-history
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_stage_history_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_stage_history_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
new=$base/builds/rar_effective_ssp_${RAR_HISTORY_BUILD_ID:-1275333}
old=$base/builds/rar_effective_ssp_1275327
out=$base/runs/rar_stage_history_${SLURM_JOB_ID:?}
test -f "$new/COMPLETE";test -f "$old/COMPLETE"
(cd "$new/source/bin" && sha256sum -c "$new/executable.sha256")
(cd "$old/source/bin" && sha256sum -c "$old/executable.sha256")
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 3221225472
printf '%s\n' 'purpose=fixed two-level static history integration and default-off regression' 'cases=3; expected_dumps_per_case=2; measured_reference_dump_bytes=91106264; estimated_total_dump_bytes=546637584; reserve_bytes=2147483648' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
set +u
source "$(cat "$new/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
unset RAR_QUMOND_COMPOSITE_EVOLVE RAR_QUMOND_GRAVITY_KICKS RAR_QUMOND_EVOLVE_DT
for mode in old_off new_off new_on; do
 case_dir=$out/$mode
 mkdir "$case_dir"
 if [ "$mode" = old_off ];then build=$old;else build=$new;fi
 cp "$build/source/bin/ramses_final3d" "$case_dir/ramses_final3d"
 cp "$build/source_commit.txt" "$case_dir/source_commit.txt"
 python3 - "$project/configs/qumond_amr_gate/qumond.nml" "$case_dir/effective.nml" <<'PY'
import os,re,sys
from pathlib import Path
s=Path(sys.argv[1]).read_text()
lev=int(os.environ.get('RAR_HISTORY_LEVELMIN','6'))
assert lev>=2
for name,value in [('levelmin',str(lev)),('levelmax',str(lev+1)),('void_refine_min_level',str(lev+1)),('use_fftw','.true.'),('gpu_fft','.false.'),('nstepmax','2'),('nsubcycle','1,1')]:
 s,n=re.subn(r'(?im)^\s*'+name+r'\s*=.*$',f' {name}={value}',s);assert n==1
for name,value in [('noutput','1'),('aout','1.1'),('foutput','1'),('fbackup','1000000000'),('nsubcycle','1,1'),('nremap','0')]:
 assert re.search(r'(?im)^\s*'+name+r'\s*=\s*'+re.escape(value)+r'\s*$',s)
Path(sys.argv[2]).write_text(s)
PY
 printf '%s\n' "effective_namelist=$case_dir/effective.nml; noutput=1,aout=1.1,tout=1e100,foutput=1,fbackup=1e9,nstepmax=2,nsubcycle=1,1; expected_full_dumps=2" > "$case_dir/output_policy.txt"
 (cd "$case_dir" && sha256sum effective.nml ramses_final3d source_commit.txt > inputs.sha256)
 if [ "$mode" = new_on ];then export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_OPERATOR=${RAR_COMPOSITE_OPERATOR:-0} RAR_QUMOND_COMPOSITE_SOLVER=${RAR_COMPOSITE_SOLVER:-0} RAR_QUMOND_COMPOSITE_PUBLISH=${RAR_COMPOSITE_PUBLISH:-0} RAR_QUMOND_COMPOSITE_EXPORT=${RAR_COMPOSITE_EXPORT:-0} RAR_QUMOND_COMPOSITE_CONTINUUM=${RAR_COMPOSITE_CONTINUUM:-0};else unset RAR_QUMOND_STAGE_HISTORY RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_PUBLISH RAR_QUMOND_COMPOSITE_EXPORT RAR_QUMOND_COMPOSITE_CONTINUUM;fi
 (cd "$case_dir" && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml > run.log 2>&1)
 grep -Fq 'Run completed' "$case_dir/run.log"
 if grep -Eiq 'failed to converge|NaN/Inf detected|segmentation fault|FATAL:|QUMOND stage history:|QUMOND composite operator:|QUMOND composite solver:|QUMOND composite publication' "$case_dir/run.log";then exit 1;fi
 count=$(find "$case_dir" -name 'data_*.h5' | wc -l)
 test "$count" -eq 2
 if [ "$mode" = new_on ];then
  grep -Fq 'RAR_QUMOND_HISTORY_V1' "$case_dir/run.log"
  if [ "${RAR_COMPOSITE_OPERATOR:-0}" = 1 ];then grep -Fq 'RAR_QUMOND_COMPOSITE_V1' "$case_dir/run.log";fi
  if [ "${RAR_COMPOSITE_SOLVER:-0}" = 1 ];then grep -Fq 'RAR_QUMOND_SOLVER_REJECTION_V1' "$case_dir/run.log";fi
  if [ "${RAR_COMPOSITE_PUBLISH:-0}" = 1 ];then grep -Fq 'RAR_QUMOND_PUBLICATION_V1' "$case_dir/run.log";fi
  if [ "${RAR_COMPOSITE_CONTINUUM:-0}" = 1 ];then grep -Fq 'RAR_QUMOND_CONTINUUM_V1' "$case_dir/run.log";fi
 fi
 (cd "$case_dir" && sha256sum run.log output_*/data_*.h5 > outputs.sha256)
done
printf '%s\n' 'EXECUTION_COMPLETE_HISTORY_AND_FORCE_EVALUATION_PENDING' > "$out/COMPLETE"
