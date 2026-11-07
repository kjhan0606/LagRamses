#!/bin/bash
#SBATCH --job-name=rar-native-disc
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_disc_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_disc_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
build=$base/builds/rar_effective_ssp_${RAR_DISC_BUILD_ID:?}
out=$base/runs/rar_native_disc_${SLURM_JOB_ID:?}
test -f "$build/COMPLETE"
(cd "$build/source/bin" && sha256sum -c "$build/executable.sha256")
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 3221225472
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' 'purpose=prescribed matched disc force comparison; not evolving galaxy formation' 'cases=newtonian_baryon,newtonian_halo,qumond_baryon; noutput=1;aout=1.1;tout=1e100;foutput=1;fbackup=1e9;nstepmax=1;expected_dumps_per_case=1;reference_dump_bytes=91106264;estimated_snapshot_bytes=273318792;reserve_including_exports_bytes=2147483648' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EXPORT=1
unset RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM
for mode in newtonian_baryon newtonian_halo qumond_baryon; do
 d=$out/$mode
 mkdir "$d"
 cp "$build/source/bin/ramses_final3d" "$d/ramses_final3d"
 cp "$build/source_commit.txt" "$d/source_commit.txt"
 python3 - "$project/configs/qumond_disc_pilot/$mode.nml" "$project/configs/qumond_amr_gate/qumond.nml" "$d/effective.nml" "$mode" <<'PY'
import os,re,sys
from pathlib import Path
source,ref,out=map(Path,sys.argv[1:4]);mode=sys.argv[4];s=source.read_text();lev=int(os.environ['RAR_DISC_LEVELMIN']);assert lev in (5,6)
for key,value in [('levelmin',str(lev)),('levelmax',str(lev+1)),('ngridtot','400000'),('use_fftw','.true.'),('gpu_fft','.false.'),('nsubcycle','1,1')]:
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1
refine=re.search(r'(?is)&REFINE_PARAMS\b.*?/',ref.read_text())[0]
refine=re.sub(r'(?im)^\s*void_refine_min_level\s*=.*$',f' void_refine_min_level={lev+1}',refine)
s+='\n'+refine+'\n'
expected={'cosmo':'.false.','pic':'.false.','hydro':'.true.','static':'.true.','sink':'.false.','poisson':'.true.','nstepmax':'1','noutput':'1','aout':'1.1','tout':'1.0d100','foutput':'1','fbackup':'1000000000','nsubcycle':'1,1','gravity_type':'-3' if mode=='newtonian_halo' else '-2','use_mond':'.true.' if mode=='qumond_baryon' else '.false.'}
for key,value in expected.items():
 matches=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s);assert len(matches)==1 and matches[0].strip().lower()==value,(key,matches,value)
params=[float(v) for v in re.search(r'(?im)^\s*gravity_params\s*=([^\n]+)',s)[1].split(',')]
assert params[:6]==[.06,.055,.025,.02,.085,.035]
assert params[6:]==([.4,.04,.35] if mode=='newtonian_halo' else [0.,.04,.35])
Path(out).write_text(s)
PY
 printf '%s\n' "effective_namelist=$d/effective.nml" > "$d/effective_namelist_path.txt"
 (cd "$d" && sha256sum effective.nml ramses_final3d source_commit.txt > inputs.sha256 && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1)
 grep -Fq 'Run completed' "$d/run.log"
 grep -Fq 'RAR_COMPOSITE_MODEL_V1' "$d/run.log"
 grep -Fq 'RAR_QUMOND_PUBLICATION_V1' "$d/run.log"
 if grep -Eiq 'failed to converge|NaN/Inf detected|segmentation fault|FATAL:|ERROR:|QUMOND stage history:|QUMOND composite publication' "$d/run.log";then exit 1;fi
 test "$(find "$d" -name 'data_*.h5' | wc -l)" -eq 1
 (cd "$d" && sha256sum run.log composite_leaf_state.tsv composite_faces.tsv output_*/data_*.h5 > outputs.sha256)
done
printf '%s\n' 'MATCHED_NATIVE_DISC_FORCES_EXECUTED_ANALYSIS_PENDING' > "$out/COMPLETE"
