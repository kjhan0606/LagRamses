#!/bin/bash
#SBATCH --job-name=rar-evolving-gas
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_evolving_gas_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_evolving_gas_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
build=$base/builds/rar_effective_ssp_${RAR_EVOLVE_BUILD_ID:?}
out=$base/runs/rar_evolving_gas_${SLURM_JOB_ID:?}
test -f "$build/COMPLETE"
(cd "$build/source/bin" && sha256sum -c "$build/executable.sha256")
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 1610612736
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' 'purpose=bounded evolving gas, native AMR common-time gravity; no particles/subgrid/dust/RT' 'cases=newtonian,qumond,qumond_half; steps=6,6,12; dt_caps=0.0002,0.0002,0.0001; expected_common_endpoint=0.0012; noutput=1;aout=1.1;tout=1e100;foutput=3,3,6;fbackup=1e9;expected_full_dumps_per_case=2;upper_dump_bytes=16777216;expected_total_dump_bytes=100663296;reserve_including_exports_and_kick_records=536870912' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EVOLVE=1 RAR_QUMOND_COMPOSITE_EXPORT=1 RAR_QUMOND_GRAVITY_KICKS=1
unset RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM
for mode in newtonian qumond qumond_half;do
 d=$out/$mode
 mkdir "$d"
 cp "$build/source/bin/ramses_final3d" "$d/ramses_final3d"
 cp "$build/source_commit.txt" "$d/source_commit.txt"
 python3 - "$project/configs/qumond_amr_gate/qumond.nml" "$d/effective.nml" "$mode" <<'PY'
import re,sys
from pathlib import Path
s=Path(sys.argv[1]).read_text();mode=sys.argv[3]
changes={'static':'.false.','use_mond':'.false.' if mode=='newtonian' else '.true.','levelmin':'4','levelmax':'5','void_refine_min_level':'5','r_refine':'0.25','ngridtot':'30000','nstepmax':'12' if mode=='qumond_half' else '6','nsubcycle':'1,1','use_fftw':'.true.','gpu_fft':'.false.','gravity_type':'0','d_region':'1.0','p_region':'1.0','courant_factor':'0.1','foutput':'6' if mode=='qumond_half' else '3'}
for key,value in changes.items():
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1
s=s.replace('&RUN_PARAMS','&RUN_PARAMS\n remap_thresh=0.0')
s=s.replace('&INIT_PARAMS',"&INIT_PARAMS\n hydro_initial_condition='rar_periodic_wave'\n u_region=0.2\n v_region=0.1\n w_region=0.05")
for key,value in {'cosmo':'.false.','pic':'.false.','hydro':'.true.','poisson':'.true.','sink':'.false.','noutput':'1','aout':'1.1','tout':'1.0d100','fbackup':'1000000000','nremap':'0','remap_thresh':'0.0'}.items():
 vals=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s);assert len(vals)==1 and vals[0].strip().lower()==value
Path(sys.argv[2]).write_text(s)
PY
 if [ "$mode" = qumond_half ];then export RAR_QUMOND_EVOLVE_DT=0.0001;else export RAR_QUMOND_EVOLVE_DT=0.0002;fi
 printf '%s\n' "effective_namelist=$d/effective.nml" "dt_cap=$RAR_QUMOND_EVOLVE_DT" > "$d/effective_input_path.txt"
 (cd "$d" && sha256sum effective.nml ramses_final3d source_commit.txt >inputs.sha256 && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1)
 grep -Fq 'Run completed' "$d/run.log"
 grep -Fq 'RAR_QUMOND_EVOLVE_END_V1' "$d/run.log"
 if grep -Eiq 'failed to converge|NaN/Inf detected|segmentation fault|FATAL:|ERROR:|QUMOND stage history:|QUMOND composite publication|QUMOND evolving gas' "$d/run.log";then exit 1;fi
 test "$(find "$d" -name 'data_*.h5' | wc -l)" -eq 2
 (cd "$d" && sha256sum run.log composite_*.tsv gravity_kicks_*.tsv output_*/data_*.h5 >outputs.sha256)
done
printf '%s\n' 'EVOLVING_GAS_EXECUTION_COMPLETE_PHYSICS_ANALYSIS_PENDING' > "$out/COMPLETE"
