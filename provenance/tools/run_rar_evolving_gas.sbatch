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
if [ "${RAR_EVOLVE_COOLING:-0}" = 1 ];then test "$free_bytes" -gt 2147483648;else test "$free_bytes" -gt 1610612736;fi
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
printf '%s\n' 'purpose=bounded evolving gas, native AMR common-time gravity; optional controlled diagnostic particles; no subgrid/dust/RT' 'cases=newtonian,qumond,qumond_half; steps=6,6,12; dt_caps=0.0002,0.0002,0.0001; expected_common_endpoint=0.0012; noutput=1;aout=1.1;tout=1e100;foutput=3,3,6;fbackup=1e9;expected_full_dumps_per_case=2;upper_dump_bytes=16777216;expected_total_dump_bytes=100663296;reserve_including_exports_and_kick_records=536870912' "free_bytes_before=$free_bytes" > "$out/output_policy.txt"
if [ "${RAR_EVOLVE_PARTICLES:-0}" = 1 ];then
 export RAR_NATIVE_PARTICLE_ROLE=diagnostic
 printf '%s\n' 'particles=4; fixed masses=0.002,0.001,0.002,0.001; fixed-cloud width=0.0625; identical mass for Newtonian/QUMOND numerical controls; not Lambda-CDM versus baryon-only QUMOND science' >> "$out/output_policy.txt"
else
 unset RAR_NATIVE_PARTICLE_ROLE
fi
if [ "${RAR_EVOLVE_COOLING:-0}" = 1 ];then
 printf '%s\n' 'thermal_profile_overrides=5 cases plus optional precision replay; newtonian/qumond/qumond_half cooling on, qumond_off cooling off, qumond_heating cool/heating on at lower initial pressure; at most 12 snapshots; upper_total_dump_bytes=201326592; reserve_all_outputs_bytes=1073741824' 'cooling=original; UV=code haardt_madau flag/default Courty model at fixed aexp1/z0; metal=true initial Z0.002; shielding=false; Jeans/eEOS/isothermal/delayed/neq/RT/dust off; T2max=1e10; source table T2_min_fix=1e-2' "cooling_varmax=${RAR_EVOLVE_COOLING_VARMAX:-4}; precision_replay=${RAR_EVOLVE_THERMAL_ACCURACY:-0}" >> "$out/output_policy.txt"
fi
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EVOLVE=1 RAR_QUMOND_COMPOSITE_EXPORT=1 RAR_QUMOND_GRAVITY_KICKS=1
unset RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM
modes='newtonian qumond qumond_half'
if [ "${RAR_EVOLVE_COOLING:-0}" = 1 ];then
 modes="$modes qumond_off qumond_heating"
 if [ "${RAR_EVOLVE_THERMAL_ACCURACY:-0}" = 1 ];then modes="$modes qumond_accuracy";fi
 export RAR_COOLING_COMMIT_LEDGER=1 RAR_GLOBAL_INVENTORY=1
else
 unset RAR_COOLING_COMMIT_LEDGER RAR_GLOBAL_INVENTORY
fi
for mode in $modes;do
 d=$out/$mode
 mkdir "$d"
 cp "$build/source/bin/ramses_final3d" "$d/ramses_final3d"
 cp "$build/source_commit.txt" "$d/source_commit.txt"
 python3 - "$project/configs/qumond_amr_gate/qumond.nml" "$d/effective.nml" "$mode" <<'PY'
import os,re,sys
from pathlib import Path
s=Path(sys.argv[1]).read_text();mode=sys.argv[3]
changes={'static':'.false.','use_mond':'.false.' if mode=='newtonian' else '.true.','levelmin':'4','levelmax':'5','void_refine_min_level':'5','r_refine':'0.25','ngridtot':'30000','nstepmax':'12' if mode=='qumond_half' else '6','nsubcycle':'1,1','use_fftw':'.true.','gpu_fft':'.false.','gravity_type':'0','d_region':'1.0','p_region':'1.0','courant_factor':'0.1','foutput':'6' if mode=='qumond_half' else '3'}
particles=os.environ.get('RAR_EVOLVE_PARTICLES','0')=='1'
if particles:changes.update(pic='.true.',nparttot='1000')
for key,value in changes.items():
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1
s=s.replace('&RUN_PARAMS','&RUN_PARAMS\n remap_thresh=0.0')
s=s.replace('&INIT_PARAMS',"&INIT_PARAMS\n hydro_initial_condition='rar_periodic_wave'\n u_region=0.2\n v_region=0.1\n w_region=0.05")
if particles:
 s=s.replace('&INIT_PARAMS',"&INIT_PARAMS\n filetype='ascii'\n initfile(1)='ic'")
 # RAMSES namelist initfile(1) names the first input level, mapped to levelmin.
 ic=Path(sys.argv[2]).parent/'ic';ic.mkdir()
 # Measured native interfaces at these transverse coordinates are x=.6875
 # and .3125; the refinement input radius is not the leaf-boundary radius.
 points=[(.6873,.51,.52,1.,0.,0.,.002),(.9998,.3,.3,1.,0.,0.,.001),(.3127,.49,.48,-1.,0.,0.,.002),(.0002,.7,.7,-1.,0.,0.,.001)]
 # The ASCII reader adds boxlen/2 to each input position.
 (ic/'ic_part').write_text(''.join(' '.join(f'{v:.17g}' for v in (*[v-.5 for v in row[:3]],*row[3:]))+'\n' for row in points))
if os.environ.get('RAR_EVOLVE_COOLING','0')=='1':
 thermal='.false.' if mode=='qumond_off' else '.true.'
 s=s.replace('&PHYSICS_PARAMS',f"&PHYSICS_PARAMS\n cooling={thermal}\n metal=.true.\n haardt_madau=.true.\n J21=0.0\n z_ave=0.1\n n_star=0.1\n T2_star=0.0\n jeans_ncells=0.0\n isothermal=.false.\n delayed_cooling=.false.\n neq_chem=.false.\n self_shielding=.false.\n cooling_method='original'\n T2max=1.0d10")
 s=s.replace('&INIT_PARAMS','&INIT_PARAMS\n var_region=0.002')
 if 'RAR_EVOLVE_COOLING_VARMAX' in os.environ:
  accuracy=float(os.environ['RAR_EVOLVE_COOLING_VARMAX']);assert accuracy>=4
  if mode=='qumond_accuracy':accuracy*=2
  s=s.replace('&PHYSICS_PARAMS',f'&PHYSICS_PARAMS\n cooling_varmax={accuracy:.17g}')
 if mode=='qumond_heating':
  s,n=re.subn(r'(?im)^\s*p_region\s*=.*$',' p_region=0.001',s);assert n==1
for key,value in {'cosmo':'.false.','pic':'.true.' if particles else '.false.','hydro':'.true.','poisson':'.true.','sink':'.false.','noutput':'1','aout':'1.1','tout':'1.0d100','fbackup':'1000000000','nremap':'0','remap_thresh':'0.0'}.items():
 vals=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s);assert len(vals)==1 and vals[0].strip().lower()==value
Path(sys.argv[2]).write_text(s)
PY
 if [ "$mode" = qumond_half ];then export RAR_QUMOND_EVOLVE_DT=0.0001;else export RAR_QUMOND_EVOLVE_DT=0.0002;fi
 printf '%s\n' "effective_namelist=$d/effective.nml" "dt_cap=$RAR_QUMOND_EVOLVE_DT" > "$d/effective_input_path.txt"
 if [ "${RAR_EVOLVE_PARTICLES:-0}" = 1 ];then
  (cd "$d" && sha256sum ic/ic_part >particle_input.sha256)
 fi
 (cd "$d" && sha256sum effective.nml ramses_final3d source_commit.txt >inputs.sha256 && srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1)
 grep -Fq 'Run completed' "$d/run.log"
 grep -Fq 'RAR_QUMOND_EVOLVE_END_V1' "$d/run.log"
 if grep -Eiq 'failed to converge|NaN/Inf detected|segmentation fault|FATAL:|ERROR:|QUMOND stage history:|QUMOND composite publication|QUMOND evolving gas|Native composite particles failed' "$d/run.log";then exit 1;fi
 test "$(find "$d" -name 'data_*.h5' | wc -l)" -eq 2
 (cd "$d" && sha256sum run.log composite_*.tsv gravity_kicks_*.tsv output_*/data_*.h5 >outputs.sha256)
 if [ "${RAR_EVOLVE_PARTICLES:-0}" = 1 ];then
  grep -Fq 'RAR_NATIVE_PARTICLE_END_V1' "$d/run.log"
  (cd "$d" && sha256sum particle_step_*.tsv particles_epoch_*.tsv >>outputs.sha256)
 fi
done
printf '%s\n' 'EVOLVING_GAS_EXECUTION_COMPLETE_PHYSICS_ANALYSIS_PENDING' > "$out/COMPLETE"
