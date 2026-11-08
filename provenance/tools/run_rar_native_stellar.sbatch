#!/bin/bash
# Bounded synthetic depletion/source-coupling check, not a galaxy calibration.
#SBATCH --job-name=rar-native-stellar
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_stellar_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_native_stellar_%j.err
set -Eeuo pipefail
umask 077
base=/gpfs/kjhan/RAR_QUMOND
project=/home/kjhan/BackUP/rar
build=$base/builds/rar_effective_ssp_${RAR_NATIVE_STELLAR_BUILD_ID:?}
old=$base/runs/rar_q_ssp_hydro_1274983
out=$base/runs/rar_native_stellar_${SLURM_JOB_ID:?}
test -f "$build/COMPLETE"
(cd "$build/source/bin" && sha256sum -c "$build/executable.sha256")
free_bytes=$(df -B1 "$base" | awk 'END {print $4}')
test "$free_bytes" -gt 2147483648
mkdir "$out"
cp "${BASH_SOURCE[0]}" "$out/job.sbatch"
cp "$project/configs/qumond_amr_gate/qumond.nml" "$out/template.nml"
cp "$old/effective.nml" "$out/stellar_template.nml"
cp "$build/source/bin/ramses_final3d" "$out/ramses_final3d"
cp "$build/source_commit.txt" "$out/source_commit.txt"
for name in phase0_yields.dat legacy_yield_table_parser_only.asc history.nml snia_effective_ssp_runtime.nml;do
 cp -p "$old/$name" "$out/$name"
done
mkdir "$out/ic"
: > "$out/ic/ic_part"
python3 - "$out" "$old" <<'PY'
import re,sys
from pathlib import Path
out,old=map(Path,sys.argv[1:]);s=(out/'template.nml').read_text()
changes={'static':'.false.','pic':'.true.','levelmin':'3','levelmax':'4',
 'void_refine_min_level':'4','r_refine':'0.25','ngridtot':'30000','nparttot':'20000',
 'nstepmax':'4','nsubcycle':'1,1','use_fftw':'.true.','gpu_fft':'.false.',
 'gravity_type':'0','d_region':'1.0','p_region':'0.001','courant_factor':'0.1','foutput':'2'}
for key,value in changes.items():
 s,n=re.subn(r'(?im)^\s*'+key+r'\s*=.*$',f' {key}={value}',s);assert n==1,(key,n)
s=s.replace('&RUN_PARAMS','&RUN_PARAMS\n remap_thresh=0.0\n rt=.false.\n radiation_transport=\'none\'')
s=s.replace('&INIT_PARAMS',"&INIT_PARAMS\n filetype='ascii'\n initfile(1)='ic'\n hydro_initial_condition='rar_periodic_wave'\n u_region=0.2\n v_region=0.1\n w_region=0.05\n var_region=0.002")
s=s.replace('&REFINE_PARAMS','&REFINE_PARAMS\n mass_sph=1.0d-5\n m_refine=2*1.0d30')
s=s.replace('&PHYSICS_PARAMS',f'''&PHYSICS_PARAMS
 cooling=.true.
 metal=.true.
 haardt_madau=.true.
 J21=0.0
 z_ave=0.1
 n_star=0.43
 del_star=0.0
 t_star=0.0
 eps_star=1000.0
 m_star=1.0
 f_w=0.0
 sf_virial=.false.
 sf_birth_properties=.false.
 T2_star=0.0
 jeans_ncells=0.0
 isothermal=.false.
 delayed_cooling=.false.
 neq_chem=.false.
 self_shielding=.false.
 cooling_method='original'
 cooling_varmax=4.0
 T2max=1.0d10
 dust_mass_enabled=.false.
 yieldtablefilename='{out}/legacy_yield_table_parser_only.asc' ''')
stellar=re.search(r'(?ims)^&stellar_enrichment_params\b.*?^/',(out/'stellar_template.nml').read_text())[0]
stellar=stellar.replace(str(old),str(out))
s,n=re.subn(r'(?ims)^&STELLAR_ENRICHMENT_PARAMS\b.*?^/',lambda m:stellar,s);assert n==1
required={'cosmo':'.false.','pic':'.true.','hydro':'.true.','poisson':'.true.',
 'sink':'.false.','use_mond':'.true.','noutput':'1','aout':'1.1','tout':'1.0d100',
 'foutput':'2','fbackup':'1000000000','nremap':'0','remap_thresh':'0.0',
 'levelmin':'3','levelmax':'4','nsubcycle':'1,1','f_w':'0.0'}
for key,value in required.items():
 vals=re.findall(r'(?im)^\s*'+key+r'\s*=\s*([^!\n]+)',s)
 assert len(vals)==1 and vals[0].strip().lower()==value,(key,vals)
(out/'effective.nml').write_text(s)
PY
printf '%s\n' 'purpose=short noncosmological native stellar endpoint integration; synthetic cap stress, not calibrated star formation' \
 'model=QUMOND; fixed levels3,4; periodic; baryons only; initially no particles; original cooling/heating and channel-resolved SSP feedback; no sinks/AGN/dust/RT' \
 'steps=4;dt_cap=0.002;maximum_end_time=0.008;eps_star=1000 deliberately stresses the accepted 50-percent depletion cap;star_mass_quantum=1e-5 code units' \
 'noutput=1;aout=1.1;tout=1e100;foutput=2;fbackup=1e9;expected_full_dumps=2;upper_dump_bytes=134217728;expected_total_dump_bytes=268435456;reserve_including_copied_SSP_tables_exports_and_ledgers=1073741824' \
 "free_bytes_before=$free_bytes" "effective_namelist=$out/effective.nml" > "$out/output_policy.txt"
set +u
source "$(cat "$build/environment_setup_path.txt")" > "$out/compiler_environment.log" 2>&1
set -u
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE I_MPI_FABRICS=shm I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export RAR_QUMOND_STAGE_HISTORY=1 RAR_QUMOND_COMPOSITE_PUBLISH=1 RAR_QUMOND_COMPOSITE_EVOLVE=1
export RAR_QUMOND_STELLAR_SOURCES=1 RAR_QUMOND_COMPOSITE_EXPORT=1 RAR_QUMOND_GRAVITY_KICKS=1
export RAR_NATIVE_PARTICLE_ROLE=baryon RAR_QUMOND_EVOLVE_DT=0.002
export RAR_STELLAR_BIRTH_LEDGER=1 RAR_STELLAR_COMMIT_LEDGER=1 RAR_COOLING_COMMIT_LEDGER=1 RAR_GLOBAL_INVENTORY=1
export PHASE0_YIELD_TABLE=$out/phase0_yields.dat PHASE0_STELLAR_ENRICHMENT=1
export PHASE0_SNIA_RUNTIME_CONTRACT=$out/snia_effective_ssp_runtime.nml
unset RAR_QUMOND_SINK_SOURCES RAR_SINK_ACCRETION_LEDGER RAR_QUMOND_COMPOSITE_OPERATOR RAR_QUMOND_COMPOSITE_SOLVER RAR_QUMOND_COMPOSITE_CONTINUUM RAR_NATIVE_PARTICLE_CLOUD_WIDTH
cd "$out"
cat output_policy.txt
cat effective.nml
sha256sum effective.nml ramses_final3d source_commit.txt phase0_yields.dat history.nml snia_effective_ssp_runtime.nml legacy_yield_table_parser_only.asc ic/ic_part >inputs.sha256
srun --mpi=pmi2 --ntasks=1 ./ramses_final3d effective.nml >run.log 2>&1
grep -Fq 'Run completed' run.log
for marker in RAR_NATIVE_STELLAR_ENDPOINT_V1 RAR_STELLAR_BIRTH_V1 RAR_STELLAR_COMMIT_V1 RAR_COOLING_COMMIT_V2;do grep -Fq "$marker" run.log;done
if grep -Eiq 'failed to converge|NaN/Inf detected|segmentation fault|FATAL:|ERROR:|QUMOND stage history:|QUMOND composite publication|QUMOND evolving gas|Native composite particles failed|Phase 0 .*failed' run.log;then exit 1;fi
test "$(find . -name 'data_*.h5' | wc -l)" -eq 2
sha256sum -c inputs.sha256
sha256sum run.log composite_*.tsv gravity_kicks_*.tsv particle_step_*.tsv particles_epoch_*.tsv output_*/data_*.h5 >outputs.sha256
printf '%s\n' 'NATIVE_STELLAR_EXECUTION_COMPLETE_NUMERICAL_AUDIT_PENDING' > COMPLETE
