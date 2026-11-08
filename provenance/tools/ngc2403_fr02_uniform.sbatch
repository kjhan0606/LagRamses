#!/bin/bash
#SBATCH --job-name=rar-fr02-uniform
#SBATCH --partition=a100_pcie
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_uniform_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/ngc2403_fr02_uniform_%j.err
set -o pipefail
for setup in /opt/ohpc/pub/intel/oneapi/setvars.sh /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh; do
  if [ -f "$setup" ]; then source "$setup" >/dev/null 2>&1; setup_rc=$?; printf "oneapi_setup_status=%s\n" "$setup_rc"; break; fi
done
set -Eeuo pipefail
trap 'printf "FR02 runner failure line=%s status=%s\n" "$LINENO" "$?" >&2' ERR
umask 077
project=/home/kjhan/BackUP/rar
build=/gpfs/kjhan/RAR_QUMOND/builds/syntax_ngc2403_gpu_448235
source_root=/gpfs/kjhan/RAR_QUMOND/sources/ngc2403_fr02_lagr_input_20261107_v2
analysis_repo=$project/build/rar_qumond_amr_history_20261107
analysis_revision=${RAR_ANALYSIS_REVISION:?}
nmesh=${RAR_FR02_NMESH:-512}
case "$nmesh" in
  128) level=7; capacity=400000; dump_limit=1073741824; total_limit=2147483648; free_limit=4294967296;;
  512) level=9; capacity=25600000; dump_limit=53687091200; total_limit=107374182400; free_limit=161061273600;;
  *) exit 1;;
esac
root=/gpfs/kjhan/RAR_QUMOND/runs/ngc2403_fr02_uniform${nmesh}_${SLURM_JOB_ID:?}
source_revision=65eeb24ddb43e22e6b6acc48b6758c386d506413
expected_executable=7baf42cc884f33e9babb09b21aca5c46c380b3c32186cd562d20971c6414061c
test -f "$build/COMPLETE" && test -f "$source_root/COMPLETE" && test ! -e "$root"
grep -Fq "$source_revision" "$build/source.txt"
(cd "$source_root" && sha256sum -c artifacts.sha256)
test "$(cat "$source_root/slurm_job_id.txt")" = 1274676
IFS= read -r build_executable <"$build/executable_path.txt"
recorded_executable=$(awk 'NR==1 {print $1}' "$build/executable.sha256")
if [ "$recorded_executable" != "$expected_executable" ]; then
  printf 'Build manifest hash mismatch: expected %s, recorded %s\n' \
    "$expected_executable" "$recorded_executable" >&2
  exit 1
fi
actual_executable=$(sha256sum "$build_executable" | awk '{print $1}')
if [ "$actual_executable" != "$expected_executable" ]; then
  printf 'Executable hash mismatch: expected %s, found %s\n' \
    "$expected_executable" "$actual_executable" >&2
  exit 1
fi
free_bytes=$(df -B1 /gpfs | awk 'END {print $4}')
# Reserve 50 GiB for each of two full 512^3 force snapshots.
test "$free_bytes" -ge "$free_limit"
gpu_model=$(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n '1p')
case "$gpu_model" in
  *A100*) ;;
  *) printf 'Expected an A100 GPU, found: %s\n' "$gpu_model" >&2; exit 1 ;;
esac
mkdir "$root"
cp -p "$source_root/generated/star_profile_code.dat" "$source_root/generated/hi_profile_code.dat" \
  "$source_root/generated/source_contract.json" "$root/"
cp -p "${BASH_SOURCE[0]}" "$root/job.sbatch"
printf '%s\n' "$analysis_revision" > "$root/analysis_source_commit.txt"
for name in analyze_ngc2403_fr02_uniform.py analyze_qumond_spherical_reference.py analyze_qumond_source_diagnostic.py; do
  git -C "$analysis_repo" show "$analysis_revision:provenance/tools/$name" > "$root/$name"
done
for model in newtonian qumond; do
  git -C "$analysis_repo" show "$analysis_revision:provenance/inputs/ngc2403_fr02/$model.nml" > "$root/$model.template.nml"
done
(cd "$root" && sha256sum *.py *.template.nml > frozen_inputs.sha256)
cp -p "$build/executable.sha256" "$root/executable.sha256"
cp -p "$build/source.txt" "$build/compiler.txt" "$build/build_options.txt" \
  "$build/gpu_devices.txt" "$root/"
cp -p "$build_executable" "$root/ramses_final3d"
executable="$root/ramses_final3d"
printf 'slurm_job_id=%s\nsource_revision=%s\nbinary_sha256=%s\nbinary_path=%s\n' \
  "$SLURM_JOB_ID" "$source_revision" "$expected_executable" "$executable" \
  > "$root/run_identity.txt"
printf '%s\n' "$gpu_model" > "$root/runtime_gpu.txt"
(cd "$root" && sha256sum -c executable.sha256 --status)
printf '%s\n' "$free_bytes" > "$root/free_bytes_at_start.txt"
cat > "$root/output_policy.json" <<JSON
{"purpose":"compiled FR02 fixed-source QUMOND/Newtonian comparison at declared resolution",
 "source":"FR02 Figure 9 HI candidate with helium 1.33; SPARC stellar M/L=0.5",
 "grid":"${nmesh}^3", "box_kpc":64.0, "nmesh":$nmesh,
 "cases":2, "nstepmax":1, "noutput":1, "aout":[1.1],
 "tout":["1.0d100"], "foutput":1, "fbackup":1000000000,
 "expected_full_dumps_per_case":1,
 "measured_comparable_bytes_per_128_dump":448291224,
 "measured_comparable_bytes_per_512_dump":28684370544,
 "reserved_bytes_per_dump":$dump_limit,"reserved_bytes_total":$total_limit,
 "free_bytes_at_start":$free_bytes,"source_revision":"$source_revision",
 "binary_sha256":"$expected_executable","analysis_revision":"$analysis_revision",
 "source_preparation_job":1274676,
 "limits":"graph-derived gas; periodic fixed sources; no live halo or cosmology; not an exact SPARC reconstruction"}
JSON
export OMP_NUM_THREADS=1 I_MPI_FABRICS=shm
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
for model in newtonian qumond; do
  case_root=$root/$model
  mkdir "$case_root"
  cp -p "$root/$model.template.nml" "$case_root/effective.nml"
  env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 - "$case_root/effective.nml" "$level" "$capacity" <<'PY_GRID'
from pathlib import Path
import re,sys
p=Path(sys.argv[1]);s=p.read_text()
for key,old,new in [('levelmin','7',sys.argv[2]),('levelmax','7',sys.argv[2]),('ngridtot','400000',sys.argv[3])]:
 s,n=re.subn(rf'(?im)^(\s*{key}\s*=\s*){old}(\s*)$',lambda m:m[1]+new+m[2],s)
 if n!=1:raise ValueError((key,n))
p.write_text(s)
PY_GRID
  env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 - "$case_root/effective.nml" "$root/source_contract.json" <<'PY_SOURCE'
import json,re,sys
from pathlib import Path
p=Path(sys.argv[1]); meta=json.loads(Path(sys.argv[2]).read_text()); t=p.read_text()
value=','.join(format(meta[k],'.17g') for k in ('star_exp_scaleheight_code','gas_exp_scaleheight_code','gas_taper_start_code','gas_taper_end_code'))
t,n=re.subn(r'(?im)^(\s*gravity_params\s*=\s*)[^\n]+',lambda m:m[1]+value,t)
if n!=1: raise SystemExit('gravity_params is not unique')
p.write_text(t)
PY_SOURCE
  cp -p "$root/star_profile_code.dat" "$root/hi_profile_code.dat" "$case_root/"
  printf '%s\n' "$case_root/effective.nml" > "$case_root/effective_namelist_path.txt"
  sha256sum "$case_root/effective.nml" "$case_root/star_profile_code.dat" "$case_root/hi_profile_code.dat" > "$case_root/inputs.sha256"
  (cd "$case_root" && sha256sum -c inputs.sha256)
  : > "$case_root/jobcontrol.txt"
  env -u LD_LIBRARY_PATH /home/kjhan/miniconda3/bin/python3 - "$case_root/effective.nml" "$model" "$root/source_contract.json" "$level" "$capacity" <<'PY'
import json,math,re,sys
from pathlib import Path
t=Path(sys.argv[1]).read_text();model=sys.argv[2]
def v(name):
 m=re.findall(rf'(?im)^\s*{name}\s*=\s*([^!\n]+)',t)
 if len(m)!=1:raise SystemExit(f'{name} not unique')
 return [x.strip().lower() for x in m[0].split(',')]
expected={'cosmo':['.false.'],'pic':['.false.'],'poisson':['.true.'],'hydro':['.true.'],'static':['.true.'],
 'sink':['.false.'],'use_mond':['.true.' if model=='qumond' else '.false.'],
 'gravity_type':['-5'],
 'units_length':['1.974833652154475d23'],'units_density':['1.0d-24'],'units_time':['3.8707679655184575e15'],
 'use_fftw':['.false.'],'gpu_fft':['.true.'],'gpu_hydro':['.false.'],'gpu_poisson':['.false.'],
 'gpu_sink':['.false.'],'gpu_scalar':['.false.'],'gpu_particle':['.false.'],'gpu_auto_tune':['.false.'],
 'levelmin':[sys.argv[4]],'levelmax':[sys.argv[4]],'ngridtot':[sys.argv[5]],'nstepmax':['1'],
 'noutput':['1'],'aout':['1.1'],'tout':['1.0d100'],'foutput':['1'],'fbackup':['1000000000'],
 'mond_type':['1'],'mond_mu_type':['1'],'a0_mond':['1.2d-8'],'g_ext_mond':['0.0','0.0','0.0'],
 'gamma':['1.6666666666666667'],'nrestart':['0'],'nremap':['0'],
 'abort_on_mg_nonconvergence':['.true.']}
for k,expected_value in expected.items():
 if v(k)!=expected_value:raise SystemExit(f'{k}: {v(k)} != {expected_value}')
meta=json.loads(Path(sys.argv[3]).read_text())
params=list(map(float,v('gravity_params')))
for actual,key in zip(params,('star_exp_scaleheight_code','gas_exp_scaleheight_code','gas_taper_start_code','gas_taper_end_code')):
 if not math.isclose(actual,meta[key],rel_tol=1e-14): raise SystemExit(f'{key} mismatch')
if len(params)!=4 or not math.isclose(float(v('units_length')[0].replace('d','e')),meta['unit_length_cm'],rel_tol=1e-14): raise SystemExit('source unit mismatch')
print(f'{model}: effective namelist audit passed')
PY
  if [ "$model" = qumond ]; then export RAMSES_QUMOND_SOURCE_DIAG=1; else unset RAMSES_QUMOND_SOURCE_DIAG || true; fi
  (cd "$case_root" && mpirun -np 1 "$executable" effective.nml >run.log 2>&1)
  grep -q 'Run completed' "$case_root/run.log"
  grep -q 'cuFFT setup:' "$case_root/run.log"
  grep -q 'FFT direct solve DONE' "$case_root/run.log"
  grep -Eq 'Main step=[[:space:]]*1' "$case_root/run.log"
  if grep -Eiq 'ERROR:|failed to converge|NaN/Inf detected|Too many errors|Aborting' "$case_root/run.log"; then exit 1; fi
  if [ "$model" = qumond ]; then grep -q 'full QUMOND with phantom density' "$case_root/run.log"; fi
  test "$(find "$case_root" -maxdepth 1 -type d -name 'output_[0-9]*' | wc -l)" -eq 1
  snapshot=$case_root/output_00001/data_00001.h5
  test -f "$snapshot" && test "$(stat -c %s "$snapshot")" -le "$dump_limit"
  (cd "$case_root" && sha256sum -c inputs.sha256)
done
env -u LD_LIBRARY_PATH PYTHONPATH="$root" \
  /home/kjhan/miniconda3/bin/python3 "$root/analyze_ngc2403_fr02_uniform.py" "$root" "$nmesh" > "$root/analysis.log" 2>&1
sha256sum "$root/analysis.json" > "$root/analysis.sha256"
(cd "$root" && sha256sum -c analysis.sha256)
test "$(du -sb "$root" | awk '{print $1}')" -le "$total_limit"
(cd "$root" && find . -type f ! -name artifacts.sha256 ! -name COMPLETE \
  -print0 | sort -z | xargs -0 sha256sum > artifacts.sha256 && sha256sum -c artifacts.sha256)
(cd "$root" && sha256sum -c frozen_inputs.sha256)
touch "$root/COMPLETE"
