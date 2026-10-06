#!/bin/bash
#SBATCH --job-name=chimes-batch-overlap
#SBATCH --partition=a100
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=0-04:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/chimes-batch-overlap-%j.log

# Matched actual 128^3 comparison: CPU reference, GPU RHS batches of 64/128,
# then bounded Nsight windows for each GPU batch size. No physics/tolerance
# changes; preserve outputs until all snapshots and parity checks pass.
set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly WORK=$ROOT/.cvode-gpu/chimes-batch-overlap-$SLURM_JOB_ID
readonly BUILD_ID=${CHIMES_BATCH_REUSE_BUILD_ID:-$SLURM_JOB_ID}
readonly BUILD=$ROOT/.chimes-real-ramses-128-build-$BUILD_ID
readonly BASE_INPUT=$ROOT/.cvode-gpu/chimes-broker-thread0-412763/cpu
readonly NML=$WORK/effective.nml
readonly ENV_FILE=$WORK/runtime.env.sh
readonly SUN=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
readonly HDF5_ROOT=/home/kjhan/local/hdf5
readonly PYTHON=/home/kjhan/miniconda3/bin/python3
mkdir "$WORK"

printf 'PURPOSE=actual_128^3_batch_size_64_vs_128_and_stream_overlap\n'
printf 'JOB_ID=%s NODE=%s PARTITION=%s\n' "$SLURM_JOB_ID" "$(hostname)" "${SLURM_JOB_PARTITION:-unknown}"
df -h "$ROOT"

# Build once if this is not a retry of the immediately preceding compilation.
if [[ -z "${CHIMES_BATCH_REUSE_BUILD_ID:-}" ]]; then
  CHIMES_REAL_BUILD_ONLY=1 bash "$ROOT/.cvode-gpu/chimes_real_ramses_128_compare.sbatch"
else
  [[ -x "$BUILD/ramses_chimes3d" ]]
  [[ -f "$ROOT/.cvode-gpu/chimes-real-ramses-128-$BUILD_ID/build.sha256" ]]
  (cd "$BUILD" && sha256sum -c "$ROOT/.cvode-gpu/chimes-real-ramses-128-$BUILD_ID/build.sha256")
fi

module purge
module load gnu13
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
readonly BINARY=$BUILD/ramses_chimes3d
readonly REPORTER=$WORK/report_chimes_counts.so
[[ -x "$BINARY" && -f "$BASE_INPUT/effective.nml" && -f "$BASE_INPUT/runtime.env.sh" ]]
cp -- "$BASE_INPUT/effective.nml" "$NML"
cp -- "$BASE_INPUT/runtime.env.sh" "$ENV_FILE"

gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n '1p')
case "$gpu_name" in *H100*|*A100*) ;; *) echo "Expected H100 or A100, got: $gpu_name" >&2; exit 2;; esac
printf 'GPU=%s\n' "$gpu_name"

# Four independent brokers require four actual streams. Preserve every other
# effective namelist field from the qualified 128^3 profile.
"$PYTHON" - "$NML" <<'PY'
import re, sys
from pathlib import Path
p = Path(sys.argv[1])
s = p.read_text()
matches = list(re.finditer(r'\bn_cuda_streams\s*=\s*\d+', s))
if len(matches) != 1:
    raise SystemExit(f'expected one n_cuda_streams setting, found {len(matches)}')
s = re.sub(r'\bn_cuda_streams\s*=\s*\d+', 'n_cuda_streams=4', s, count=1)
if 'nstepmax=2' not in s or 'noutput=1' not in s or 'foutput=2' not in s:
    raise SystemExit('unexpected effective profile; refusing to launch')
p.write_text(s)
PY

gcc -std=c11 -O2 -fPIC -shared \
  "$ROOT/.cvode-gpu/report_chimes_counts.c" -ldl -o "$REPORTER"
export LD_LIBRARY_PATH="$BUILD/chimes-install:$SUN/lib64:$SUN/lib:$HDF5_ROOT/lib:/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm
unset LD_PRELOAD SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND \
  SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_RHS_DIAGNOSTIC SNRT_CHIMES_CVODE_DIAGNOSTICS \
  SNRT_CHIMES_BUNDLE_PROFILE SNRT_CHIMES_TEST_REFERENCE SNRT_CHIMES_TEST_STATE_FILE

sha256sum "$BINARY" "$NML" "$ENV_FILE" "$REPORTER" > "$WORK/build-input.sha256"
printf 'BINARY=%s\nEFFECTIVE_NAMELIST=%s\n' "$BINARY" "$NML"
printf 'GRID=128^3 NSTEP=2 MPI=1 OMP=8 GPU_BROKERS=4 CUDA_STREAMS=4\n'
printf 'ONLY_VARIANT=SNRT_CHIMES_GPU_BATCH_CELLS 64 versus 128; compact results\n'

for arm in cpu cuda64 cuda128; do
  mkdir "$WORK/$arm"
  cp -- "$NML" "$WORK/$arm/effective.nml"
  cp -- "$ENV_FILE" "$WORK/$arm/runtime.env.sh"
done

run_arm() {
  local arm=$1 backend=$2 batch=$3 rc=0 wall
  printf 'RUN_BEGIN arm=%s backend=%s batch_cells=%s time=%s\n' "$arm" "$backend" "$batch" "$(date -Is)"
  set +e
  (
    cd "$WORK/$arm"
    # shellcheck disable=SC1091
    source runtime.env.sh
    export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
    export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm
    unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND
    if [[ "$batch" == 0 ]]; then
      unset SNRT_CHIMES_GPU_BATCH_CELLS SNRT_CHIMES_CELL_SCHEDULER \
        SNRT_CHIMES_GPU_BROKERS SNRT_CHIMES_GPU_WORKER SNRT_CHIMES_GPU_RESULT
    else
      export SNRT_CHIMES_GPU_BATCH_CELLS=$batch
      export SNRT_CHIMES_CELL_SCHEDULER=level_queue SNRT_CHIMES_GPU_BROKERS=4
      export SNRT_CHIMES_GPU_WORKER=any SNRT_CHIMES_GPU_RESULT=compact
    fi
    TIMEFORMAT='ELAPSED_SEC=%R LAUNCHER_USER_SEC=%U LAUNCHER_SYS_SEC=%S'
    exec 3>&2
    { time mpirun -genv SNRT_CHIMES_COMPUTE_BACKEND "$backend" \
        -genv LD_PRELOAD "$REPORTER" -np 1 "$BINARY" effective.nml 2>&3; } 2>timing.txt
  ) > "$WORK/$arm/ramses.log" 2>&1
  rc=$?
  set -e
  wall=$(sed -nE 's/^ELAPSED_SEC=([^ ]+).*/\1/p' "$WORK/$arm/timing.txt" | tail -n 1)
  printf 'RUN_END arm=%s rc=%s wall_s=%s time=%s\n' "$arm" "$rc" "${wall:-unknown}" "$(date -Is)"
  rg 'CHIMES GPU batch policy|CHIMES_LEVEL_QUEUE|CHIMES_CVODE_COUNTS|CHIMES_RHS_COUNTS|CHIMES_GPU_WORKER_MASK|CHIMES_RESULT_TRANSFER|FFT direct solve DONE|SNRT dark material level=|ERROR|FATAL|nonconvergence|NONFINITE' \
    "$WORK/$arm/ramses.log" || true
  [[ "$rc" -eq 0 && -n "$wall" ]] || return 1
  rg -q 'FFT direct solve DONE' "$WORK/$arm/ramses.log"
  rg -q 'SNRT dark material begin' "$WORK/$arm/ramses.log"
  ! rg -q 'SNRT RT enabled|ERROR|FATAL|nonconvergence|NONFINITE' "$WORK/$arm/ramses.log"
  if [[ "$batch" != 0 ]]; then
    rg -q "CHIMES GPU batch policy: cells=$batch maximum=128" "$WORK/$arm/ramses.log"
    rg -q 'CHIMES_GPU_WORKER_MASK value=15' "$WORK/$arm/ramses.log"
    rg -q 'CHIMES_RHS_COUNTS gpu_calls=[1-9][0-9]* .*gpu_batches=[1-9][0-9]* .*device_errors=0 .*streams_used=4' "$WORK/$arm/ramses.log"
  fi
}

run_arm cpu cpu 0
run_arm cuda64 cuda_integrated 64
run_arm cuda128 cuda_integrated 128

shopt -s nullglob
cpu_out=("$WORK/cpu"/output_[0-9][0-9][0-9][0-9][0-9])
gpu64_out=("$WORK/cuda64"/output_[0-9][0-9][0-9][0-9][0-9])
gpu128_out=("$WORK/cuda128"/output_[0-9][0-9][0-9][0-9][0-9])
[[ ${#cpu_out[@]} -eq 1 && ${#gpu64_out[@]} -eq 1 && ${#gpu128_out[@]} -eq 1 ]]

# Compare every HDF5 field against the same CPU reference; keep the report
# compact while making the numerical admission threshold explicit.
"$PYTHON" - "${cpu_out[0]}" "${gpu64_out[0]}" "${gpu128_out[0]}" "$WORK/field-comparison.tsv" <<'PY'
import sys
from pathlib import Path
import h5py
import numpy as np

cpu, gpu64, gpu128, report = map(Path, sys.argv[1:])
def h5files(d):
    return {p.name:p for p in d.glob('*.h5') if h5py.is_hdf5(p)}
def datasets(g, prefix=''):
    for name, obj in g.items():
        key=f'{prefix}/{name}'
        if isinstance(obj, h5py.Dataset): yield key,obj
        elif isinstance(obj,h5py.Group): yield from datasets(obj,key)

rows=[]
for label, other in (('cpu_vs_batch64',gpu64),('cpu_vs_batch128',gpu128)):
    af,bf=h5files(cpu),h5files(other)
    if not af or af.keys()!=bf.keys(): raise SystemExit(f'{label}: HDF5 file mismatch')
    worst=0.0; nonfinite=0; discrete=0; count=0
    for filename in sorted(af):
        with h5py.File(af[filename],'r') as a,h5py.File(bf[filename],'r') as b:
            ad=dict(datasets(a)); bd=dict(datasets(b))
            if ad.keys()!=bd.keys(): raise SystemExit(f'{label}: dataset mismatch in {filename}')
            for key in sorted(ad):
                xobj,yobj=ad[key],bd[key]
                if xobj.shape!=yobj.shape or xobj.dtype.kind!=yobj.dtype.kind:
                    raise SystemExit(f'{label}: layout mismatch {filename}:{key}')
                x=np.asarray(xobj[...]); y=np.asarray(yobj[...]); count+=1
                if x.dtype.kind=='f':
                    good=np.isfinite(x)&np.isfinite(y)
                    nonfinite+=int(good.size-np.count_nonzero(good))
                    if np.any(good):
                        delta=np.abs(x[good]-y[good])
                        peak=max(float(np.max(np.abs(x[good]))),float(np.max(np.abs(y[good]))))
                        worst=max(worst,float(np.max(delta))/max(peak,1e-300))
                elif x.dtype.kind in 'iub':
                    discrete+=int(np.count_nonzero(x!=y))
                elif not np.array_equal(x,y):
                    discrete+=1
    rows.append((label,count,nonfinite,discrete,worst))
    print(f'FIELD_PARITY {label} datasets={count} nonfinite={nonfinite} discrete_mismatches={discrete} max_field_peak_rel={worst:.12g} limit=1e-8')
    if nonfinite or discrete or worst>1e-8: raise SystemExit(f'{label}: numerical parity failed')
with report.open('w') as f:
    f.write('comparison\tdatasets\tnonfinite_pairs\tdiscrete_mismatches\tmax_field_peak_relative\tlimit\n')
    for row in rows: f.write('\t'.join(map(str,row))+'\t1e-8\n')
PY

"$PYTHON" - "$WORK" <<'PY'
import re,sys
from pathlib import Path
root=Path(sys.argv[1]); rows=[]
for arm in ('cpu','cuda64','cuda128'):
    text=(root/arm/'timing.txt').read_text()
    m=re.search(r'ELAPSED_SEC=([0-9.]+)',text)
    if not m: raise SystemExit(f'missing timing for {arm}')
    rows.append((arm,float(m.group(1))))
base=rows[0][1]
with (root/'wall-comparison.tsv').open('w') as f:
    f.write('arm\twall_s\tspeedup_vs_cpu\n')
    for arm,wall in rows:
        print(f'PERFORMANCE {arm} wall_s={wall:.3f} speedup_vs_cpu={base/wall:.6f}')
        f.write(f'{arm}\t{wall:.6f}\t{base/wall:.9f}\n')
PY

# A matched 20-second Nsight window per batch size establishes real cross-
# stream kernel overlap; profiled runs are excluded from wall-clock timing.
for batch in 64 128; do
  trace=$WORK/timeline-$batch
  mkdir "$trace"
  cp -- "$NML" "$trace/effective.nml"
  cp -- "$ENV_FILE" "$trace/runtime.env.sh"
  "$PYTHON" - "$trace/effective.nml" <<'PY'
import re,sys
from pathlib import Path
p=Path(sys.argv[1]); s=p.read_text()
if len(re.findall(r'\bfoutput\s*=\s*2\b',s))!=1:
    raise SystemExit('expected exactly one periodic output clock')
p.write_text(re.sub(r'\bfoutput\s*=\s*2\b','foutput=1000000',s,count=1))
PY
  (
    cd "$trace"
    source runtime.env.sh
    export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
    export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm
    export SNRT_CHIMES_CELL_SCHEDULER=level_queue SNRT_CHIMES_GPU_BROKERS=4
    export SNRT_CHIMES_GPU_WORKER=any SNRT_CHIMES_GPU_RESULT=compact
    export SNRT_CHIMES_GPU_BATCH_CELLS=$batch SNRT_CHIMES_TIMELINE=1
    timeout 900 mpirun -genv SNRT_CHIMES_COMPUTE_BACKEND cuda_integrated -np 1 \
      nsys profile --trace=cuda,nvtx --sample=none --cpuctxsw=none \
      --capture-range=cudaProfilerApi --capture-range-end=stop-shutdown --kill=sigterm \
      --force-overwrite=false --output="$trace/chimes-$batch" "$BINARY" effective.nml \
      > "$trace/profile-run.log" 2>&1 || true
  )
  rg -q 'CHIMES_TIMELINE_BEGIN' "$trace/profile-run.log"
  rg -q 'CHIMES_TIMELINE_END' "$trace/profile-run.log"
  [[ -f "$trace/chimes-$batch.nsys-rep" ]]
  nsys export --type sqlite --output="$trace/chimes-$batch.sqlite" "$trace/chimes-$batch.nsys-rep" >/dev/null
  nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum \
    "$trace/chimes-$batch.nsys-rep" > "$trace/summary.txt" 2>&1
  "$PYTHON" - "$trace/chimes-$batch.sqlite" "$trace/overlap.tsv" "$batch" <<'PY'
import sqlite3,sys
db=sqlite3.connect('file:'+sys.argv[1]+'?mode=ro',uri=True)
q='''WITH events AS (
 SELECT start AS t,1 AS d FROM CUPTI_ACTIVITY_KIND_KERNEL
 UNION ALL SELECT end AS t,-1 AS d FROM CUPTI_ACTIVITY_KIND_KERNEL
), grouped AS (SELECT t,SUM(d) AS d FROM events GROUP BY t), running AS (
 SELECT t,LEAD(t) OVER (ORDER BY t) AS t_next,
        SUM(d) OVER (ORDER BY t ROWS UNBOUNDED PRECEDING) AS active FROM grouped
) SELECT (SELECT COUNT(*) FROM CUPTI_ACTIVITY_KIND_KERNEL),
         (SELECT COUNT(DISTINCT streamId) FROM CUPTI_ACTIVITY_KIND_KERNEL),
         COALESCE(MAX(active),0),
         COALESCE(SUM(CASE WHEN active>1 AND t_next>t THEN t_next-t ELSE 0 END),0)
  FROM running WHERE t_next IS NOT NULL'''
kernels,streams,peak,overlap=db.execute(q).fetchone()
with open(sys.argv[2],'w') as f:
    f.write('kernel_instances\tstreams\tpeak_overlapping_kernels\tkernel_overlap_ns\n')
    f.write(f'{kernels}\t{streams}\t{peak}\t{overlap}\n')
print(f'GPU_TIMELINE batch={sys.argv[3]} kernels={kernels} streams={streams} peak_overlapping_kernels={peak} overlapped_kernel_ns={overlap}')
if not kernels or streams<2 or peak<2 or overlap<=0:
    raise SystemExit('no cross-stream kernel overlap observed')
PY
  raw=("$trace"/output_[0-9][0-9][0-9][0-9][0-9])
  [[ ${#raw[@]} -eq 0 ]] || { echo "Unexpected profiler output under $trace"; exit 3; }
done

printf 'batch\tkernels\tstreams\tpeak_overlapping_kernels\tkernel_overlap_ns\n' > "$WORK/timeline-overlap.tsv"
for batch in 64 128; do
  tail -n 1 "$WORK/timeline-$batch/overlap.tsv" | awk -v b="$batch" '{print b"\t"$0}' >> "$WORK/timeline-overlap.tsv"
done

# Stage evaluation passed: remove only the three exact temporary snapshots;
# retain effective inputs, logs, comparison tables, build hashes and traces.
for arm in cpu cuda64 cuda128; do
  raw=("$WORK/$arm"/output_[0-9][0-9][0-9][0-9][0-9])
  [[ ${#raw[@]} -eq 1 && -d "${raw[0]}" ]]
  printf 'CLEANUP_RAW_OUTPUT path=%s\n' "${raw[0]}"
  rm -r -- "${raw[0]}"
done
printf 'STAGE_EVALUATION=PASS\n' > "$WORK/cleanup-manifest.txt"
printf 'Retained: effective namelist/env, binary/input hashes, per-arm logs, field parity, wall comparison, 64/128 timeline summaries and Nsight reports.\n' >> "$WORK/cleanup-manifest.txt"
echo "CHIMES_BATCH_OVERLAP_QUALIFICATION=PASS results=$WORK"
