#!/bin/bash
#SBATCH --job-name=chimes-timeline
#SBATCH --partition=h200,h100,a100
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.cvode-gpu/chimes-timeline-%j.log
set -euo pipefail
readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly SUN=/gpfs/kjhan/LRD_JWST/.snrt-performance.jDx9Iz/sundials-release
readonly BASE=$ROOT/.cvode-gpu/chimes-broker-thread0-411811/cpu
readonly WORK=$ROOT/.cvode-gpu/chimes-timeline-$SLURM_JOB_ID
readonly BUILD=$ROOT/.chimes-real-ramses-128-build-$SLURM_JOB_ID
mkdir "$WORK"
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
command -v nsys
nvidia-smi --query-gpu=name,uuid,memory.total --format=csv
printf 'PURPOSE=one 20-second CPU/CUDA timeline from actual 128^3 chemistry, OMP8 brokers4 GPU1 streams4\n'
printf 'OUTPUT=noutput=1 aout=1.1 tout=1e100; foutput=fbackup=1000000; no scheduled dumps before intentional profile-stop termination\n'
printf 'STORAGE=bounded 20-second trace; allow 5 GiB trace budget, no simulation snapshot expected\n'
df -h "$ROOT"
SNRT_CHIMES_TILE_REGRESSION_ONLY=1 bash "$ROOT/simulation/snrt/tools/chimes_integrated_gpu_qualification.sbatch"
CHIMES_REAL_REUSE_BUILD= CHIMES_REAL_BUILD_ONLY=1 bash "$ROOT/.cvode-gpu/chimes_real_ramses_128_compare.sbatch"
module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
cp "$BASE/effective.nml" "$WORK/effective.nml"
cp "$BASE/runtime.env.sh" "$WORK/runtime.env.sh"
sed -i 's/foutput=2/foutput=1000000/' "$WORK/effective.nml"
rg -q 'n_cuda_streams=4' "$WORK/effective.nml"
rg -q 'foutput=1000000' "$WORK/effective.nml"
printf 'EFFECTIVE_NAMELIST=%s/effective.nml\n' "$WORK"
sha256sum "$BUILD/ramses_chimes3d" "$WORK/effective.nml" > "$WORK/build-input.sha256"
cd "$WORK"
source runtime.env.sh
unset LD_PRELOAD SNRT_CHIMES_RHS_VERIFY SNRT_CHIMES_RHS_DIAGNOSTIC
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=FALSE
export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm
export LD_LIBRARY_PATH="$BUILD/chimes-install:$SUN/lib64:$SUN/lib:/home/kjhan/local/hdf5/lib:/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}"
export SNRT_CHIMES_CELL_SCHEDULER=level_queue SNRT_CHIMES_GPU_BROKERS=4
export SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated SNRT_CHIMES_GPU_WORKER=any
export SNRT_CHIMES_TIMELINE=1
set +e
timeout 900 mpirun -np 1 nsys profile --trace=cuda,nvtx --sample=none --cpuctxsw=none \
  --capture-range=cudaProfilerApi --capture-range-end=stop-shutdown --kill=sigterm \
  --force-overwrite=false --output="$WORK/broker" "$BUILD/ramses_chimes3d" effective.nml \
  > "$WORK/profile-run.log" 2>&1
rc=$?
set -e
printf 'PROFILE_LAUNCH_EXIT=%s (target termination is intentional after capture)\n' "$rc"
rg 'CHIMES_TIMELINE|CUDA pool:|ERROR|FATAL' "$WORK/profile-run.log" || true
rg -q 'CHIMES_TIMELINE_BEGIN' "$WORK/profile-run.log"
rg -q 'CHIMES_TIMELINE_END' "$WORK/profile-run.log"
[[ -f "$WORK/broker.nsys-rep" ]]
nsys export --type sqlite --output="$WORK/broker.sqlite" "$WORK/broker.nsys-rep"
nsys stats --report cuda_api_sum,cuda_gpu_kern_sum,cuda_gpu_mem_time_sum,nvtx_sum \
  "$WORK/broker.nsys-rep" > "$WORK/timeline-summary.txt" 2>&1
rg -q 'reactions' "$WORK/timeline-summary.txt"
shopt -s nullglob
raw=("$WORK"/output_[0-9][0-9][0-9][0-9][0-9])
[[ ${#raw[@]} == 0 ]] || { echo 'Unexpected dump: retain for diagnosis'; exit 3; }
printf 'NO_RAW_OUTPUT_CREATED\n' > "$WORK/cleanup.txt"
echo "TIMELINE_CAPTURE_COMPLETE report=$WORK/timeline-summary.txt"
