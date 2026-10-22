#!/usr/bin/env bash
# Prospective qualification of the approved 0.125/0.0625 Courant pair.
# No dust, transported RT, CUDA chemistry, or new resolution campaign.
#SBATCH --job-name=grackle-accuracy
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.grackle-accuracy-%j.log
set -euo pipefail
root=/gpfs/kjhan/chimes-cvode-gpu-20261001
inputs=${1:?preflight input directory required}
run=$root/.grackle-neq-accuracy-${SLURM_JOB_ID}
build=$run/build
test ! -e "$run"
test -r "$inputs/candidate/effective.nml" && test -r "$inputs/reference/effective.nml"
(cd "$root"; sha256sum --check "$inputs/assets.sha256")
mkdir -p "$build/bin" "$build/run" "$run/candidate" "$run/reference"
cp "$inputs/candidate/effective.nml" "$run/candidate/effective.nml"
cp "$inputs/reference/effective.nml" "$run/reference/effective.nml"
cp "$inputs/analyzer.py" "$inputs/temperature.c" "$inputs/Makefile" "$inputs/assets.sha256" "$run/"
cp "$run/Makefile" "$build/bin/Makefile"
cp "$inputs/submitted.sbatch" "$run/submitted.sbatch"
cp "$inputs/preflight.md" "$run/preflight.md"
for folder in patch hydro mhd pm amr poisson aton utils; do
    ln -s "$root/$folder" "$build/$folder"
done
module purge
module load gnu13 intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 \
    intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=8 OMP_STACKSIZE=64M OMP_PROC_BIND=close
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
uvdata=/home/kjhan/BACKUP/Eunha.A1/Prerequisites/Cooling_Grackle/input/CloudyData_UVB=HM2012.h5
python=/home/kjhan/miniconda3/bin/python3
cd "$build/bin"
make -j1 GRACKLE_NEQ=1 HDF5=1 USE_CUDA=0 SNRT=0 CHIMES=0 DUST_LIVE=0 \
    EXEC=ramses_grackle USE_FFTW=0 NVECTOR=32 > "$run/build.log" 2>&1
binary=$build/bin/ramses_grackle3d
test -x "$binary"
(cd "$root"; sha256sum --check "$inputs/assets.sha256") > "$run/source_identity_check.txt"
mpicc -O2 -fopenmp -I/home/kjhan/local/include "$run/temperature.c" \
    -L/home/kjhan/local/lib64 -lgrackle -L/home/kjhan/local/lib -lhdf5 -lz \
    -lgfortran -lquadmath -lm -Wl,-rpath,/home/kjhan/local/lib64 \
    -Wl,-rpath,/home/kjhan/local/lib -o "$run/temperature_api"
sha256sum "$binary" "$run/temperature_api" "$run/candidate/effective.nml" \
    "$run/reference/effective.nml" > "$run/executables_inputs.sha256"
"$python" "$run/analyzer.py" --accuracy-controls "$run/temperature_api" "$uvdata" \
    > "$run/analyzer_controls.log" 2>&1
grep -q ACCURACY_ANALYZER_CONTROLS_PASS "$run/analyzer_controls.log"
for name in candidate reference; do
    (
        cd "$run/$name"
        TIMEFORMAT='wall=%3R user=%3U sys=%3S'
        { time timeout 4800s srun --ntasks=1 --cpus-per-task=8 "$binary" effective.nml; } \
            > run.log 2> wall.txt
    )
    grep -q 'Run completed' "$run/$name/run.log"
    if grep -Ei 'ERROR:|Aborting|Too many errors|Fine multigrid Poisson failed to converge|NaN_CHK.*=[[:space:]]*[1-9][0-9]*' \
        "$run/$name/run.log"; then exit 1; fi
done
"$python" "$run/analyzer.py" --accuracy-policy "$run" "$run/temperature_api" "$uvdata" \
    > "$run/analysis.log" 2>&1
printf 'GRACKLE_APPROVED_ACCURACY_PASS %s\n' "$run"
# Preserve both snapshots for the requested final Astra review and any follow-up.
# Remove only after the evaluation is finished, with a cleanup manifest.
