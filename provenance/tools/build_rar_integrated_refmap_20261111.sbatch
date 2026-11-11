#!/bin/bash
#SBATCH --job-name=rar-integrated-map-build
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_integrated_map_build_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_integrated_map_build_%j.err
set -Eeuo pipefail

project=/home/kjhan/BackUP/rar
build_source="$project/build/rar_integrated_refmap_20261111"
expected_revision=246958ad0a0f4b4ef14eb27237f190245b46f2f1
base=/gpfs/kjhan/RAR_QUMOND
record="$base/builds/rar_integrated_refmap_20261111"

module purge
module load intel/tbb/2022.3
module load intel/compiler-rt/2025.3.0
module load intel/umf/1.0.2
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
command -v mpiifx
test "$(git -C "$build_source" rev-parse HEAD)" = "$expected_revision"
test -z "$(git -C "$build_source" status --porcelain)"
test ! -e "$record"
mkdir "$record"

hostname > "$record/hostname.txt"
module list > "$record/modules.txt" 2>&1
mpiifx --version > "$record/compiler.txt" 2>&1
git -C "$build_source" log -1 --format=fuller > "$record/source.txt"
git -C "$build_source" remote -v > "$record/remotes.txt"
cp -p "${BASH_SOURCE[0]}" "$record/job.sbatch"

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
cd "$build_source/bin"

printf '%s\n' 'NVECTOR=32 HDF5=1 USE_FFTW=1 USE_CUDA=0 SNRT=0' > "$record/build_options.txt"
link_helpers='dust_pah_atomization.o dust_pah_hydrogen.o dust_pah_live_model.o dust_pah_radiation.o dust_stochastic.o snrt_dust_contract.o snrt_dust_ir.o snrt_moving_scatter.o'
printf '%s\n' "$link_helpers" > "$record/additional_link_objects.txt"
make clean > "$record/make.log" 2>&1
make NVECTOR=32 HDF5=1 USE_FFTW=1 USE_CUDA=0 SNRT=0 SNRT_MODOBJ="$link_helpers" >> "$record/make.log" 2>&1

test -x ramses_final3d
cp -p ramses_final3d "$record/ramses_final3d"
sha256sum "$record/ramses_final3d" > "$record/executable.sha256"
sha256sum \
    "$build_source/patch/lagRamses/amr_parameters.jaehyun.f90" \
    "$build_source/patch/lagRamses/flag_utils.kjhan.f90" \
    "$build_source/patch/lagRamses/read_hydro_params.f90" \
    "$build_source/patch/lagRamses/rho_fine.kjhan.f90" \
    "$record/ramses_final3d" "$record/job.sbatch" > "$record/PAYLOAD.sha256"
sha256sum -c "$record/PAYLOAD.sha256" --status
printf '%s\n' 'BUILD_SUCCESS' > "$record/COMPLETE"
printf 'Build complete. Record: %s\n' "$record"
