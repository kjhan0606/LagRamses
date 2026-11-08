#!/usr/bin/env bash
# Isolated CPU build for the RAR exploratory no-dust/no-RT hydro source run.
# This job compiles only; it does not initialize or evolve a simulation.
#SBATCH --job-name=rar-ssp-build
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:25:00
#SBATCH --output=/gpfs/kjhan/RAR_QUMOND/logs/rar_ssp_build_%j.out
#SBATCH --error=/gpfs/kjhan/RAR_QUMOND/logs/rar_ssp_build_%j.err

set -Eeuo pipefail
umask 077

project=/home/kjhan/BackUP/rar
source_repo=${RAR_BUILD_SOURCE_REPO:-$project/code/lagRamses}
source_revision=${RAR_BUILD_REVISION:-69c85e4b4f49f7c57d4019c8487c04fa2b70bba8}
base=/gpfs/kjhan/RAR_QUMOND
record="$base/builds/rar_effective_ssp_${SLURM_JOB_ID:?}"

test "$(git -C "$source_repo" rev-parse HEAD)" = "$source_revision"
test -z "$(git -C "$source_repo" status --porcelain)"
mkdir -m 700 "$record"
mkdir -m 700 "$record/source"
hostname >"$record/hostname.txt"
git -C "$source_repo" log -1 --format=fuller >"$record/source_revision.txt"
git -C "$source_repo" status --short --branch >"$record/source_status.txt"
printf '%s\n' "$source_revision" >"$record/source_commit.txt"
cp -p "${BASH_SOURCE[0]}" "$record/job.sbatch"

git -C "$source_repo" archive --format=tar "$source_revision" | tar -xf - -C "$record/source"
test -f "$record/source/AGENTS.md"

environment_setup_rc=127
for compiler_setup in /opt/ohpc/pub/intel/oneapi/setvars.sh /opt/ohpc/pub/intel/oneapi-2025.03/setvars.sh; do
  if [ -f "$compiler_setup" ]; then
    printf '%s\n' "$compiler_setup" >"$record/environment_setup_path.txt"
    # oneAPI's Fortran setup probes an optional variable that is unset here.
    # Load it before restoring the runner's strict shell options.
    set +e
    set +u
    source "$compiler_setup" >"$record/intel_environment.log" 2>&1
    environment_setup_rc=$?
    set -euo pipefail
    break
  fi
done
printf '%s\n' "$environment_setup_rc" >"$record/environment_setup_rc.txt"
command -v mpiifx >"$record/compiler_path.txt"
mpiifx --version >"$record/compiler.txt" 2>&1

export OMP_NUM_THREADS=1
export OMP_DYNAMIC=FALSE
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}

cd "$record/source/bin"
cat >"$record/build_options.txt" <<'EOF'
NVECTOR=32
HDF5=1
USE_FFTW=1
USE_CUDA=0
SNRT=0
PHASE0_STELLAR_ENRICHMENT=1
DUST_LIVE=0
DUST_IRON=0
DUST_PAH=0
DUST_PAH_CHARGE=0
DUST_PAH_H=0
DUST_DYNAMICS=0
CHIMES=0
EOF
link_helpers='dust_pah_atomization.o dust_pah_hydrogen.o dust_pah_live_model.o dust_pah_radiation.o dust_stochastic.o snrt_dust_contract.o snrt_dust_ir.o snrt_moving_scatter.o'
printf '%s\n' "$link_helpers" >"$record/additional_link_objects.txt"

make clean >"$record/make.log" 2>&1
make NVECTOR=32 HDF5=1 USE_FFTW=1 USE_CUDA=0 SNRT=0 \
  PHASE0_STELLAR_ENRICHMENT=1 DUST_LIVE=0 DUST_IRON=0 DUST_PAH=0 \
  DUST_PAH_CHARGE=0 DUST_PAH_H=0 DUST_DYNAMICS=0 CHIMES=0 \
  SNRT_MODOBJ="$link_helpers" >>"$record/make.log" 2>&1

test -x ramses_final3d
grep -Fq ' -o ramses_final3d' "$record/make.log"
sha256sum ramses_final3d >"$record/executable.sha256"
printf '%s\n' "$record/source/bin/ramses_final3d" >"$record/executable_path.txt"
printf '%s\n' 'SNRT=0; DUST_LIVE=0; DUST_IRON=0; DUST_PAH=0; DUST_DYNAMICS=0; hydro solver; HDF5 output' >"$record/physics_build.txt"
printf '%s\n' 'BUILD_SUCCESS_COMPILE_ONLY' >"$record/COMPLETE"
printf 'Isolated build complete: %s\n' "$record"
