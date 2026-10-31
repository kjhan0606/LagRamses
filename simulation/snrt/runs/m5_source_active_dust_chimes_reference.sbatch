#!/usr/bin/env bash
# Reference-only source-active integration: M5 + live two-size dust + CHIMES.
# The AGN-only profile isolates AGN photons; effective_population can opt into
# a combined stellar-Sed/high-chi thermal-AGN attribution run.
#SBATCH --job-name=m5-src-dust-chem
# Keep Slurm's initial working directory inside the canonical LagRamses tree,
# even when sbatch is invoked from a different project or a home directory.
#SBATCH --chdir=/gpfs/kjhan/chimes-cvode-gpu-20261001
# GPU defaults are deliberately typed; override both options for another card:
# A10: --partition=a10 --gres=gpu:A10:1
# A40: --partition=a40 --gres=gpu:A40:1
# A100: --partition=a100 --gres=gpu:A100:1
# H200: --partition=h200 --gres=gpu:H200:1
# CPU-only CHIMES: submit on grammar with --partition=normal --gres=none and M5_CHIMES_BACKEND=cpu.
#SBATCH --partition=h100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:H100:1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-source-dust-chimes-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly ASSETS=/gpfs/kjhan/LRD_JWST
readonly PROFILE=${M5_SOURCE_PROFILE:-agn_reference}
readonly EFFECTIVE_AGN_IC_MODE=${M5_EFFECTIVE_AGN_IC_MODE:-small_bh_control}
readonly VERIFY_RESTART=${M5_VERIFY_RESTART:-0}
readonly AMR_INTERFACE=${M5_AMR_INTERFACE_CHECK:-0}
readonly REQUESTED_STELLAR_PACKAGE=${M5_STELLAR_PACKAGE:-}
# Extend the effective-SSP profile to stellar terminal events when requested.
# Final-dump capture is limited to one output for post-run state inspection.
NSTEP_MAX=3
if [[ "$PROFILE" == effective_population ]]; then NSTEP_MAX=${M5_EFFECTIVE_NSTEP:-4}; fi
readonly NSTEP_MAX
readonly CAPTURE_FINAL_DUMP=${M5_CAPTURE_FINAL_DUMP:-0}
readonly RUN_ID=${M5_RUN_ID:-${SLURM_JOB_ID:-}}
readonly OMP_THREADS=${SLURM_CPUS_PER_TASK:-${M5_OMP_THREADS:-8}}
readonly MPI_RANKS=${SLURM_NTASKS:-${M5_MPI_RANKS:-1}}
readonly CHIMES_BACKEND=${M5_CHIMES_BACKEND:-cpu}
readonly CHIMES_GPU_BROKERS=${M5_CHIMES_GPU_BROKERS:-1}
readonly CUDA_STREAMS=${M5_CUDA_STREAMS:-1}
CELL_SCHEDULER=${M5_CHIMES_CELL_SCHEDULER:-}
if [[ -z "$CELL_SCHEDULER" ]]; then
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then CELL_SCHEDULER=level_queue; else CELL_SCHEDULER=cohort; fi
fi
readonly CELL_SCHEDULER
readonly TRANSPORT_BACKEND=${M5_TRANSPORT_BACKEND:-auto}
readonly MATERIAL_BACKEND=${M5_MATERIAL_BACKEND:-auto}
case "$EFFECTIVE_AGN_IC_MODE" in
  small_bh_control|thermal_seed_reference) ;;
  *) echo 'M5_EFFECTIVE_AGN_IC_MODE must be small_bh_control or thermal_seed_reference' >&2; exit 2 ;;
esac
if [[ "$PROFILE" != effective_population && "$EFFECTIVE_AGN_IC_MODE" != small_bh_control ]]; then
  echo 'M5_EFFECTIVE_AGN_IC_MODE applies only to the effective_population profile' >&2
  exit 2
fi
case "$CELL_SCHEDULER" in
  cohort|level_queue) ;;
  *) echo 'M5_CHIMES_CELL_SCHEDULER must be cohort or level_queue' >&2; exit 2 ;;
esac
test -n "${SLURM_JOB_ID:-}" || { echo 'This runner must be submitted through Slurm' >&2; exit 2; }
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  case "${SLURM_JOB_PARTITION:-}" in
    h200|h100|a100|a40|a10) ;;
    *) echo "Hybrid CHIMES requires an approved GPU partition, got: ${SLURM_JOB_PARTITION:-unknown}" >&2; exit 2 ;;
  esac
else
  [[ "${SLURM_JOB_PARTITION:-}" == normal ]] || {
    echo "CPU-only CHIMES must use grammar's normal queue, got: ${SLURM_JOB_PARTITION:-unknown}" >&2
    exit 2
  }
fi
[[ "$OMP_THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'OMP thread count must be a positive integer' >&2; exit 2; }
if [[ "$CHIMES_BACKEND" == hybrid && "$PROFILE" == effective_population && "$AMR_INTERFACE" == 1 ]] &&
   ((OMP_THREADS<4)); then
  echo 'This RT-on per-cell CHIMES hybrid path needs at least 4 OpenMP workers; refusing silent CPU fallback' >&2
  exit 2
fi
case "$CHIMES_BACKEND" in
  cpu|hybrid) ;;
  *) echo 'M5_CHIMES_BACKEND must be cpu or hybrid' >&2; exit 2 ;;
esac
[[ "$CHIMES_GPU_BROKERS" =~ ^[1-9][0-9]*$ ]] || {
  echo 'M5_CHIMES_GPU_BROKERS must be a positive integer' >&2; exit 2;
}
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  ((CHIMES_GPU_BROKERS < OMP_THREADS)) || {
    echo 'Hybrid CHIMES needs at least one OpenMP worker outside the GPU broker set' >&2; exit 2;
  }
  [[ "$CUDA_STREAMS" =~ ^[1-9][0-9]*$ ]] && ((CUDA_STREAMS <= 16)) || {
    echo 'M5_CUDA_STREAMS must be an integer from 1 through 16' >&2; exit 2;
  }
elif [[ "$CHIMES_GPU_BROKERS" != 1 || "$CUDA_STREAMS" != 1 ]]; then
  echo 'GPU broker/stream settings apply only to the hybrid CHIMES backend' >&2; exit 2
fi
if [[ "$CHIMES_BACKEND" == cpu && "$CELL_SCHEDULER" != cohort ]]; then
  echo 'CPU CHIMES uses cohort scheduling; level_queue requires the CUDA hybrid backend' >&2
  exit 2
fi
if [[ "$CHIMES_BACKEND" == hybrid && "$CELL_SCHEDULER" != level_queue ]]; then
  echo 'Hybrid CHIMES requires level_queue scheduling' >&2
  exit 2
fi
case "$TRANSPORT_BACKEND:$MATERIAL_BACKEND" in
  auto:auto|openmp:auto|openmp:openmp|auto:openmp) ;;
  *) echo 'Invalid M5 diagnostic backend selection' >&2; exit 2 ;;
esac
case "$VERIFY_RESTART" in
  0) ;;
  1) test "$PROFILE" = effective_population ;;
  *) echo 'M5_VERIFY_RESTART must be 0 or 1' >&2; exit 2 ;;
esac
case "$AMR_INTERFACE" in
  0) ;;
  1) test "$PROFILE" = effective_population; test "$VERIFY_RESTART" = 0 ;;
  *) echo 'M5_AMR_INTERFACE_CHECK must be 0 or 1' >&2; exit 2 ;;
esac
if [[ "$PROFILE" == effective_population ]]; then
  [[ "$NSTEP_MAX" =~ ^[1-9][0-9]{0,2}$ ]] && ((NSTEP_MAX<=512)) || {
    echo 'M5_EFFECTIVE_NSTEP must be an integer from 1 through 512' >&2; exit 2;
  }
else
  [[ "$CAPTURE_FINAL_DUMP" == 0 ]] || {
    echo 'M5_CAPTURE_FINAL_DUMP requires the effective_population profile' >&2; exit 2;
  }
fi
case "$CAPTURE_FINAL_DUMP" in
  0|1) ;;
  *) echo 'M5_CAPTURE_FINAL_DUMP must be 0 or 1' >&2; exit 2 ;;
esac
if [[ "$VERIFY_RESTART" == 1 || "$AMR_INTERFACE" == 1 ]]; then
  [[ "$NSTEP_MAX" == 4 && "$CAPTURE_FINAL_DUMP" == 0 ]] || {
    echo 'Restart and AMR-interface contracts require the original four-step, no-final-dump profile' >&2
    exit 2
  }
fi
if [[ "$CAPTURE_FINAL_DUMP" == 1 && "$NSTEP_MAX" -le 4 ]]; then
  echo 'Final-dump capture is reserved for an extended effective_population run' >&2
  exit 2
fi
case "$PROFILE" in
  agn_reference)
    STELLAR_PACKAGE=''
    STELLAR_INPUT=$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5
    STELLAR_HISTORY=$STELLAR_INPUT/history.nml
    STELLAR_YIELDS=$STELLAR_INPUT/yields.dat
    STELLAR_NODES=''
    STELLAR_BASE_SOURCE=''
    STELLAR_BINDING=''
    TEMPLATE=$ROOT/simulation/snrt/config/snrt_stellar_agn_dust_reference_smoke.nml
    WORK=$ROOT/.m5-source-dust-chimes-$RUN_ID
    SOURCE_IC=$ROOT/simulation/snrt/config/snrt_agn_accretion_smoke.ic_sink
    # Isolate AGN photons. Stellar feedback settings in the namelist are
    # unchanged; only the independent stellar radiation source is omitted.
    SED=''
    ;;
  effective_population)
    TEMPLATE=$ASSETS/.effective-coupled.Q9Kw29/short/run.nml
    WORK=$ROOT/.m5-effective-dust-chimes-$RUN_ID
    if [[ "$VERIFY_RESTART" == 1 ]]; then
      WORK=$ROOT/.m5-effective-restart-$RUN_ID
    fi
    if [[ "$AMR_INTERFACE" == 1 ]]; then
      WORK=$ROOT/.m5-effective-amr-$RUN_ID
    fi
    SOURCE_IC=$ASSETS/.effective-coupled.Q9Kw29/short/ic_sink
    expected_source_ic_sha=c9c1d969915a7987d24c23ac9facd39855d53bc7bb94f5ad5aff90ea6c40ea37
    if [[ "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]]; then
      # Reuse the already exercised high-chi AGN sink seed while retaining the
      # effective stellar SED and all other profile inputs. This is a
      # nonproduction source-attribution control, not a physical seed choice.
      SOURCE_IC=$ROOT/simulation/snrt/config/snrt_agn_accretion_smoke.ic_sink
      expected_source_ic_sha=971d2f0f2b0c4731fd2b9801ebc5c363dd9cd51268660eb3cf801e7a34d88a54
    fi
    : "${REQUESTED_STELLAR_PACKAGE:?set M5_STELLAR_PACKAGE to the generated source-package directory}"
    STELLAR_PACKAGE=$(realpath "$REQUESTED_STELLAR_PACKAGE")
    STELLAR_INPUT=$STELLAR_PACKAGE/input-v5
    STELLAR_HISTORY=$STELLAR_INPUT/history.nml
    STELLAR_YIELDS=$STELLAR_INPUT/yields.dat
    STELLAR_NODES=$STELLAR_INPUT/nodes.dat
    STELLAR_BASE_SOURCE=$STELLAR_INPUT/source.nml
    STELLAR_BINDING=$STELLAR_PACKAGE/effective-ssp-binding.json
    SED=$STELLAR_PACKAGE/source-effective-ssp.nml
    # Reuse the evaluated effective profile, not a newly chosen stellar model.
    printf '%s  %s\n' \
      1cc6534578fcaffa84c02d0dbc037973d32691aa059efc738bc9249dfa2cd358 "$TEMPLATE" \
      "$expected_source_ic_sha" "$SOURCE_IC" \
      | sha256sum --check --strict
    ;;
  *) echo "Unknown M5_SOURCE_PROFILE: $PROFILE" >&2; exit 2 ;;
esac
readonly TEMPLATE WORK SOURCE_IC SED
readonly SNIA_CONTRACT=$ROOT/simulation/snrt/config/fp2_snia_effective_population_runtime_v1.nml
readonly RANK_BIND=$ROOT/simulation/snrt/runs/slurm_rank_bind.sh
readonly TABLE=$ASSETS/.dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5
readonly GROUP_DIR=$ASSETS/.dust-extension.AOz7mU/groups-v7
readonly GROUP_CONTRACT=$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
readonly SECONDARY=$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml
readonly DUST_CONTRACT=$ROOT/simulation/snrt/config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml
readonly CHIMES_REPORTER_SOURCE=$ROOT/simulation/snrt/tools/chimes_runtime_counts_reporter.c
readonly HDF5_ROOT=/home/kjhan/local/hdf5
readonly SUN=$ASSETS/.snrt-performance.jDx9Iz/sundials-release

: "${M5_CHIMES_BINARY:?set M5_CHIMES_BINARY to a verified CHIMES-enabled M5 executable}"
readonly BINARY=$(realpath "$M5_CHIMES_BINARY")
readonly CHIMES_INSTALL=$(realpath "$(dirname "$BINARY")/../chimes-install")
readonly CHIMES_BUILD_RESULT=$(dirname "$CHIMES_INSTALL")/build-result.txt
test -f "$CHIMES_BUILD_RESULT" || { echo 'Missing CHIMES build provenance' >&2; exit 3; }
readonly CHIMES_LIBRARY_MODE=$(awk -F= '$1=="CHIMES_LIBRARY_MODE" {mode=$2} END {print mode}' "$CHIMES_BUILD_RESULT")
case "$CHIMES_LIBRARY_MODE" in
  reference|trace) ;;
  *) echo "Build provenance has no recognized CHIMES library mode: $CHIMES_LIBRARY_MODE" >&2; exit 3 ;;
esac
if [[ -n "${M5_CHIMES_LIBRARY_MODE:-}" && "$M5_CHIMES_LIBRARY_MODE" != "$CHIMES_LIBRARY_MODE" ]]; then
  echo "Requested CHIMES library mode $M5_CHIMES_LIBRARY_MODE does not match build mode $CHIMES_LIBRARY_MODE" >&2
  exit 3
fi
test -x "$BINARY"
test -f "$CHIMES_INSTALL/libchimes.so"
grep -aFq 'SNRT_MN_AGN_SOURCE accepted_packets=' "$BINARY" || {
  echo 'Binary lacks the committed M5 AGN source receipt; rebuild from a source revision containing f260266' >&2
  exit 3
}
if [[ "$PROFILE" == effective_population && "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]] && \
  ! grep -aFq 'SNRT_MN_STELLAR_SOURCE accepted_stars=' "$BINARY"; then
  echo 'Combined source attribution requires the M5 stellar-source receipt; rebuild from the current source revision' >&2
  exit 3
fi
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  grep -aFq 'hybrid (host CVODE/LU, asynchronous device RHS/cooling)' "$BINARY" || {
    echo 'Hybrid CHIMES requested with a binary lacking the CHIMES_CUDA backend' >&2
    exit 3
  }
else
  if ldd "$BINARY" | grep -E 'lib(cudart|cuda)\.so'; then
    echo 'CPU-only CHIMES requested with a binary linked to CUDA runtime libraries' >&2
    exit 3
  fi
fi
test ! -e "$WORK"
mkdir "$WORK"
cp "$RANK_BIND" "$WORK/rank-bind.sh"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17
if [[ "$CHIMES_BACKEND" == hybrid ]]; then module load cuda/13.0.2; fi
readonly CHIMES_REPORTER=$WORK/chimes_counts_reporter.so
gcc -std=c11 -O2 -fPIC -shared "$CHIMES_REPORTER_SOURCE" -ldl -o "$CHIMES_REPORTER"

if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n '1p')
  case "$gpu_name" in
    *A10*|*A40*|*A100*|*H100*|*H200*) ;;
    *) echo "Refusing GPU outside the approved Slurm partitions: $gpu_name" >&2; exit 2 ;;
  esac
  printf 'GPU=%s\n' "$gpu_name" | tee "$WORK/launch.txt"
else
  printf 'GPU=none (CPU-only)\n' | tee "$WORK/launch.txt"
fi
printf 'RUN_ID=%s\nBINARY=%s\n' "$RUN_ID" "$BINARY" | tee -a "$WORK/launch.txt"
printf 'CHIMES_LIBRARY_MODE=%s\n' "$CHIMES_LIBRARY_MODE" | tee -a "$WORK/launch.txt"
printf 'SLURM_JOB_ID=%s PARTITION=%s NODELIST=%s GPU_COUNT=%s GPU_IDS=%s GPU_GRES=%s\n' \
  "$SLURM_JOB_ID" "${SLURM_JOB_PARTITION:-unknown}" "${SLURM_JOB_NODELIST:-unknown}" \
  "${SLURM_GPUS_ON_NODE:-unknown}" "${SLURM_JOB_GPUS:-none}" "${SLURM_JOB_GRES:-none}" \
  | tee -a "$WORK/launch.txt"

test -f "$TABLE" -a -f "$GROUP_DIR/manifest.json"
test -f "$SUN/include/cvode/cvode.h" -a -f "$HDF5_ROOT/include/hdf5.h"
test -f "$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5"
test -f "$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5"
if [[ "$PROFILE" == agn_reference ]]; then
  test -f "$STELLAR_HISTORY" -a -f "$STELLAR_YIELDS"
fi
if [[ "$PROFILE" == effective_population ]]; then
  test -d "$STELLAR_INPUT"
  test -f "$STELLAR_INPUT/manifest.json" -a -f "$STELLAR_HISTORY"
  test -f "$STELLAR_YIELDS" -a -f "$STELLAR_NODES" -a -f "$STELLAR_BASE_SOURCE"
  test -f "$SED" -a -f "$STELLAR_BINDING"
  python3 -B - "$STELLAR_PACKAGE" "$SNIA_CONTRACT" <<'PY'
import hashlib
import json
from pathlib import Path
import re
import sys

package, snia_contract = map(Path, sys.argv[1:])
inputs = package / "input-v5"
manifest_path = inputs / "manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
binding_path = package / "effective-ssp-binding.json"
binding = json.loads(binding_path.read_text(encoding="utf-8"))

def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

if manifest.get("model") != "parsec_mixed_lowmass_truncated_v1":
    raise SystemExit("stellar source package model identity mismatch")
if manifest.get("imf_id") != 1 or manifest.get("imf_denominator") != [0.08, 600]:
    raise SystemExit("stellar source package is not normalized for Kroupa IMF id 1")
for name, expected in manifest.get("output_sha256", {}).items():
    if sha(inputs / name) != expected:
        raise SystemExit(f"stellar source package hash mismatch: {name}")

base_path = inputs / "source.nml"
effective_path = package / "source-effective-ssp.nml"
base = base_path.read_text(encoding="ascii")
effective = effective_path.read_text(encoding="ascii")
if not re.search(r"\bimf_id=1, population_id=0, binary_fraction=0\b", base):
    raise SystemExit("base SED does not identify Kroupa single-star evolution")
if re.sub(r"\bpopulation_id=0, binary_fraction=0\b",
          "population_id=2, binary_fraction=0", base, count=1) != effective:
    raise SystemExit("effective-SSP SED is not the declared population-only overlay")
history = (inputs / "history.nml").read_text(encoding="ascii")
if not re.search(r"\binput_imf_id=1\b", history):
    raise SystemExit("stellar history IMF does not match the Kroupa package")
source_hash = re.search(r"source_sha256='([0-9a-f]{64})'", base)
if source_hash is None or sha(inputs / "nodes.dat") != source_hash.group(1):
    raise SystemExit("SED source hash does not match its node history")

checks = {
    "source_package_manifest_sha256": sha(manifest_path),
    "source_package_model": manifest["model"],
    "source_package_imf_id": manifest["imf_id"],
    "source_package_node_history_sha256": sha(inputs / "nodes.dat"),
    "source_namelist_sha256": sha(base_path),
    "effective_ssp_namelist_sha256": sha(effective_path),
    "snia_runtime_contract_sha256": sha(snia_contract),
}
for key, value in checks.items():
    if binding.get(key) != value:
        raise SystemExit(f"stellar runtime binding mismatch: {key}")
if binding.get("population_id_overlay") != {"from": 0, "to": 2}:
    raise SystemExit("effective-SSP population overlay mismatch")
if binding.get("status") != "reference_control_not_production_approval":
    raise SystemExit("stellar source package lost its reference-only status")
print("M5_STELLAR_PACKAGE_BINDING_PASS model=PARSEC_mixed_lowmass imf=Kroupa population=effective_SSP status=reference_control")
PY
fi
for i in {1..9}; do test -f "$GROUP_DIR/group_$(printf '%02d' "$i").hdf5"; done
cp -- "$TEMPLATE" "$WORK/effective.nml"
cp -- "$SOURCE_IC" "$WORK/ic_sink"
test "$(grep -Ec '^[[:space:]]*n_cuda_streams=1$' "$WORK/effective.nml")" -eq 1
sed -i "s/^[[:space:]]*n_cuda_streams=1$/  n_cuda_streams=$CUDA_STREAMS/" "$WORK/effective.nml"
grep -Fqx "  n_cuda_streams=$CUDA_STREAMS" "$WORK/effective.nml"
if [[ "$PROFILE" == effective_population ]]; then
  test "$(grep -Ec '^[[:space:]]*nstepmax=4$' "$WORK/effective.nml")" -eq 1
  sed -i "s/^[[:space:]]*nstepmax=4$/  nstepmax=$NSTEP_MAX/" "$WORK/effective.nml"
  test "$(grep -Ec '^[[:space:]]*rt=\.false\.$' "$WORK/effective.nml")" -eq 1
  test "$(grep -Ec '^[[:space:]]*foutput=2$' "$WORK/effective.nml")" -eq 1
  sed -i \
    -e "s/^[[:space:]]*rt=\.false\.$/  radiation_transport='snrt_mn', snrt_moment_order=5/" \
    -e 's/^[[:space:]]*foutput=2$/  foutput=1000000/' \
    "$WORK/effective.nml"
  grep -Fq "population_model='effective_ssp'" "$WORK/effective.nml"
  grep -Fq 'use_snia=.true.' "$WORK/effective.nml"
  grep -Eq "^[[:space:]]*nstepmax=$NSTEP_MAX$" "$WORK/effective.nml"
  history_path_count=$(grep -Ec '^[[:space:]]*high_mass_history_path=' "$WORK/effective.nml")
  [[ "$history_path_count" == 1 ]] || {
    echo "Expected one high_mass_history_path in effective template, found $history_path_count" >&2
    exit 3
  }
  sed -i "s|^[[:space:]]*high_mass_history_path=.*$|  high_mass_history_path='$STELLAR_HISTORY'|" \
    "$WORK/effective.nml"
  grep -Fqx "  high_mass_history_path='$STELLAR_HISTORY'" "$WORK/effective.nml"
  test -f "$SNIA_CONTRACT"
else
test "$(grep -Fc "radiation_transport='snrt_sn'" "$WORK/effective.nml")" -eq 1
test "$(grep -Ec '^[[:space:]]*foutput=1$' "$WORK/effective.nml")" -eq 1
test "$(grep -Fc "cooling_method='original'" "$WORK/effective.nml")" -eq 1
sed -i \
  -e "s/radiation_transport='snrt_sn'/radiation_transport='snrt_mn', snrt_moment_order=5/" \
  -e 's/^[[:space:]]*foutput=1$/  foutput=1000000/' \
  -e 's/gamma=1.6666666667d0/gamma=1.6666666666666667d0/' \
  -e 's/cooling=.false./cooling=.true./' \
  -e "s/feedback_mode='legacy'/feedback_mode='channel_resolved'/" \
  -e 's/imf_mass_max_msun=120d0/imf_mass_max_msun=600d0/' \
  -e 's/channel_mass_min_msun=0.8d0,1d0,8d0,3d0,140d0/channel_mass_min_msun=2d0,2d0,2d0,3d0,2d0/' \
  -e 's/channel_mass_max_msun=120d0,8d0,120d0,8d0,260d0/channel_mass_max_msun=600d0,600d0,600d0,8d0,600d0/' \
  -e 's/use_wind=.false./use_wind=.true./' \
  -e 's/use_agb=.false./use_agb=.true./' \
  -e 's/use_snii=.false./use_snii=.true./' \
  -e 's/use_pisn=.false./use_pisn=.true./' \
  -e 's/var_region(1,12)=0.01d0/var_region(1,12)=0.00098066685346035322d0/' \
  "$WORK/effective.nml"
dust_settings=$(printf '%s\n' \
  "  dust_mass_enabled=.true., dust_cooling='chimes_neq_v1'," \
  "  dust_mass_model='carbon_olivine_2size_v1', dust_material_model='dl01_composition_v1'," \
  "  dust_optics_model='d03_transport_v1', dust_size_radius_cm=1d-6,1d-5," \
  "  dust_size_density=2.2d0,3.8d0, dust_growth=.true., dust_sputtering=.true.," \
  "  dust_coagulation=.true., dust_shattering=.true., dust_sn_shocks=.false.")
stellar_settings=$(printf '%s\n' \
  "  fate_policy='user_selected_model_v1'" \
  "  high_mass_preset='source_consistent'" \
  "  high_mass_history_path='$STELLAR_HISTORY'")
awk -v dust="$dust_settings" -v stellar="$stellar_settings" \
  'index($0,"cooling_method") {print; print dust; next}
   index($0,"feedback_mode") {print; print stellar; next}
   {print}' \
  "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
dust_ic_settings=$(printf '%s\n' \
  '  var_region(1,4)=0.0069879518072289157d0' \
  '  var_region(1,6)=0.0011207621182404036d0' \
  '  var_region(1,8)=0.00042028579434015139d0' \
  '  var_region(1,9)=0.00049033342673017661d0' \
  '  var_region(1,15)=0.003012048192771084d0' \
  '  var_region(1,16)=0.003012048192771084d0' \
  '  var_region(1,17)=0.0022590361445783136d0' \
  '  var_region(1,18)=0.000753012048192771d0' \
  '  var_region(1,19)=0.0022590361445783136d0' \
  '  var_region(1,20)=0.000753012048192771d0')
awk -v dust_ic="$dust_ic_settings" \
  'index($0,"var_region(1,14)") {print; print dust_ic; next} {print}' \
  "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
grep -Fq "population_model='single_star_ssp'" "$WORK/effective.nml"
if [[ "$PROFILE" == effective_population ]]; then
  # This profile is bound to the Kroupa id=1 source package above.
  grep -Eq '^[[:space:]]*imf_id=1$' "$WORK/effective.nml"
  grep -Eq '^[[:space:]]*input_imf_id=1([,[:space:]]|$)' "$STELLAR_HISTORY"
else
  # The AGN-only profile has no stellar SED and retains its separate Chabrier
  # id=2 feedback history. Do not validate it against the effective-profile IMF.
  grep -Eq '^[[:space:]]*imf_id=2$' "$WORK/effective.nml"
  grep -Eq '^[[:space:]]*input_imf_id=2([,[:space:]]|$)' "$STELLAR_HISTORY"
fi
grep -Fq 'imf_mass_max_msun=600d0' "$WORK/effective.nml"
grep -Fq 'channel_mass_min_msun=2d0,2d0,2d0,3d0,2d0' "$WORK/effective.nml"
grep -Fq 'channel_mass_max_msun=600d0,600d0,600d0,8d0,600d0' "$WORK/effective.nml"
grep -Fq "fate_policy='user_selected_model_v1'" "$WORK/effective.nml"
grep -Fq 'use_pisn=.true.' "$WORK/effective.nml"
grep -Fq 'gamma=1.6666666666666667d0' "$WORK/effective.nml"
grep -Fq 'dust_mass_enabled=.true.' "$WORK/effective.nml"
grep -Fq 'cooling=.true.' "$WORK/effective.nml"
grep -Fq 'var_region(1,20)=0.000753012048192771d0' "$WORK/effective.nml"
fi
if [[ "$CAPTURE_FINAL_DUMP" == 1 ]]; then
  test "$(grep -Ec '^[[:space:]]*foutput=1000000$' "$WORK/effective.nml")" -eq 1
  sed -i "s/^[[:space:]]*foutput=1000000$/  foutput=$NSTEP_MAX/" "$WORK/effective.nml"
fi
grep -Fq "radiation_transport='snrt_mn', snrt_moment_order=5" "$WORK/effective.nml"
grep -Fq "dust_cooling='chimes_neq_v1'" "$WORK/effective.nml"
grep -Fq "feedback_mode='channel_resolved'" "$WORK/effective.nml"
test "$(grep -Ec '^[[:space:]]*noutput=1$' "$WORK/effective.nml")" -eq 1
test "$(grep -Ec '^[[:space:]]*aout=2\.0d0$' "$WORK/effective.nml")" -eq 1
test "$(grep -Ec '^[[:space:]]*tout=1\.0d30$' "$WORK/effective.nml")" -eq 1
if [[ "$CAPTURE_FINAL_DUMP" == 1 ]]; then
  grep -Eq "^[[:space:]]*foutput=$NSTEP_MAX$" "$WORK/effective.nml"
else
  grep -Eq '^[[:space:]]*foutput=1000000$' "$WORK/effective.nml"
fi
grep -Eq '^[[:space:]]*fbackup=1000000$' "$WORK/effective.nml"
case "$PROFILE" in
  effective_population) grep -Eq "^[[:space:]]*nstepmax=$NSTEP_MAX$" "$WORK/effective.nml" ;;
  agn_reference) grep -Eq '^[[:space:]]*nstepmax=3$' "$WORK/effective.nml" ;;
esac
if [[ "$AMR_INTERFACE" == 1 ]]; then
  # Exercise the changed coarse/fine transport consumer, not a new galaxy
  # model: eight central level2 cells refine to level3, leaving coarse leaves.
  # Keep all source/material settings and the no-dump policy unchanged.
  sed -e 's/^[[:space:]]*levelmax=4$/  levelmax=3/' \
      -e 's/^[[:space:]]*nexpand=1$/  nexpand=0/' \
      "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
  mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
  awk '/^&REFINE_PARAMS/ {
         print; print "  void_refine=.true., void_refine_min_level=3"
         print "  x_refine=.5d0, y_refine=.5d0, z_refine=.5d0"
         print "  r_refine=.3d0, exp_refine=10d0"
         print "  interpol_var=1, interpol_type=1"; next
       } {print}' "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
  mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
  grep -Eq '^[[:space:]]*levelmax=3$' "$WORK/effective.nml"
  grep -Fq 'void_refine=.true., void_refine_min_level=3' "$WORK/effective.nml"
fi
sha256sum "$TEMPLATE" "$WORK/effective.nml" "$WORK/ic_sink" "$WORK/rank-bind.sh" \
  "$ROOT/simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch" "$BINARY" \
  "$CHIMES_REPORTER_SOURCE" "$CHIMES_REPORTER" \
  "$CHIMES_INSTALL/libchimes.so" "$CHIMES_BUILD_RESULT" "$TABLE" "$GROUP_DIR/manifest.json" \
  "$GROUP_DIR"/group_*.hdf5 "$GROUP_CONTRACT" "$SECONDARY" "$DUST_CONTRACT" \
  "$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5" \
  "$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5" \
  > "$WORK/inputs.sha256"
if [[ -n "$SED" ]]; then sha256sum "$SED" >> "$WORK/inputs.sha256"; fi
if [[ "$PROFILE" == effective_population ]]; then
  sha256sum "$SNIA_CONTRACT" "$STELLAR_INPUT/manifest.json" "$STELLAR_HISTORY" \
    "$STELLAR_YIELDS" "$STELLAR_NODES" "$STELLAR_BASE_SOURCE" "$STELLAR_BINDING" \
    >> "$WORK/inputs.sha256"
elif [[ "$PROFILE" == agn_reference ]]; then
  sha256sum "$STELLAR_HISTORY" "$STELLAR_YIELDS" >> "$WORK/inputs.sha256"
fi

runtime_libs="$CHIMES_INSTALL:$SUN/lib64:$SUN/lib:$HDF5_ROOT/lib"
if [[ "$CHIMES_BACKEND" == hybrid ]]; then runtime_libs="$runtime_libs:/opt/ohpc/pub/cuda/13.0.2/lib64"; fi
export LD_LIBRARY_PATH="$runtime_libs:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS="$OMP_THREADS" OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
export OMP_STACKSIZE=512M KMP_STACKSIZE=512M I_MPI_FABRICS=shm I_MPI_PIN=0
export SNRT_BACKEND="$TRANSPORT_BACKEND" SNRT_DUST_BACKEND="$MATERIAL_BACKEND" SNRT_AGN_MODEL=partition_reference_v1
export SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 SNRT_ALLOW_REFERENCE_CONTROL=1
export SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY"
if [[ -n "$SED" ]]; then export SNRT_STELLAR_SED="$SED"; else unset SNRT_STELLAR_SED; fi
export SNRT_DUST_CONTRACT="$DUST_CONTRACT"
export SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1
export SNRT_CHIMES_BINARY="$BINARY" SNRT_CHIMES_MAIN_DATA="$TABLE" SNRT_CHIMES_GROUP_DIR="$GROUP_DIR"
export SNRT_CHIMES_BAND_TABLE=$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5
export SNRT_CHIMES_MOLECULAR_TABLE=$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5
export SNRT_CHIMES_COMPUTE_BACKEND="$CHIMES_BACKEND" SNRT_PERF_DIAGNOSTIC=1
if [[ "$PROFILE" == effective_population ]]; then
  export PHASE0_YIELD_TABLE="$STELLAR_YIELDS"
else
  export PHASE0_YIELD_TABLE="$STELLAR_YIELDS"
fi
unset SNRT_DRIVER_TEST_SEED_SOURCE SNRT_RT_ENABLE SNRT_RT_TX_TEST_FAIL_STAGE SNRT_RT_TX_TEST_FAIL_LEAF
unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_CELL_SCHEDULER
unset SNRT_CHIMES_GPU_BROKERS SNRT_CHIMES_GPU_WORKER SNRT_CHIMES_GPU_BATCH_CELLS
unset SNRT_CHIMES_GPU_RESULT SNRT_CHIMES_GPU_COOLING
export SNRT_CHIMES_CELL_SCHEDULER="$CELL_SCHEDULER"
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  export SNRT_CHIMES_GPU_WORKER=any
  export SNRT_CHIMES_GPU_BROKERS="$CHIMES_GPU_BROKERS" SNRT_CHIMES_GPU_BATCH_CELLS=64
  export SNRT_CHIMES_GPU_RESULT=compact SNRT_CHIMES_GPU_COOLING=1
fi
unset PHASE0_SNIA_RUNTIME_CONTRACT
if [[ "$PROFILE" == effective_population ]]; then
  export PHASE0_SNIA_RUNTIME_CONTRACT="$SNIA_CONTRACT"
  unset SNRT_RT_LEVEL
fi

printf 'PROFILE=%s_M5_live_dust_CHIMES_%s_NONPRODUCTION\n' "$PROFILE" "$CHIMES_LIBRARY_MODE" | tee -a "$WORK/launch.txt"
if [[ "$PROFILE" == effective_population ]]; then
  printf 'EFFECTIVE_AGN_IC_MODE=%s\nEFFECTIVE_AGN_IC_SHA256=%s\n' \
    "$EFFECTIVE_AGN_IC_MODE" "$(sha256sum "$SOURCE_IC" | awk '{print $1}')" \
    | tee -a "$WORK/launch.txt"
fi
stellar_sed_status=enabled
if [[ -z "$SED" ]]; then stellar_sed_status=disabled; fi
printf 'MPI=%s OMP=%s; real_BH_accretion_sink; stellar_SED=%s; AGN/dust contracts=reference_control; CHIMES_backend=%s; VERIFY_RESTART=%s\n' \
  "$MPI_RANKS" "$OMP_NUM_THREADS" "$stellar_sed_status" "$CHIMES_BACKEND" "$VERIFY_RESTART" \
  | tee -a "$WORK/launch.txt"
printf 'TRANSPORT_BACKEND=%s MATERIAL_BACKEND=%s\n' "$TRANSPORT_BACKEND" "$MATERIAL_BACKEND" \
  | tee -a "$WORK/launch.txt"
printf 'CHIMES_cell_scheduler=%s\n' "$CELL_SCHEDULER" | tee -a "$WORK/launch.txt"
printf 'CHIMES_GPU_BROKERS=%s\n' "$CHIMES_GPU_BROKERS" | tee -a "$WORK/launch.txt"
printf 'CUDA_STREAMS=%s\n' "$CUDA_STREAMS" | tee -a "$WORK/launch.txt"
printf 'AMR_INTERFACE_CHECK=%s\n' "$AMR_INTERFACE" | tee -a "$WORK/launch.txt"
printf 'NSTEP_MAX=%s CAPTURE_FINAL_DUMP=%s\n' "$NSTEP_MAX" "$CAPTURE_FINAL_DUMP" | tee -a "$WORK/launch.txt"
# Explicit OMP_PLACES is needed: generic "cores" can be reconstructed over
# the whole shared allocation during MPI/OMP initialization even when srun
# initially binds its tasks. Reuse the verified per-local-rank CPU slices.
run_m5_case() {
  local directory=$1
  (
    cd "$directory"
    srun --export=ALL,LD_PRELOAD="$CHIMES_REPORTER" --mpi=pmi2 \
      --ntasks="$SLURM_NTASKS" --cpus-per-task="$SLURM_CPUS_PER_TASK" --cpu-bind=none \
      bash "$WORK/rank-bind.sh" "$BINARY" effective.nml > run.log 2>&1
  )
}
check_hybrid_process_reports() {
  local logfile=$1 minimum_gpu_reports=$2
  local report_ranks gpu_count_ranks gpu_mask_ranks
  grep -Fq 'CHIMES integrated compute backend: hybrid (host CVODE/LU, asynchronous device RHS/cooling)' "$logfile"
  report_ranks=$(grep -Ec 'CHIMES_RHS_COUNTS gpu_calls=[0-9]+ .*device_errors=0 .*streams_used=[0-9]+' "$logfile" || true)
  if [[ "$report_ranks" -ne "$MPI_RANKS" ]]; then
    echo "FAIL: CHIMES RHS counters from $report_ranks/$MPI_RANKS ranks in $logfile" >&2
    return 6
  fi
  gpu_mask_ranks=$(grep -Ec 'CHIMES_GPU_WORKER_MASK value=[0-9]+' "$logfile" || true)
  if [[ "$gpu_mask_ranks" -ne "$MPI_RANKS" ]]; then
    echo "FAIL: CHIMES GPU worker mask reported from $gpu_mask_ranks/$MPI_RANKS ranks in $logfile" >&2
    return 6
  fi
  gpu_count_ranks=$(grep -Ec 'CHIMES_RHS_COUNTS gpu_calls=[1-9][0-9]* .*device_errors=0 .*streams_used=[1-9][0-9]*' "$logfile" || true)
  if [[ "$gpu_count_ranks" -lt "$minimum_gpu_reports" ]]; then
    echo "FAIL: successful CHIMES GPU RHS reports from $gpu_count_ranks ranks; need at least $minimum_gpu_reports in $logfile" >&2
    return 6
  fi
}
if [[ "$VERIFY_RESTART" == 1 ]]; then
  mkdir "$WORK/continuous" "$WORK/resumed"
  cp "$WORK/ic_sink" "$WORK/continuous/ic_sink"
  cp "$WORK/ic_sink" "$WORK/resumed/ic_sink"
  sed 's/^[[:space:]]*foutput=1000000$/  foutput=2/' "$WORK/effective.nml" > "$WORK/continuous/effective.nml"
  sed -e 's/^[[:space:]]*foutput=1000000$/  foutput=2/' \
      -e 's/^[[:space:]]*nrestart=0$/  nrestart=1/' \
      "$WORK/effective.nml" > "$WORK/resumed/effective.nml"
  grep -Eq '^[[:space:]]*nrestart=0$' "$WORK/continuous/effective.nml"
  grep -Eq '^[[:space:]]*nrestart=1$' "$WORK/resumed/effective.nml"
  grep -Eq '^[[:space:]]*foutput=2$' "$WORK/continuous/effective.nml"
  grep -Eq '^[[:space:]]*foutput=2$' "$WORK/resumed/effective.nml"
  sha256sum "$WORK/continuous/effective.nml" "$WORK/resumed/effective.nml" >> "$WORK/inputs.sha256"
  printf 'OUTPUT_POLICY=noutput1 aout2.0d0(future) tout1.0d30 nstepmax4 foutput2 fbackup1000000; 3 raw dumps <=312MiB estimated (measured ~97MiB/dump), retained until evaluation\nCONTINUOUS_NML=%s\nRESTART_NML=%s\n' \
    "$WORK/continuous/effective.nml" "$WORK/resumed/effective.nml" | tee -a "$WORK/launch.txt"
  run_m5_case "$WORK/continuous"
  grep -q 'Run completed' "$WORK/continuous/run.log"
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/continuous/run.log" "$MPI_RANKS"; fi
  test -f "$WORK/continuous/output_00001/COMPLETE"
  test -f "$WORK/continuous/output_00002/COMPLETE"
  ln -s "$WORK/continuous/output_00001" "$WORK/resumed/output_00001"
  run_m5_case "$WORK/resumed"
  grep -q 'Run completed' "$WORK/resumed/run.log"
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/resumed/run.log" "$MPI_RANKS"; fi
  grep -Fq 'SNRT_HDF5_RADIATION_RESTORE_PASS' "$WORK/resumed/run.log"
  if grep -Eiq 'ERROR:|Aborting|forrtl:|failed to converge|nonfinite|NaN_CHK.*=[[:space:]]*[1-9]|M_N.*(rejected|failed)|CHIMES cold-stage rejection|SNRT_CHIMES_COLD_RECEIVER_REJECT|SNRT DUST_LIVE pre-commit validation failed' "$WORK/resumed/run.log"; then
    echo 'FAIL: restarted runtime rejected; all evidence retained' >&2; exit 3
  fi
  cd "$WORK/continuous"
else
  if [[ "$CAPTURE_FINAL_DUMP" == 1 ]]; then
    printf 'OUTPUT_POLICY=one final full dump at coarse step %s; noutput1 aout2.0d0(future) tout1.0d30 fbackup1000000\n' \
      "$NSTEP_MAX" | tee -a "$WORK/launch.txt"
  else
    printf 'OUTPUT_POLICY=raw_outputs0; future dummy clock; foutput/fbackup1000000\n' | tee -a "$WORK/launch.txt"
  fi
  run_m5_case "$WORK"
  cd "$WORK"
fi

grep -q 'Run completed' run.log
actual_main_steps=$(awk '$1=="Main" && $2=="step=" {value=$3} END {print value+0}' run.log)
if [[ "$actual_main_steps" != "$NSTEP_MAX" ]]; then
  echo "FAIL: reached main step $actual_main_steps; expected $NSTEP_MAX" >&2
  exit 3
fi
# A final-dump run exits below, so enforce real device work before that branch.
# Restart runs validate both continuous and resumed logs in their own path.
if [[ "$CHIMES_BACKEND" == hybrid && "$VERIFY_RESTART" != 1 ]]; then
  check_hybrid_process_reports run.log "$MPI_RANKS"
fi
grep -Fq 'Dust mass/thermal evolution (dust_mass_enabled)=T' run.log
grep -Fq 'DUST_MASS model=carbon_olivine_2size_v1' run.log
grep -Fq 'DUST_SIZE: four masses' run.log
grep -Fq 'DUST_MATERIAL: local DL01 composition' run.log
grep -Fq 'CHIMES built in double precision.' run.log
grep -Fq 'SNRT M_N live AMR transport' run.log
grep -Eq 'SNRT_MN_COMMIT level/order/substeps=.*[[:space:]]5[[:space:]]' run.log
awk '$1=="SNRT_MN_LEDGER" && $3+0>0 && $4+0>0 && $9+0>0 {found=1} END {exit !found}' run.log
if [[ "$CELL_SCHEDULER" == level_queue ]]; then
  queue_gpu_required=0
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then queue_gpu_required=1; fi
  # Level-queue continuations deliberately suppress per-cell dark timers:
  # fiber suspension makes those thread-local timings ambiguous. Require the
  # authoritative M5 material queue report instead, including full cell
  # completion; hybrid mode must also show a GPU broker.
  awk -v need_gpu="$queue_gpu_required" '
    $1=="CHIMES_LEVEL_QUEUE" {
      cells=0; completed=0; brokers=0
      for(i=2;i<=NF;i++) {
        split($i,field,"=")
        if(field[1]=="cells") cells=field[2]+0
        else if(field[1]=="completed") completed=field[2]+0
        else if(field[1]=="gpu_brokers") brokers=field[2]+0
      }
      if(cells>0 && completed==cells && (!need_gpu || brokers>0)) found=1
    }
    END {exit !found}
  ' run.log
else
  grep -Eq 'SNRT PERF DARK' run.log
fi
if [[ "$PROFILE" == effective_population ]]; then
  grep -Eq 'Stellar population model id:[[:space:]]*2' run.log
  grep -Fq 'SNRT stellar population binding: common_imf_mixed_evolution' run.log
  grep -Fq 'SNIa mass accounting: effective_ssp' run.log
  grep -Eq 'New star=[[:space:]]*[1-9]' run.log
fi
momentum_assert_cells=$(awk '
  $1=="SNRT" && $2=="PERF" && $3=="TILE" && $5~/^first=/ && $6~/^last=/ {
    sub(/^first=/,"",$5); sub(/^last=/,"",$6)
    if($6+0 >= $5+0) cells += $6-$5+1
  }
  END {print cells+0}
' run.log)
if ((momentum_assert_cells<=0)); then
  echo 'FAIL: no successfully committed material-cell momentum-assertion coverage' >&2
  exit 5
fi
if grep -Eiq 'ERROR:|Aborting|forrtl:|failed to converge|nonfinite|NaN_CHK.*=[[:space:]]*[1-9]|M_N (radiation momentum finite check failed|RAMSES (dust )?momentum closure failed)|CHIMES cold-stage rejection|SNRT_CHIMES_COLD_RECEIVER_REJECT|SNRT DUST_LIVE pre-commit validation failed' run.log; then
  echo 'FAIL: runtime error/nonfinite/nonconvergence marker; preserving evidence' >&2
  exit 3
fi
if [[ "$PROFILE" == agn_reference ]] && grep -Eq 'SNRT stellar population binding:|SNRT radiation population \(' run.log; then
  echo 'FAIL: AGN-only source run unexpectedly loaded a stellar radiation population' >&2
  exit 5
fi
shopt -s nullglob
raw_outputs=(output_[0-9][0-9][0-9][0-9][0-9])
if [[ "$VERIFY_RESTART" == 1 ]]; then
  test "${raw_outputs[*]}" = 'output_00001 output_00002'
  # Record exact raw inputs before comparison: a numerical failure must still
  # leave enough provenance to identify and reproduce the three compared dumps.
  sha256sum "$WORK/continuous/output_00001/data_00001.h5" \
    "$WORK/continuous/output_00002/data_00002.h5" \
    "$WORK/resumed/output_00002/data_00002.h5" > "$WORK/raw.sha256"
  stat -c '%s %n' "$WORK/continuous/output_00001/data_00001.h5" \
    "$WORK/continuous/output_00002/data_00002.h5" \
    "$WORK/resumed/output_00002/data_00002.h5" > "$WORK/raw-bytes.txt"
  python3 -B "$ROOT/simulation/snrt/tests/p0_hdf5_stellar_restart.py" \
    --compare-m5-coupled "$WORK/continuous/output_00002/data_00002.h5" \
      "$WORK/resumed/output_00002/data_00002.h5" \
    --source-checkpoint "$WORK/continuous/output_00001/data_00001.h5" \
    > "$WORK/restart-comparison.json"
  printf 'M5_COUPLED_RESTART_EXECUTION_PASS profile=effective_population CHIMES_backend=%s raw_dumps_retained=3 pending_driver_evaluation=true\n' "$CHIMES_BACKEND" \
    | tee "$WORK/result.txt"
  exit 0
fi
agn_evidence=not_requested
agn_packets=0
agn_photons=0
agn_energy_erg=0
if [[ "$PROFILE" == effective_population || "$PROFILE" == agn_reference ]]; then
  agn_summary=$(python3 -B - "$PROFILE" "$WORK/agn_coarse_state_v1.jsonl" run.log <<'PY'
import json
import math
import re
import sys
from pathlib import Path

profile, state_path, log_path = sys.argv[1:]
states = [
    json.loads(line)
    for line in Path(state_path).read_text(encoding="utf-8").splitlines()
    if line.strip()
]
active = []
for row in states:
    if row.get("record_type") != "agn_coarse_state":
        continue
    try:
        mass = float(row.get("coarse_accreted_bh_mass_code", 0.0))
        luminosity = float(row.get("bolometric_luminosity_erg_s", 0.0))
    except (TypeError, ValueError):
        continue
    if math.isfinite(mass) and mass > 0.0 and math.isfinite(luminosity) and luminosity > 0.0:
        active.append(row)

log_text = Path(log_path).read_text(encoding="utf-8", errors="replace")
receipt_re = re.compile(
    r"^\s*SNRT_MN_AGN_SOURCE\s+accepted_packets=(\d+)\s+"
    r"photons/energy_erg=\s*([^\s]+)\s+([^\s]+)\s*$"
)
packet_count = 0
photon_count = 0.0
energy_erg = 0.0
receipt_lines = 0
for line in log_text.splitlines():
    if "SNRT_MN_AGN_SOURCE" not in line:
        continue
    match = receipt_re.match(line)
    if not match:
        raise SystemExit("FAIL: malformed SNRT_MN_AGN_SOURCE receipt")
    packets = int(match.group(1))
    photons, energy = map(float, match.group(2, 3))
    if packets < 0 or not math.isfinite(photons) or not math.isfinite(energy) or photons < 0.0 or energy < 0.0:
        raise SystemExit("FAIL: invalid SNRT_MN_AGN_SOURCE receipt values")
    packet_count += packets
    photon_count += photons
    energy_erg += energy
    receipt_lines += 1
if receipt_lines == 0:
    raise SystemExit("FAIL: runtime emitted no M5 AGN source receipts")
if packet_count > 0 and (photon_count <= 0.0 or energy_erg <= 0.0):
    raise SystemExit("FAIL: accepted AGN packets have no positive photon/energy receipt")
if packet_count == 0 and (photon_count != 0.0 or energy_erg != 0.0):
    raise SystemExit("FAIL: zero accepted AGN packets have nonzero photon/energy receipt")

modes = {str(row.get("feedback_mode", "")).strip().upper() for row in active}
if any(mode not in {"JET", "THERMAL"} for mode in modes):
    raise SystemExit("FAIL: active AGN sidecar has an unknown feedback mode")
if profile == "agn_reference":
    if "THERMAL" not in modes:
        raise SystemExit("FAIL: AGN reference has no luminous thermal feedback state")
    if packet_count <= 0:
        raise SystemExit("FAIL: thermal AGN state produced no committed M5 photon packet")
    state = "THERMAL_radiative_packets_committed"
elif not active:
    if packet_count != 0:
        raise SystemExit("FAIL: M5 AGN packets exist without an active accretion state")
    state = "no_active_accretion"
elif packet_count > 0:
    if "THERMAL" in modes:
        state = "THERMAL_radiative_packets_committed"
    else:
        state = "JET_state_with_radiative_packets"
elif "THERMAL" in modes:
    raise SystemExit("FAIL: luminous thermal AGN state has no committed M5 photon packet")
elif modes == {"JET"}:
    state = "JET_mode_zero_radiative_packets"
else:
    state = "inactive_no_radiative_packets"

print(f"{state} {packet_count} {photon_count:.8e} {energy_erg:.8e}")
PY
  )
  read -r agn_evidence agn_packets agn_photons agn_energy_erg <<< "$agn_summary"
fi
if [[ "$PROFILE" == effective_population && "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]]; then
  grep -Eq '^[[:space:]]*haardt_madau=\.false\.$' effective.nml || {
    echo 'Combined source attribution requires the external UV background to be disabled' >&2
    exit 5
  }
  python3 -B - run.log <<'PY' | tee "$WORK/combined-source-attribution.txt"
import math
import re
import sys

path = sys.argv[1]
agn_re = re.compile(
    r"^\s*SNRT_MN_AGN_SOURCE\s+accepted_packets=(\d+)\s+"
    r"photons/energy_erg=\s*([^\s]+)\s+([^\s]+)\s*$"
)
stellar_re = re.compile(
    r"^\s*SNRT_MN_STELLAR_SOURCE\s+accepted_stars=(\d+)\s+"
    r"photons/energy_erg=\s*([^\s]+)\s+([^\s]+)\s*$"
)
commit_count = 0
agn_receipts = []
stellar_receipts = []
with open(path, encoding="utf-8", errors="replace") as stream:
    for line in stream:
        if "SNRT_MN_LEDGER source/absorbed/escaped/projection/residual/gas_heat/dust_heat/N=" in line:
            commit_count += 1
        elif "SNRT_MN_AGN_SOURCE" in line:
            match = agn_re.match(line)
            if match is None:
                raise SystemExit("FAIL: malformed M5 AGN source receipt")
            packets = int(match.group(1))
            photons, energy = (
                float(value.replace("D", "E").replace("d", "E"))
                for value in match.group(2, 3)
            )
            if packets < 0 or not math.isfinite(photons) or not math.isfinite(energy):
                raise SystemExit("FAIL: invalid M5 AGN source receipt")
            if photons < 0.0 or energy < 0.0:
                raise SystemExit("FAIL: negative M5 AGN source receipt")
            agn_receipts.append((packets, photons, energy))
        elif "SNRT_MN_STELLAR_SOURCE" in line:
            match = stellar_re.match(line)
            if match is None:
                raise SystemExit("FAIL: malformed M5 stellar source receipt")
            stars = int(match.group(1))
            photons, energy = (
                float(value.replace("D", "E").replace("d", "E"))
                for value in match.group(2, 3)
            )
            if stars < 0 or not math.isfinite(photons) or not math.isfinite(energy):
                raise SystemExit("FAIL: invalid M5 stellar source receipt")
            if photons < 0.0 or energy < 0.0:
                raise SystemExit("FAIL: negative M5 stellar source receipt")
            stellar_receipts.append((stars, photons, energy))
if not commit_count or len(agn_receipts) != commit_count or len(stellar_receipts) != commit_count:
    raise SystemExit("FAIL: M5 source receipts are missing or misaligned")

stellar_total = 0.0
agn_total = 0.0
overlap_commits = 0
for (packets, agn_photons, agn_energy), (stars, stellar_photons, stellar_energy) in zip(
    agn_receipts, stellar_receipts
):
    if packets > 0 and agn_photons > 0.0 and agn_energy <= 0.0:
        raise SystemExit("FAIL: accepted AGN photons have no positive energy receipt")
    if stars > 0 and stellar_photons > 0.0 and stellar_energy <= 0.0:
        raise SystemExit("FAIL: accepted stellar photons have no positive energy receipt")
    if packets > 0 and agn_photons > 0.0 and stars > 0 and stellar_photons > 0.0:
        overlap_commits += 1
        stellar_total += stellar_photons
        agn_total += agn_photons
if overlap_commits == 0:
    raise SystemExit("FAIL: no M5 commit contains both stellar and thermal-AGN photons")
print(
    "M5_COMBINED_SOURCE_ATTRIBUTION_PASS "
    f"overlap_commits={overlap_commits} "
    f"stellar_photons={stellar_total:.8e} "
    f"thermal_agn_photons={agn_total:.8e}"
)
PY
fi
if [[ "$CAPTURE_FINAL_DUMP" == 1 ]]; then
  test "${#raw_outputs[@]}" -eq 1
  test "${raw_outputs[0]}" = output_00001
  test -f output_00001/COMPLETE
  test -s output_00001/data_00001.h5
  raw_output_bytes=$(stat -c '%s' output_00001/data_00001.h5)
  printf 'M5_SOURCE_ACTIVE_DUST_CHIMES_%s_PASS profile=%s M5_order=5 source_scope=effective_population AGN_state=%s AGN_packets=%s AGN_photons=%s AGN_photon_energy_erg=%s M5_total_source_ledger=nonzero dust_mass=enabled CHIMES_backend=%s momentum_assert_cells=%s reached_main_step=%s final_dump_bytes=%s raw_outputs=1 pending_evaluation=true reference_controls=NONPRODUCTION\n' \
    "${CHIMES_LIBRARY_MODE^^}" "$PROFILE" "$agn_evidence" "$agn_packets" "$agn_photons" "$agn_energy_erg" \
    "$CHIMES_BACKEND" "$momentum_assert_cells" "$actual_main_steps" "$raw_output_bytes" | tee "$WORK/result.txt"
  exit 0
fi
if ((${#raw_outputs[@]})); then
  printf 'FAIL: unexpected raw output retained for inspection: %s\n' "${raw_outputs[*]}" >&2
  exit 4
fi
if [[ "$AMR_INTERFACE" == 1 ]]; then
  grep -Eq 'Level[[:space:]]+3 has[[:space:]]+8 grids' run.log
  grep -Eq 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]+2[[:space:]]+5' run.log
  grep -Eq 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]+3[[:space:]]+5' run.log
  printf 'M5_COUPLED_AMR_INTERFACE_PASS coarse_leaf_cells=56 fine_leaf_cells=64 MPI=%s OMP=%s raw_outputs=0 reference_controls=NONPRODUCTION\n' \
    "$MPI_RANKS" "$OMP_NUM_THREADS" | tee "$WORK/result.txt"
  exit 0
fi
source_scope=effective_population
if [[ "$PROFILE" == agn_reference ]]; then source_scope=AGN_only_point_source; fi
if [[ "$PROFILE" == effective_population && "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]]; then
  source_scope=stellar_plus_thermal_AGN_reference
fi
printf 'M5_SOURCE_ACTIVE_DUST_CHIMES_%s_PASS profile=%s M5_order=5 source_scope=%s AGN_state=%s AGN_packets=%s AGN_photons=%s AGN_photon_energy_erg=%s M5_total_source_ledger=nonzero dust_mass=enabled CHIMES_backend=%s momentum_assert_cells=%s reference_controls=NONPRODUCTION raw_outputs=0\n' \
  "${CHIMES_LIBRARY_MODE^^}" "$PROFILE" "$source_scope" "$agn_evidence" "$agn_packets" "$agn_photons" "$agn_energy_erg" "$CHIMES_BACKEND" "$momentum_assert_cells" \
  | tee "$WORK/result.txt"
