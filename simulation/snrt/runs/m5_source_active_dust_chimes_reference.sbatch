#!/usr/bin/env bash
# Reference-only source-active integration: real BH accretion + M5 + live
# two-size dust + CHIMES. The AGN profile disables stellar photons so its M5
# source ledger is attributable to the AGN; physical tables remain references.
#SBATCH --job-name=m5-src-dust-chem
#SBATCH --partition=a10,a40,a100,h100
# For H200, override at submission: --partition=h200 --gres=gpu:H200:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-source-dust-chimes-%j.log

set -euo pipefail

readonly ROOT=/gpfs/kjhan/chimes-cvode-gpu-20261001
readonly ASSETS=/gpfs/kjhan/LRD_JWST
readonly PROFILE=${M5_SOURCE_PROFILE:-agn_reference}
readonly VERIFY_RESTART=${M5_VERIFY_RESTART:-0}
readonly AMR_INTERFACE=${M5_AMR_INTERFACE_CHECK:-0}
readonly RUN_ID=${M5_RUN_ID:-${SLURM_JOB_ID:-}}
readonly OMP_THREADS=${SLURM_CPUS_PER_TASK:-${M5_OMP_THREADS:-8}}
readonly MPI_RANKS=${SLURM_NTASKS:-${M5_MPI_RANKS:-1}}
readonly CHIMES_BACKEND=${M5_CHIMES_BACKEND:-cpu}
readonly TRANSPORT_BACKEND=${M5_TRANSPORT_BACKEND:-auto}
readonly MATERIAL_BACKEND=${M5_MATERIAL_BACKEND:-auto}
test -n "${SLURM_JOB_ID:-}" || { echo 'This GPU runner must be submitted through Slurm' >&2; exit 2; }
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
case "$PROFILE" in
  agn_reference)
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
    SED=$ASSETS/.effective-population.gQodLz/inputs/source.nml
    # Reuse the evaluated effective profile, not a newly chosen stellar model.
    printf '%s  %s\n' \
      1cc6534578fcaffa84c02d0dbc037973d32691aa059efc738bc9249dfa2cd358 "$TEMPLATE" \
      c9c1d969915a7987d24c23ac9facd39855d53bc7bb94f5ad5aff90ea6c40ea37 "$SOURCE_IC" \
      19504d2845e4ea8faf46ddced0c2e01e561dda3dfa804449221bb559ede2aea9 "$SED" \
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
test -x "$BINARY"
test -f "$CHIMES_INSTALL/libchimes.so"
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  grep -aFq 'hybrid (host CVODE/LU, asynchronous device RHS/cooling)' "$BINARY" || {
    echo 'Hybrid CHIMES requested with a binary lacking the CHIMES_CUDA backend' >&2
    exit 3
  }
fi
test ! -e "$WORK"
mkdir "$WORK"
cp "$RANK_BIND" "$WORK/rank-bind.sh"

module purge
module load intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
module load intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2
readonly CHIMES_REPORTER=$WORK/chimes_counts_reporter.so
gcc -std=c11 -O2 -fPIC -shared "$CHIMES_REPORTER_SOURCE" -ldl -o "$CHIMES_REPORTER"

gpu_name=$(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n '1p')
case "$gpu_name" in
  *A10*|*A40*|*A100*|*H100*|*H200*) ;;
  *) echo "Refusing GPU outside the approved Slurm partitions: $gpu_name" >&2; exit 2 ;;
esac
printf 'GPU=%s\nRUN_ID=%s\nBINARY=%s\n' "$gpu_name" "$RUN_ID" "$BINARY" \
  | tee "$WORK/launch.txt"
printf 'SLURM_JOB_ID=%s PARTITION=%s NODELIST=%s GPUS=%s\n' \
  "$SLURM_JOB_ID" "${SLURM_JOB_PARTITION:-unknown}" "${SLURM_JOB_NODELIST:-unknown}" \
  "${SLURM_JOB_GPUS:-${SLURM_JOB_GRES:-unknown}}" | tee -a "$WORK/launch.txt"

test -f "$TABLE" -a -f "$GROUP_DIR/manifest.json"
test -f "$SUN/include/cvode/cvode.h" -a -f "$HDF5_ROOT/include/hdf5.h"
test -f "$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5"
test -f "$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5"
test -f "$ASSETS/.effective-population.gQodLz/inputs/history.nml"
test -f "$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat"
for i in {1..9}; do test -f "$GROUP_DIR/group_$(printf '%02d' "$i").hdf5"; done
cp -- "$TEMPLATE" "$WORK/effective.nml"
cp -- "$SOURCE_IC" "$WORK/ic_sink"
if [[ "$PROFILE" == effective_population ]]; then
  test "$(grep -Ec '^[[:space:]]*rt=\.false\.$' "$WORK/effective.nml")" -eq 1
  test "$(grep -Ec '^[[:space:]]*foutput=2$' "$WORK/effective.nml")" -eq 1
  sed -i \
    -e "s/^[[:space:]]*rt=\.false\.$/  radiation_transport='snrt_mn', snrt_moment_order=5/" \
    -e 's/^[[:space:]]*foutput=2$/  foutput=1000000/' \
    "$WORK/effective.nml"
  grep -Fq "population_model='effective_ssp'" "$WORK/effective.nml"
  grep -Fq 'use_snia=.true.' "$WORK/effective.nml"
  grep -Fq 'nstepmax=4' "$WORK/effective.nml"
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
  -e 's/imf_id=2/imf_id=1/' \
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
  "  high_mass_history_path='$ASSETS/.effective-population.gQodLz/inputs/history.nml'")
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
grep -Fq "radiation_transport='snrt_mn', snrt_moment_order=5" "$WORK/effective.nml"
grep -Fq "dust_cooling='chimes_neq_v1'" "$WORK/effective.nml"
test "$(grep -Ec '^[[:space:]]*noutput=1$' "$WORK/effective.nml")" -eq 1
grep -Eq '^[[:space:]]*aout=2\.0d0$' "$WORK/effective.nml"
grep -Eq '^[[:space:]]*tout=1\.0d30$' "$WORK/effective.nml"
case "$PROFILE" in
  effective_population) grep -Eq '^[[:space:]]*nstepmax=4$' "$WORK/effective.nml" ;;
  agn_reference) grep -Eq '^[[:space:]]*nstepmax=3$' "$WORK/effective.nml" ;;
esac
grep -Eq '^[[:space:]]*foutput=1000000$' "$WORK/effective.nml"
grep -Fq "feedback_mode='channel_resolved'" "$WORK/effective.nml"
grep -Fq "population_model='single_star_ssp'" "$WORK/effective.nml"
grep -Fq 'imf_id=1' "$WORK/effective.nml"
grep -Fq 'imf_mass_max_msun=600d0' "$WORK/effective.nml"
grep -Fq 'channel_mass_min_msun=2d0,2d0,2d0,3d0,2d0' "$WORK/effective.nml"
grep -Fq 'channel_mass_max_msun=600d0,600d0,600d0,8d0,600d0' "$WORK/effective.nml"
grep -Fq "fate_policy='user_selected_model_v1'" "$WORK/effective.nml"
grep -Fq 'use_pisn=.true.' "$WORK/effective.nml"
grep -Fq 'gamma=1.6666666666666667d0' "$WORK/effective.nml"
grep -Fq 'dust_mass_enabled=.true.' "$WORK/effective.nml"
grep -Fq 'cooling=.true.' "$WORK/effective.nml"
grep -Fq 'nstepmax=3' "$WORK/effective.nml"
grep -Fq 'var_region(1,20)=0.000753012048192771d0' "$WORK/effective.nml"
fi
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
grep -Fq "radiation_transport='snrt_mn', snrt_moment_order=5" "$WORK/effective.nml"
grep -Fq "dust_cooling='chimes_neq_v1'" "$WORK/effective.nml"
grep -Eq '^[[:space:]]*foutput=1000000$' "$WORK/effective.nml"
grep -Fq 'fbackup=1000000' "$WORK/effective.nml"
sha256sum "$TEMPLATE" "$WORK/effective.nml" "$WORK/ic_sink" "$WORK/rank-bind.sh" \
  "$ROOT/simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch" "$BINARY" \
  "$CHIMES_REPORTER_SOURCE" "$CHIMES_REPORTER" \
  "$CHIMES_INSTALL/libchimes.so" "$TABLE" "$GROUP_DIR/manifest.json" \
  "$GROUP_DIR"/group_*.hdf5 "$GROUP_CONTRACT" "$SECONDARY" "$DUST_CONTRACT" \
  "$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5" \
  "$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5" \
  "$ASSETS/.effective-population.gQodLz/inputs/history.nml" \
  "$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat" \
  > "$WORK/inputs.sha256"
if [[ -n "$SED" ]]; then sha256sum "$SED" >> "$WORK/inputs.sha256"; fi
if [[ "$PROFILE" == effective_population ]]; then
  sha256sum "$SNIA_CONTRACT" "$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/nodes.dat" \
    >> "$WORK/inputs.sha256"
fi

export LD_LIBRARY_PATH="$CHIMES_INSTALL:$SUN/lib64:$SUN/lib:$HDF5_ROOT/lib:/opt/ohpc/pub/cuda/13.0.2/lib64:${LD_LIBRARY_PATH:-}"
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
export PHASE0_YIELD_TABLE=$ASSETS/.mixed-lowmass.Rbb2fZ/input-v5/yields.dat
unset SNRT_DRIVER_TEST_SEED_SOURCE SNRT_RT_ENABLE SNRT_RT_TX_TEST_FAIL_STAGE SNRT_RT_TX_TEST_FAIL_LEAF
unset SNRT_CHIMES_CVODE_BACKEND SNRT_CHIMES_RHS_BACKEND SNRT_CHIMES_CELL_SCHEDULER
unset SNRT_CHIMES_GPU_BROKERS SNRT_CHIMES_GPU_WORKER SNRT_CHIMES_GPU_BATCH_CELLS
unset SNRT_CHIMES_GPU_RESULT SNRT_CHIMES_GPU_COOLING
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  export SNRT_CHIMES_CELL_SCHEDULER=level_queue SNRT_CHIMES_GPU_WORKER=any
  export SNRT_CHIMES_GPU_BROKERS=1 SNRT_CHIMES_GPU_BATCH_CELLS=64
  export SNRT_CHIMES_GPU_RESULT=compact SNRT_CHIMES_GPU_COOLING=1
fi
unset PHASE0_SNIA_RUNTIME_CONTRACT
if [[ "$PROFILE" == effective_population ]]; then
  export PHASE0_SNIA_RUNTIME_CONTRACT="$SNIA_CONTRACT"
  unset SNRT_RT_LEVEL
fi

printf 'PROFILE=%s_M5_live_dust_CHIMES_reference_NONPRODUCTION\n' "$PROFILE" | tee -a "$WORK/launch.txt"
stellar_sed_status=enabled
if [[ -z "$SED" ]]; then stellar_sed_status=disabled; fi
printf 'MPI=%s OMP=%s; real_BH_accretion_sink; stellar_SED=%s; AGN/dust contracts=reference_control; CHIMES_backend=%s; VERIFY_RESTART=%s\n' \
  "$MPI_RANKS" "$OMP_NUM_THREADS" "$stellar_sed_status" "$CHIMES_BACKEND" "$VERIFY_RESTART" \
  | tee -a "$WORK/launch.txt"
printf 'TRANSPORT_BACKEND=%s MATERIAL_BACKEND=%s\n' "$TRANSPORT_BACKEND" "$MATERIAL_BACKEND" \
  | tee -a "$WORK/launch.txt"
printf 'AMR_INTERFACE_CHECK=%s\n' "$AMR_INTERFACE" | tee -a "$WORK/launch.txt"
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
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/continuous/run.log" 1; fi
  test -f "$WORK/continuous/output_00001/COMPLETE"
  test -f "$WORK/continuous/output_00002/COMPLETE"
  ln -s "$WORK/continuous/output_00001" "$WORK/resumed/output_00001"
  run_m5_case "$WORK/resumed"
  grep -q 'Run completed' "$WORK/resumed/run.log"
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/resumed/run.log" 0; fi
  grep -Fq 'SNRT_HDF5_RADIATION_RESTORE_PASS' "$WORK/resumed/run.log"
  if grep -Eiq 'ERROR:|Aborting|forrtl:|failed to converge|nonfinite|NaN_CHK.*=[[:space:]]*[1-9]|M_N.*(rejected|failed)|CHIMES cold-stage rejection|SNRT_CHIMES_COLD_RECEIVER_REJECT|SNRT DUST_LIVE pre-commit validation failed' "$WORK/resumed/run.log"; then
    echo 'FAIL: restarted runtime rejected; all evidence retained' >&2; exit 3
  fi
  cd "$WORK/continuous"
else
  printf 'OUTPUT_POLICY=raw_outputs0; future dummy clock; foutput/fbackup1000000\n' | tee -a "$WORK/launch.txt"
  run_m5_case "$WORK"
  cd "$WORK"
fi

grep -q 'Run completed' run.log
grep -Fq 'Dust mass/thermal evolution (dust_mass_enabled)=T' run.log
grep -Fq 'DUST_MASS model=carbon_olivine_2size_v1' run.log
grep -Fq 'DUST_SIZE: four masses' run.log
grep -Fq 'DUST_MATERIAL: local DL01 composition' run.log
grep -Fq 'CHIMES built in double precision.' run.log
grep -Fq 'SNRT M_N live AMR transport' run.log
grep -Eq 'SNRT_MN_COMMIT level/order/substeps=.*[[:space:]]5[[:space:]]' run.log
awk '$1=="SNRT_MN_LEDGER" && $3+0>0 && $4+0>0 && $9+0>0 {found=1} END {exit !found}' run.log
grep -Eq 'SNRT PERF DARK' run.log
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
  python3 -B "$ROOT/simulation/snrt/tests/p0_hdf5_stellar_restart.py" \
    --compare-m5-coupled "$WORK/continuous/output_00002/data_00002.h5" \
      "$WORK/resumed/output_00002/data_00002.h5" \
    --source-checkpoint "$WORK/continuous/output_00001/data_00001.h5" \
    > "$WORK/restart-comparison.json"
  # Hash/size collection belongs on the compute node. Keep raw data until
  # the driver has inspected the result, including ordinary-source debit.
  sha256sum "$WORK/continuous/output_00001/data_00001.h5" \
    "$WORK/continuous/output_00002/data_00002.h5" \
    "$WORK/resumed/output_00002/data_00002.h5" > "$WORK/raw.sha256"
  stat -c '%s %n' "$WORK/continuous/output_00001/data_00001.h5" \
    "$WORK/continuous/output_00002/data_00002.h5" \
    "$WORK/resumed/output_00002/data_00002.h5" > "$WORK/raw-bytes.txt"
  printf 'M5_COUPLED_RESTART_EXECUTION_PASS profile=effective_population CHIMES_backend=%s raw_dumps_retained=3 pending_driver_evaluation=true\n' "$CHIMES_BACKEND" \
    | tee "$WORK/result.txt"
  exit 0
fi
if ((${#raw_outputs[@]})); then
  printf 'FAIL: unexpected raw output retained for inspection: %s\n' "${raw_outputs[*]}" >&2
  exit 4
fi
if [[ "$AMR_INTERFACE" == 1 ]]; then
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then
    check_hybrid_process_reports run.log "$MPI_RANKS"
  fi
  grep -Eq 'Level[[:space:]]+3 has[[:space:]]+8 grids' run.log
  grep -Eq 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]+2[[:space:]]+5' run.log
  grep -Eq 'SNRT_MN_COMMIT level/order/substeps=[[:space:]]+3[[:space:]]+5' run.log
  printf 'M5_COUPLED_AMR_INTERFACE_PASS coarse_leaf_cells=56 fine_leaf_cells=64 MPI=%s OMP=%s raw_outputs=0 reference_controls=NONPRODUCTION\n' \
    "$MPI_RANKS" "$OMP_NUM_THREADS" | tee "$WORK/result.txt"
  exit 0
fi
agn_evidence=not_requested
if [[ "$PROFILE" == agn_reference ]]; then
  python3 -B - "$WORK/agn_coarse_state_v1.jsonl" <<'PY'
import json
import math
import sys

records = [
    json.loads(line)
    for line in open(sys.argv[1], encoding="utf-8")
    if line.strip()
]
if not any(
    row.get("record_type") == "agn_coarse_state"
    and row.get("feedback_mode") == "THERMAL"
    and math.isfinite(float(row.get("bolometric_luminosity_erg_s", 0.0)))
    and float(row["bolometric_luminosity_erg_s"]) > 0.0
    and math.isfinite(float(row.get("coarse_accreted_bh_mass_code", 0.0)))
    and float(row["coarse_accreted_bh_mass_code"]) > 0.0
    for row in records
):
    raise SystemExit("FAIL: AGN profile has no accreting, thermally luminous coarse state")
PY
  agn_evidence=thermal_luminous_accretion_state
fi
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  check_hybrid_process_reports run.log 1
fi
source_scope=effective_population
if [[ "$PROFILE" == agn_reference ]]; then source_scope=AGN_only_point_source; fi
printf 'M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS profile=%s M5_order=5 source_scope=%s AGN_state=%s M5_total_source_ledger=nonzero dust_mass=enabled CHIMES_backend=%s momentum_assert_cells=%s reference_controls=NONPRODUCTION raw_outputs=0\n' \
  "$PROFILE" "$source_scope" "$agn_evidence" "$CHIMES_BACKEND" "$momentum_assert_cells" \
  | tee "$WORK/result.txt"
