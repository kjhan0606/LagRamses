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
# A100 PCIe: --partition=a100_pcie --gres=gpu:A100:1
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
readonly CANDIDATE_25G_DIAGNOSTIC=${M5_25G_AGN_DIAGNOSTIC:-0}
readonly CANDIDATE_25G_STELLAR_AGN=${M5_25G_STELLAR_AGN_DIAGNOSTIC:-0}
readonly CANDIDATE_25G_RECEIVER=${M5_25G_RECEIVER:-hot_atomic_d03}
readonly CANDIDATE_RESTART_ORDERING=${M5_25G_RESTART_ORDERING:-}
readonly AGN_RESTART_TRACE=${M5_25G_AGN_RESTART_TRACE:-0}
readonly AGN_STENCIL_TRACE=${M5_25G_AGN_STENCIL_TRACE:-0}
readonly AGN_DONOR_KEYS=${M5_25G_AGN_DONOR_KEYS:-}
readonly AGN_DONOR_SNAPSHOT=${M5_25G_AGN_DONOR_SNAPSHOT:-0}
readonly AGN_DONOR_OCT=${M5_25G_AGN_DONOR_OCT:-}
readonly AGN_LEVEL_ENTRY_TRACE=${M5_25G_AGN_LEVEL_ENTRY_TRACE:-0}
readonly EFFECTIVE_AGN_IC_MODE=${M5_EFFECTIVE_AGN_IC_MODE:-small_bh_control}
readonly EFFECTIVE_AGN_ENABLED=${M5_EFFECTIVE_AGN_ENABLED:-1}
readonly VERIFY_RESTART=${M5_VERIFY_RESTART:-0}
readonly AMR_INTERFACE=${M5_AMR_INTERFACE_CHECK:-0}
readonly REQUESTED_STELLAR_PACKAGE=${M5_STELLAR_PACKAGE:-}
# Extend the effective-SSP profile to stellar terminal events when requested.
# Final-dump capture is limited to one output for post-run state inspection.
NSTEP_MAX=3
if [[ "$PROFILE" == effective_population ]]; then NSTEP_MAX=${M5_EFFECTIVE_NSTEP:-4}; fi
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 ]]; then NSTEP_MAX=4; fi
readonly NSTEP_MAX
readonly CAPTURE_FINAL_DUMP=${M5_CAPTURE_FINAL_DUMP:-0}
readonly RUN_ID=${M5_RUN_ID:-${SLURM_JOB_ID:-}}
readonly OMP_THREADS=${SLURM_CPUS_PER_TASK:-${M5_OMP_THREADS:-8}}
readonly MPI_RANKS=${SLURM_NTASKS:-${M5_MPI_RANKS:-1}}
test -n "${SLURM_JOB_ID:-}" || { echo 'This runner must be submitted through Slurm' >&2; exit 2; }
if [[ -n "${M5_CHIMES_BACKEND:-}" ]]; then
  CHIMES_BACKEND=$M5_CHIMES_BACKEND
else
  case "${SLURM_JOB_PARTITION:-}" in
    h200|h100|a100|a100_pcie|a40|a10) CHIMES_BACKEND=hybrid ;;
    normal) CHIMES_BACKEND=cpu ;;
    *)
      echo "Set M5_CHIMES_BACKEND explicitly for Slurm partition ${SLURM_JOB_PARTITION:-unknown}" >&2
      exit 2
      ;;
  esac
fi
readonly CHIMES_BACKEND
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
  small_bh_control|thermal_seed_reference|agn25_m8_candidate) ;;
  *) echo 'M5_EFFECTIVE_AGN_IC_MODE must be small_bh_control, thermal_seed_reference, or agn25_m8_candidate' >&2; exit 2 ;;
esac
case "$EFFECTIVE_AGN_ENABLED" in
  0|1) ;;
  *) echo 'M5_EFFECTIVE_AGN_ENABLED must be 0 or 1' >&2; exit 2 ;;
esac
if [[ "$EFFECTIVE_AGN_ENABLED" == 0 && \
      ( "$PROFILE" != effective_population || "$CANDIDATE_25G_DIAGNOSTIC" == 1 ) ]]; then
  echo 'M5_EFFECTIVE_AGN_ENABLED=0 is limited to the effective_population reference profile' >&2
  exit 2
fi
case "$CANDIDATE_25G_STELLAR_AGN" in
  0|1) ;;
  *) echo 'M5_25G_STELLAR_AGN_DIAGNOSTIC must be 0 or 1' >&2; exit 2 ;;
esac
if [[ "$PROFILE" != effective_population && "$EFFECTIVE_AGN_IC_MODE" != small_bh_control ]]; then
  echo 'M5_EFFECTIVE_AGN_IC_MODE applies only to the effective_population profile' >&2
  exit 2
fi
if [[ "$EFFECTIVE_AGN_IC_MODE" == agn25_m8_candidate && "$CANDIDATE_25G_STELLAR_AGN" != 1 ]]; then
  echo 'agn25_m8_candidate is reserved for the bounded 25-group stellar+AGN diagnostic' >&2
  exit 2
fi
if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 && \
      ( "$PROFILE" != effective_population || "$EFFECTIVE_AGN_IC_MODE" != agn25_m8_candidate ) ]]; then
  echo '25-group stellar+AGN diagnostics require effective_population with agn25_m8_candidate' >&2
  exit 2
fi
if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 && "${M5_25G_COURANT_FACTOR:-}" != '5.0d-8' ]]; then
  echo 'The bounded 25-group stellar+AGN diagnostic requires M5_25G_COURANT_FACTOR=5.0d-8' >&2
  exit 2
fi
case "$CELL_SCHEDULER" in
  cohort|level_queue) ;;
  *) echo 'M5_CHIMES_CELL_SCHEDULER must be cohort or level_queue' >&2; exit 2 ;;
esac
if [[ "$CHIMES_BACKEND" == hybrid ]]; then
  case "${SLURM_JOB_PARTITION:-}" in
    h200|h100|a100|a100_pcie|a40|a10) ;;
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
  1)
    [[ "$PROFILE" == effective_population || \
       ( "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$PROFILE" == agn_reference && \
         "$CANDIDATE_25G_RECEIVER" == cold_transition_d03 ) ]] || {
      echo 'M5_VERIFY_RESTART is limited to effective_population or the cold 25-group diagnostic' >&2
      exit 2
    }
    ;;
  *) echo 'M5_VERIFY_RESTART must be 0 or 1' >&2; exit 2 ;;
esac
case "$CANDIDATE_RESTART_ORDERING" in
  '') ;;
  ksection)
    [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 ]] || {
      echo 'M5_25G_RESTART_ORDERING is restricted to the 25-group candidate restart diagnostic' >&2
      exit 2
    }
    ;;
  *) echo 'M5_25G_RESTART_ORDERING must be empty or ksection' >&2; exit 2 ;;
esac
case "$AGN_RESTART_TRACE" in
  0) ;;
  1)
    [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 && \
       -z "$CANDIDATE_RESTART_ORDERING" ]] || {
      echo 'M5_25G_AGN_RESTART_TRACE is restricted to the Hilbert candidate restart diagnostic' >&2
      exit 2
    }
    ;;
  *) echo 'M5_25G_AGN_RESTART_TRACE must be 0 or 1' >&2; exit 2 ;;
esac
case "$AGN_DONOR_SNAPSHOT" in
  0)
    [[ -z "$AGN_DONOR_OCT" ]] || {
      echo 'M5_25G_AGN_DONOR_OCT requires donor snapshot mode' >&2
      exit 2
    }
    ;;
  1)
    [[ "$AGN_RESTART_TRACE" == 1 && "$CANDIDATE_25G_DIAGNOSTIC" == 1 && \
       "$VERIFY_RESTART" == 1 && -z "$CANDIDATE_RESTART_ORDERING" && \
       -n "$AGN_DONOR_OCT" ]] || {
      echo 'M5_25G_AGN_DONOR_SNAPSHOT is restricted to a targeted Hilbert candidate restart diagnostic with an oct key' >&2
      exit 2
    }
    ;;
  *) echo 'M5_25G_AGN_DONOR_SNAPSHOT must be 0 or 1' >&2; exit 2 ;;
esac
case "$AGN_LEVEL_ENTRY_TRACE" in
  0) ;;
  1)
    [[ "$AGN_RESTART_TRACE" == 1 && "$AGN_DONOR_SNAPSHOT" == 1 && \
       "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 && \
       -z "$CANDIDATE_RESTART_ORDERING" ]] || {
      echo 'M5_25G_AGN_LEVEL_ENTRY_TRACE requires the targeted Hilbert restart diagnostic' >&2
      exit 2
    }
    ;;
  *) echo 'M5_25G_AGN_LEVEL_ENTRY_TRACE must be 0 or 1' >&2; exit 2 ;;
esac
case "$AGN_STENCIL_TRACE" in
  0) ;;
  1)
    [[ "$AGN_RESTART_TRACE" == 1 && "$CANDIDATE_25G_DIAGNOSTIC" == 1 && \
       "$VERIFY_RESTART" == 1 && -z "$CANDIDATE_RESTART_ORDERING" ]] || {
      echo 'M5_25G_AGN_STENCIL_TRACE is restricted to the traced Hilbert candidate restart diagnostic' >&2
      exit 2
    }
    ;;
  *) echo 'M5_25G_AGN_STENCIL_TRACE must be 0 or 1' >&2; exit 2 ;;
esac
if [[ -n "$AGN_DONOR_KEYS" ]]; then
  [[ ("$AGN_STENCIL_TRACE" == 1 || "$AGN_DONOR_SNAPSHOT" == 1) && \
     -f "$AGN_DONOR_KEYS" ]] || {
    echo 'M5_25G_AGN_DONOR_KEYS requires stencil or targeted donor-snapshot tracing and an existing key file' >&2
    exit 2
  }
fi
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
case "$CANDIDATE_25G_DIAGNOSTIC" in
  0)
    [[ -z "${M5_25G_RECEIVER:-}" ]] || {
      echo 'M5_25G_RECEIVER is valid only with M5_25G_AGN_DIAGNOSTIC=1' >&2
      exit 2
    }
    [[ "$CANDIDATE_25G_STELLAR_AGN" == 0 ]] || {
      echo 'M5_25G_STELLAR_AGN_DIAGNOSTIC also requires M5_25G_AGN_DIAGNOSTIC=1' >&2
      exit 2
    }
    ;;
  1)
    if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 ]]; then
      [[ "$PROFILE" == effective_population && "$AMR_INTERFACE" == 0 && \
         "$CAPTURE_FINAL_DUMP" == 0 && "$VERIFY_RESTART" == 0 && \
         "$NSTEP_MAX" == 4 && \
         "$CANDIDATE_25G_RECEIVER" == cold_transition_d03 && \
         "$CANDIDATE_RESTART_ORDERING" == '' && "$AGN_RESTART_TRACE" == 0 && \
         "$AGN_STENCIL_TRACE" == 0 && "$AGN_DONOR_SNAPSHOT" == 0 ]] || {
        echo 'The joint 25-group fixture is restricted to effective_population, cold_transition_d03, and a non-restart bounded run' >&2
        exit 2
      }
    else
      [[ "$PROFILE" == agn_reference && "$AMR_INTERFACE" == 0 && \
         "$CAPTURE_FINAL_DUMP" == 0 ]] || {
        echo 'The 25-group AGN-only candidate requires the AGN diagnostic fixture' >&2
        exit 2
      }
    fi
    case "$CANDIDATE_25G_RECEIVER" in
      hot_atomic_d03|cold_transition_d03) ;;
      *) echo 'M5_25G_RECEIVER must be hot_atomic_d03 or cold_transition_d03' >&2; exit 2 ;;
    esac
    if [[ "$VERIFY_RESTART" == 1 ]]; then
      [[ "$CANDIDATE_25G_RECEIVER" == cold_transition_d03 && "$NSTEP_MAX" == 4 ]] || {
        echo 'Candidate restart is limited to the four-step cold-transition receiver' >&2
        exit 2
      }
    fi
    if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 && "$VERIFY_RESTART" == 1 ]]; then
      echo '25-group stellar+AGN diagnostic does not yet claim restart equivalence' >&2
      exit 2
    fi
    ;;
  *) echo 'M5_25G_AGN_DIAGNOSTIC must be 0 or 1' >&2; exit 2 ;;
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
    if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
      # Match the fixed Su M8, mdot=-1 source state without a forced packet.
      SOURCE_IC=$ROOT/simulation/snrt/config/snrt_agn25_m8_mdotm1_smoke.ic_sink
    else
      SOURCE_IC=$ROOT/simulation/snrt/config/snrt_agn_accretion_smoke.ic_sink
    fi
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
    if [[ "$EFFECTIVE_AGN_IC_MODE" == agn25_m8_candidate ]]; then
      SOURCE_IC=$ROOT/simulation/snrt/config/snrt_agn25_m8_mdotm1_smoke.ic_sink
      expected_source_ic_sha=ad2b21307ac2217ece39cf610d9071aea3758045dabc78c562eac3aa949c24cf
    elif [[ "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]]; then
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
GROUP_DIR=$ASSETS/.dust-extension.AOz7mU/groups-v7
GROUP_CONTRACT=$ROOT/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
readonly SECONDARY=$ROOT/simulation/snrt/config/snrt_secondary_table_contract_v1.nml
DUST_CONTRACT=$ROOT/simulation/snrt/config/dust_dl01_bulk_030_scattering_exchange_reference_v4.nml
CHIMES_ATOMIC_TABLE=$ASSETS/.chimes-spectral.TrQxAG/bank-final/atomic_shells.h5
CHIMES_MOLECULAR_TABLE=$ASSETS/.chimes-molecular.U6ZQCr/bank/molecular_nodes.h5
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  GROUP_DIR=${M5_25G_CHIMES_GROUP_DIR:-$ROOT/simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_chimes_groups_live25}
  GROUP_CONTRACT=$ROOT/simulation/snrt/config/snrt_group_contract_agn25_candidate_diagnostic_v1.nml
  DUST_CONTRACT=$ROOT/simulation/snrt/config/dust_agn25_d03_live_ir_candidate_v1.nml
  CHIMES_ATOMIC_TABLE=$ROOT/simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_chimes_bank128/atomic_shells.h5
  CHIMES_MOLECULAR_TABLE=$ROOT/simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_chimes_bank128/molecular_nodes.h5
fi
readonly GROUP_DIR GROUP_CONTRACT DUST_CONTRACT CHIMES_ATOMIC_TABLE CHIMES_MOLECULAR_TABLE
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
grep -aFq 'SNRT_MN_AGN_ENERGY input/injected/unrepresented_erg=' "$BINARY" || {
  echo 'Binary lacks the M5 AGN energy-partition receipt; rebuild from the current source' >&2
  exit 3
}
if [[ "$PROFILE" == effective_population && \
      ( "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference || "$CANDIDATE_25G_STELLAR_AGN" == 1 ) ]] && \
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
test -f "$CHIMES_ATOMIC_TABLE"
test -f "$CHIMES_MOLECULAR_TABLE"
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  python3 -B - "$ROOT/simulation/snrt/data/snrt_agn25_joint_diag_20261102_v1_asset_manifest.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path('/gpfs/kjhan/chimes-cvode-gpu-20261001')
manifest=json.loads(Path(sys.argv[1]).read_text())
if manifest.get('group_profile',{}).get('runtime_admissible') is not False:
    raise SystemExit('25-group candidate asset manifest must remain runtime-inadmissible')
for name,digest in manifest['files_sha256'].items():
    if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
        raise SystemExit(f'candidate asset SHA256 mismatch: {name}')
print('SNRT_AGN25_CANDIDATE_ASSET_HASHES_OK')
PY
  python3 -B - "$GROUP_DIR/manifest.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
p=Path(sys.argv[1]); d=json.loads(p.read_text())
expected='439a3f1317ebdec5aa0054026b82d9ee0d059104a2ba461ca1aeb44d790ed58e'
if d.get('group_set_id')!='snrt_agn25_joint_diag_20261102_v1' or d.get('group_edges_sha256')!=expected:
    raise SystemExit('25-group CHIMES initialization tables do not match the compiled candidate profile')
groups=d.get('groups',[])
if len(groups)!=25: raise SystemExit(f'expected 25 CHIMES initialization tables, found {len(groups)}')
for i,g in enumerate(groups,1):
    f=p.parent/g['path']
    if g.get('index')!=i or hashlib.sha256(f.read_bytes()).hexdigest()!=g.get('sha256'):
        raise SystemExit(f'25-group CHIMES initialization table hash mismatch: {f}')
print('SNRT_AGN25_CHIMES_GROUP_TABLES_OK count=25')
PY
fi
if [[ "$PROFILE" == agn_reference ]]; then
  test -f "$STELLAR_HISTORY" -a -f "$STELLAR_YIELDS"
fi
if [[ "$PROFILE" == effective_population ]]; then
  test -d "$STELLAR_INPUT"
  test -f "$STELLAR_INPUT/manifest.json" -a -f "$STELLAR_HISTORY"
  test -f "$STELLAR_YIELDS" -a -f "$STELLAR_NODES" -a -f "$STELLAR_BASE_SOURCE"
  test -f "$SED" -a -f "$STELLAR_BINDING"
  python3 -B - "$STELLAR_PACKAGE" "$SNIA_CONTRACT" "$CANDIDATE_25G_STELLAR_AGN" <<'PY'
import hashlib
import json
from pathlib import Path
import re
import sys

package, snia_contract = map(Path, sys.argv[1:3])
joint_25g = sys.argv[3] == "1"
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
if joint_25g:
    expected_profile = "snrt_agn25_joint_diag_20261102_v1"
    expected_edges = "439a3f1317ebdec5aa0054026b82d9ee0d059104a2ba461ca1aeb44d790ed58e"
    for label, text in (("base", base), ("effective", effective)):
        if not re.search(rf"\bgroup_set_id='{expected_profile}'", text):
            raise SystemExit(f"{label} stellar SED is not bound to the 25-group profile")
        if not re.search(rf"\bedges_sha256='{expected_edges}'", text):
            raise SystemExit(f"{label} stellar SED group-edge identity mismatch")
        if not re.search(r"\btransport_sha256='658583a70de61273c05b8bac039366cec3c351ed4baafb52e120dd6aefa8ba70'", text):
            raise SystemExit(f"{label} stellar SED transport-ledger identity mismatch")
    if manifest.get("group_set_id") != expected_profile or manifest.get("group_edges_sha256") != expected_edges:
        raise SystemExit("stellar source package manifest is not bound to the 25-group candidate")
print("M5_STELLAR_PACKAGE_BINDING_PASS model=PARSEC_mixed_lowmass imf=Kroupa population=effective_SSP status=reference_control")
PY
fi
GROUP_COUNT=9
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then GROUP_COUNT=25; fi
for ((i=1; i<=GROUP_COUNT; i++)); do test -f "$GROUP_DIR/group_$(printf '%02d' "$i").hdf5"; done
cp -- "$TEMPLATE" "$WORK/effective.nml"
cp -- "$SOURCE_IC" "$WORK/ic_sink"
if [[ -n "$CANDIDATE_RESTART_ORDERING" ]]; then
  python3 -B - "$WORK/effective.nml" "$CANDIDATE_RESTART_ORDERING" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
ordering = sys.argv[2]
text = path.read_text(encoding="ascii")
if ordering != "ksection":
    raise SystemExit("unsupported candidate restart ordering")
if len(re.findall(r"(?im)^&RUN_PARAMS\s*$", text)) != 1:
    raise SystemExit("effective namelist must contain exactly one RUN_PARAMS block")
if re.search(r"(?im)^\s*ordering\s*=", text):
    raise SystemExit("effective namelist already defines ordering")
text, count = re.subn(
    r"(?im)^(&RUN_PARAMS\s*)$",
    r"\1\n  ordering='ksection'",
    text,
    count=1,
)
if count != 1:
    raise SystemExit("failed to set ksection ordering")
path.write_text(text, encoding="ascii")
PY
  grep -Fqx "  ordering='ksection'" "$WORK/effective.nml"
fi
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
  if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 ]]; then
    # Keep the M8 stellar+AGN capture fixture physically normalized. The
    # non-cosmological Bondi implementation assumes G_code=1, so choose the
    # time unit from G_cgs*rho_unit*time_unit^2=1. Preserve the former gas
    # temperature and set a 91.687 km/s gas-sink relative speed to recover the
    # target lambda_Edd=0.1 after restoring the physical G factor. This is a
    # controlled source-coupling fixture, not a calibrated galaxy model.
    grep -Eq '^[[:space:]]*d_region\(1\)=1\.0d-3$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*p_region\(1\)=1\.0d-8$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*units_time=3\.1556926d16$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*units_length=3\.0857d24$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*boxlen=1\.0d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*x_center\(1\)=0\.5d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*y_center\(1\)=0\.5d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*z_center\(1\)=0\.5d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*length_x\(1\)=1\.0d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*length_y\(1\)=1\.0d0$' "$WORK/effective.nml"
    grep -Eq '^[[:space:]]*length_z\(1\)=1\.0d0$' "$WORK/effective.nml"
    python3 -B - "$WORK/ic_sink" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="ascii")
text, count = re.subn(
    r"(?m)^(\s*)4\.07708582420825d0(\s)",
    r"\g<1>4.07708582420825d-9\g<2>",
    text,
    count=1,
)
if count != 1:
    raise SystemExit("candidate IC did not contain the expected original M8 mass token")
path.write_text(text, encoding="ascii")
print("M5_CANDIDATE_IC_CODE_MASS_RESCALED factor=1e-9")
PY
    sed -i \
      -e 's/d_region(1)=1.0d-3/d_region(1)=2.9168505392559529d0/' \
      -e 's/p_region(1)=1.0d-8/p_region(1)=2.6437799050783008d-7/' \
      -e 's/u_region(1)=0.0d0/u_region(1)=8.927020847868951d-3/' \
      -e 's/units_time=3.1556926d16/units_time=3.0043481067814525d15/' \
      -e 's/boxlen=1.0d0/boxlen=0.001d0/' \
      -e 's/x_center(1)=0.5d0/x_center(1)=0.0005d0/' \
      -e 's/y_center(1)=0.5d0/y_center(1)=0.0005d0/' \
      -e 's/z_center(1)=0.5d0/z_center(1)=0.0005d0/' \
      -e 's/length_x(1)=1.0d0/length_x(1)=0.001d0/' \
      -e 's/length_y(1)=1.0d0/length_y(1)=0.001d0/' \
      -e 's/length_z(1)=1.0d0/length_z(1)=0.001d0/' \
      -e 's/dust_condensation=.1d0,.1d0,.1d0/dust_condensation=0d0,0d0,0d0/' \
      -e 's/dust_growth=.true./dust_growth=.false./' \
      -e 's/dust_sputtering=.true./dust_sputtering=.false./' \
      -e 's/dust_coagulation=.true./dust_coagulation=.false./' \
      -e 's/dust_shattering=.true./dust_shattering=.false./' \
      "$WORK/effective.nml"
    if ! grep -q '^[[:space:]]*dust_cooling=' "$WORK/effective.nml"; then
      awk 'index($0,"cooling_method=") {
             print
             print "  dust_cooling=\047chimes_neq_v1\047"
             print "  dust_iron_model=\047none\047"
             print "  dust_pah_model=\047none\047"
             next
           } {print}' \
        "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
      mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
    else
      grep -Fq "dust_cooling='chimes_neq_v1'" "$WORK/effective.nml"
      if ! grep -q '^[[:space:]]*dust_iron_model=' "$WORK/effective.nml" && \
         ! grep -q '^[[:space:]]*dust_pah_model=' "$WORK/effective.nml"; then
        awk 'index($0,"dust_cooling=") {
               print
               print "  dust_iron_model=\047none\047"
               print "  dust_pah_model=\047none\047"
               next
             } {print}' \
          "$WORK/effective.nml" > "$WORK/effective.nml.tmp"
        mv -- "$WORK/effective.nml.tmp" "$WORK/effective.nml"
      fi
    fi
    if [[ -n "${M5_25G_COURANT_FACTOR:-}" ]]; then
      [[ "$M5_25G_COURANT_FACTOR" == '5.0d-8' ]] || {
        echo 'Only the reviewed candidate interval factor 5.0d-8 is supported' >&2; exit 2;
      }
      test "$(grep -Ec '^[[:space:]]*courant_factor=0\.000003d0$' "$WORK/effective.nml")" -eq 1
      sed -i 's/^[[:space:]]*courant_factor=0\.000003d0$/  courant_factor=5.0d-8/' "$WORK/effective.nml"
    fi
    grep -Fqx '  d_region(1)=2.9168505392559529d0' "$WORK/effective.nml"
    grep -Fqx '  p_region(1)=2.6437799050783008d-7' "$WORK/effective.nml"
    grep -Fqx '  u_region(1)=8.927020847868951d-3' "$WORK/effective.nml"
    grep -Fqx '  units_time=3.0043481067814525d15' "$WORK/effective.nml"
    grep -Fqx '  units_length=3.0857d24' "$WORK/effective.nml"
    grep -Fqx '  boxlen=0.001d0' "$WORK/effective.nml"
    grep -Fqx '  x_center(1)=0.0005d0' "$WORK/effective.nml"
    grep -Fqx '  y_center(1)=0.0005d0' "$WORK/effective.nml"
    grep -Fqx '  z_center(1)=0.0005d0' "$WORK/effective.nml"
    grep -Fqx '  length_x(1)=0.001d0' "$WORK/effective.nml"
    grep -Fqx '  length_y(1)=0.001d0' "$WORK/effective.nml"
    grep -Fqx '  length_z(1)=0.001d0' "$WORK/effective.nml"
    grep -Fqx '  dust_condensation=0d0,0d0,0d0' "$WORK/effective.nml"
    grep -Fqx '  dust_growth=.false.' "$WORK/effective.nml"
    grep -Fqx '  dust_sputtering=.false.' "$WORK/effective.nml"
    grep -Fqx '  dust_coagulation=.false.' "$WORK/effective.nml"
    grep -Fqx '  dust_shattering=.false.' "$WORK/effective.nml"
    grep -Fq "dust_cooling='chimes_neq_v1'" "$WORK/effective.nml"
    grep -Fq "dust_iron_model='none'" "$WORK/effective.nml"
    grep -Fq "dust_pah_model='none'" "$WORK/effective.nml"
    python3 -B - "$WORK/effective.nml" "$WORK/ic_sink" <<'PY'
import math
import re
import sys
from pathlib import Path

namelist = Path(sys.argv[1]).read_text(encoding="ascii")
records = Path(sys.argv[2]).read_text(encoding="ascii").splitlines()
sink_record = next(line.split("!", 1)[0].strip() for line in records
                   if line.split("!", 1)[0].strip())
mass_code = float(sink_record.split()[0].replace("D", "E").replace("d", "e"))

def nml_real(name):
    match = re.search(rf"(?im)^\s*{re.escape(name)}\s*=\s*([^,\s/]+)", namelist)
    if match is None:
        raise SystemExit(f"candidate M8 mass check: missing {name}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))

mass_msun = mass_code * nml_real("units_density") * nml_real("units_length")**3 / 1.98847e33
if not math.isfinite(mass_msun) or abs(mass_msun / 1.0e8 - 1.0) > 0.01:
    raise SystemExit(f"candidate sink IC is not M8 in effective units: {mass_msun:.8e} Msun")
print(f"M5_CANDIDATE_BH_MASS_UNIT_PASS mass_msun={mass_msun:.8e}")
g_code = 6.67408e-8 * nml_real("units_density") * nml_real("units_time")**2
if not math.isfinite(g_code) or abs(g_code - 1.0) > 1e-6:
    raise SystemExit(f"candidate Bondi gravity units are inconsistent: G_code={g_code:.12g}")
velocity_kms = nml_real("u_region(1)") * nml_real("units_length") / nml_real("units_time") / 1e5
if not math.isfinite(velocity_kms) or abs(velocity_kms / 91.68747179493549 - 1.0) > 1e-6:
    raise SystemExit(f"candidate gas-sink velocity is inconsistent: {velocity_kms:.8g} km/s")
print(f"M5_CANDIDATE_BONDI_UNIT_PASS G_code={g_code:.12g} relative_velocity_kms={velocity_kms:.8g}")
PY
  fi
else
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 ]]; then
  test "$(grep -Ec '^[[:space:]]*nstepmax=3$' "$WORK/effective.nml")" -eq 1
  sed -i 's/^[[:space:]]*nstepmax=3$/  nstepmax=4/' "$WORK/effective.nml"
fi
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  dust_settings=$(printf '%s\n' \
    "  dust_mass_enabled=.true., dust_cooling='chimes_neq_v1'," \
    "  dust_mass_model='carbon_olivine_2size_v1', dust_material_model='dl01_composition_v1'," \
    "  dust_optics_model='d03_transport_v1', dust_size_radius_cm=1d-6,1d-5," \
    "  dust_iron_model='none', dust_pah_model='none'," \
    "  dust_size_density=2.2d0,3.8d0, dust_growth=.false., dust_sputtering=.false.," \
    "  dust_coagulation=.false., dust_shattering=.false., dust_sn_shocks=.false.," \
    "  dust_condensation=0d0,0d0,0d0")
else
  dust_settings=$(printf '%s\n' \
    "  dust_mass_enabled=.true., dust_cooling='chimes_neq_v1'," \
    "  dust_mass_model='carbon_olivine_2size_v1', dust_material_model='dl01_composition_v1'," \
    "  dust_optics_model='d03_transport_v1', dust_size_radius_cm=1d-6,1d-5," \
    "  dust_size_density=2.2d0,3.8d0, dust_growth=.true., dust_sputtering=.true.," \
    "  dust_coagulation=.true., dust_shattering=.true., dust_sn_shocks=.false.")
fi
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
if [[ -n "${M5_25G_COURANT_FACTOR:-}" && "$CANDIDATE_25G_DIAGNOSTIC" != 1 ]]; then
  echo 'M5_25G_COURANT_FACTOR is restricted to the fixed 25-group candidate diagnostic' >&2
  exit 2
fi
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  # The old AGN-only fixture violated the non-cosmological Bondi contract
  # G_code=G_cgs*rho_unit*time_unit^2=1. Use the independently exercised,
  # physically normalized M8/lambda_Edd=0.1 gas/BH state. This deliberately
  # changes the invalid fixture's physical initial conditions; it does not
  # claim equivalence to the earlier runs. Keep the 1 kpc physical box while
  # expressing it in Mpc code-length units.
  test "$(grep -Ec '^[[:space:]]*d_region\(1\)=1\.0d-3$' "$WORK/effective.nml")" -eq 1
  test "$(grep -Ec '^[[:space:]]*units_time=3\.1556926d13$' "$WORK/effective.nml")" -eq 1
  test "$(grep -Ec '^[[:space:]]*units_length=3\.0857d21$' "$WORK/effective.nml")" -eq 1
  sed -i \
    -e 's/d_region(1)=1.0d-3/d_region(1)=2.9168505392559529d0/' \
    -e 's/p_region(1)=1.0d-8/p_region(1)=2.6437799050783008d-7/' \
    -e 's/u_region(1)=0.0d0/u_region(1)=8.927020847868951d-3/' \
    -e 's/units_time=3.1556926d13/units_time=3.0043481067814525d15/' \
    -e 's/units_length=3.0857d21/units_length=3.0857d24/' \
    -e 's/boxlen=1.0d0/boxlen=0.001d0/' \
    -e 's/x_center(1)=0.5d0/x_center(1)=0.0005d0/' \
    -e 's/y_center(1)=0.5d0/y_center(1)=0.0005d0/' \
    -e 's/z_center(1)=0.5d0/z_center(1)=0.0005d0/' \
    -e 's/length_x(1)=1.0d0/length_x(1)=0.001d0/' \
    -e 's/length_y(1)=1.0d0/length_y(1)=0.001d0/' \
    -e 's/length_z(1)=1.0d0/length_z(1)=0.001d0/' \
    -e 's/eps_star=0.1d0/eps_star=0.0d0/' \
    "$WORK/effective.nml"
  python3 -B - "$WORK/ic_sink" <<'PY'
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="ascii")
text, count = re.subn(
    r"(?m)^(\s*)4\.07708582420825d0(\s)",
    r"\g<1>4.07708582420825d-9\g<2>",
    text,
    count=1,
)
if count != 1:
    raise SystemExit("candidate IC did not contain the expected original M8 mass token")
path.write_text(text, encoding="ascii")
print("M5_CANDIDATE_IC_CODE_MASS_RESCALED factor=1e-9")
PY
  if [[ -n "${M5_25G_COURANT_FACTOR:-}" ]]; then
    # This candidate-only short interval changes the coupled accretion/RT dt;
    # it never rescales accepted source energy or the Su spectral partition.
    [[ "$M5_25G_COURANT_FACTOR" == '5.0d-8' || \
       "$M5_25G_COURANT_FACTOR" == '2.5d-8' ]] || {
      echo 'Only reviewed candidate interval factors 5.0d-8 or 2.5d-8 are supported' >&2; exit 2;
    }
    test "$(grep -Ec '^[[:space:]]*courant_factor=0\.005d0$' "$WORK/effective.nml")" -eq 1
    sed -i "s/^[[:space:]]*courant_factor=0\.005d0$/  courant_factor=$M5_25G_COURANT_FACTOR/" \
      "$WORK/effective.nml"
    grep -Fqx "  courant_factor=$M5_25G_COURANT_FACTOR" "$WORK/effective.nml"
  fi
  grep -Fq 'dust_growth=.false., dust_sputtering=.false.' "$WORK/effective.nml"
  grep -Fq 'dust_coagulation=.false., dust_shattering=.false.' "$WORK/effective.nml"
  grep -Fq 'dust_condensation=0d0,0d0,0d0' "$WORK/effective.nml"
  grep -Fq "dust_iron_model='none', dust_pah_model='none'" "$WORK/effective.nml"
  grep -Fq 'd_region(1)=2.9168505392559529d0' "$WORK/effective.nml"
  grep -Fq 'p_region(1)=2.6437799050783008d-7' "$WORK/effective.nml"
  grep -Fq 'u_region(1)=8.927020847868951d-3' "$WORK/effective.nml"
  grep -Fq 'units_time=3.0043481067814525d15' "$WORK/effective.nml"
  grep -Fq 'units_length=3.0857d24' "$WORK/effective.nml"
  grep -Fq 'boxlen=0.001d0' "$WORK/effective.nml"
  grep -Fq 'x_center(1)=0.0005d0' "$WORK/effective.nml"
  grep -Fq 'y_center(1)=0.0005d0' "$WORK/effective.nml"
  grep -Fq 'z_center(1)=0.0005d0' "$WORK/effective.nml"
  grep -Fq 'length_x(1)=0.001d0' "$WORK/effective.nml"
  grep -Fq 'length_y(1)=0.001d0' "$WORK/effective.nml"
  grep -Fq 'length_z(1)=0.001d0' "$WORK/effective.nml"
  grep -Fq 'eps_star=0.0d0' "$WORK/effective.nml"
  printf 'AGN_RESTART_PHYSICAL_FIXTURE=normalized_M8_lambda0.1_from_448663; prior_invalid_fixture_not_physical_baseline\n' \
    | tee -a "$WORK/launch.txt"
  python3 -B - "$WORK/effective.nml" "$WORK/ic_sink" <<'PY'
import math
import re
import sys
from pathlib import Path

namelist = Path(sys.argv[1]).read_text(encoding="ascii")
sink_record = next(line.split("!", 1)[0].strip()
                   for line in Path(sys.argv[2]).read_text(encoding="ascii").splitlines()
                   if line.split("!", 1)[0].strip())
mass_code = float(sink_record.split()[0].replace("D", "E").replace("d", "e"))

def nml_real(name):
    match = re.search(rf"(?im)^\s*{re.escape(name)}\s*=\s*([^,\s/]+)", namelist)
    if match is None:
        raise SystemExit(f"candidate M8 mass/unit check: missing {name}")
    return float(match.group(1).replace("D", "E").replace("d", "e"))

mass_msun = mass_code * nml_real("units_density") * nml_real("units_length")**3 / 1.98847e33
if not math.isfinite(mass_msun) or abs(mass_msun / 1.0e8 - 1.0) > 0.01:
    raise SystemExit(f"candidate sink IC is not M8 in normalized units: {mass_msun:.8e} Msun")
g_code = 6.67408e-8 * nml_real("units_density") * nml_real("units_time")**2
if not math.isfinite(g_code) or abs(g_code - 1.0) > 1e-6:
    raise SystemExit(f"candidate Bondi gravity units are inconsistent: G_code={g_code:.12g}")
velocity_kms = nml_real("u_region(1)") * nml_real("units_length") / nml_real("units_time") / 1e5
if not math.isfinite(velocity_kms) or abs(velocity_kms / 91.68747179493549 - 1.0) > 1e-6:
    raise SystemExit(f"candidate gas-sink velocity is inconsistent: {velocity_kms:.8g} km/s")
if nml_real("eps_star") != 0.0:
    raise SystemExit("AGN-only restart diagnostic must suppress stellar births")
print(f"M5_AGN25_NORMALIZED_INPUT_PASS mass_msun={mass_msun:.8e} "
      f"G_code={g_code:.12g} relative_velocity_kms={velocity_kms:.8g} "
      "stellar_birth_efficiency=0")
PY
fi
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$PROFILE" == agn_reference ]]; then
  # var_region(1,14) is specific dust energy. Preserve the reference
  # composition's 10 K dust initialization when changing the velocity unit.
  test "$(grep -Ec '^[[:space:]]*var_region\(1,14\)=6\.78681493041457d-14$' "$WORK/effective.nml")" -eq 1
  sed -i 's/var_region(1,14)=6.78681493041457d-14/var_region(1,14)=6.151444748722844d-16/' \
    "$WORK/effective.nml"
  grep -Fq 'var_region(1,14)=6.151444748722844d-16' "$WORK/effective.nml"
fi
if [[ "$PROFILE" == effective_population ]]; then
  # This profile is bound to the Kroupa id=1 source package above.
  grep -Fq "population_model='effective_ssp'" "$WORK/effective.nml"
  grep -Eq '^[[:space:]]*imf_id=1$' "$WORK/effective.nml"
  grep -Eq '^[[:space:]]*input_imf_id=1([,[:space:]]|$)' "$STELLAR_HISTORY"
else
  # The AGN-only profile has no stellar SED and retains its separate Chabrier
  # id=2 feedback history. Do not validate it against the effective-profile IMF.
  grep -Fq "population_model='single_star_ssp'" "$WORK/effective.nml"
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
  agn_reference)
    if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 ]]; then
      grep -Eq '^[[:space:]]*nstepmax=4$' "$WORK/effective.nml"
    else
      grep -Eq '^[[:space:]]*nstepmax=3$' "$WORK/effective.nml"
    fi
    ;;
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
if [[ "$PROFILE" == effective_population && "$EFFECTIVE_AGN_ENABLED" == 0 ]]; then
  agn_true_count=$(grep -Eic '^[[:space:]]*agn[[:space:]]*=[[:space:]]*\.true\.' \
    "$WORK/effective.nml" || true)
  sink_agn_true_count=$(grep -Eic '^[[:space:]]*sink_AGN[[:space:]]*=[[:space:]]*\.true\.' \
    "$WORK/effective.nml" || true)
  if ((agn_true_count == 0 || sink_agn_true_count == 0)); then
    echo 'Cannot construct AGN-off effective_population input: expected AGN flags are missing' >&2
    exit 2
  fi
  sed -i -E \
    -e 's/^([[:space:]]*agn[[:space:]]*=[[:space:]]*)\.true\./\1.false./I' \
    -e 's/^([[:space:]]*sink_AGN[[:space:]]*=[[:space:]]*)\.true\./\1.false./I' \
    "$WORK/effective.nml"
  if grep -Eiq '^[[:space:]]*(agn|sink_AGN)[[:space:]]*=[[:space:]]*\.true\.' \
       "$WORK/effective.nml" || \
     [[ $(grep -Eic '^[[:space:]]*agn[[:space:]]*=[[:space:]]*\.false\.' "$WORK/effective.nml" || true) -eq 0 ]] || \
     [[ $(grep -Eic '^[[:space:]]*sink_AGN[[:space:]]*=[[:space:]]*\.false\.' "$WORK/effective.nml" || true) -eq 0 ]]; then
    echo 'Failed to disable both AGN and sink_AGN namelist flags' >&2
    exit 2
  fi
fi
sha256sum "$TEMPLATE" "$WORK/effective.nml" "$WORK/ic_sink" "$WORK/rank-bind.sh" \
  "$ROOT/simulation/snrt/runs/m5_source_active_dust_chimes_reference.sbatch" "$BINARY" \
  "$CHIMES_REPORTER_SOURCE" "$CHIMES_REPORTER" \
  "$CHIMES_INSTALL/libchimes.so" "$CHIMES_BUILD_RESULT" "$TABLE" "$GROUP_DIR/manifest.json" \
  "$GROUP_DIR"/group_*.hdf5 "$GROUP_CONTRACT" "$SECONDARY" "$DUST_CONTRACT" \
  "$CHIMES_ATOMIC_TABLE" "$CHIMES_MOLECULAR_TABLE" \
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
if [[ "$AGN_RESTART_TRACE" == 1 ]]; then export SNRT_AGN_RESTART_TRACE=1; else unset SNRT_AGN_RESTART_TRACE || true; fi
if [[ "$AGN_STENCIL_TRACE" == 1 ]]; then export SNRT_AGN_RESTART_STENCIL_TRACE=1; else unset SNRT_AGN_RESTART_STENCIL_TRACE || true; fi
if [[ -n "$AGN_DONOR_KEYS" ]]; then
  export SNRT_AGN_RESTART_DONOR_KEYS="$AGN_DONOR_KEYS"
  sha256sum "$AGN_DONOR_KEYS" >> "$WORK/inputs.sha256"
  printf 'AGN_DONOR_KEYS=%s\n' "$AGN_DONOR_KEYS" | tee -a "$WORK/launch.txt"
else
  unset SNRT_AGN_RESTART_DONOR_KEYS || true
fi
if [[ "$AGN_DONOR_SNAPSHOT" == 1 ]]; then
  export SNRT_AGN_RESTART_DONOR_SNAPSHOT=1
  export SNRT_AGN_RESTART_DONOR_OCT="$AGN_DONOR_OCT"
else
  unset SNRT_AGN_RESTART_DONOR_SNAPSHOT SNRT_AGN_RESTART_DONOR_OCT || true
fi
if [[ "$AGN_LEVEL_ENTRY_TRACE" == 1 ]]; then
  export SNRT_AGN_RESTART_LEVEL_ENTRY_TRACE=1
else
  unset SNRT_AGN_RESTART_LEVEL_ENTRY_TRACE || true
fi
printf 'AGN_RESTART_TRACE=%s\n' "$AGN_RESTART_TRACE" | tee -a "$WORK/launch.txt"
printf 'AGN_STENCIL_TRACE=%s\n' "$AGN_STENCIL_TRACE" | tee -a "$WORK/launch.txt"
printf 'AGN_DONOR_SNAPSHOT=%s OCT=%s\n' "$AGN_DONOR_SNAPSHOT" "$AGN_DONOR_OCT" | tee -a "$WORK/launch.txt"
printf 'AGN_LEVEL_ENTRY_TRACE=%s\n' "$AGN_LEVEL_ENTRY_TRACE" | tee -a "$WORK/launch.txt"
export SNRT_BACKEND="$TRANSPORT_BACKEND" SNRT_DUST_BACKEND="$MATERIAL_BACKEND" SNRT_AGN_MODEL=partition_reference_v1
export SNRT_REDUCED_C=0.01 SNRT_RT_LEVEL=3 SNRT_ALLOW_REFERENCE_CONTROL=1
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  case "$CANDIDATE_25G_RECEIVER" in
    hot_atomic_d03)
      export SNRT_SPECTRAL_MODEL=chimes_hot_atomic_d03_absorption_diag_v1
      unset SNRT_CHIMES_CAPTURE_ONE
      ;;
    cold_transition_d03)
      export SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1
      export SNRT_CHIMES_CAPTURE_ONE=1
      ;;
  esac
  export SNRT_ALLOW_25G_AGN_DIAGNOSTIC=1
  if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 ]]; then
    export SNRT_REQUIRE_25G_STELLAR_AGN=1
  else
    unset SNRT_REQUIRE_25G_STELLAR_AGN
  fi
  if [[ "$VERIFY_RESTART" == 1 ]]; then
    export SNRT_ALLOW_25G_AGN_DIAGNOSTIC_RESTART=1
    unset SNRT_HYDRO_ENTRY_STOP
  else
    unset SNRT_ALLOW_25G_AGN_DIAGNOSTIC_RESTART
    export SNRT_HYDRO_ENTRY_STOP=2
  fi
else
  export SNRT_SPECTRAL_MODEL=chimes_transition_d03_maxent128_fs2010_v1
  unset SNRT_ALLOW_25G_AGN_DIAGNOSTIC SNRT_ALLOW_25G_AGN_DIAGNOSTIC_RESTART \
    SNRT_REQUIRE_25G_STELLAR_AGN SNRT_HYDRO_ENTRY_STOP SNRT_CHIMES_CAPTURE_ONE
fi
export SNRT_GROUP_CONTRACT="$GROUP_CONTRACT" SNRT_SECONDARY_TABLE_CONTRACT="$SECONDARY"
if [[ -n "$SED" ]]; then export SNRT_STELLAR_SED="$SED"; else unset SNRT_STELLAR_SED; fi
export SNRT_DUST_CONTRACT="$DUST_CONTRACT"
export SNRT_CHIMES_BINARY="$BINARY" SNRT_CHIMES_MAIN_DATA="$TABLE" SNRT_CHIMES_GROUP_DIR="$GROUP_DIR"
export SNRT_CHIMES_BAND_TABLE="$CHIMES_ATOMIC_TABLE"
export SNRT_CHIMES_MOLECULAR_TABLE="$CHIMES_MOLECULAR_TABLE"
export SNRT_CHIMES_COMPUTE_BACKEND="$CHIMES_BACKEND"
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  export SNRT_PERF_DIAGNOSTIC=0
else
  export SNRT_PERF_DIAGNOSTIC=1
fi
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
printf 'EFFECTIVE_AGN_ENABLED=%s\n' "$EFFECTIVE_AGN_ENABLED" | tee -a "$WORK/launch.txt"
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  printf '25G_CANDIDATE_RECEIVER=%s; metallic_Fe_grains=disabled; profile_runtime_admissible=false\n' \
    "$CANDIDATE_25G_RECEIVER" | tee -a "$WORK/launch.txt"
  printf '25G_SOURCE_MODE=%s\n' \
    "$([[ "$CANDIDATE_25G_STELLAR_AGN" == 1 ]] && printf 'PARSEC_stellar_plus_Su_M8_AGN' || printf 'Su_M8_AGN_only')" \
    | tee -a "$WORK/launch.txt"
fi
printf 'TRANSPORT_BACKEND=%s MATERIAL_BACKEND=%s\n' "$TRANSPORT_BACKEND" "$MATERIAL_BACKEND" \
  | tee -a "$WORK/launch.txt"
printf 'CHIMES_cell_scheduler=%s\n' "$CELL_SCHEDULER" | tee -a "$WORK/launch.txt"
printf 'CHIMES_GPU_BROKERS=%s\n' "$CHIMES_GPU_BROKERS" | tee -a "$WORK/launch.txt"
printf 'CUDA_STREAMS=%s\n' "$CUDA_STREAMS" | tee -a "$WORK/launch.txt"
printf 'AMR_INTERFACE_CHECK=%s\n' "$AMR_INTERFACE" | tee -a "$WORK/launch.txt"
printf 'NSTEP_MAX=%s CAPTURE_FINAL_DUMP=%s\n' "$NSTEP_MAX" "$CAPTURE_FINAL_DUMP" | tee -a "$WORK/launch.txt"
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  printf 'CANDIDATE_COURANT_FACTOR=%s\n' \
    "$(awk -F= '/^[[:space:]]*courant_factor=/ {gsub(/[[:space:]]/,"",$2); print $2; exit}' "$WORK/effective.nml")" \
    | tee -a "$WORK/launch.txt"
fi
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
check_candidate_bondi_state() {
  local ledger=$1
  python3 -B - "$ledger" <<'PY'
import json
import math
import sys
from pathlib import Path

path = Path(sys.argv[1])
if not path.is_file():
    raise SystemExit(f"candidate restart branch is missing AGN state ledger: {path}")
rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
active = [row for row in rows if int(row.get("nstep_coarse", 0)) > 0 and
          float(row.get("eddington_ratio", 0.0)) > 0.0]
if not active:
    raise SystemExit(f"candidate restart branch has no positive AGN accretion sample: {path}")
state = min(active, key=lambda row: int(row["nstep_coarse"]))
mass = float(state.get("mass_msun", math.nan))
lambda_edd = float(state.get("eddington_ratio", math.nan))
if not math.isfinite(mass) or abs(mass / 1.0e8 - 1.0) > 0.01:
    raise SystemExit(f"candidate restart branch is not M8: mass_msun={mass:.8e}")
if not math.isfinite(lambda_edd) or abs(lambda_edd / 0.1 - 1.0) > 0.1:
    raise SystemExit(f"candidate restart branch left the lambda_Edd=0.1 target: {lambda_edd:.8e}")
print(f"M5_AGN25_RESTART_BONDI_STATE_PASS branch={path.parent.name} "
      f"nstep={state['nstep_coarse']} mass_msun={mass:.8e} lambda_edd={lambda_edd:.8e}")
PY
}
check_candidate_no_stellar_births() {
  local logfile=$1
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$PROFILE" == agn_reference ]] && \
     grep -Eq 'New star=[[:space:]]*[1-9]' "$logfile"; then
    echo "FAIL: AGN-only diagnostic formed a stellar particle; inspect $logfile" >&2
    return 1
  fi
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
check_oct_gravity_trace() {
  local logfile=$1 report=$2
  local stage_before=pre_gravity_sync stage_after=post_gravity_sync
  local trace_mode=gravity_sync preserve_storage=1
  if [[ "$AGN_LEVEL_ENTRY_TRACE" == 1 ]]; then
    stage_before=pre_level3_entry_i
    stage_after=pre_level3_negative_sync
    trace_mode=level3_operator_boundaries
    preserve_storage=0
  fi
  python3 -B - "$logfile" "$report" "$stage_before" "$stage_after" \
    "$trace_mode" "$preserve_storage" <<'PY'
import json, math, re, sys
from pathlib import Path

log_path, report_path = map(Path, sys.argv[1:3])
stage_before, stage_after, trace_mode = sys.argv[3:6]
preserve_storage = sys.argv[6] == '1'
level_entry_trace = trace_mode == 'level3_operator_boundaries'
stages = (stage_before, stage_after)
entry_stages = set()
observed_stages = set()
boundary_stages = (
    'post_level3_repartition', 'post_level3_output',
    'pre_level3_feedback_operators', 'post_level3_kinetic_feedback',
    'post_level3_agn_feedback', 'post_level3_create_sink',
    'post_level3_eeos', 'post_level3_feedback',
    'post_level3_particle_sort',
) if level_entry_trace else ()
meta = re.compile(
    r'^SNRT_AGN_RESTART_DONOR stage=(\S+) rank=(\d+) step=(\d+) '
    r'consumer_mask=(\d+) level=(\d+) cpu_map_grid=(\d+) copy_kind=(\d+) '
    r'grid=(\d+) cell_addr=(\d+) child=(\d+) '
    r'grid_key\(ix,iy,iz\)=\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)')
records = {}
for line_number, line in enumerate(log_path.open(encoding='utf-8', errors='replace'), 1):
    if 'SNRT_AGN_RESTART_DONOR missing' in line:
        raise SystemExit(f'{log_path}:{line_number}: donor snapshot reported a missing key')
    if not line.startswith('SNRT_AGN_RESTART_DONOR stage='):
        continue
    m = meta.match(line)
    if not m:
        raise SystemExit(f'{log_path}:{line_number}: malformed donor snapshot record')
    stage = m.group(1)
    if level_entry_trace:
        if stage.startswith(stage_before):
            if re.fullmatch(re.escape(stage_before) + r'[1-9][0-9]*', stage) is None:
                raise SystemExit(f'{log_path}:{line_number}: malformed level-entry index in stage {stage}')
            entry_stages.add(stage)
        elif stage in boundary_stages or stage == stage_after:
            observed_stages.add(stage)
        else:
            # The same bounded runner retains its established restart-trace
            # stages; this report selects only the requested AMR boundaries.
            continue
    elif stage not in stages:
        raise SystemExit(f'{log_path}:{line_number}: unexpected targeted stage {stage}')
    rank, step, mask, level, cpu_map_grid, copy_kind, grid, cell_addr, child, ix, iy, iz = map(int, m.groups()[1:])
    if (step, mask, level, ix, iy, iz) != (2, 3, 3, 3, 7, 3):
        raise SystemExit(f'{log_path}:{line_number}: donor key escaped the requested step/oct')
    if rank not in (1, 2) or copy_kind not in (0, 1) or child not in range(1, 9):
        raise SystemExit(f'{log_path}:{line_number}: invalid rank, storage kind, or oct child')
    try:
        state_txt, rest = line.split(' rho/mx/my/mz/E/eint/c2=', 1)[1].split(' f_xyz=', 1)
        force_txt, dt_txt = rest.split(' dtnew=', 1)
        state = [float(x.replace('D', 'E').replace('d', 'E')) for x in state_txt.split()]
        force = [float(x.replace('D', 'E').replace('d', 'E')) for x in force_txt.split()]
        dt = float(dt_txt.split()[0].replace('D', 'E').replace('d', 'E'))
    except (ValueError, IndexError) as exc:
        raise SystemExit(f'{log_path}:{line_number}: malformed hydro/force/dt values: {exc}')
    if len(state) != 7 or len(force) != 3 or not all(math.isfinite(x) for x in state + force + [dt]):
        raise SystemExit(f'{log_path}:{line_number}: nonfinite or incomplete snapshot values')
    if state[0] <= 0.0 or dt <= 0.0 or grid <= 0 or cell_addr <= 0:
        raise SystemExit(f'{log_path}:{line_number}: invalid density, timestep, or storage address')
    key = (stage, rank, child)
    if key in records:
        raise SystemExit(f'{log_path}:{line_number}: duplicate physical child {key}')
    records[key] = {
        'copy_kind': copy_kind, 'grid': grid, 'cell_addr': cell_addr,
        'cpu_map_grid': cpu_map_grid, 'rho_mom_E_eint_c2': state,
        'force_xyz': force, 'dtnew': dt,
    }
if level_entry_trace:
    if not entry_stages:
        raise SystemExit(f'{log_path}: no level-3 entry donor snapshots were emitted')
    stages = tuple(sorted(entry_stages,
                          key=lambda s: int(s[len(stage_before):])) +
                   list(boundary_stages) + [stage_after])
    if observed_stages != set(boundary_stages + (stage_after,)):
        raise SystemExit(f'{log_path}: incomplete operator-boundary trace; '
                         f'missing={sorted(set(boundary_stages + (stage_after,))-observed_stages)}; '
                         f'extra={sorted(observed_stages-set(boundary_stages + (stage_after,)))}')
expected = {(stage, rank, child) for stage in stages
            for rank in (1, 2) for child in range(1, 9)}
if set(records) != expected:
    raise SystemExit(f'{log_path}: incomplete oct trace; missing={sorted(expected-set(records))}; extra={sorted(set(records)-expected)}')
if preserve_storage:
    for rank in (1, 2):
        for child in range(1, 9):
            before = records[(stage_before, rank, child)]
            after = records[(stage_after, rank, child)]
            if any(before[k] != after[k] for k in ('copy_kind', 'grid', 'cell_addr')):
                raise SystemExit(f'{log_path}: storage identity changed across sync rank={rank} child={child}')
summary = {
    'status': 'PASS', 'log': str(log_path), 'rows': len(records),
    'target_oct': [3, 7, 3], 'coarse_step': 2,
    'trace_mode': trace_mode, 'stages': list(stages),
    'ranks': [1, 2], 'children_per_rank_per_stage': 8,
    'records': [
        {'stage': stage, 'rank': rank, 'child': child, **records[(stage, rank, child)]}
        for stage in stages
        for rank in (1, 2) for child in range(1, 9)
    ],
}
report_path.write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
print(f'SNRT_AGN_RESTART_OCT_TRACE_PASS rows={len(records)} mode={trace_mode} report={report_path}')
PY
}
gpu_reports_minimum=$MPI_RANKS
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" == 1 ]]; then
  # Restart equivalence does not certify CUDA throughput. This bounded fixture
  # can legitimately make no CHIMES GPU calls; record that separately.
  gpu_reports_minimum=0
  printf 'GPU_CHIMES_EXECUTION_GATE=not_required_for_restart_equivalence_only; GPU performance unverified\n' \
    | tee -a "$WORK/launch.txt"
fi
if [[ "$VERIFY_RESTART" == 1 ]]; then
  oct_trace_status=0
  hilbert_trace_status=0
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
  printf 'OUTPUT_POLICY=noutput1 aout2.0d0(future) tout1.0d30 nstepmax4 foutput2 fbackup1000000; 3 raw dumps <=330MiB budget (measured 109.4MiB/dump), checksummed then deleted after comparator\nRESTART_ORDERING_OVERRIDE=%s\nCONTINUOUS_NML=%s\nRESTART_NML=%s\n' \
    "${CANDIDATE_RESTART_ORDERING:-default}" \
    "$WORK/continuous/effective.nml" "$WORK/resumed/effective.nml" | tee -a "$WORK/launch.txt"
  run_m5_case "$WORK/continuous"
  grep -q 'Run completed' "$WORK/continuous/run.log"
  check_candidate_no_stellar_births "$WORK/continuous/run.log"
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$PROFILE" == agn_reference ]]; then
    check_candidate_bondi_state "$WORK/continuous/agn_coarse_state_v1.jsonl"
  fi
  if [[ "$AGN_DONOR_SNAPSHOT" == 1 ]] && ! check_oct_gravity_trace \
       "$WORK/continuous/run.log" "$WORK/oct-trace-continuous.json"; then
    oct_trace_status=1
  fi
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/continuous/run.log" "$gpu_reports_minimum"; fi
  test -f "$WORK/continuous/output_00001/COMPLETE"
  test -f "$WORK/continuous/output_00002/COMPLETE"
  ln -s "$WORK/continuous/output_00001" "$WORK/resumed/output_00001"
  run_m5_case "$WORK/resumed"
  grep -q 'Run completed' "$WORK/resumed/run.log"
  check_candidate_no_stellar_births "$WORK/resumed/run.log"
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$PROFILE" == agn_reference ]]; then
    check_candidate_bondi_state "$WORK/resumed/agn_coarse_state_v1.jsonl"
  fi
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
    resumed_steps=$(awk '$1=="Main" && $2=="step=" {value=$3} END {print value+0}' "$WORK/resumed/run.log")
    if [[ "$resumed_steps" != "$NSTEP_MAX" ]]; then
      echo "FAIL: candidate resumed run reached $resumed_steps steps; expected $NSTEP_MAX" >&2
      exit 3
    fi
  fi
  if [[ "$CHIMES_BACKEND" == hybrid ]]; then check_hybrid_process_reports "$WORK/resumed/run.log" "$gpu_reports_minimum"; fi
  grep -Fq 'SNRT_HDF5_RADIATION_RESTORE_PASS' "$WORK/resumed/run.log"
  if [[ "$CANDIDATE_RESTART_ORDERING" == ksection ]]; then
    grep -Fq 'HDF5: ksection tree restored' "$WORK/resumed/run.log"
    if grep -Fq 'VARIABLE NCPU / CROSS-ORDERING RESTART' "$WORK/resumed/run.log"; then
      echo 'FAIL: ksection restart unexpectedly entered the repartition path' >&2
      exit 3
    fi
  elif [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$MPI_RANKS" == 2 ]]; then
    grep -Eiq 'File ncpu[[:space:]]*=[[:space:]]*2.*ordering[[:space:]]*=[[:space:]]*hilbert' "$WORK/resumed/run.log" || hilbert_trace_status=1
    grep -Fq 'Load balance will be forced on first coarse step' "$WORK/resumed/run.log" || hilbert_trace_status=1
    grep -Fq 'HDF5: grid-balanced Hilbert bound_key rebuilt' "$WORK/resumed/run.log" || hilbert_trace_status=1
  fi
  if [[ "$AGN_DONOR_SNAPSHOT" == 1 ]] && ! check_oct_gravity_trace \
       "$WORK/resumed/run.log" "$WORK/oct-trace-resumed.json"; then
    oct_trace_status=1
  fi
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

if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 && "$VERIFY_RESTART" != 1 ]]; then
  # This fixture deliberately exits after its bounded source-bearing M5 step.
  grep -Fq 'SNRT bounded 25-group candidate source-bearing M5 transaction completed' run.log
else
  grep -q 'Run completed' run.log
fi
actual_main_steps=$(awk '$1=="Main" && $2=="step=" {value=$3} END {print value+0}' run.log)
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  ((actual_main_steps>=1)) || { echo 'FAIL: candidate fixture stopped before a coarse step' >&2; exit 3; }
  if [[ "$VERIFY_RESTART" == 1 && "$actual_main_steps" != "$NSTEP_MAX" ]]; then
    echo "FAIL: candidate continuous run reached $actual_main_steps steps; expected $NSTEP_MAX" >&2
    exit 3
  fi
  if [[ "$VERIFY_RESTART" != 1 ]]; then
    grep -Fq 'SNRT bounded 25-group candidate source-bearing M5 transaction completed' run.log
  fi
  grep -Fq 'SNRT spectral contract loaded: groups=25 status=candidate_explicit_sed' run.log
  grep -Fq 'snrt_agn25_joint_diag_20261102_v1_su_m8_mdotm1_fesc1' run.log
  if [[ "$CANDIDATE_25G_RECEIVER" == cold_transition_d03 ]]; then
    grep -Fq 'SNRT spectral absorption model=chimes_transition_d03_maxent128_fs2010_v1' run.log
    # CHIMES_CAPTURE records are emitted only when the solver rejects a cell.
    # Successful runs instead prove coupling through the live M_N/CHIMES
    # capture receipts checked below.
  fi
else
  if [[ "$actual_main_steps" != "$NSTEP_MAX" ]]; then
    echo "FAIL: reached main step $actual_main_steps; expected $NSTEP_MAX" >&2
    exit 3
  fi
fi
# A final-dump run exits below, so enforce real device work before that branch.
# Restart runs validate both continuous and resumed logs in their own path.
if [[ "$CHIMES_BACKEND" == hybrid && "$VERIFY_RESTART" != 1 && "$CANDIDATE_25G_DIAGNOSTIC" != 1 ]]; then
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  python3 -B - run.log <<'PY' | tee "$WORK/candidate-capture-check.txt"
import json,math,re,sys
from pathlib import Path
text=Path(sys.argv[1]).read_text(encoding='utf-8',errors='replace')
def rows(prefix,n):
    out=[]
    for line in text.splitlines():
        if prefix not in line: continue
        tail=line.split('=',1)[-1]
        vals=[float(x.replace('D','E').replace('d','E')) for x in tail.split()]
        if len(vals)!=n or not all(math.isfinite(x) for x in vals):
            raise SystemExit(f'malformed {prefix} receipt')
        out.append(vals)
    return out
ledger=rows('SNRT_MN_LEDGER source/absorbed/escaped/projection/residual/gas_heat/dust_heat/N=',8)
energy=rows('SNRT_MN_ENERGY absorbed/escaped/projection/residual/scale=',5)
capture=rows('SNRT_MN_CANDIDATE_CAPTURE gas_N/dust_N/gas_E_erg/dust_E_erg/chem_state_change_code_integral=',5)
source_re=re.compile(r'^\s*SNRT_MN_AGN_SOURCE\s+accepted_packets=(\d+)\s+photons/energy_erg=\s*([^\s]+)\s+([^\s]+)\s*$')
sources=[]
for line in text.splitlines():
    if 'SNRT_MN_AGN_SOURCE' not in line: continue
    m=source_re.match(line)
    if not m: raise SystemExit('malformed AGN N/E receipt')
    sources.append((int(m.group(1)),float(m.group(2).replace('D','E').replace('d','E')),
                    float(m.group(3).replace('D','E').replace('d','E'))))
if not ledger or not energy or not capture or not sources:
    raise SystemExit('missing candidate M5/CHIMES capture or closure receipts')
if not any(p>0 and n>0 and e>0 for p,n,e in sources):
    raise SystemExit('candidate source emitted no positive AGN photon/energy receipt')
max_rel=0.0
for x in ledger:
    source_n,absorbed_n,escaped_n,projection_n,residual_n,gas_heat,dust_heat,final_n=x
    scale=max(abs(source_n),abs(final_n),abs(absorbed_n)+abs(escaped_n),1e-300)
    max_rel=max(max_rel,abs(residual_n)/scale)
    if abs(residual_n)>1e-9*scale:
        raise SystemExit('candidate photon-number closure exceeds 1e-9 relative')
for x in energy:
    absorbed,escaped,projection,residual,scale=x
    # A dark pre-source M5 transaction has no closure scale. Admit it only
    # when every radiation-energy ledger term is exactly zero; source-bearing
    # transactions still require a positive scale and the strict residual gate.
    if scale==0:
        if any(v!=0 for v in (absorbed,escaped,projection,residual)):
            raise SystemExit('candidate dark transaction has a zero scale but nonzero energy ledger')
        continue
    if scale<0: raise SystemExit('candidate solver reported a negative energy closure scale')
    max_rel=max(max_rel,abs(residual)/scale)
    if abs(residual)>1e-9*scale:
        raise SystemExit('candidate radiation-energy closure exceeds 1e-9 relative')
gas_n=dust_n=gas_e=dust_e=chem_l1=0.0
for x in capture:
    gas_n+=x[0];dust_n+=x[1];gas_e+=x[2];dust_e+=x[3];chem_l1+=x[4]
if min(gas_n,dust_n,gas_e,dust_e,chem_l1)<=0:
    raise SystemExit('candidate did not produce positive gas+dust capture and a CHIMES state change')
coarse_path=Path('agn_coarse_state_v1.jsonl')
if not coarse_path.is_file():
    raise SystemExit('candidate is missing the physical AGN coarse-state ledger')
states=[json.loads(line) for line in coarse_path.read_text().splitlines() if line.strip()]
active=[row for row in states if int(row.get('nstep_coarse',0))>0 and
        float(row.get('eddington_ratio',0.0))>0.0]
if not active:
    raise SystemExit('candidate has no positive post-initialization AGN accretion state')
state=min(active,key=lambda row:int(row['nstep_coarse']))
lambda_edd=float(state['eddington_ratio'])
mass_msun=float(state['mass_msun'])
if not math.isfinite(lambda_edd) or abs(lambda_edd/0.1-1.0)>0.1:
    raise SystemExit(f'candidate does not reproduce M8 log10(mdot)=-1: lambda_Edd={lambda_edd:.8e}')
if not math.isfinite(mass_msun) or abs(mass_msun/1.0e8-1.0)>0.01:
    raise SystemExit(f'candidate AGN ledger is not M8: mass_msun={mass_msun:.8e}')
print('SNRT_AGN25_LIVE_CAPTURE_SMOKE_PASS '
      f'gas_photons={gas_n:.8e} dust_photons={dust_n:.8e} '
      f'gas_energy_erg={gas_e:.8e} dust_energy_erg={dust_e:.8e} '
      f'chem_state_change_code_integral={chem_l1:.8e} '
      f'lambda_edd={lambda_edd:.8e} mass_msun={mass_msun:.8e} '
      f'radiation_closure_relative={max_rel:.3e}')
PY
elif [[ "$CELL_SCHEDULER" == level_queue ]]; then
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" != 1 ]] && ((momentum_assert_cells<=0)); then
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
  comparison_status=0
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
    python3 -B "$ROOT/simulation/snrt/tests/p0_hdf5_stellar_restart.py" \
      --compare-m5-candidate-restart "$WORK/continuous/output_00002/data_00002.h5" \
        "$WORK/resumed/output_00002/data_00002.h5" \
      --source-checkpoint "$WORK/continuous/output_00001/data_00001.h5" \
      --chimes-proto-header "$CHIMES_INSTALL/src/chimes_proto.h" \
      --candidate-restart-rtol 1.0e-10 \
      > "$WORK/restart-comparison.json" || comparison_status=$?
  else
    python3 -B "$ROOT/simulation/snrt/tests/p0_hdf5_stellar_restart.py" \
      --compare-m5-coupled "$WORK/continuous/output_00002/data_00002.h5" \
        "$WORK/resumed/output_00002/data_00002.h5" \
      --source-checkpoint "$WORK/continuous/output_00001/data_00001.h5" \
      > "$WORK/restart-comparison.json" || comparison_status=$?
  fi
  # Keep logs, comparator JSON, and hashes for audit; remove the three large
  # simulation dumps immediately after comparison, including on a failed gate.
  rm -rf -- "$WORK/continuous/output_00001" "$WORK/continuous/output_00002" \
    "$WORK/resumed/output_00001" "$WORK/resumed/output_00002"
  if ((oct_trace_status != 0)); then
    echo 'FAIL: targeted oct force/timestep trace did not cover all required records; raw HDF5 dumps were removed after evaluation' >&2
  fi
  if ((hilbert_trace_status != 0)); then
    echo 'FAIL: resumed log did not positively confirm the expected Hilbert reconstruction; raw HDF5 dumps were removed after evaluation' >&2
  fi
  if ((comparison_status != 0)); then
    echo "FAIL: restart comparator exited $comparison_status; raw HDF5 dumps were removed after evaluation" >&2
    exit "$comparison_status"
  fi
  if ((oct_trace_status != 0)); then
    exit 6
  fi
  if ((hilbert_trace_status != 0)); then
    exit 7
  fi
  if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
    printf 'M5_AGN25_CANDIDATE_HDF5_RESTART_PASS receiver=%s CHIMES_backend=%s restart_rtol=1.0e-10 raw_dumps_retained=0 raw_dump_hashes_retained=3 runtime_admissible=false pending_driver_evaluation=true\n' \
      "$CANDIDATE_25G_RECEIVER" "$CHIMES_BACKEND" | tee "$WORK/result.txt"
  else
    printf 'M5_COUPLED_RESTART_EXECUTION_PASS profile=effective_population CHIMES_backend=%s raw_dumps_retained=0 raw_dump_hashes_retained=3 pending_driver_evaluation=true\n' \
      "$CHIMES_BACKEND" | tee "$WORK/result.txt"
  fi
  exit 0
fi
agn_evidence=not_requested
agn_packets=0
agn_photons=0
agn_energy_erg=0
# The bounded 25-group fixture injects a fixed, auditable AGN packet to test
# spectral capture; it does not require the live accretion/thermal-feedback
# sidecar that full AGN reference runs must produce. Its dedicated capture
# gate below still requires a positive packet receipt and closes M5 ledgers.
if [[ "$CANDIDATE_25G_DIAGNOSTIC" != 1 && ( "$PROFILE" == effective_population || "$PROFILE" == agn_reference ) ]]; then
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
if [[ "$PROFILE" == effective_population && \
      ( "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference || "$CANDIDATE_25G_STELLAR_AGN" == 1 ) ]]; then
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
agn_energy_re = re.compile(
    r"^\s*SNRT_MN_AGN_ENERGY\s+input/injected/unrepresented_erg=\s*"
    r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s*$"
)
commit_count = 0
agn_receipts = []
stellar_receipts = []
agn_energy_receipts = []
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
        elif "SNRT_MN_AGN_ENERGY" in line:
            match = agn_energy_re.match(line)
            if match is None:
                raise SystemExit("FAIL: malformed M5 AGN energy partition receipt")
            values = tuple(
                float(value.replace("D", "E").replace("d", "E"))
                for value in match.groups()
            )
            if any(not math.isfinite(value) or value < 0.0 for value in values):
                raise SystemExit("FAIL: invalid M5 AGN energy partition receipt")
            agn_energy_receipts.append(values)
if (not commit_count or len(agn_receipts) != commit_count or
        len(stellar_receipts) != commit_count or len(agn_energy_receipts) != commit_count):
    raise SystemExit("FAIL: M5 source receipts are missing or misaligned")

stellar_total = 0.0
agn_total = 0.0
overlap_commits = 0
max_agn_energy_partition_rel = 0.0
for (packets, agn_photons, agn_energy), (stars, stellar_photons, stellar_energy), (
    agn_input_energy, agn_injected_energy, agn_unrepresented_energy
) in zip(
    agn_receipts, stellar_receipts, agn_energy_receipts
):
    energy_partition_error = abs(
        agn_input_energy - agn_injected_energy - agn_unrepresented_energy
    )
    energy_partition_scale = max(agn_input_energy, agn_injected_energy, agn_unrepresented_energy)
    if energy_partition_error > 1.0e-10 * max(energy_partition_scale, 1.0e-300):
        raise SystemExit("FAIL: AGN input != represented + unrepresented energy")
    if energy_partition_scale > 0.0:
        max_agn_energy_partition_rel = max(
            max_agn_energy_partition_rel, energy_partition_error / energy_partition_scale
        )
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
    f"thermal_agn_photons={agn_total:.8e} "
    f"agn_energy_partition_max_rel={max_agn_energy_partition_rel:.3e}"
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
if [[ "$CANDIDATE_25G_DIAGNOSTIC" == 1 ]]; then
  if [[ "$CANDIDATE_25G_STELLAR_AGN" == 1 ]]; then
    printf 'M5_AGN25_CANDIDATE_CAPTURE_PASS profile=%s source_scope=PARSEC_stellar_plus_Su_M8_AGN receiver=%s M5_order=5 capture_check=candidate-capture-check.txt combined_source_attribution=combined-source-attribution.txt raw_outputs=0 runtime_admissible=false\n' \
      "$PROFILE" "$CANDIDATE_25G_RECEIVER" | tee "$WORK/result.txt"
  else
    printf 'M5_AGN25_CANDIDATE_CAPTURE_PASS profile=%s receiver=%s M5_order=5 capture_check=candidate-capture-check.txt source_receipts=run.log raw_outputs=0 runtime_admissible=false\n' \
      "$PROFILE" "$CANDIDATE_25G_RECEIVER" | tee "$WORK/result.txt"
  fi
  exit 0
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
if [[ "$PROFILE" == effective_population && \
      "$EFFECTIVE_AGN_IC_MODE" == thermal_seed_reference ]]; then
  source_scope=stellar_plus_thermal_AGN_reference
fi
printf 'M5_SOURCE_ACTIVE_DUST_CHIMES_%s_PASS profile=%s M5_order=5 source_scope=%s AGN_state=%s AGN_packets=%s AGN_photons=%s AGN_photon_energy_erg=%s M5_total_source_ledger=nonzero dust_mass=enabled CHIMES_backend=%s momentum_assert_cells=%s reference_controls=NONPRODUCTION raw_outputs=0\n' \
  "${CHIMES_LIBRARY_MODE^^}" "$PROFILE" "$source_scope" "$agn_evidence" "$agn_packets" "$agn_photons" "$agn_energy_erg" "$CHIMES_BACKEND" "$momentum_assert_cells" \
  | tee "$WORK/result.txt"
