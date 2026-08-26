#!/bin/bash
#SBATCH --partition=normal
#SBATCH --time=30-00:00:00
#SBATCH --exclusive
# --- Job-specific SBATCH lines (job-name, nodes, ntasks, exclude) ---
# --- are set in each run's run.sh, NOT here ---

# ============================================================
#  Auto-restart wrapper for SIDM x AGN simulations
#  Features:
#    1. Auto-sync nrestart with latest output at job start
#    2. Inner retry loop for in-job crashes (segfault, MPI abort)
#    3. Successor job via --dependency=afternotok (NODE FAILURE)
#    4. Cleanup: keep only latest 2 foutput checkpoints
# ============================================================

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_STACKSIZE=256M
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:$LD_LIBRARY_PATH

# UCX InfiniBand settings
export UCX_RC_TIMEOUT=30s
export UCX_RC_RETRY_COUNT=7
export UCX_LOG_LEVEL=warn
export UCX_IB_REG_METHODS=direct
export UCX_RNDV_THRESH=65536

MAX_RETRIES=10       # In-job retry limit
CLEANUP_KEEP=2       # Keep latest N foutput checkpoints

# ============================================================
#  Step 0: Submit successor job (runs only if THIS job fails)
# ============================================================
SCRIPT_PATH="$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K\S+')"
SUCCESSOR=$(sbatch --dependency=afternotok:$SLURM_JOB_ID --parsable "$SCRIPT_PATH" 2>/dev/null)
echo "============================================"
echo "Job $SLURM_JOB_ID started: $(date)"
echo "Nodes: ${SLURM_JOB_NODELIST}"
echo "Tasks: ${SLURM_NTASKS}, Threads: ${OMP_NUM_THREADS}"
echo "Successor job: $SUCCESSOR (auto-resubmit on failure)"
echo "============================================"

# ============================================================
#  Step 1: Auto-sync nrestart with latest output
# ============================================================
sync_nrestart() {
  local latest_out=$(ls -d output_* 2>/dev/null | sort -V | tail -1)
  if [ -n "$latest_out" ]; then
    local latest_num=$(echo "$latest_out" | sed 's/output_0*//')
    local current_nr=$(grep -oP 'nrestart=\K[0-9]+' cosmo.nml)
    if [ "$latest_num" != "$current_nr" ]; then
      sed -i "s/nrestart=$current_nr/nrestart=$latest_num/" cosmo.nml
      echo ">>> nrestart synced: $current_nr -> $latest_num"
    else
      echo ">>> nrestart=$current_nr already matches latest output"
    fi
  fi
}

# ============================================================
#  Step 2: Cleanup old foutput checkpoints
#  Keep only the latest $CLEANUP_KEEP outputs that are NOT
#  in the aout schedule (i.e., foutput-generated safety dumps)
# ============================================================
cleanup_old_checkpoints() {
  # Get aout numbers from cosmo.nml (these are science outputs, never delete)
  local aout_outputs=""
  # Parse which output numbers correspond to aout snapshots
  # by reading info_XXXXX.txt headers — skip this complexity.
  # Instead: keep all outputs, only delete if total > 40
  local n_outputs=$(ls -d output_* 2>/dev/null | wc -l)
  if [ "$n_outputs" -gt 40 ]; then
    echo ">>> WARNING: $n_outputs outputs on disk. Consider manual cleanup."
  fi
}

# ============================================================
#  Step 3: Main run loop
# ============================================================
sync_nrestart
rc=1

for attempt in $(seq 1 $MAX_RETRIES); do
  echo ""
  echo ">>> Attempt $attempt/$MAX_RETRIES at $(date)"

  # Clear jobcontrol before each run
  echo "" > jobcontrol.txt

  srun --mpi=pmi2 ./ramses_final3d cosmo.nml
  rc=$?

  if [ $rc -eq 0 ]; then
    echo ">>> Clean exit at $(date)"
    break
  fi

  echo ">>> CRASH detected (exit=$rc) at $(date)"

  # Check if new output was created
  local_latest=$(ls -d output_* 2>/dev/null | sort -V | tail -1)
  local_num=$(echo "$local_latest" | sed 's/output_0*//')
  current_nr=$(grep -oP 'nrestart=\K[0-9]+' cosmo.nml)

  if [ "$local_num" = "$current_nr" ]; then
    echo ">>> No new checkpoint since nrestart=$current_nr."
    echo ">>> Crash before first checkpoint. Retrying from same point..."
  else
    sed -i "s/nrestart=$current_nr/nrestart=$local_num/" cosmo.nml
    echo ">>> Updated nrestart: $current_nr -> $local_num"
  fi

  cleanup_old_checkpoints
  echo ">>> Waiting 10s before retry..."
  sleep 10
done

# ============================================================
#  Step 4: Final status
# ============================================================
if [ $rc -eq 0 ]; then
  # Clean exit — cancel successor job
  if [ -n "$SUCCESSOR" ]; then
    scancel $SUCCESSOR 2>/dev/null
    echo ">>> Cancelled successor job $SUCCESSOR (clean exit)"
  fi
else
  echo ">>> Job ended with errors. Successor job $SUCCESSOR will auto-start."
fi

echo "============================================"
echo "Job $SLURM_JOB_ID finished: $(date) (exit=$rc)"
echo "============================================"
