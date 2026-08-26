#!/bin/bash
#SBATCH --job-name=sidm5_128
#SBATCH --partition=normal
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=2
#SBATCH --time=30-00:00:00
#SBATCH --output=run_%j.log
#SBATCH --error=run_%j.err
#SBATCH --exclusive
#SBATCH --exclude=grammar[007,022,026-027,039,066,069,084,092,095,100]

# ============================================================
#  Auto-restart wrapper for SIDM x AGN simulations
#  Features:
#    1. Auto-sync nrestart with latest output at job start
#    2. Inner retry loop for in-job crashes (segfault, MPI abort)
#    3. Successor job via --dependency=afternotok (NODE FAILURE)
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

MAX_RETRIES=10

# --- Step 0: Submit successor job (runs only if THIS job fails) ---
SCRIPT_DIR=$(scontrol show job $SLURM_JOB_ID | grep -oP 'WorkDir=\K\S+')
SUCCESSOR=$(sbatch --dependency=afternotok:$SLURM_JOB_ID --parsable "$SCRIPT_DIR/run.sh" 2>/dev/null)

echo "============================================"
echo "Job $SLURM_JOB_ID started: $(date)"
echo "Nodes: ${SLURM_JOB_NODELIST}"
echo "Tasks: ${SLURM_NTASKS}, Threads: ${OMP_NUM_THREADS}"
echo "Successor job: $SUCCESSOR (auto-resubmit on failure)"
echo "============================================"

# --- Step 1: Auto-sync nrestart with latest output ---
latest_out=$(ls -d output_* 2>/dev/null | sort -V | tail -1)
if [ -n "$latest_out" ]; then
  latest_num=$(echo "$latest_out" | sed 's/output_0*//; s/\..*//')
  current_nr=$(grep -oP 'nrestart=\K[0-9]+' cosmo.nml)
  if [ "$latest_num" != "$current_nr" ]; then
    sed -i "s/nrestart=$current_nr/nrestart=$latest_num/" cosmo.nml
    echo ">>> nrestart synced: $current_nr -> $latest_num"
  else
    echo ">>> nrestart=$current_nr OK"
  fi
fi

# --- Step 2: Main run loop ---
rc=1
for attempt in $(seq 1 $MAX_RETRIES); do
  echo ""
  echo ">>> Attempt $attempt/$MAX_RETRIES at $(date)"
  echo "" > jobcontrol.txt

  srun --mpi=pmi2 ./ramses_final3d cosmo.nml
  rc=$?

  [ $rc -eq 0 ] && echo ">>> Clean exit at $(date)" && break

  echo ">>> CRASH (exit=$rc) at $(date)"

  # Update nrestart from latest output
  new_out=$(ls -d output_* 2>/dev/null | sort -V | tail -1)
  new_num=$(echo "$new_out" | sed 's/output_0*//; s/\..*//')
  old_nr=$(grep -oP 'nrestart=\K[0-9]+' cosmo.nml)

  if [ "$new_num" = "$old_nr" ]; then
    echo ">>> No new checkpoint. Retrying from nrestart=$old_nr..."
  else
    sed -i "s/nrestart=$old_nr/nrestart=$new_num/" cosmo.nml
    echo ">>> nrestart: $old_nr -> $new_num"
  fi

  sleep 10
done

# --- Step 3: Cancel successor if clean exit ---
if [ $rc -eq 0 ]; then
  [ -n "$SUCCESSOR" ] && scancel $SUCCESSOR 2>/dev/null
  echo ">>> Cancelled successor $SUCCESSOR (clean exit)"
else
  echo ">>> Successor $SUCCESSOR will auto-start"
fi

echo "============================================"
echo "Job $SLURM_JOB_ID finished: $(date) (rc=$rc)"
echo "============================================"
