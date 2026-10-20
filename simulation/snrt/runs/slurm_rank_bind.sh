#!/usr/bin/env bash
# Same explicit rank/OMP binding used by the verified Grackle replay.
# Invoke with srun --cpu-bind=none: select disjoint slices of its allocation.
set -euo pipefail
allowed=$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)
IFS=, read -ra ranges <<< "$allowed"
cpus=()
for item in "${ranges[@]}"; do
    if [[ $item =~ ^([0-9]+)-([0-9]+)$ ]]; then
        for ((cpu=${BASH_REMATCH[1]}; cpu<=${BASH_REMATCH[2]}; cpu++)); do
            cpus+=("$cpu")
        done
    elif [[ $item =~ ^[0-9]+$ ]]; then
        cpus+=("$item")
    else
        echo "Invalid allowed CPU list: $allowed" >&2; exit 1
    fi
done
width=${SLURM_CPUS_PER_TASK:?}
local_rank=${SLURM_LOCALID:?}
[[ $width =~ ^[1-9][0-9]*$ && $local_rank =~ ^[0-9]+$ ]]
start=$((local_rank*width))
((start+width<=${#cpus[@]})) || { echo 'Insufficient allowed CPUs for disjoint rank binding' >&2; exit 1; }
selected=()
places=()
for ((i=start; i<start+width; i++)); do
    selected+=("${cpus[i]}")
    places+=("{${cpus[i]}}")
done
cpu_list=$(IFS=,; echo "${selected[*]}")
export OMP_PLACES=$(IFS=,; echo "${places[*]}")
export I_MPI_PIN=0
unset KMP_AFFINITY
echo "RANK_BIND local_rank=$local_rank cpus=$cpu_list places=$OMP_PLACES" >&2
exec taskset -c "$cpu_list" "$@"
