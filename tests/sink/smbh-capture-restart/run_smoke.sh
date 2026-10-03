#!/usr/bin/env bash
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  echo "usage: $0 --run-dir NEW_DIR --binary RAMSES --bytes-per-output BYTES [--reserve-bytes BYTES] --execute" >&2
  exit 2
}

run_dir=
binary=
bytes_per_output=
reserve_bytes=1073741824
execute=false
while (($#)); do
  case "$1" in
    --run-dir) run_dir=${2:?}; shift 2 ;;
    --binary) binary=${2:?}; shift 2 ;;
    --bytes-per-output) bytes_per_output=${2:?}; shift 2 ;;
    --reserve-bytes) reserve_bytes=${2:?}; shift 2 ;;
    --execute) execute=true; shift ;;
    *) usage ;;
  esac
done

[[ -n $run_dir && -n $binary && -n $bytes_per_output ]] || usage
[[ $execute == true ]] || { echo 'CAPTURE-RESTART: --execute is required; nothing launched' >&2; exit 2; }
[[ -n ${SLURM_JOB_ID:-} ]] || {
  echo 'CAPTURE-RESTART: a live Slurm allocation (SLURM_JOB_ID) is required; refusing login-node execution' >&2
  exit 2
}
command -v srun >/dev/null 2>&1 || { echo 'CAPTURE-RESTART: srun is unavailable' >&2; exit 1; }
binary=$(realpath -- "$binary")
[[ -x $binary ]] || { echo "CAPTURE-RESTART: binary is not executable: $binary" >&2; exit 1; }
python3 "$here/prepare_and_validate.py" prepare "$run_dir" \
  --bytes-per-output "$bytes_per_output" --reserve-bytes "$reserve_bytes"
run_dir=$(realpath -- "$run_dir")

source_root=$(git -C "$here" rev-parse --show-toplevel)
{
  echo "slurm_job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "binary=$binary"
  echo "binary_sha256=$(sha256sum "$binary" | awk '{print $1}')"
  echo "source_root=$source_root"
  echo "source_commit=$(git -C "$source_root" rev-parse HEAD)"
  echo "source_branch=$(git -C "$source_root" branch --show-current)"
  echo 'source_status_begin'
  git -C "$source_root" status --porcelain --untracked-files=normal
  echo 'source_status_end'
  sha256sum "$run_dir/stage1.nml" "$run_dir/stage2.nml" "$run_dir/ic_sink"
} | tee "$run_dir/preflight_provenance.txt"

(
  cd "$run_dir"
  srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    "$binary" stage1.nml >stage1.log 2>&1
  [[ -f output_00001/COMPLETE ]]
  [[ -f output_00001/SMBH_CAPTURE_LINEAGE ]]
  srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
    "$binary" stage2.nml >stage2.log 2>&1
)
python3 "$here/prepare_and_validate.py" postcheck "$run_dir"
