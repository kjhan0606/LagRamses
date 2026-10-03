#!/usr/bin/env bash
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  echo "prepare: $0 --prepare --run-dir NEW_DIR --binary RAMSES --source-tree CLEAN_TREE --expected-source-commit COMMIT --bytes-per-output BYTES [--reserve-bytes BYTES]" >&2
  echo "execute: $0 --execute --run-dir PREPARED_DIR --binary RAMSES --manifest-sha256 SHA256" >&2
  echo "manual:  $0 --manual-lageunha --run-dir PREPARED_DIR --binary RAMSES --manifest-sha256 SHA256" >&2
  exit 2
}

run_dir=
binary=
bytes_per_output=
reserve_bytes=1073741824
execute=false
prepare=false
manual_lageunha=false
expected_source_commit=
source_tree=
manifest_sha256=
while (($#)); do
  case "$1" in
    --run-dir) run_dir=${2:?}; shift 2 ;;
    --binary) binary=${2:?}; shift 2 ;;
    --bytes-per-output) bytes_per_output=${2:?}; shift 2 ;;
    --reserve-bytes) reserve_bytes=${2:?}; shift 2 ;;
    --expected-source-commit) expected_source_commit=${2:?}; shift 2 ;;
    --source-tree) source_tree=${2:?}; shift 2 ;;
    --manifest-sha256) manifest_sha256=${2:?}; shift 2 ;;
    --prepare) prepare=true; shift ;;
    --execute) execute=true; shift ;;
    --manual-lageunha) manual_lageunha=true; shift ;;
    *) usage ;;
  esac
done

[[ -n $run_dir && -n $binary ]] || usage
mode_count=0
[[ $prepare == true ]] && ((mode_count+=1))
[[ $execute == true ]] && ((mode_count+=1))
[[ $manual_lageunha == true ]] && ((mode_count+=1))
[[ $mode_count == 1 ]] || { echo 'CAPTURE-RESTART: choose exactly one execution mode' >&2; exit 2; }
binary=$(realpath -- "$binary")
[[ -x $binary ]] || { echo "CAPTURE-RESTART: binary is not executable: $binary" >&2; exit 1; }
if [[ $prepare == true ]]; then
  [[ -n $bytes_per_output && -n $expected_source_commit && -n $source_tree ]] || usage
  python3 "$here/prepare_and_validate.py" prepare "$run_dir" \
    --binary "$binary" --expected-source-commit "$expected_source_commit" \
    --source-tree "$source_tree" \
    --bytes-per-output "$bytes_per_output" --reserve-bytes "$reserve_bytes"
  exit 0
fi
host=$(hostname -s)
if [[ $execute == true ]]; then
  python3 "$here/prepare_and_validate.py" execution-context --execution-mode slurm \
    --hostname "$host" --slurm-job-id "${SLURM_JOB_ID:-}"
  command -v srun >/dev/null 2>&1 || { echo 'CAPTURE-RESTART: srun is unavailable' >&2; exit 1; }
  execution_mode=slurm
else
  python3 "$here/prepare_and_validate.py" execution-context --execution-mode manual-lageunha \
    --hostname "$host"
  execution_mode=manual-lageunha
fi
[[ -n $manifest_sha256 ]] || usage
run_dir=$(realpath -- "$run_dir")
python3 "$here/prepare_and_validate.py" launch-check "$run_dir" --binary "$binary" \
  --manifest-sha256 "$manifest_sha256"

source_root=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_tree"])' \
  "$run_dir/preflight_manifest.json")
source_commit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["expected_source_commit"])' \
  "$run_dir/preflight_manifest.json")
{
  echo "execution_mode=$execution_mode"
  echo "slurm_job_id=${SLURM_JOB_ID:-none}"
  echo "host=$(hostname)"
  echo "binary=$binary"
  echo "binary_sha256=$(sha256sum "$binary" | awk '{print $1}')"
  echo "manifest=$run_dir/preflight_manifest.json"
  echo "manifest_sha256=$manifest_sha256"
  echo "source_root=$source_root"
  echo "source_commit=$source_commit"
  echo "source_branch=$(git -C "$source_root" branch --show-current)"
  echo 'source_status_begin'
  git -C "$source_root" status --porcelain --untracked-files=normal
  echo 'source_status_end'
  sha256sum "$run_dir/stage1.nml" "$run_dir/stage2.nml" "$run_dir/ic_sink"
  sha256sum "$run_dir/preflight_manifest.json"
} | tee "$run_dir/preflight_provenance.txt"

(
  cd "$run_dir"
  if [[ $execution_mode == slurm ]]; then
    srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
      "$binary" stage1.nml >stage1.log 2>&1
  else
    "$binary" stage1.nml >stage1.log 2>&1
  fi
  [[ -f output_00001/COMPLETE ]]
  [[ -f output_00001/SMBH_CAPTURE_LINEAGE ]]
  if [[ $execution_mode == slurm ]]; then
    srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
      "$binary" stage2.nml >stage2.log 2>&1
  else
    "$binary" stage2.nml >stage2.log 2>&1
  fi
)
python3 "$here/prepare_and_validate.py" postcheck "$run_dir"
