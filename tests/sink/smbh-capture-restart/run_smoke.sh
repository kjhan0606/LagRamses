#!/usr/bin/env bash
set -euo pipefail

here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
  echo "prepare: $0 --prepare --run-dir NEW_DIR --binary RAMSES --source-tree CLEAN_TREE --expected-source-commit COMMIT --yield-table-source FILE --bytes-per-output BYTES [--reserve-bytes BYTES] [--mpi-ranks 1|2]" >&2
  echo "execute: $0 --execute --run-dir PREPARED_DIR --binary RAMSES --manifest-sha256 SHA256" >&2
  echo "manual:  $0 --manual-lageunha --run-dir PREPARED_DIR --binary RAMSES --manifest-sha256 SHA256" >&2
  echo "mpi2:    $0 --manual-lageunha-mpi2 --run-dir PREPARED_DIR --binary RAMSES --manifest-sha256 SHA256" >&2
  exit 2
}

run_dir=
binary=
bytes_per_output=
reserve_bytes=1073741824
execute=false
prepare=false
manual_lageunha=false
manual_lageunha_mpi2=false
mpi_ranks=1
mpi_ranks_specified=false
expected_source_commit=
source_tree=
manifest_sha256=
yield_table_source=
while (($#)); do
  case "$1" in
    --run-dir) run_dir=${2:?}; shift 2 ;;
    --binary) binary=${2:?}; shift 2 ;;
    --bytes-per-output) bytes_per_output=${2:?}; shift 2 ;;
    --reserve-bytes) reserve_bytes=${2:?}; shift 2 ;;
    --expected-source-commit) expected_source_commit=${2:?}; shift 2 ;;
    --source-tree) source_tree=${2:?}; shift 2 ;;
    --manifest-sha256) manifest_sha256=${2:?}; shift 2 ;;
    --yield-table-source) yield_table_source=${2:?}; shift 2 ;;
    --mpi-ranks) mpi_ranks=${2:?}; mpi_ranks_specified=true; shift 2 ;;
    --prepare) prepare=true; shift ;;
    --execute) execute=true; shift ;;
    --manual-lageunha) manual_lageunha=true; shift ;;
    --manual-lageunha-mpi2) manual_lageunha_mpi2=true; shift ;;
    *) usage ;;
  esac
done

[[ -n $run_dir && -n $binary ]] || usage
mode_count=0
[[ $prepare == true ]] && ((mode_count+=1))
[[ $execute == true ]] && ((mode_count+=1))
[[ $manual_lageunha == true ]] && ((mode_count+=1))
[[ $manual_lageunha_mpi2 == true ]] && ((mode_count+=1))
[[ $mode_count == 1 ]] || { echo 'CAPTURE-RESTART: choose exactly one execution mode' >&2; exit 2; }
binary=$(realpath -- "$binary")
[[ -x $binary ]] || { echo "CAPTURE-RESTART: binary is not executable: $binary" >&2; exit 1; }
if [[ $prepare == true ]]; then
  [[ $mpi_ranks == 1 || $mpi_ranks == 2 ]] || usage
  [[ -n $bytes_per_output && -n $expected_source_commit && -n $source_tree && -n $yield_table_source ]] || usage
  python3 "$here/prepare_and_validate.py" prepare "$run_dir" \
    --binary "$binary" --expected-source-commit "$expected_source_commit" \
    --source-tree "$source_tree" \
    --yield-table-source "$yield_table_source" \
    --bytes-per-output "$bytes_per_output" --reserve-bytes "$reserve_bytes" \
    --mpi-ranks "$mpi_ranks"
  exit 0
fi
[[ $mpi_ranks_specified == false ]] || usage
host=$(hostname -s)
if [[ $execute == true ]]; then
  python3 "$here/prepare_and_validate.py" execution-context --execution-mode slurm \
    --hostname "$host" --slurm-job-id "${SLURM_JOB_ID:-}"
  command -v srun >/dev/null 2>&1 || { echo 'CAPTURE-RESTART: srun is unavailable' >&2; exit 1; }
  execution_mode=slurm
elif [[ $manual_lageunha_mpi2 == true ]]; then
  python3 "$here/prepare_and_validate.py" execution-context --execution-mode manual-lageunha-mpi2 \
    --hostname "$host"
  command -v mpirun >/dev/null 2>&1 || { echo 'CAPTURE-RESTART: mpirun is unavailable' >&2; exit 1; }
  mpi_runner=$(realpath -- "$(command -v mpirun)")
  mpi_prefix=$(dirname -- "$(dirname -- "$mpi_runner")")
  [[ $mpi_runner == /opt/ohpc/pub/intel/oneapi/mpi/*/bin/mpirun ]] || {
    echo "CAPTURE-RESTART: expected the Intel MPI launcher, found $mpi_runner" >&2; exit 1;
  }
  mpi_libraries=$(ldd "$binary")
  [[ $mpi_libraries == *"libmpi.so.12 => $mpi_prefix/lib/"* &&
     $mpi_libraries == *"libmpifort.so.12 => $mpi_prefix/lib/"* &&
     $mpi_libraries != *"not found"* ]] || {
    echo 'CAPTURE-RESTART: binary MPI libraries do not match the Intel MPI launcher' >&2; exit 1;
  }
  mpi_version=$(mpirun --version | head -n 1)
  [[ $mpi_version == Intel*MPI* ]] || {
    echo 'CAPTURE-RESTART: MPI launcher version is not Intel MPI' >&2; exit 1;
  }
  execution_mode=manual-lageunha-mpi2
else
  python3 "$here/prepare_and_validate.py" execution-context --execution-mode manual-lageunha \
    --hostname "$host"
  execution_mode=manual-lageunha
fi
[[ -n $manifest_sha256 ]] || usage
run_dir=$(realpath -- "$run_dir")
expected_mpi_ranks=1
[[ $execution_mode == manual-lageunha-mpi2 ]] && expected_mpi_ranks=2
python3 "$here/prepare_and_validate.py" launch-check "$run_dir" --binary "$binary" \
  --manifest-sha256 "$manifest_sha256" --expected-mpi-ranks "$expected_mpi_ranks"

source_root=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_tree"])' \
  "$run_dir/preflight_manifest.json")
source_commit=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["expected_source_commit"])' \
  "$run_dir/preflight_manifest.json")
{
  echo "execution_mode=$execution_mode"
  echo "expected_mpi_ranks=$expected_mpi_ranks"
  if [[ $execution_mode == manual-lageunha-mpi2 ]]; then
    echo "mpi_runner=$mpi_runner"
    echo "mpi_runner_sha256=$(sha256sum "$mpi_runner" | awk '{print $1}')"
    echo "mpi_version=$mpi_version"
    echo "mpi_prefix=$mpi_prefix"
  fi
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
  sha256sum "$run_dir/stage1.nml" "$run_dir/stage2.nml" "$run_dir/ic_sink" \
    "$run_dir/yield_table.asc"
  sha256sum "$run_dir/preflight_manifest.json"
} | tee "$run_dir/preflight_provenance.txt"

(
  cd "$run_dir"
  if [[ $execution_mode == slurm ]]; then
    srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
      "$binary" stage1.nml >stage1.log 2>&1
  elif [[ $execution_mode == manual-lageunha-mpi2 ]]; then
    OMP_NUM_THREADS=1 mpirun -n 2 "$binary" stage1.nml >stage1.log 2>&1
  else
    "$binary" stage1.nml >stage1.log 2>&1
  fi
  [[ -f output_00001/COMPLETE ]]
  [[ -f output_00001/SMBH_CAPTURE_LINEAGE ]]
  if [[ $execution_mode == slurm ]]; then
    srun --kill-on-bad-exit=1 --nodes=1 --ntasks=1 --cpus-per-task=1 \
      "$binary" stage2.nml >stage2.log 2>&1
  elif [[ $execution_mode == manual-lageunha-mpi2 ]]; then
    OMP_NUM_THREADS=1 mpirun -n 2 "$binary" stage2.nml >stage2.log 2>&1
  else
    "$binary" stage2.nml >stage2.log 2>&1
  fi
)
python3 "$here/prepare_and_validate.py" postcheck "$run_dir" \
  --expected-mpi-ranks "$expected_mpi_ranks"
