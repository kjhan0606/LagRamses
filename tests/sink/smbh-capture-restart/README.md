# SMBH capture-ledger restart smoke

This synthetic smoke starts five SMBHs in a two-member group and a
transitive three-member chain. Stage 1 writes `output_00001`; stage 2 restarts
that COMPLETE checkpoint and writes `output_00002`. Postcheck requires one
committed `BINARY`, one committed `MULTIPLE`, member counts 2+3, both lineage
sidecars, and a validator-selected active restart branch. The original
one-rank path remains available; a separately sealed two-rank manual path
checks MPI rank count in both outputs.

The fixture is periodic: `BOUNDARY_PARAMS` is intentionally absent, selecting
RAMSES's default `nboundary=0` topology.  This isolates the capture/restart
contract for five central sinks; it is a fixture workaround, not a repair or
regression diagnosis for the six-boundary `make_grid_fine` failure preserved
in the `21d62e4` Lageunha run.  Nonperiodic AMR initialization remains
unvalidated and requires a separate minimal reproducer and source-level
investigation before any repair.  The existing Bondi physics test remains
unchanged.  The preflight rejects any added `BOUNDARY_PARAMS` block so a
boundary-topology change cannot be mistaken for a capture-ledger test change.

The workflow is deliberately two phase. `--prepare` is safe on a login node:
it creates only two effective namelists, `ic_sink`, a run-local
`yield_table.asc`, and a JSON audit manifest;
it never builds or runs RAMSES. The run directory must not exist. `--execute`
works only inside Slurm and only with that existing directory. It revalidates
the exact input hashes and policies, clean source worktree and HEAD, binary path
and SHA-256, pristine output state, and current free space before `srun`.
Execution also requires the manifest SHA-256 printed by preparation, so the
approved manifest cannot silently redefine its own hashes or space estimate.

The required `--yield-table-source` is copied rather than referenced in place.
Both effective namelists contain the exact absolute run-local path, and the
copy's size and SHA-256 are sealed into the manifest and revalidated before
launch. This table satisfies the current unconditional sink initialization
open; all enrichment channels and AGN remain disabled, so its presence is not
evidence that chemistry or feedback is active.

The manifest records complete parsed stage-1/stage-2 namelists and absolute
paths, the output schedule (`noutput=1`, no `aout`, `tout=1.0d100`,
`foutput=1`, `fbackup=1000000`), two projected periodic outputs, measured or
conservative bytes per output, total projection, reserve, clean source commit,
and binary identity. The binary must live inside that source worktree. This
strongly binds the launch to an unchanged build artifact and clean checkout,
but SHA-256 plus location cannot independently prove compiler inputs; retain
the build log if stronger toolchain provenance is required.

## 1. Build a fresh CPU binary from the eventual committed HEAD

First commit the runner changes that are to be tested. Do not use the dirty hub
checkout as either the runner or build source. On the actual `LagEunha` CPU
host, the commands below resolve the current HEAD and make one clean detached
worktree containing the runner, validator, source, and resulting binary. They
refuse an existing build path and pin the compiler/MPI module versions:

```bash
set -euo pipefail
[[ $(hostname -s | tr '[:upper:]' '[:lower:]') == lageunha ]] || exit 1
module purge
module load intel/tbb/2022.3 intel/umf/1.0.2 \
  intel/compiler-rt/2025.3.0 intel/compiler/2025.3.0 intel/mpi/2021.17
repo=/home/kjhan/BACKUP/lagRamses
commit=$(git -C "$repo" rev-parse HEAD)
short_commit=$(git -C "$repo" rev-parse --short=8 HEAD)
build_tree=/home/kjhan/BACKUP/lagRamses-build-capture-smoke-$short_commit
[[ ! -e $build_tree ]] || { echo "refusing existing build path: $build_tree" >&2; exit 1; }
git -C "$repo" worktree add --detach "$build_tree" "$commit"
test -z "$(git -C "$build_tree" status --porcelain --untracked-files=normal)"
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
make -C "$build_tree/bin" clean
make -C "$build_tree/bin" -j1 HDF5=1 USE_FFTW=1
test -x "$build_tree/bin/ramses_final3d"
test -z "$(git -C "$build_tree" status --porcelain --untracked-files=normal)"
echo "source_commit=$commit"
sha256sum "$build_tree/bin/ramses_final3d"
```

These are manual Lageunha CPU-build instructions; no GPU or Slurm allocation is
needed. No build is performed by the runner. Retain the terminal/build log and
keep the worktree and binary unchanged through execution.

Use `-j1` explicitly. On Lageunha, a clean `-j4` build exposed a Fortran module
dependency race (`amr_commons` was compiled before `amr_parameters.mod`), while
the serial command above succeeded. The clean test build at commit `35a307d`
produced `bin/ramses_final3d` with SHA-256
`8a89cc49466e0bbf748e7e172b9216a2f83aa73955f58af3bb009b96a2930421`.
That digest is evidence for that exact build only; the dynamic committed-HEAD
workflow must use and approve the digest printed by its own clean build.

## 2. Prepare and inspect on Lageunha

On the same Lageunha server that will execute the manual smoke, choose a
never-before-used `/home` run path and a measured/conservative output size:

```bash
build_tree=/home/kjhan/BACKUP/lagRamses-build-capture-smoke-COMMIT8
commit=$(git -C "$build_tree" rev-parse HEAD)
short_commit=$(git -C "$build_tree" rev-parse --short=8 HEAD)
binary="$build_tree/bin/ramses_final3d"
yield_table_source=/gpfs/kjhan/Run_JWST/opt_run/yield_table.asc
run_dir=/home/kjhan/BACKUP/smbh-capture-restart-$short_commit
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --prepare --run-dir "$run_dir" --binary "$binary" \
  --source-tree "$build_tree" --expected-source-commit "$commit" \
  --yield-table-source "$yield_table_source" \
  --bytes-per-output 500000000
python3 -m json.tool "$run_dir/preflight_manifest.json" | less
manifest_sha256=$(sha256sum "$run_dir/preflight_manifest.json" | awk '{print $1}')
# Required manual evidence on Lageunha; retain this output with the manifest.
/usr/lpp/mmfs/bin/mmlsquota -u "$USER" aicpuhome
df -h /home
```

Audit both effective namelists, output schedule, projected bytes, reserve/free
space, Lustre quota evidence, source SHA, and binary SHA before submission.
For `/scratch`, `lfs quota` values of zero mean no explicit numeric user ceiling
was reported; they do not create an invented "free quota" value. The runner
records that evidence and separately gates on filesystem free bytes. If either
the quota or hard-limit field is nonzero, it also requires the projected bytes
plus reserve to fit below the stricter nonzero ceiling at prepare and launch.
For the manual `/home` path, the runner still gates on `/home` filesystem free
bytes. The operator must also capture the GPFS `mmlsquota` report above. The
observed `quota=0, limit=0` means no explicit numeric user ceiling was reported,
not infinite invented free quota.

## 3. Execute manually on Lageunha (CPU-only)

For this one-rank CPU smoke, the normal path is manual execution on the actual
`LagEunha` host. The runner checks the hostname case-insensitively, performs the
same sealed launch revalidation, and invokes the binary directly and serially:

```bash
build_tree=/home/kjhan/BACKUP/lagRamses-build-capture-smoke-COMMIT8
binary="$build_tree/bin/ramses_final3d"
short_commit=$(git -C "$build_tree" rev-parse --short=8 HEAD)
run_dir=/home/kjhan/BACKUP/smbh-capture-restart-$short_commit
manifest_sha256=REPLACE_WITH_APPROVED_64_HEX_DIGEST
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --manual-lageunha --run-dir "$run_dir" --binary "$binary" \
  --manifest-sha256 "$manifest_sha256"
```

`--manual-lageunha` refuses other hosts and never uses `srun`. It is not an
authorization to run a GPU workload manually.

## 4. Two-rank manual MPI regression on LagEunha

For a two-rank MPI regression on LagEunha, prepare a new, never-used run
directory with `--mpi-ranks 2` and a separately audited output-size estimate.
Do not reuse the one-rank manifest or output directory. The manifest binds the
rank count, and `--manual-lageunha-mpi2` refuses a manifest prepared for one
rank. The runner also requires the Intel MPI launcher and the same Intel MPI
library prefix in the binary. Old schema-2 serial manifests are intentionally
not accepted by the new schema-3 runner; prepare a fresh run rather than
rewriting an approved manifest. Load the build-matched modules, then run:

```bash
module purge
module load intel/tbb/2022.3 intel/umf/1.0.2 \
  intel/compiler-rt/2025.3.0 intel/compiler/2025.3.0 intel/mpi/2021.17
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
build_tree=/home/kjhan/BACKUP/lagRamses-build-capture-smoke-COMMIT8
binary="$build_tree/bin/ramses_final3d"
commit=$(git -C "$build_tree" rev-parse HEAD)
short_commit=$(git -C "$build_tree" rev-parse --short=8 HEAD)
run_dir=/home/kjhan/BACKUP/smbh-capture-restart-mpi2-$short_commit
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --prepare --run-dir "$run_dir" --binary "$binary" \
  --source-tree "$build_tree" --expected-source-commit "$commit" \
  --yield-table-source /gpfs/kjhan/Run_JWST/opt_run/yield_table.asc \
  --bytes-per-output 1000000000 --mpi-ranks 2
# Inspect the effective namelists, output projection, quota/free-space report,
# and source/binary identities before continuing.
/usr/lpp/mmfs/bin/mmlsquota -u "$USER" aicpuhome
df -h /home
```

Stop here. Inspect the effective files, projected storage, quota/free-space
evidence, clean source and binary hash. Only after that review, use the
manifest digest printed by preparation in a separate command:

```bash
manifest_sha256=$(sha256sum "$run_dir/preflight_manifest.json" | awk '{print $1}')
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --manual-lageunha-mpi2 --run-dir "$run_dir" --binary "$binary" \
  --manifest-sha256 "$manifest_sha256"
```

The runner uses exactly two MPI ranks with one OpenMP thread each, checks
`ncpu=2` in both RAMSES output headers, then applies the same mass, momentum,
`BINARY`/`MULTIPLE`, lineage, and restart-continuity postcheck as the serial
smoke. A failed MPI run must be retained for diagnosis, never relabelled as a
passed one-rank run.

## 5. Optional execution under Slurm

Slurm may use `/scratch`, but it needs a distinct preparation there (including
the same explicit `--yield-table-source` argument); a sealed
manual `/home` manifest cannot be moved. On the Slurm login host where
`/scratch/$USER` exists, invoke the same clean-worktree `--prepare` command with
the unique `run_dir=/scratch/$USER/smbh-capture-restart-$short_commit-slurm`.
Inspect its filesystem/Lustre report and approve that manifest SHA. Then use
those exact paths:

```bash
sbatch <<'SBATCH'
#!/usr/bin/env bash
#SBATCH --job-name=capture_restart_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
set -euo pipefail
build_tree=/home/kjhan/BACKUP/lagRamses-build-capture-smoke-COMMIT8
binary="$build_tree/bin/ramses_final3d"
short_commit=$(git -C "$build_tree" rev-parse --short=8 HEAD)
run_dir=/scratch/$USER/smbh-capture-restart-$short_commit-slurm
manifest_sha256=REPLACE_WITH_APPROVED_64_HEX_DIGEST
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --execute --run-dir "$run_dir" --binary "$binary" \
  --manifest-sha256 "$manifest_sha256"
SBATCH
```

The requested GPU satisfies cluster allocation policy only; the smoke runs one
CPU MPI rank. Allocation-time provenance records job/host, source state, binary
and manifest identities, and all effective input hashes. Any input edit,
binary replacement, unexpected output/log/ledger, source drift, or inadequate
space aborts before RAMSES. The fixture is not an MPI scaling or scientific
binary-formation test.
