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

The first two-rank execution on LagEunha passed on 2026-10-04 using clean
source commit `cb924745be266bf15fe5c7b28261124d1d48a2c3` and binary
SHA-256 `c6a534f87c82174ab58391779d200b9f2dd9ab597cb8ebd54350443c90c4425d`.
The prepared manifest SHA-256 was
`646c49feba2e2e4e259f700bea4e3f1e2cc3d25138fc9a6f2778474c0355107a`.
The preserved run is
`/home/kjhan/BACKUP/smbh-capture-restart-mpi2-cb92474-20261004/`;
its ledger SHA-256 is
`914f8f6dbd4b666bc71933ceacfaa0c69271862c0edd94639b224e36b505744f`.
Both outputs report `ncpu=2`; postcheck found one committed `BINARY`
(two members) and one committed `MULTIPLE` (three members), no duplicate or
incomplete batches, mass conservation, and sink-state continuity across the
restart. An independent read-only review confirmed the ledger and lineage.
`committed_batch_seq` is per attempt: output 1 has one committed batch, while
the new restart attempt has none and correctly writes zero in output 2.
This smoke uses zero-velocity sinks, so its momentum check is algebraically
correct but trivial; it does not validate sink motion or a new capture after
restart. It also does not validate nonperiodic AMR boundaries or a physical
galaxy-merger capture rate.

### Latest clean-HEAD two-rank rerun — 2026-10-12

The smoke was rerun from clean detached source commit
`56bcd0d9c9715e7d0419aaf61df7ee3490105947`, using the pinned Intel MPI
2021.17 toolchain and binary SHA-256
`3eb8b58beed20783387da1773d5bef76dec0667214e9bef0317ec21f997d2c53`.
The schema-3 preflight manifest SHA-256 was
`94bc870eb7d797eaa6d78cc8c7772e5be740d231f65a8bddca0f8f28229fb604`.
The preserved run is
`/home/kjhan/BACKUP/smbh-capture-restart-mpi2-56bcd0d0/`; its ledger SHA-256
is `7270ed5699e72c3f86f6fc7b1c2b0ea5fa3ea65c831196b2a241010cd2c1f0bd`.
Both outputs have `COMPLETE`, report `ncpu=2`, and advance coarse steps 1 then
2. The postcheck passed with one committed batch, two unique complete events
(one `BINARY`, one `MULTIPLE` with 2+3 members), zero duplicate or incomplete
batches/events, zero invalid JSON lines, and a validator-selected active
restart branch. Mass and sink-state restart-continuity checks passed. This is
the latest clean-HEAD two-rank smoke; the zero-velocity momentum and
nonperiodic-boundary limitations above still apply.

### Moving-pair fixture — runtime verified 2026-10-12

The working-tree `ic_sink` now gives each capture group nonzero internal
tangential velocities while keeping each group's initial centre-of-mass
momentum zero. The preparation gate checks those two properties, and the
ledger postcheck now requires every recorded pair to have nonzero relative
speed and angular momentum. Post-compaction mass, periodic COM position, COM
velocity, and momentum are checked against the pre-compaction event; the
restart output is then checked for sink-state continuity. The run uses the
same clean-HEAD binary as the two-rank test because the Fortran source is
unchanged, with the fixture variant separately hash-bound in its manifest.

The first Slurm submission, job `420719`, remained pending under
`QOSMaxGRESPerUser` and was cancelled before starting; it touched no solver
outputs. A first serial attempt on lageunha then stopped before initialization
because `libimf.so` was missing; that failed directory is preserved at
`/gpfs/kjhan/lagramses_capture_restart_smoke/moving-pair-56bcd0d0-20261012/`
and contains no `output_*` directories. Loading the pinned Intel compiler
runtime, TBB, UMF, compiler, and MPI modules resolved the dynamic-library
dependency. The successful serial run used a fresh directory:
`/gpfs/kjhan/lagramses_capture_restart_smoke/moving-pair-56bcd0d0-20261012-r2/`.

The binary is from clean source commit
`56bcd0d9c9715e7d0419aaf61df7ee3490105947`, SHA-256
`3eb8b58beed20783387da1773d5bef76dec0667214e9bef0317ec21f997d2c53`. The
prepared `ic_sink` SHA-256 is
`553cbbe7952c6f53feed39560f38a55c1a51e272e0bdffa26798d3f734445d90`; the r2
manifest SHA-256 is
`e03d6fedfe206e10ac54fc21e8984e69b7ff74f14975cf388e86ec64d497893c`.
The effective schedule wrote only outputs at coarse steps 1 and 2
(`noutput=1`, `tout=1.0d100`, `foutput=1`, `fbackup=1000000`). The preflight
reserved 2 GB for outputs plus 1 GiB reserve.

Both outputs have `COMPLETE`. The ledger and restart postchecks pass: one
committed batch, one `BINARY`, one 3-member `MULTIPLE`, two unique complete
events, no duplicate/incomplete records or invalid JSON, nonzero relative
speed and angular momentum for every pair, conserved post-compaction mass,
periodic COM position, velocity and momentum, and sink-state continuity across
restart. The first working-tree postcheck exposed output quantization: RAMSES
sink CSV coordinates are printed to ten decimal places while the ledger stores
full precision. The comparison now applies the existing relative tolerance to
the periodic-equivalent coordinate, with a regression reproducing the Fortran
CSV precision. Fixture tests pass (50 tests and 4 subtests).

### Two-rank moving-pair fixture — runtime verified 2026-10-29

The moving-pair fixture was rerun with two MPI ranks from clean source commit
`227a7497c6b80620ce865e541d3af4321d8bc935`. The binary SHA-256 is
`49d3d84ae2b5db6d08cde140b228eaf91678065fe9d5b86c79dd92f3e3133c7d`, and
the sealed schema-3 manifest SHA-256 is
`b683043635a7d9e5644b3c910cc781fdffe55616ee4020df0a1f4fc0b16a1400`. The
manifest's `ic_sink` SHA-256 is
`553cbbe7952c6f53feed39560f38a55c1a51e272e0bdffa26798d3f734445d90`, the
nonzero-motion fixture hash. The preserved run is
`/home/kjhan/BACKUP/smbh-capture-restart-mpi2-moving-227a7490-20261029-r1/`;
its capture ledger SHA-256 is
`01b62de80800fa5d281faf28624fb71dea83805d9721b807b60e053c5eca5771a`.

The clean-worktree unit/preflight suite passed 44 tests plus 4 subtests, and
the same worktree's postcheck independently passed on the runtime output. Both
outputs report `ncpu=2` and are attached to the active restart branch. The
schema-v2 ledger has one committed batch, one `BINARY`, one three-member
`MULTIPLE`, two unique complete events, and zero duplicate/incomplete events
or batches and invalid JSON lines. Every recorded pair has nonzero relative
speed and angular momentum. The postcheck also verifies post-compaction mass,
periodic COM position, COM velocity and momentum against the pre-compaction
ledger, then checks sink-state continuity across restart. This is a synthetic
restart/ledger regression, not a capture-rate or galaxy-merger physics
validation. It does not close the nonperiodic-AMR-boundary or FOF-backend
equivalence gaps. The Fortran source was unchanged by this fixture commit.

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

## 6. Differential FoF backend unit test

The active sink merge selects the sequential FoF implementation below 300
sinks and the OMP+MPI tree implementation at or above 300. The moving-pair
restart fixture has only five sinks, so it does not exercise that tree path.
The standalone test-only driver calls both backends on an identical nine-sink
fixture and compares component membership, independent of group labels. The
fixture covers periodic wrapping, transitive links, an isolated sink, and
components crossing both OpenMP-thread and MPI-rank partitions.

Submit from the lagRamses repository root:

```bash
sbatch tests/sink/smbh-capture-restart/run_fof_backend_equivalence.slurm
```

The job compiles the current `patch/lagRamses/kjhan.f90` with `_KJHAN_TEST`,
then runs two MPI ranks with two OpenMP threads each. Intel MPI is bound to
Slurm's PMI2 plugin using `I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so` and
`srun --mpi=pmi2`. Logs are written to `/gpfs/kjhan/fof-equivalence-<job>.out`
and `.err`.

Job `1234486` completed with exit 0. Its output reports that the sequential
and OMP+MPI partitions both match the expected four components. Source hash:
`28e2e7e228457810dce2d3911f6980e9dda90c08d2776cd05603a8a9d293d357`; driver
hash: `e3d3818f7bc347fb157561c766eb94dc2f6d1a109e9add07c650e18fae400b07`.
This is a direct backend-level grouping regression, not a full `merge_sink`
dispatch test: `_KJHAN_TEST` omits the binding-energy `vrel_merge` filter, and
the production threshold and canonicalization remain unchanged.
