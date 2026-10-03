# SMBH capture-ledger restart smoke

This serial synthetic smoke starts five SMBHs in a two-member group and a
transitive three-member chain. Stage 1 writes `output_00001`; stage 2 restarts
that COMPLETE checkpoint and writes `output_00002`. Postcheck requires one
committed `BINARY`, one committed `MULTIPLE`, member counts 2+3, both lineage
sidecars, and a validator-selected active restart branch.

The workflow is deliberately two phase. `--prepare` is safe on a login node:
it creates only two effective namelists, `ic_sink`, and a JSON audit manifest;
it never builds or runs RAMSES. The run directory must not exist. `--execute`
works only inside Slurm and only with that existing directory. It revalidates
the exact input hashes and policies, clean source worktree and HEAD, binary path
and SHA-256, pristine output state, and current free space before `srun`.
Execution also requires the manifest SHA-256 printed by preparation, so the
approved manifest cannot silently redefine its own hashes or space estimate.

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
checkout as either the runner or build source. The job below resolves the then-
current HEAD and makes one clean detached worktree containing the runner,
validator, source, and resulting binary.

Submit the compilation itself to Slurm; do not compile on the login node. This
job refuses an existing build path, creates a detached worktree, pins the
compiler/MPI module versions, and leaves an external build log:

```bash
repo=/home/kjhan/BACKUP/lagRamses
commit=$(git -C "$repo" rev-parse HEAD)
sbatch --export=ALL,LAGRAMSES_REPO="$repo",LAGRAMSES_COMMIT="$commit" <<'SBATCH'
#!/usr/bin/env bash
#SBATCH --job-name=build_capture_smoke
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/%u/lagRamses-build-capture-%j.log
set -euo pipefail
module purge
module load intel/compiler/2025.3.0
module load intel/mpi/2021.17
repo=$LAGRAMSES_REPO
commit=$LAGRAMSES_COMMIT
build_tree=/scratch/$USER/lagRamses-build-capture-smoke
[[ ! -e $build_tree ]] || { echo "refusing existing build path: $build_tree" >&2; exit 1; }
git -C "$repo" worktree add --detach "$build_tree" "$commit"
test -z "$(git -C "$build_tree" status --porcelain --untracked-files=normal)"
export LD_LIBRARY_PATH=/home/kjhan/local/hdf5/lib:/home/kjhan/local/lib:${LD_LIBRARY_PATH:-}
make -C "$build_tree/bin" clean
make -C "$build_tree/bin" HDF5=1 USE_FFTW=1
test -x "$build_tree/bin/ramses_final3d"
test -z "$(git -C "$build_tree" status --porcelain --untracked-files=normal)"
echo "source_commit=$commit"
sha256sum "$build_tree/bin/ramses_final3d"
SBATCH
```

Wait for a successful build job and retain its Slurm log. These are instructions
only; no build is performed by the runner. Keep the worktree and binary
unchanged through execution.

## 2. Prepare and inspect on the login node

Choose a never-before-used run path and a measured/conservative output size:

```bash
build_tree=/scratch/$USER/lagRamses-build-capture-smoke
commit=$(git -C "$build_tree" rev-parse HEAD)
binary="$build_tree/bin/ramses_final3d"
run_dir=/scratch/$USER/smbh-capture-restart-prepared
"$build_tree/tests/sink/smbh-capture-restart/run_smoke.sh" \
  --prepare --run-dir "$run_dir" --binary "$binary" \
  --source-tree "$build_tree" --expected-source-commit "$commit" \
  --bytes-per-output 500000000
python3 -m json.tool "$run_dir/preflight_manifest.json" | less
manifest_sha256=$(sha256sum "$run_dir/preflight_manifest.json" | awk '{print $1}')
```

Audit both effective namelists, output schedule, projected bytes, reserve/free
space, Lustre quota evidence, source SHA, and binary SHA before submission.
For `/scratch`, `lfs quota` values of zero mean no explicit numeric user ceiling
was reported; they do not create an invented "free quota" value. The runner
records that evidence and separately gates on filesystem free bytes. If either
the quota or hard-limit field is nonzero, it also requires the projected bytes
plus reserve to fit below the stricter nonzero ceiling at prepare and launch.

## 3. Execute that prepared directory under Slurm

Use the identical absolute paths. Do not derive a new run directory from the
job ID: preparation has already fixed and sealed it.

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
build_tree=/scratch/$USER/lagRamses-build-capture-smoke
binary="$build_tree/bin/ramses_final3d"
run_dir=/scratch/$USER/smbh-capture-restart-prepared
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
