# SMBH capture-ledger restart smoke

This serial, synthetic smoke starts five SMBHs in two separated FoF components:
a two-member group and a transitive three-member chain.  Stage 1 advances to
coarse step 1 and must create `output_00001`; stage 2 restarts that COMPLETE
checkpoint, advances to step 2, and must create `output_00002`.  The postcheck
requires one committed `BINARY`, one committed `MULTIPLE`, member counts 2+3,
both COMPLETE lineage sidecars, and a validator-selected active restart branch.
`vrel_merge=.false.` is explicit so this fixture isolates FoF capture geometry;
`nsinkmax=5` holds the complete initial population and compaction only reduces it.
The run is explicitly noncosmological with hydro, particles, Poisson, and sinks
enabled. Postcheck also requires lineage coarse steps 1 and 2, exactly the two
expected output directories, a stage-1 log proving five loaded seed sinks, and
no warning, fatal, or standalone NaN diagnostic in either stage log.

The runner is deliberately fail-closed: the run directory must not exist, the
binary must be executable, `--execute` is mandatory, and the caller must give a
measured or conservative byte estimate for each of the two outputs.  Preflight
prints the exact effective namelist paths, output mechanisms, expected storage,
reserve, and current free space before invoking RAMSES.

Submit from the repository root; the job script supplies the required compute
allocation and the runner uses one-task `srun` steps:

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
cd /absolute/path/to/lagRamses
tests/sink/smbh-capture-restart/run_smoke.sh \
  --run-dir /scratch/$USER/smbh-capture-restart-${SLURM_JOB_ID} \
  --binary "$PWD/bin/ramses3d" \
  --bytes-per-output 500000000 \
  --execute
SBATCH
```

The runner refuses execution unless `SLURM_JOB_ID` is present. It records the
binary SHA-256, source commit/branch/status, and hashes of both effective
namelists and `ic_sink` in `preflight_provenance.txt` before the first `srun`.
The GPU is reserved solely to satisfy cluster allocation policy; the smoke
still executes one CPU rank and does not exercise GPU code. Use a freshly built
executable from the intended capture-ledger commit. The recorded source status
and binary hash preserve identities independently, but do not prove that an
arbitrary pre-existing executable was built from the recorded source revision.

This fixture assumes the active `patch/lagRamses` build and its version-2
capture lineage support.  It is intentionally a one-rank smoke; MPI behavior
and scientific binary formation are outside its claim.
