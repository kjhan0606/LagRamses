# Claude Opus 5 end audit — SNRT AGN source-to-RT transaction atomicity

Date: 2026-09-07 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Auditor: Claude Opus 5, invoked through `claude -p --model opus` without
`--bare`
Mode: read-only; no repository files, jobs, commits, or pushes were changed
by the auditor

## Verdict

**CONDITIONAL PASS** for the declared serial, fresh-start SNRT AGN
source-to-RT/chemistry/dust ordering scope.

The ordering repair is internally correct: the transaction snapshot precedes
source deposition, multigroup deposition is validated before mutation, and
pending AGN energy is cleared only after the coupled global commit. No scoped
blocking Fortran defect was found.

## Verified

- The active `bin/Makefile` selects `../patch/lagRamses` first in `VPATH` and
  includes the active SNRT driver object.
- The source writes through the same state-slot mapping captured by the
  transaction snapshot, so the admitted serial path can restore injected
  photons.
- All post-allocation rollback and normal-exit paths release the source mask;
  DUST_LIVE cleanup is correctly guarded.
- The clear operation is after RT/chemistry and live-dust commit, and all
  failure returns examined by the audit precede it.

## Conditions and dispositions

1. **HIGH — retry wording overclaims runtime behavior.** `clean_stop` calls
   `MPI_ABORT`; there is no in-process retry loop or checkpoint persistence for
   `agn_pending_erg`. “Retryable” must therefore be qualified as
   “retained in memory until terminal failure,” not durable restart retry.
2. **MEDIUM — driver-level failure coverage is incomplete.** The added native
   source smoke reconstructs the ordering but does not invoke
   `snrt_ramses_advance_level` or exercise the production mask, locator, RT
   non-convergence, DUST_LIVE precommit, and commit-rejection branches. This
   remains a follow-up verification item before unconditional PASS.
3. **MEDIUM — MPI replicated-ledger behavior is pre-existing.** Live AGN is
   currently restricted to the serial admitted domain. A two-phase global
   commit and replicated pending-energy policy are required before MPI AGN is
   admitted.
4. **MEDIUM — zeroing contract needs clearer documentation.**
   `snrt_agn_source_commit` clears the accepted event rather than debiting a
   partial amount; the new comment should preserve that invariant while the
   clear window spans the coupled solve.
5. **LOW — provenance wording is inaccurate.** The evidence refers to an
   older `patch/cuRamses` SNRT copy, but no such `snrt_*` mirror exists in the
   active tree. The source-selection conclusion remains correct.
6. **LOW — dead pre-source `incoming_intensity` copy should be removed or
   explicitly justified; it is overwritten before transport consumes it.**
7. **INFO — transaction atomicity is leaf-slot scoped.** The coarse-flux
   correction is not restored by the leaf-slot restore routine; this boundary
   must remain explicit.

## Evidence boundary

The implementation evidence records successful native source/transaction
smokes and active-driver compile checks. The Opus audit itself was read-only
and did not rerun writers. Existing linked binaries were not accepted as
HEAD-qualified evidence for `210b611` because the recorded binary predates the
commit and embeds an older dirty identity. A fresh full-link record and a
driver-faithful failure harness remain verification work.

After the audit, the following bounded rechecks were run from the active
`/gpfs/kjhan/LRD_JWST` tree:

- `make -C bin SNRT=1 USE_CUDA=1 DUST_LIVE=0 snrt_ramses_driver.o` passed,
  selecting `../patch/lagRamses/snrt_ramses_driver.f90` through the active
  VPATH.
- `make -B -C bin SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1
  snrt_ramses_driver.o` passed for the conditional live-dust path.
- `simulation/snrt/tests/run_fp15_agn_efficiency.sh` passed, including
  `SNRT_AGN_SOURCE_COUPLED_ROLLBACK_PENDING_PASS`,
  `SNRT_AGN_SOURCE_COUPLED_COMMIT_PENDING_CLEAR_PASS`, and
  `SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK`.
- Fresh full links also passed without overwriting the default executable:
  `make -C bin -B -j4 SNRT=1 DUST_LIVE=0 USE_CUDA=1 HDF5=1
  EXEC=/gpfs/kjhan/LRD_JWST/.snrt-atomicity-link.se08tt/ramses_snrt_atomicity_cpu
  ramses` produced `ramses_snrt_atomicity_cpu3d`, and the corresponding
  `DUST_LIVE=1` command produced
  `/gpfs/kjhan/LRD_JWST/.snrt-atomicity-dust-link.PIQ814/ramses_snrt_atomicity_dust3d`.
  Both link lines include the active `patch/lagRamses` SNRT driver, HDF5,
  CUDA, and FFTW libraries.
- These links embed `210b611542903147a9197836d0e3f2f755506cca-dirty`: the
  source-tree commit is identified exactly, but the worktree also contains
  the pending audit-closure edits and an unrelated pre-existing tracked
  change. They are therefore fresh working-tree links, not clean-tree
  revision artifacts.

These rechecks do not close the two remaining evidence conditions above: they
are not a production-driver failure harness and do not constitute a clean
exact-revision artifact.

The conditional result does not approve physical AGN SED/MAD choices,
restart/migration persistence, simultaneous legacy-plus-SNRT ownership,
MPI AGN, or initialized live RAMSES evolution.
