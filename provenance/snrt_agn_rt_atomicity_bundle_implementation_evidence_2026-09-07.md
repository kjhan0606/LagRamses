# SNRT AGN source-to-RT transaction atomicity — implementation evidence

Date: 2026-09-07 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Parent: AGN accepted-fuel/overlapping-deposition and F-P2.6 native RT transaction
Status: bounded engineering repair implemented; physical AGN SED and production
live qualification remain unapproved.

## Purpose

The previous source path consumed `agn_pending_erg` immediately after a
successful multigroup deposition, while the RT transaction snapshot was made
after that deposition. A later coupled rollback could therefore restore the
photon state without restoring the accepted-event fuel marker. This repair
keeps the accepted event pending in memory until the complete source →
RT/chemistry → dust commit succeeds. It does not provide durable restart
retry after RAMSES reaches its terminal failure path.

## Implementation

- `patch/lagRamses/snrt_ramses_driver.f90` now begins the native RT transaction
  before the AGN source loop, so its snapshot contains the pre-source photon,
  H/He, neutral-H and thermal state.
- A per-call `source_transaction_ok` mask records only source deposits whose
  complete spectral-group validation and atomic deposition succeeded. The
  source loop no longer clears `agn_pending_erg`.
- Every coupled rollback path releases the mask and restores the pre-source
  persistent state. The pending energy consequently remains uncleared until
  the process reaches a successful coupled commit.
- After global RT commit and the live-dust commit, the driver clears pending
  energy only for the source entries recorded in that mask. The helper
  `snrt_agn_source_commit` now documents this caller contract.
- The active source remains the `patch/lagRamses` implementation. With the
  production `bin/Makefile` VPATH, `patch/lagRamses/snrt_ramses_driver.f90`
  is selected first; no duplicate SNRT driver mirror exists in
  `patch/cuRamses`, and no mirror was edited.

## Native evidence

The existing Intel native AGN runner was extended to compile
`snrt_rt_transaction.f90` and to exercise the coupled ordering:

```
SNRT_AGN_SOURCE_COUPLED_ROLLBACK_PENDING_PASS
SNRT_AGN_SOURCE_COUPLED_COMMIT_PENDING_CLEAR_PASS
SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK helper=compiled source_api=compiled runtime=disabled
```

The rollback case snapshots an empty photon state, deposits both source
groups, restores the transaction, and verifies exact photon restoration and
unchanged pending energy. The retry case deposits again, commits the coupled
state, then clears pending energy exactly once.

Production-source compile checks also passed:

```
make -C bin SNRT=1 USE_CUDA=1 DUST_LIVE=0 snrt_ramses_driver.o
mpiifx ... -DSNRT -DDUST_LIVE -DNVAR=30 -c \
  ../patch/lagRamses/snrt_ramses_driver.f90
```

The first command rebuilt the active SNRT module dependencies and compiled the
driver. The second explicitly compiled the same active driver with the
`DUST_LIVE` conditional path. Existing unrelated dirty files and scratch
directories were not staged.

After the bundle-end review, the active tree was rechecked with fresh full
links in separate `/gpfs` scratch directories:

```
make -C bin -B -j4 SNRT=1 DUST_LIVE=0 USE_CUDA=1 HDF5=1 \
  EXEC=/gpfs/kjhan/LRD_JWST/.snrt-atomicity-link.se08tt/ramses_snrt_atomicity_cpu \
  ramses
make -C bin -B -j4 SNRT=1 DUST_LIVE=1 USE_CUDA=1 HDF5=1 \
  EXEC=/gpfs/kjhan/LRD_JWST/.snrt-atomicity-dust-link.PIQ814/ramses_snrt_atomicity_dust \
  ramses
```

Both full links passed, producing separate `*3d` executables without
overwriting `bin/ramses_final3d`. The embedded source identity was
`210b611542903147a9197836d0e3f2f755506cca-dirty`; the commit component is
exact, while the worktree also contained pending audit-closure edits and an
unrelated tracked change. These are fresh working-tree links, not clean-tree
revision artifacts.

The exact-revision link condition was then verified in the clean detached
worktree `/gpfs/kjhan/LRD_JWST/.snrt-clean-b5b20e1` at commit
`b5b20e12b66b17efd943be9562224be79be1a7b8`. Serial (`-j1`) full links passed
for both `DUST_LIVE=0` and `DUST_LIVE=1`, producing
`ramses_snrt_clean_cpu3d` and `ramses_snrt_clean_dust3d`, respectively. The
link logs embed the exact clean commit identity. Serial build is used because
the initial module dependency graph is not safe for a cold `-j4` build.

## Boundary and remaining limits

This closes source-ledger consumption ordering for the current serial,
fresh-start SNRT AGN path. On a coupled failure the marker remains uncleared
until terminal stop, but it is not persisted for a later restart. It does not
approve physical AGN spectra, MAD partition choices, restart/migration
persistence, simultaneous legacy plus SNRT ownership, or an initialized live
RAMSES evolution. The existing transaction commit API and those broader
runtime qualifications remain under their governing bundle records.

## Post-audit bounded closure

Following the Opus read-only audit, the active driver was tightened without
changing the admitted serial/fresh-start scope. The accepted-energy ledger is
validated before `snrt_transaction_begin`, and the final clear now requires
both allocated arrays to cover `nsink`; this closes the audit's invalid-input
window and malformed-ledger guard conditions. The source smoke no longer
contains the stale pre-repair snapshot comment or its no-op check.

Post-edit evidence:

- `git diff --check` passed.
- `simulation/snrt/tests/run_fp15_agn_efficiency.sh` passed, including the
  coupled rollback/pending and coupled commit/clear markers.
- Active CPU (`DUST_LIVE=0`) and DUST (`DUST_LIVE=1`, `HDF5=1`) driver object
  compiles passed through the `../patch/lagRamses` VPATH entry.

These are operator rechecks after the Opus verdict, not a replacement for the
still-open driver-faithful call through `snrt_ramses_advance_level`; the next
bundle-end audit must review both the closure edits and that harness.

The older
`provenance/agn_accepted_fuel_overlap_bundle_evidence_2026-09-06.md` remains a
historical record of the pre-repair implementation. This record is the current
source-ordering evidence for the repair.
