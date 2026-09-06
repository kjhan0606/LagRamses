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
keeps one accepted event retryable until the complete source → RT/chemistry →
dust commit succeeds.

## Implementation

- `patch/lagRamses/snrt_ramses_driver.f90` now begins the native RT transaction
  before the AGN source loop, so its snapshot contains the pre-source photon,
  H/He, neutral-H and thermal state.
- A per-call `source_transaction_ok` mask records only source deposits whose
  complete spectral-group validation and atomic deposition succeeded. The
  source loop no longer clears `agn_pending_erg`.
- Every coupled rollback path releases the mask and restores the pre-source
  persistent state. The pending energy is consequently retained for a retry.
- After global RT commit and the live-dust commit, the driver clears pending
  energy only for the source entries recorded in that mask. The helper
  `snrt_agn_source_commit` now documents this caller contract.
- The active source remains the `patch/lagRamses` implementation. With the
  production `bin/Makefile` VPATH, `patch/lagRamses/snrt_ramses_driver.f90`
  is selected before the older `patch/cuRamses` copy. No mirror was edited.

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

## Boundary and remaining limits

This closes source-ledger consumption ordering for the current serial,
fresh-start SNRT AGN path. It does not approve physical AGN spectra, MAD
partition choices, restart/migration persistence, simultaneous legacy plus
SNRT ownership, or an initialized live RAMSES evolution. The existing
transaction commit API and those broader runtime qualifications remain under
their governing bundle records.

The older
`provenance/agn_accepted_fuel_overlap_bundle_evidence_2026-09-06.md` remains a
historical record of the pre-repair implementation. This record is the current
source-ordering evidence for the repair.
