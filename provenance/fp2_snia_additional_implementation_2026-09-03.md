# F-P2 common/Fable findings — additional implementation

Date: 2026-09-03

Status: **physical baseline approved; runtime activation remains gated**.

This record covers the implementation after the three-model F-P2 bundle audit.
It records the explicit physical baseline approval without enabling SNIa
runtime deposition.

## Implemented

- SNIa event budgets now carry actual tracked ejecta and net-yield arrays.
  Returned mass plus terminal remnant is required to equal the WD debit; a
  tracked-ejecta sum above returned mass is rejected.
- WD reservoir capacity is owned by the configured WD-producing channel
  (AGB in the default contract), not by the aggregate remnant mass.  A SNIa
  debit removes mass from that owning channel, records the SNIa terminal
  remnant in its own channel, updates tracked/untracked return ledgers, and
  recomputes living mass transactionally.
- Cell deposition carries element and total-metal densities.  Bulk kinetic
  energy includes the event-momentum cross term and event-momentum kinetic
  term; nonzero event momentum requires an explicit coupling declaration.
- The production RAMSES bridge has a row-major `unew(local_cell,variable)`
  adapter with target bounds, uniqueness, local MPI-owner checks, scratch
  transaction, weighted scatter, and chemical-field propagation.
- Legacy prompt-table SNIa is disabled by default and requires an explicit
  legacy-only namelist opt-in; channel-resolved mode rejects that opt-in.
- The DTD near-`alpha=-1` handoff uses a higher-order series branch and has
  tests on both sides of the handoff.
- The population realization JSON contract is mirrored into a complete
  Fortran namelist loader.  Binary-fraction semantics, metallicity-factor
  source identity, expectation-versus-seeded-Poisson policy, immutable
  revision format, and approval id are explicit; unapproved or incomplete
  realizations fail closed.
- The production source identity contract now includes every SNIa object that
  the `/gpfs` `bin/Makefile` actually links.  The Fable historical
  reproduction test now distinguishes resolved F3/F4/F7/F8 from still-open
  findings.
- The approved F-P2 baseline is Maoz, Mannucci & Brandt (2012)'s field
  power-law DTD (`alpha=-1.07`, `0.04--13.7 Gyr`, `1.3e-3 events/Msun`) on
  the project Kroupa-like basis, with binary fraction `0.5` recorded as
  metadata and not multiplied into the empirical rate.
- HESMA record `yysd4-xap92`, model `n100`, is the explicit event source.
  Its integrated stable-element mass is used for returned mass and WD debit,
  its documented inner-zero/outer-half-bin profile estimate supplies
  `1.5063100005966762e51 erg`, and unresolved isotropic ejecta use a zero net
  source-frame momentum vector.  The deterministic promotion tool emits the
  checksum-bound approved event asset.
- Approval id `FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ` binds the source,
  population, conversion code, and source staging commit.  The approval is
  production-ready as a physical baseline but not publication-ready for net
  nucleosynthesis or metallicity trends; runtime remains disabled until the
  tested bridge is connected to the real AMR/MPI caller.

## Verification

- `run_fp2_snia_dtd_contract.sh`: PASS, including native and production
  contract tests, HESMA admission negatives, and mirror compilation.
- `run_g2_population_ledger.sh`: PASS.
- `run_stellar_feedback_policy_unit.sh`: PASS.
- `run_stellar_residual_deposition_unit.sh`: PASS, including the row-major
  RAMSES adapter and nonlocal-owner transactional rejection.
- `run_g1_native_contract.sh`: PASS; JAX 0.11.1 CPU differential has six
  zero-difference queries.
- Production compile targets `stellar_ramses_runtime.o`,
  `stellar_ramses_bridge.o`, and `feedback.kjhan3.o`: PASS.
- `/gpfs` production-linked RAMSES build and source parity: PASS; the
  resulting binary remains a smoke/build artifact, not a SNIa physics run.
- `test_fable_sn_agn_reproduction.py`: PASS; current disposition is 11
  reproduced, 2 partial, and 4 resolved historical findings.

## Still blocked

The approved source is not yet called by an AMR/MPI neighbour-selection
routine, so the runtime still refuses SNIa activation.  Full net-yield and
metallicity-sensitive population extensions remain open; no large RAMSES run
is authorized by this change.
