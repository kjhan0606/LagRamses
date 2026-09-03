# Claude Opus 5 audit — P0.2 time and cumulative-interval semantics

Date: 2026-09-02  
Repository: `/gpfs/kjhan/LRD_JWST`  
Scope: P0.2 / G1.1 age units, RAMSES time conversion, cumulative interval
semantics, and production progress transaction.  The audit was read-only.

## Verdict: CONDITIONAL PASS

The implementation and test evidence satisfy the requested P0.2 algorithmic
contract, but the gate is not closed because the two newly promoted
production modules are still untracked in the current worktree.  The linked
binary is valid evidence for this worktree, but it is not reproducible from a
commit until the modified patch set is committed and the production evidence
is recorded against that commit.

## Confirmed by the audit

- `stellar_yield_tables.f90` converts on-disk `age_yr` to the in-memory
  `age_gyr` field exactly once.  The embedded generator emits the same Gyr
  convention.
- `stellar_native_units.f90` centralizes the RAMSES convention and retains an
  explicit `aexp**2` denominator.  Runtime code forms `texp-tpp`, seeds the
  progress state from `indtab`, and converts both endpoints.
- Production and native source increments evaluate the cumulative source at
  `previous_age_gyr` and `current_age_gyr`, then apply
  `C(current_age_gyr)-C(previous_age_gyr)`.  No forward `age+timestep` call
  remains.
- The native G1 test covers the non-uniform `0→0.25→1` telescope, zero-width
  interval, restart/repeat, stale-age rejection, and abort/retry behavior.
  `G1_NATIVE_CONTRACT_RUN_OK` and `G1_NATIVE_JAX_DIFFERENTIAL_OK` passed.
- The production `-B` build, link, smoke, and parity evidence passed as
  `P0_PRODUCTION_LINKED_CONTRACT_OK`; progress aborts on runtime failure paths
  and exports `indtab` only after gas and particle writes.

## Required closure action

Commit the two new production modules
`patch/lagRamses/stellar_native_units.f90` and
`patch/lagRamses/stellar_progress_contract.f90` together with the modified
source/config/Makefile set, then re-record the production-linked evidence.
No simulation launch is authorized by this audit.

## Later-gate caveat

The current runtime re-derives both endpoint ages using the current `aexp`.
For a cosmological run, exact cross-step physical-age telescoping requires
the restart/time-history contract to define how the prior endpoint's
expansion-factor history is persisted.  This is retained as an explicit
follow-up for the HDF5 restart/time-history gate, along with yield physics,
MPI/AMR restart, and deposition qualification; it is not evidence that the
per-call interval subtraction itself is forward-shifted.
