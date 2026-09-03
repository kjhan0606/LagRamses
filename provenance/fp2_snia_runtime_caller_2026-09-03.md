# F-P2 runtime caller handoff — 2026-09-03

Status: **caller connected; SNIa activation remains disabled**.

## Scope

The approved F-P2 Maoz/HESMA baseline is now connected to the production
stellar runtime under `/gpfs/kjhan/LRD_JWST`.  This change does not relax the
independent `production_source_model_supported()` gate and does not authorize
a RAMSES evolution run.

## Runtime path

1. `phase0_initialize` requires the explicit
   `PHASE0_SNIA_RUNTIME_CONTRACT` handoff and reads the ordered population,
   physical-event, and thermal-coupling groups.
2. All three groups must validate and share the approved source commit and
   named approval id `FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ`.
3. `deposit_one_star` evaluates the interval-integrated DTD, debits the
   channel-owned WD reservoir transactionally, and checks combined generic +
   SNIa mass closure.
4. `locate_star_cell` resolves the actual AMR leaf cell through RAMSES's
   `get3cubefather`/`ICELL_OF` path.  The current runtime policy is one-cell
   NGP with owner rank `myid-1`; it is not yet a neighbour stencil or MPI
   exchange.
5. `deposit_snia_budget_to_unew` validates local ownership and bounds,
   constructs variable-major scratch increments, then atomically scatters to
   RAMSES's row-major `unew(local_cell,variable)` layout.

## Restart and idempotence boundary

For the approved baseline, `terminal_remnant=0` and
`returned_mass=WD_debit`.  A repeated call at an already committed age is a
no-op through `stellar_progress_contract`.  On a normal restart, the prior
SNIa return is reconstructed from persisted particle mass and the generic
cumulative ledger; no second in-memory cursor is assumed.  This invariant is
rejected for a future nonzero-remnant model until a versioned per-particle
debit payload is added.

This is normal retry/restart idempotence, not a hard-crash exactly-once claim:
the existing RAMSES hydro-array write and `indtab` checkpoint commit are not
one atomic transaction.  A pending-event journal is retained as a separate
production-hardening task.

## Evidence

- F-P2 native contract runner: ordered handoff, DTD, physical contract, cell
  deposition, event ledger, source admission, and production mirrors pass.
- Production `/gpfs` Makefile: `stellar_ramses_runtime.o`,
  `stellar_ramses_bridge.o`, and `feedback.kjhan3.o` compile and link.
- Source parity: native and production physical/cell/driver mirrors match.
- Runtime handoff audit: contract, caller, and bridge are present and linked;
  activation remains explicitly false.
- Existing bridge unit evidence covers weighted two-cell conservation,
  scale conversion, duplicate-target rejection, and nonlocal-owner rejection.

## Qualification result

The bundled qualification is now recorded by the automated evidence:

- normal retry/restart reconstruction passes in
  `fp2_snia_runtime_accounting_test.f90`;
- weighted multi-cell conservation, duplicate-target rejection, and
  non-local-owner rejection pass in the RAMSES bridge unit test;
- the linked production binary rejects missing/mismatched runtime contracts
  and reaches the independent production gate for a valid handoff, without
  entering evolution.

The runtime remains disabled.  A hard-crash exactly-once guarantee is not
claimed because RAMSES `unew` writes and the existing `indtab` checkpoint
commit are not one atomic transaction.  Full net-yield and
metallicity-sensitive population extensions remain publication work.
