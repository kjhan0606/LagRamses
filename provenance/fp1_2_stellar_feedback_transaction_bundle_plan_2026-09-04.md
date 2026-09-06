# F-P1.2 stellar feedback transaction bundle plan — 2026-09-04

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Status: **Claude Opus 5 final bundle re-audit returned PASS; F-P1.2 engineering
bundle closed; deferred follow-up remains outside this bundle**

## Purpose and current gap

The next high-level feedback task is to make the stellar source-to-cell update
transactional.  `stellar_ramses_runtime:deposit_one_star` currently stages
the SNIa bridge separately, then updates the generic stellar mass, momentum,
total-energy, delayed-cooling, metallicity, and element fields directly in
`unew`; in particular, the SNIa path can mutate real `unew` before later
generic-source validation.  A later error or progress failure can therefore
leave a partially applied source, and the generic source does not have the
same explicit row-major scratch/scatter boundary as the SNIa path.  The
caller is also OpenMP-parallel (`feedback.kjhan3.f90`), so the commit must
provide a named synchronization boundary rather than rely on a serial-local
ownership assumption.

This bundle closes that software correctness boundary without selecting a new
yield source or claiming physical closure.  It is directly part of the
production-ready RT/stellar-AGN-feedback/dust objective because thermal and
momentum feedback must be applied exactly once and with a declared energy
owner before live coupling is eligible.

## Scope

### A. One source transaction for generic stellar and SNIa channels

1. Build a complete target-cell delta in scratch storage with the RAMSES
   row-major `unew(local_cell,variable)` convention.  The transaction state
   must include the staged `unew` delta and the proposed `mp_after` and
   `indtab_after`; no real state is changed during preparation.
2. Build the SNIa contribution through a non-mutating cell-increment builder
   or equivalent scratch path.  The mutating SNIa scatter routine must not be
   called while preparing this mixed generic+SNIa transaction.
3. Validate field indices and non-overlap, source/budget finiteness,
   non-negative returned and tracked mass, source-cell volume, and unit
   conversion.  The target must resolve through the existing `son`/
   `ICELL_OF` path to an addressable row in `unew` and appear exactly once in
   the target list; virtual/reception rows are legal RAMSES targets and must
   remain eligible for the existing `make_virtual_reverse_dp` reconciliation.
   Reject unresolved, duplicate, and out-of-range targets.  Do not reject a
   virtual row merely because it is not owned by the current rank; cross-rank
   atomicity is explicitly deferred.  Validate the full post-update field map
   (density, total energy, momentum, delayed cooling, total metal, and every
   active element) before the first shared mutation.
4. Prepare progress completion locally, including the proposed particle and
   progress values.  Once shared mutation begins, the commit path contains no
   fallible call or error exit.  The first shared mutation must occur under the
   selected synchronization mechanism: the hashed `omp_lock_t` keyed by target
   cell is the implemented and tested mechanism.  A named
   `stellar_feedback_transaction` OpenMP critical fallback is not implemented
   or tested in this bundle and is deferred; no fallback behavior is claimed.
   Re-read
   the current `unew` row under the lock and validate increment/local-row
   invariants against that row; positivity checks apply only to owned rows,
   because a virtual row can contain a partial ghost contribution.  Apply the
   complete row plus `mp`/`indtab` exactly once.  No new distributed or
   cross-rank atomic transaction is claimed.
5. Preserve the current no-journal limitation: a hard process crash between
   the hydro write and persisted particle/progress state is not closed by this
   bundle and remains explicit negative evidence.

### B. Thermal and momentum ownership contract

1. Define `source%energy` as the source event/internal energy contribution in
   the declared gas/source frame, separate from the kinetic energy of each
   returned-mass component and its source momentum.
2. Deposit momentum once in that declared frame.  Preserve the already
   correct per-channel bulk kinetic, cross, and source-momentum kinetic terms
   for each independent returned-mass component before summing staged deltas.
   Do not compute kinetic energy from a merged generic+SNIa net momentum,
   which loses counter-streaming energy.  Reject a nonzero source momentum
   with zero corresponding returned mass, or fail closed for any unsupported
   nonzero generic-momentum case; B2 is a preservation requirement, not a
   claim of a new physical correction.
3. Deposit the resulting total hydro energy exactly once.  Do not add the
   same event energy to a second thermal or delayed-cooling receiver.
4. Retain the current delayed-cooling field convention as an explicitly
   channel-owned SNII returned-mass tracer used by the existing threshold and
   decay prescription, not as a blast-energy reservoir.  Winds, AGB, SNIa,
   and PISN must not populate it.  A new delayed-cooling energy model or its
   physical calibration is outside this bundle and will not be inferred from
   a yield table.
5. Retain the reduced-chemistry rule: tracked ejecta enter their mapped
   element fields, while the untracked residual enters only the generic
   metallicity field.  `net_yield` remains diagnostic and is never deposited
   as gas mass.
6. Make the field layout explicit in the implementation and evidence:
   `unew(cell,variable)` is row-major, and the hydro total-energy slot is
   `ndim+2`, not a passive `inener` slot.  The complete map, including
   `idelay`, must be checked for valid non-overlapping indices.  Port the
   existing delayed-cooling/non-overlap implementation from
   `simulation/snrt/native/phase0/stellar_ramses_field_map.f90` into the patch
   path, connect it to runtime validation, and add
   `stellar_ramses_field_map.o` to the runtime prerequisite list.  The stale
   `energy_index=inener` copy is a named implementation trap and is forbidden.
7. Resolve the dimensional contract in one way: this bundle asserts
   `ndim==3` at the runtime boundary, matching the `NDIM==3` feedback caller
   and the three-component SNIa bridge.  A future generalized `ndim` path is
   outside this bundle; no mixed three-component energy with lower-dimensional
   momentum may be introduced.

### C. Evidence and source/build wiring

The bundle adds or updates:

- a native Fortran transaction/smoke test with successful mixed-source
  staging, invalid-field or non-finite rejection, builder zero-delta failure
  behavior plus the architectural no-production-state mutation boundary,
  mass/energy/momentum closure, opposed-momentum and zero-mass checks,
  channel-owned delayed-reservoir checks, a pinned multi-thread full-row
  same-cell case, and the `ndim==3` contract;
- a runtime/source-order audit proving that all checks precede the one commit,
  that the commit path has no fallible post-write call, and that `mp`/`indtab`
  progress follows the commit;
- a Python/source-contract regression for row-major orientation, full
  field-map non-overlap, failure-injection byte identity in the independent
  transaction model (with production-state non-mutation enforced by the
  builder signatures), virtual-row acceptance, and two same-cell sources
  retaining both increments;
- the `/gpfs` Makefile prerequisite/link graph and focused syntax/compile
  evidence where the local compiler permits it;
- a negative evidence record for the unclosed hard-crash journal and for any
  physical source/energy prescription not supplied by the approved source
  package.

## Explicit non-goals

- no G2 physical yield/source selection, canonical physical nodes, or
  40–120 M☉ fate promotion;
- no new stellar SED, AGN SED, obscuration, dust, scattering, IR, or
  radiation-pressure physics;
- no change to the runtime activation gate and no RAMSES production run;
- no global AMR/MPI neighbour transaction or generic checkpoint/HDF5 repair;
- no hard-crash exactly-once journal claim and no invention of SNII delayed
  cooling energy, momentum, or fallback values.

The OpenMP critical/cell-lock scope closes concurrent in-process same-cell
updates only.  Virtual/reception rows remain legal and are reconciled by
RAMSES's existing reverse virtual-cell exchange, but this bundle does not
claim atomicity across MPI ranks, a process crash, or distributed neighbour
deposition; those remain explicit long-term work.

## Acceptance gate

The bundle is complete only when the source transaction, field-map, and
channel-ownership tests pass; the production source shows one staged commit
boundary; the `/gpfs` source/build evidence is refreshed; runtime activation
remains fail-closed; and Claude Opus 5 completes a read-only bundle-end audit.

The result will be reported as engineering/transactional closure only. It
does not make the stellar physical-source, dust, live-coupling, or publication
gates green.  The amended plan incorporates the Opus primary verdict and the
Sol backup conditions; the implementation may now begin.  If a later Opus
audit cannot issue a bundle-end verdict, the governed Sol fallback rule still
applies and the technical failure is recorded.

GPT-5.6-Sol is not an automatic independent auditor. Under the current
governance it is called only if Opus 5 does not perform the audit/cannot issue
a verdict, or if the operator explicitly requests a separate confirmation.
