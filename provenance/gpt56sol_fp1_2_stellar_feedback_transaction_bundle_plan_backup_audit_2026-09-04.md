# GPT-5.6-Sol backup plan audit — F-P1.2 stellar feedback transaction

- Date: 2026-09-04 (KST)
- Project: `kjhan0606/LagRamses`
- Workspace: `/gpfs/kjhan/LRD_JWST`
- Role: **backup adjudication only**, invoked because Claude Opus 5
  performed repository review but did not issue a resumable final verdict.
  This is not an independent parallel audit.
- Plan audited:
  `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`

## Decision

**APPROVE WITH CHANGES.** Sol judged the bundle necessary, bounded, and
aligned with the high-level feedback objective, but required the transaction
boundary and concurrency semantics to be made explicit before implementation.

## Evidence

- The change is necessary: the current runtime sends the SNIa contribution to
  `unew` before generic-source validation/deposition, leaving a partial update
  possible (`patch/lagRamses/stellar_ramses_runtime.f90`, around line 626).
- The caller is OpenMP-parallel (`patch/lagRamses/feedback.kjhan3.f90`, around
  lines 75 and 88), so a simple serial-ownership statement is insufficient.
- Existing SNIa cell-increment construction is a useful non-mutating staging
  primitive (`patch/lagRamses/stellar_snia_cell_deposition.f90`, around line
  121).
- The delayed-cooling field is a decaying SN-tagged mass-density tracer, not a
  blast-energy reservoir.

## Mandatory amendments incorporated into the plan

1. Stage a complete `unew` delta plus `mp_after` and `indtab_after`; no
   mutating SNIa scatter may occur during preparation.
2. Use an explicit named OpenMP critical section or equivalent cell lock. The
   post-update row must be checked against the current `unew` value inside the
   synchronization boundary so concurrent stars cannot lose updates.
3. Require a valid local target and exactly one MPI owner; reject nonlocal,
   duplicate, out-of-range, and unresolved targets.
4. Prepare progress completion before shared writes. Once the first shared
   mutation occurs, the commit path must contain no fallible calls or error
   exits.
5. Define the momentum frame unambiguously and calculate kinetic/cross terms
   per independently returned-mass component, or fail closed for unsupported
   nonzero generic momentum. Do not calculate kinetic energy from a merged net
   momentum when counter-streaming components would be lost.
6. Validate the full field map, including delayed cooling, total metal, active
   elements, density, momentum, and total energy. Use the explicit layout
   `unew(cell,variable)` and the hydro total-energy index `ndim+2`.
7. Describe delayed cooling only as an SNII returned-mass tracer used by the
   existing threshold/decay prescription; do not claim this bundle validates
   its physical calibration.

## Required acceptance evidence

- Failure injection at every preparation check leaves `unew`, `mp`, and
  `indtab` byte-identical.
- Mixed generic+SNIa tests close mass, metals, momentum, and energy, including
  opposed momenta and zero-mass rejection.
- Two OpenMP sources targeting one cell retain both increments; repeating a
  committed age is an exact no-op.
- Asymmetric sentinels prove row-major orientation and complete field-map
  non-overlap.
- Source-order evidence shows one synchronized, non-failing commit boundary
  and exactly-once progress advancement.

## Deferred limitations

Physical source approval, the 40–120 M☉ fate gap, delayed-cooling calibration,
stellar/AGN SEDs, dust, radiation pressure, distributed-neighbour deposition,
persistent crash journaling, runtime production, and publication approval
remain closed/deferred.
