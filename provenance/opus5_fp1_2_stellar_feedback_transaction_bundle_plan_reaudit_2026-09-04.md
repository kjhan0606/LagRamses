# Claude Opus 5 plan re-audit — F-P1.2 stellar feedback transaction

- Date: 2026-09-04 (KST)
- Project: `kjhan0606/LagRamses`
- Workspace: `/gpfs/kjhan/LRD_JWST`
- Role: **primary auditor**; this re-audit follows the recorded GPT-5.6-Sol
  backup adjudication after the first Opus review stopped without a verdict.
- Audited plan:
  `provenance/fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md`

## Decision

**APPROVE WITH CHANGES.** The bundle is the correct next high-level feedback
step, its premise is real, and its scope is honest. The following mandatory
changes must be integrated before implementation.

## Mandatory changes

1. **Correct the MPI target rule.** Do not require a target to have exactly one
   local MPI owner or reject every nonlocal/virtual row. RAMSES intentionally
   deposits into virtual/reception rows and later reconciles them through
   `make_virtual_reverse_dp`; rejecting those rows would break valid boundary
   deposits and can reach the existing `MPI_ABORT` path. Restate target
   validation as resolved `son`/`ICELL_OF` index bounds and target-list
   uniqueness. Virtual/reception rows are legal; cross-rank atomicity remains
   deferred. The in-lock check must therefore validate the increment/local-row
   invariant, not positivity of a partial virtual row. Positivity checks are
   restricted to owned rows.

2. **Connect and complete the field map.** The runtime currently has no full
   field-map validation and still hardcodes energy `5` and momentum indices
   `(/2,3,4/)`. Port the existing
   `simulation/snrt/native/phase0/stellar_ramses_field_map.f90` delayed-cooling
   and non-overlap logic into the patch implementation, add
   `stellar_ramses_field_map.o` to the runtime prerequisite list, and explicitly
   record the stale-copy trap where `energy_index=inener`. The production
   energy slot must be `ndim+2`, with `inener` excluded.

3. **Specify synchronization.** Choose a concrete mechanism rather than
   leaving “critical or equivalent” open. Preferred implementation: hashed
   `omp_lock_t` keyed by target cell, with one named global `critical` fallback
   if the local build cannot support the lock array. Add a one-thread versus
   N-thread comparison over the same particle set, equal to round-off.

4. **Resolve dimensional consistency.** The current path writes momentum only
   through `1..ndim`, sums all three components in kinetic energy, and the SNIa
   bridge writes three components. Assert `ndim==3` for this bundle or
   generalize momentum and kinetic-energy expressions together; add this to the
   acceptance gate.

## Confirmations and cautions

- The independent generic/SNIa kinetic-energy treatment is already correct in
  the current paths; the refactor must preserve it and must not introduce
  merged-net-momentum energy loss. B2 is preservation, not a new physical fix.
- `idelay` defaults/aliases to `imetal` on this metal-enabled Phase-0 path, so
  full non-overlap validation is defense-in-depth rather than a claimed live
  bug closure.
- The delayed-cooling field remains an SNII-only returned-mass tracer, and the
  reduced-chemistry residual remains mapped to the generic metallicity field.
- The process-crash journal and cross-rank atomicity limitations are accurately
  stated and must remain in the closure record.

## Required acceptance evidence

1. Complete scratch delta for `unew(cell,var)`, `mp_after`, and `indtab_after`,
   with no mutating SNIa scatter during preparation.
2. Index-bounds/uniqueness target validation accepting virtual rows; increment
   revalidation inside the lock; owned-row positivity only.
3. Connected full field map with density, `ndim+2` energy, momentum, `idelay`,
   total metal, and active elements; `inener` explicitly excluded.
4. Concrete cell-lock/fallback synchronization, one commit boundary, and no
   fallible call after the first shared write.
5. Failure injection leaves `unew`, `mp`, and `indtab` byte-identical.
6. Mixed generic+SNIa closure including opposed momenta and zero-mass rejection;
   two same-cell concurrent sources retain both increments.
7. One-thread/N-thread equivalence, repeated-age exact no-op, focused build and
   syntax evidence.
8. Explicit `ndim==3` assertion or consistently generalized expressions.
9. Negative evidence for hard-crash journaling and cross-rank atomicity.
