# Claude Opus 5 P0.3 HDF5 restart-state re-audit

Date: 2026-09-02  
Project root: `/gpfs/kjhan/LRD_JWST`  
Auditor: Claude Opus 5 (read-only independent gate audit)  
Scope: P0.3 HDF5 stellar restart state, linked RAMSES fixtures only; not a
production simulation

## Verdict

**CONDITIONAL PASS — P0.3 remains open.**

The closure actions addressed the original high-risk implementation findings,
but the evidence still does not establish continuation equivalence for actual
`PTYPE_STAR` particles. The gate must not be promoted to closed until the
remaining runtime evidence is collected.

## Confirmed closure actions

- The linked HDF5 writer now emits non-empty `tpp`, `mp0`, and `indtab`
  datasets for three active particles.
- The linked reader restores distinct nonzero values and exposes a value-level
  assertion for the restored state.
- A `ncpu_file=1 -> ncpu=2` reader run succeeds.
- Linked Fortran negative runs reject a missing schema marker and a bad
  `indtab` extent with nonzero exit status.
- Release-critical HDF5 reads now use a checked helper that validates rank,
  global extent, hyperslab bounds, read status, and close status.
- The binary and HDF5 writer fixtures compare the serialized release fields
  byte-for-byte for their common fixture.

## Remaining findings

1. **HIGH — no true `PTYPE_STAR` continuation equivalence.** The fixture
   currently uses sink-enabled allocation with `stellar=.true.` removed and
   RT disabled. It exercises the `star .or. sink` serialization allocation,
   but does not create a real star, advance it through feedback, restart it,
   and compare released mass/metals with an uninterrupted continuation.

2. **HIGH — writer nonzero fidelity is not yet demonstrated.** The linked
   writer starts with all-zero release fields; Python then injects distinct
   values into the temporary HDF5 checkpoint. This proves reader fidelity but
   not that the production writer serialized a nonzero stellar state.

3. **MEDIUM — binary/HDF5 comparison is currently zero-valued.** The
   comparison is structurally useful, but both writer payloads are zero and
   the parser assumes the final three records are `tpp`, `mp0`, and `indtab`
   without an independent record-order assertion. Repeat it with a nonzero
   payload and explicit ordering evidence.

4. **MEDIUM — multi-rank coverage is shallow.** The two-rank case restores
   three particles but does not exercise a zero-particle rank, per-rank
   release-value assertions, or non-contiguous active slots/`levelp` holes.

5. **MEDIUM — evidence metadata must distinguish injected data.** The JSON
   should record which values came from the linked writer and which were
   injected only into the temporary reader fixture. The tautological
   `release_state_bitwise_zero_fixture` field should be removed or replaced
   by explicit writer-payload assertions.

## Required next actions

1. Make the linked writer fixture produce at least one nonzero, distinct
   `tpp/mp0/indtab` tuple and assert it directly from the HDF5 output.
2. Add a real `PTYPE_STAR` continuation fixture: uninterrupted and
   checkpoint/restart branches must agree on released mass/metals and cursor
   state within the declared floating-point tolerance.
3. Re-run binary↔HDF5 serialization comparison with the nonzero payload and
   assert field record ordering.
4. Extend the MPI fixture with a zero-particle rank and explicit per-rank
   checks; add a non-contiguous active-slot checkpoint if the production
   particle writer can represent one.
5. Update the evidence schema and gate plans with these distinctions.

P0.3 therefore remains a **conditional pass**, and P0.4 must not be promoted
until the high-priority actions above are complete and re-audited.
