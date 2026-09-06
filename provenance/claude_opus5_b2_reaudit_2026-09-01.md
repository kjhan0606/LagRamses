# B2 re-audit — Claude Opus 5

Date: 2026-09-01  
Mode: read-only; source and recorded values inspected, calculations not rerun  
Scope: B2 only

## Verdict: CONDITIONAL PASS

All three mandatory actions from the first audit are closed. Two new defects
remain in the gate/document layer; neither touches the physics or invalidates
the recorded numbers.

## Original mandatory actions

1. **All-run convergence gates — CLOSED.** The validator builds the tuple
   baseline, dust, secondary-off, and secondary-on, then gates every maximum
   fixed-point residual and H-ledger L1 error; Solver B is gated separately.
   The artifact test mirrors those assertions. The worst values are
   `2.3961e-5` residual and `5.6038e-5` H-ledger L1.
2. **Helium scope — CLOSED.** The B2 document explicitly states that every
   front run has `n_He=0`, makes no coupled-He validation claim, and defers it
   to a later H+He gate.
3. **Retired-cap relabel — CLOSED.** The criterion and table now call zero cap
   activity and unit scale exact structural invariants.

The tightened A/B threshold, disclosed shared core, fixture magnitude bands,
CLI 20-iteration minimum, explicit Solver-B energy, sibling density floors,
and zero-density regression were also verified by source inspection.

## New B2 blockers

### B2-1 — Incomplete source-hash closure

The artifact records hashes for the validator, `multiphysics.py`, and
`conservative_hydrogen.py`, but not other result-determining modules such as
`implicit.py`, `transport.py`, `primordial.py`, `dust.py`, `secondary.py`,
`primordial_cooling.py`, `quadrature.py`, `shadow.py`, and `sources.py`. The
document's statement that any source edit forces a fresh B2 run is therefore
too broad, especially in an intentionally dirty worktree where file hashes are
the primary artifact-to-code anchor.

Closure: bind the full import closure, or narrow the documentation claim to the
three covered files.

### B2-2 — Headline ledger is not worst-of-five

The result table reports the baseline H-ledger L1 error `1.94513e-6` without a
qualifier, although the strengthened gate covers five runs and the worst value
is `5.60382e-5` in `secondary_200ev_on`. The fixed-point row already happens to
be worst-of-five.

Closure: report the worst-of-five value with its owner, or explicitly label
the row as baseline and show the worst value.

## Later-gate improvements, not B2 blockers

- Add one mesh/angular refinement datum for the radius ratio, whose 4.41%
  error has limited margin against the 5% threshold.
- Normalize the dust-run H ledger by gas-absorbed rather than total-absorbed
  photons.
- Assert each run dictionary records at least 20 iterations.
- Disclose that the inert-gas shadow needs only one opacity iteration.
- Avoid residual shadow-energy literals, improve payload invariant key names,
  and add zero-H conservative-primordial coverage.
- Retain the later plans for a relative near-neutral residual and reactive-gas
  shadow.

The algorithm, C2-Ray-style closure, ledger definitions, and five recorded run
families remain consistent with the source as read. Closing B2-1 and B2-2 is
sufficient for an H-only B2 PASS, not an overall production/publication PASS.
