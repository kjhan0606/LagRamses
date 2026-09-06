# Codex self end audit: DUST-10 native dust contract admission

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Auditor: Codex, independent end-of-bundle source/evidence review
- Date: 2026-09-06
- Verdict: **PASS within native-admission scope; CONDITIONAL for live/publication dust physics**

## Findings

1. The loader uses an explicit versioned namelist and bounded arrays.  It
   rejects missing/invalid versions, unknown statuses, incomplete source/hash
   identity, invalid array bounds, non-monotonic groups or thermal rows,
   non-positive reference dust mass, and non-finite/negative numerical data.
2. State publication is atomic at the contract level: reset occurs before
   read, and published arrays/metadata are filled only after all checks pass.
   The invalid fixture confirms that a failed load leaves
   `loaded=.false.`, `runtime_allowed=.false.`, and zero group count.
3. Candidate and production semantics are separated.  The valid fixture loads
   for inspection but cannot pass the runtime gate without both approved
   statuses and an approval identifier.
4. The existing consolidated build and runtime gates pass after the new object
   is added to the SNRT production graph.  Both available native compilers
   pass the focused loader smoke.
5. No unrelated RAMSES, feedback, AMR, CUDA, or live dust behavior was changed.
   The live driver remains `ZERO_SCAFFOLD`, as required because the current
   hydro state has no dedicated dust field.

## Residual conditions

- The native layer checks identity-token syntax and contract consistency, but
  upstream sidecar construction and actual file SHA-256 verification remain
  outside the RAMSES process.
- The tracked thermal values are test fixtures.  Physical Draine/WD01 mixture
  approval, temperature-dependent grain closure, dust-gas exchange, IR
  re-emission, radiation pressure, destruction/growth, and production
  convergence are not promoted.
- A future live adapter still needs a dedicated persistent RAMSES dust state,
  native sidecar-to-state binding, and restart/MPI semantics before replacing
  `ZERO_SCAFFOLD`.

The DUST-10 declared engineering objective is met.
