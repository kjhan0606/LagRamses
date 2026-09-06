# Claude Opus 5 final B2 re-audit

Date: 2026-09-01  
Mode: read-only; logic and recorded values inspected, calculations not rerun  
Scope: H-only B2

## Verdict: PASS

**B2-1 closed. B2-2 closed. No new B2 blocker.**

## B2-1 — source-hash closure

- The validator writes `provenance.snrt_core_sha256` as a sorted mapping over
  every current `snrt_core/*.py`; the canonical artifact contains 32 entries.
- The artifact test reconstructs that mapping and requires dictionary
  equality. A core edit changes a value; an addition or removal changes the
  key set. The validator in `tools/` is anchored separately.
- Every module named by the prior audit is covered: `implicit`, `transport`,
  `primordial`, `dust`, `secondary`, `primordial_cooling`, `quadrature`,
  `shadow`, and `sources`.
- No external data file determines this B2 result. The dust fixture is built
  inline and does not call the sidecar reader.
- The B2 document now claims exactly the enforced scope: validator plus all 32
  current core modules.

The audit did not independently recompute hash values; it verified coverage,
reconstruction logic, and internal consistency. Local artifact tests recompute
the values.

## B2-2 — worst-of-five headline

The result table reports the true worst H-ledger L1 relative error,
`5.60382e-5` from `secondary_200ev_on`. The artifact values are baseline
`1.945125e-6`, dust `2.302548e-6`, secondary-off `4.953354e-5`, secondary-on
`5.6038197e-5`, and Solver B `2.048173e-6`. The residual row `2.39611e-5` is
also the true maximum, tied by baseline and Solver B.

## Additional verified improvements

- Dust H-ledger normalization uses gas-absorbed photons.
- All five run dictionaries are asserted to record at least 20 iterations;
  the CLI rejects fewer than 20.
- The one-iteration inert-gas shadow is disclosed and matches the zero gas
  cross sections.
- Zero-H conservative-primordial coverage is present.

## Later-gate improvements, not B2 blockers

- Switch from non-recursive `glob` to `rglob` if `snrt_core` gains subpackages.
- Bind or separately certify the recorded environment for backend-sensitive
  float32 replay.
- Add a refinement datum for the 4.41% Strömgren-radius error.
- Record shadow opacity-iteration count directly in the artifact.
- Add a footnote on secondary-ionization count normalization.

This PASS is limited to H-only B2. Helium coupling, the physical dust
prescription, source SEDs, thermal-atlas licensing, yields, and live RAMSES
coupling remain later independent gates.
