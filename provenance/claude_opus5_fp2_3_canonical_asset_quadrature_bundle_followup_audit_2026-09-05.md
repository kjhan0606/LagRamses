# Claude Opus 5 follow-up audit — F-P2.3 — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Auditor: Claude Opus 5 (`claude-opus-5`)
Mode: read-only; file reads/content search only; no edits, shell, tests, jobs,
hash recomputation, or network
Audit session: `846f4fb1-e17a-434a-a276-15b77dd85970`

## Verdict

`PASS` — for the declared engineering/provenance scope of F-P2.3 only.

Both blocking findings from the first same-day F-P2.3 audit were genuinely
repaired. The listed medium/minor repairs introduced no new mismatch between
builder, validator, tests, and serialized artifacts. This does not grant
astrophysical, obscuration, dust-physics, production, publication, or live
RAMSES approval.

## Closed findings

- The publication-facing `simulation/snrt/AGN_NINE_GROUP_VALIDATION.md` now
  carries the current date/status, current pilot and explicit hashes, complete
  reproduction commands, the explicit engineering-control section, and the
  `synthetic_non_physical_wiring_fixture` boundary.
- The explicit validator now compares reconstructed group means and energy
  fractions with serialized metadata. Both new criteria are true in the
  regenerated explicit artifact.
- Empty-group support intervals are derived from SED support, the explicit
  limit text agrees with enforced closed-fraction behavior, the evidence
  distinguishes artifact-independent from implementation-independent checks,
  and empty-group representative energies are geometric means as documented.
- The manifest includes an expected-failure reproduction record for the
  preserved full-CFL probe, the convergence metric is described as a
  base/refined difference rather than an error bound, and pilot/explicit
  artifact tests both rerun the validator in a temporary output and compare
  the complete criteria vector.

## Independent confirmation

The auditor independently checked the flat explicit `f_E` fixture and found
the represented fraction, group rates, group means, hard/soft ratios, and the
reported `6.32e-7` group-0 refinement difference consistent with closed-form
expectations. It confirmed the algebraic equivalence of the `sed.py` and
`primordial.py` mean-energy definitions, the shared `f_E -> q_E` convention,
and the offset-grid/closed-form dust oracles.

It also confirmed that all ten manifest assets are checked for exact ID/path,
live size, and live SHA-256; both artifact tests use repository-root-relative
`simulation/snrt` attestation, current HEAD, hash sweeps, and fresh validator
runs; the dirty tree is honestly recorded; and the pilot/explicit transport
controls use the documented 32 time-averaged absorption iterations.

## Non-blocking findings carried forward

- The validation markdown records the two validation-JSON hashes without a
  checker; either add a checker or treat them as informational.
- The working-tree status digest is path-level rather than content-level;
  content freshness is supplied by the explicit SHA-256 sweeps, and a clean
  tree must be required before publication.
- The explicit artifact test is a consistency gate rather than a separate
  numerical-value pin beyond its all-criteria result.
- Explicit status expectations still derive from the luminosity-scaled rate;
  the current fixture has escape fraction one. A future physical-source
  contract should derive this from the per-normalization moment.
- Per-group dust status tokens and the pilot null-identity quadrature record
  can be tightened later; `snrt_core/sed.py` is covered through the explicit
  closure manifest.

These are non-blocking engineering follow-ups. Physical AGN/stellar SED
selection, out-of-range SED-fraction policy, mixed STAR+AGN dust, dust-to-metal
normalization, `[40,120] M_sun` yields, scattering/IR/grain evolution, live
Fortran RAMSES coupling, B3 timestep/spatial convergence, and publication
deposit identifiers remain later bundles.

## Acceptance gates

All F-P2.3 gates were judged PASS: complete asset manifest; symmetric artifact
freshness/attestation; explicit SED rates, means, fractions, and Verner closure
reconstruction; photon/Verner/dust quadrature convergence; fail-closed
negative paths; passing pilot/explicit validators; and the reported focused
tests, compilation, and diff check.

No next bundle was started by the auditor.
