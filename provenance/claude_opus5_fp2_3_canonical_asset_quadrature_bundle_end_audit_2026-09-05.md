# Claude Opus 5 bundled end audit — F-P2.3 — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Auditor: Claude Opus 5 (`claude-opus-5`)
Audit mode: read-only; no tests, jobs, edits, or external network access
Audit session: `09fbacf9-a701-4190-b236-0853dc06e8a3`

## Verdict

`CONDITIONAL PASS`

The numerical core is correct, and the auditor independently checked the
explicit SED arithmetic rather than trusting only serialized values. Q1, Q2,
and the fail-closed behavior of Q3/Q4 are closed. Unconditional PASS is
blocked by one stale publication-facing document and one explicit gate that
does not independently check the group mean-energy array consumed by
transport. These are local F-P2.3 repairs, not evidence of a failed RT
algorithm.

## Independent checks reported by the auditor

For the canonical explicit flat `f_E` fixture on `[0.01,10000] eV`, the
auditor verified:

- the represented fraction is unity;
- the group-0 rate agrees with `f_E/EV_ERG * ln(100)`;
- the group-8 rate agrees with `f_E/EV_ERG * ln(5)`;
- the group-8 mean agrees with `8000/ln(5)`; and
- the reported base/refined group-0 difference agrees with the expected
  composite-trapezoid `h^2/6` behavior (`6.320e-7`).

The auditor also confirmed that the `sed.py` and `primordial.py` mean-energy
definitions are algebraically identical under
`q_E=f_E/(E*EV_ERG)`, and that the source/dust tests provide separate analytic,
closed-form, and offset-grid checks.

## Findings

### MAJOR-1 — stale publication-facing validation document (blocking)

`simulation/snrt/AGN_NINE_GROUP_VALIDATION.md` still reports the 2026-09-02
state, including a final Opus PASS and old hashes. It is the document pointed
to by `P4_AGN_RATE_LEDGER.md` and therefore contradicts the live canonical
artifacts. The live values include validator hash prefix `79a7664d1816`,
photon metadata prefix `3865aa70bd5d`, external-manifest prefix
`769e546d99c5`, and transport-control prefix `77eae4b74b67`; the old document
records different prefixes. It also lacks the explicit-mode section,
reproduction commands, and the synthetic non-physical fixture label.

Disposition: update the document with the current pilot/explicit status,
hashes, reproduction commands, quadrature/attestation contract, and explicit
engineering-only boundary. Do not retain the stale final-PASS claim.

### MAJOR-2 — explicit gate omits reconstructed group mean energies (blocking)

`tools/validate_agn_nine_group_ledger.py` reconstructs explicit rates,
quadrature convergence, and Verner closure, but does not compare
`explicit_moments.photon_weighted_mean_energy_ev` or
`explicit_moments.group_energy_fraction_per_norm` with the serialized
`groups[*].photon_weighted_mean_energy_ev` and
`groups[*].energy_fraction_of_lbol`. Transport consumes the serialized group
means, so the existing transport criterion checks the output against an
unverified array. The present values are correct, but an in-band wrong mean
could pass.

Disposition: add explicit independent `allclose` criteria for serialized
group means and energy fractions at the same declared numerical tolerance,
then regenerate both validation artifacts and rerun their artifact tests.

### MEDIUM-3 — empty explicit support interval convention mismatch

The explicit builder writes `[low,high]` for a group even when its source
photon rate is zero, while the validator currently expects `None` for every
zero-rate explicit group. This is latent for the current fixture and fails
closed, but will affect a physical SED with zero-valued portions.

Disposition: derive `sed_supported_interval_ev` from SED support, not photon
rate. Use `None` only when the group lies outside the tabulated support; retain
the full interval for a supported but zero-valued group.

### MEDIUM-4 — explicit support policy versus out-of-range photons

The explicit reader requires the SED to cover every group edge and requires
the group energy fractions to reproduce the whole represented table fraction.
Thus a table extending above `10000 eV` is rejected rather than truncated,
while one inherited limit says photons above the group table are excluded.

Disposition: make the explicit contract and limit text agree. Either require
an exactly group-covered table and say no out-of-range photons are accepted,
or design a future physical-SED contract that records the omitted fraction.
This is a prerequisite for the physical-SED bundle, not an astrophysical
approval in F-P2.3.

### MEDIUM-5 — reconstruction is artifact-independent, not implementation-independent

The validator recomputes expectations from the same integration functions
used by the builder. This is sufficient for stale/tampered artifact and
contract detection, but cannot by itself detect a shared algorithmic error.
The independent offset-grid and closed-form checks in
`tests/source_sed_dust_closure.py` are the actual algorithm-diverse evidence.

Disposition: correct the evidence wording to distinguish artifact-independent
reconstruction from independent numerical or analytic oracle coverage.

### MINOR-6 — empty-group docstring/behavior and unused argument

`primordial._closure_once` documents a geometric-mean representative energy
for an empty group, but leaves the value at zero. The
`interpolation_convention` parameter is accepted but unused in that helper.

Disposition: return the documented in-band geometric mean for an allowed empty
group and remove or use the unused parameter.

### MINOR-7 — no reproduction entry for the preserved failed probe

The manifest has reproduction commands for passing artifacts but none for
`agn9_full_cfl_failed_probe`.

Disposition: add the exact expected-failure command or explicitly document
that the probe is a preserved, non-reproducible imported artifact with its
hash and expected gate result.

### MINOR-8 — convergence difference is not an error bound

The `5e-6` base/refined criterion is a conservative refinement difference,
not a formal Richardson error bound.

Disposition: describe it as a convergence-difference gate, not an error bound.

### MINOR-9 — escape fraction is taken from the artifact under test

The explicit expected rate reconstruction reads `escape_fraction` from the
metadata under test. The value is currently one and constrained by the ledger,
but the provenance is not fully independent.

Disposition: retain for the fixture or bind the expected escape fraction in a
separate driver-side contract in the future physical source bundle.

### MINOR-10 — fresh validator rerun asymmetry

The explicit artifact test reruns the validator in a temporary output and
compares criteria with the canonical artifact; the pilot artifact test only
checks the stored artifact against live files.

Disposition: add the same fresh temporary validator rerun and criteria equality
check to the pilot test.

## Acceptance-gate disposition

| F-P2.3 gate | Opus disposition |
| --- | --- |
| All declared assets have current path, size, and SHA-256 | PASS; all ten entries are rehashed and restatted with exact ID-set checks. |
| Pilot and explicit artifacts enforce current hashes and scoped attestation | PASS; both use the actual repository root and honestly record the dirty tree. |
| Explicit support and rates are independently reconstructed | PASS with MAJOR-2 and MEDIUM-5 attached; rates/closure are reconstructed, but means/fractions need explicit criteria and implementation-diverse wording. |
| Photon, Verner, and source-weighted dust quadrature converges | PASS; all three fail closed and current maxima are below `5e-6`. |
| Negative status/duplicate/payload/manifest paths fail closed | PASS; the reported negative cases cover these contracts. |
| Canonical pilot and explicit validators pass | PASS as recorded; both artifacts pass, with the above gate gap. |
| Focused tests, compilation, and diff check | Accepted as reported; not re-executed by the auditor. |

## Scope conclusion

The audit grants no physical AGN/stellar SED approval, no obscuration/escape
calibration, no mixed STAR+AGN dust admission, no `[40,120] M_sun` yield
resolution, no scattering/IR/grain evolution, no live Fortran RHD coupling,
and no production/publication claim. MEDIUM-3 and MEDIUM-4 should be carried
explicitly into the physical-SED bundle if not repaired here.

No next bundle was started by the auditor.
