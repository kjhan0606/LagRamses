# Claude Opus 5 end-of-bundle audit — F-P2 source SED/dust closure — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Audit mode: read-only; `Read/Grep/Glob` only; no repository files, jobs, or
builds were modified or launched by the auditor.
Auditor: Claude Opus 5 (`claude-opus-5`), session
`474c1ab3-a9a1-40e3-81bc-24da0b91ccc2`.

## Overall verdict

**CONDITIONAL PASS.** The source-SED identity architecture and the
single-source SED → photon-ledger → gas/dust closure are mathematically
coherent, and the fail-closed runner boundary is real. The bundle is not yet
closed as an engineering gate because the explicit-SED metadata contains
pilot-only claims, the source-weighted Draine closure lacks an independent
numerical test, and recorded provenance hashes are not all enforced by the
loader. No blocker was found that invalidates the bundle concept; these
findings block an unqualified bundle close and any promotion claim.

## Findings

### MAJOR-1 — explicit-SED AGN metadata is not truthful (blocks this bundle)

In `simulation/snrt/tools/p4_build_agn_photon_ledger.py`, the explicit SED
path still derives `sed_supported_interval_ev` from the pilot constant
`SED_MIN_EV=10 eV`, retains the Sazonov citation in the common `reference`
field, exposes the unused pilot `lyman_nu_lnu_fraction` normalization, and
retains pilot-only limit text. This conflicts with the explicit path's actual
full group-edge support and makes output provenance ambiguous. The existing
source-SED test checks identity/hash and energy-fraction sum but does not catch
these fields.

### MAJOR-2 — source-weighted v2 Draine numerics lack independent verification
(blocks this bundle)

`build_source_weighted_opacity_metadata` is exercised only for schema,
identity, and positivity. The test SED and Draine samples lie on group edges,
so the interpolation/union-grid calculation is not meaningfully exercised and
there is no independent recomputation of
`∫q_E κ_abs dE / ∫q_E dE` or
`∫E q_E κ_abs dE / ∫q_E κ_abs dE`. This is listed in the plan's acceptance
evidence but is not yet implemented.

### MAJOR-3 — recorded opacity provenance is not enforced (blocks the claim)

The v2 sidecar records `group_edges_sha256` and the Draine
`source_table.sha256`, but `read_dust_opacity_metadata` does not validate
either. The v2 sidecar also lacks a builder hash, and the closure arrays are
trusted after structural validation. A self-consistent hand-authored v2
sidecar could therefore pass without proving that its opacity arrays came from
the declared Draine table and group-edge file. The identity and numerical edge
matching are enforced; the source-table/builder provenance is not.

### MAJOR-4 — mixed STAR+AGN admission is incomplete (mostly next bundle)

The merger's aggregate photon/absorber/excess-energy arithmetic is correct and
component-only dust is prevented by fail-closed identity mismatch. However:

1. the aggregate identity cannot be reproduced by the existing source-bound
   dust builder, so no aggregate v2 dust sidecar can currently be supplied;
2. the merged CSV drops `aexp`, preventing the normal static-input attachment
   path;
3. the mixed identity hashes component CSVs but not the metadata JSON that
   supplies the aggregate closure arrays; and
4. a mixed result can be labeled candidate even when a component is only a
   null-identity reference control.

Items 1–2 belong to a subsequent mixed-source admission bundle. Items 3–4
should be repaired before treating the mixed metadata as a complete contract.

### MINOR findings

- The strict `2e-6` group energy-closure tolerance can reject a legitimate
  coarsely sampled SED because linear interpolation of `q_E` at an interior
  edge is not exactly linear in energy fraction.
- Loader defaults are permissive (`require_source_match=False`), and the
  runner echoes raw JSON status instead of the validated closure status.
- `escape_fraction` is folded into total photon rates but not the per-Lbol
  metadata fields, so non-unit escape fractions can confuse recomputation.
- The plan's negative-test list is missing explicit raw-SED-hash and v2
  group-edge mismatch cases.

### Notes

- The grouped closures are weighted by the emitted, unattenuated spectrum;
  within-group hardening in optically thick cells remains a later transport
  approximation.
- The path-free identity is intentionally bound to SED bytes/contract fields,
  while group edges are bound separately.
- Empty groups use zero closure and a geometric-mean representative energy;
  this is acceptable for inactive groups and is recorded by the caller.

## What is closed and what remains deferred

Closed for the candidate engineering boundary: dimensional SED conversion,
support and group-boundary handling, shared explicit-source weighting for
photon/H/He gas closure, source-bound v2 dust schema and runner identity
matching, preserved null-identity v1 controls, and deterministic mixed-source
gas aggregation with duplicate-ID rejection.

Deferred correctly: BPASS/CCSN/AGB physical source approval and the
`[40,120] M_sun` seam; AGN SED/escape/obscuration and LRD calibration; dust
scattering, grain temperature, IR re-emission, destruction/growth, full
radiation pressure; live RT–RAMSES feedback; production convergence; and all
publication claims.

The auditor did not recompute local shell hashes in its restricted session.
The local `/gpfs` verification independently reran the focused tests,
`git diff --check`, canonical AGN validation, and the listed SHA-256 commands;
those results are recorded in the implementation evidence file.

## Recommended next coherent bundle

F-P2.1 — source-closure verification and honest metadata:

1. split explicit-SED and pilot metadata construction so support, reference,
   normalization, and limits are truthful;
2. add independent Draine source-weighted opacity recomputation using
   deliberately offset SED/table/group samples plus a closed-form power-law
   case; and
3. enforce the already recorded group-edge and Draine source hashes and add a
   v2 builder hash.

Keep aggregate STAR+AGN dust admission in a later bundle after this repair.
Per project cadence, this recommendation is not started until the driver
approves the next bundle plan.
