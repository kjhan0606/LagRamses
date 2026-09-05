# F-P2.1 source-closure verification and honest metadata bundle — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent: F-P2 source SED/dust closure; driver approval received 2026-09-05.
Status: implementation/evidence complete; Claude Opus 5 end-of-bundle audit
returned CONDITIONAL PASS; next bundle approval required.

## Objective

Resolve the three implementation findings from the Claude Opus 5 F-P2
end-of-bundle audit without widening scope:

1. make explicit-tabulated-SED AGN metadata distinct from the parameterized
   Sazonov pilot metadata;
2. independently verify the source-weighted Draine closure on non-coincident
   sample/group grids and against a power-law reference; and
3. enforce the source-Sed, group-edge, Draine-table, and builder hashes already
   recorded by a v2 sidecar.

This is a repair/verification bundle for the single-source candidate path. It
does not admit a physical AGN SED, change the dust physics scope, or fix the
separate mixed STAR+AGN admission design.

## Work packages

### R1 — honest AGN metadata

- Build pilot and explicit metadata blocks separately.
- Derive explicit group support from the actual validated SED rather than the
  pilot's 10 eV constant.
- Keep intrinsic per-`L_bol` quantities and escaped emitted totals explicit so
  non-unit escape fractions remain auditable.

### R2 — independent Draine verification

- Use deliberately offset SED, Draine-table, and group-edge samples.
- Independently recompute the source-photon-weighted opacity and absorbed
  photon energy in the test, without calling the builder's integration code.
- Include a power-law SED/opacity case with closed-form group averages.
- Add negative tests for malformed/mismatched source and edge provenance.

### R3 — enforced provenance

- Add source SED input-file, Draine source-table, group-edge-file, and builder
  code hashes to the v2 contract where not already present.
- Recompute and validate those hashes in `read_dust_opacity_metadata`.
- Pass the photon-ledger group-edge hash into P4/P5 and require it for a
  source-bound v2 attachment; retain null-identity v1 reference controls.

## Acceptance gates

- Explicit AGN metadata contains no pilot-only reference, support, or limit
  claim; pilot metadata remains labeled reference control.
- Independent numerical and closed-form Draine checks pass on offset grids.
- Tampered/missing SED, Draine, edge, or builder provenance is rejected before
  a runner output is created.
- Existing P4/P5, AGN nine-group, source-ledger, and reference-control tests
  remain passing; `git diff --check` and Python compilation pass.
- No production activation, publication claim, live RAMSES run, or mixed-source
  dust admission is introduced.

After these gates, one bundled Claude Opus 5 read-only audit will decide
whether F-P2.1 is closed or conditional. No additional auditor is invoked
unless the active governance requires a backup after a missing/indeterminate
Opus verdict.

## Audit disposition — 2026-09-05

Claude Opus 5 returned `CONDITIONAL PASS`. The explicit/pilot metadata split,
independent offset-grid and closed-form Draine checks, and four-leg
source/edge/table/builder re-hashing were accepted as real repairs. The audit
found no blocker and did not call for a backup auditor. It recommends F-P2.2
for closure-code dependency manifests, sidecar-payload integrity binding, a
fixed dust-status vocabulary, propagated dust provenance attributes, canonical
explicit-SED nine-group coverage, and working-tree cleanliness attestation.
F-P2.2 is not started pending driver approval; mixed STAR+AGN admission
remains later.
