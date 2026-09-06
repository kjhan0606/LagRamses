# Claude Opus 5 G2 source-package staging final audit record — 2026-09-02

Auditor: `claude-opus-5` CLI  
Scope: read-only audit of the G2 review-only candidate source staging and
Sukhbold high-mass projection. No canonical rows, runtime deposition, commit,
push, or RAMSES integration was performed by the auditor.

## Verdict at audit time

`G2 SOURCE-PACKAGE ENGINEERING FINAL VERDICT`: **CONDITIONAL PASS**  
`G2 SOURCE-PACKAGE SCIENTIFIC FINAL VERDICT`: **BLOCK**

Opus independently confirmed the aggregate fail-closed behavior, component
identity coupling, implosion wind semantics, high-mass source exposure, source
fingerprints, and the intentional scientific blockers. It identified four
remaining hardening items:

- add an ownership mutation regression;
- scope the Z9.6 mass-budget bound and report high-mass diagnostics separately;
- reject an empty/truncated acquisition manifest with a coverage floor;
- add Limongi/NuGrid parser-mutation regressions.

## Disposition after the audit

All four items were implemented after the audit:

- `g2_sukhbold_channel_projection.py` now tests ownership mutation along with
  missing firewall/component and bad channel mutations;
- `g2_sukhbold2016_candidate.py` labels the existing mass-budget bound as
  Z9.6/9--12 M☉ scope and publishes a separate unbounded diagnostic for the
  six available high-mass yield tables;
- `audit_g2_candidate_sources.py` requires the eleven expected manifest
  candidate IDs and a minimum 65-file manifest, and propagates coverage failure
  to status/exit code;
- `g2_candidate_sources.py` tests empty/truncated manifests and malformed
  Limongi/NuGrid parser inputs;
- implosion radioactive entries retain `ejecta_msun: null`, and their
  wind-only/nonnegative claims are derived from parsed header/value data.

Post-disposition local evidence:

- `G2_CANDIDATE_SOURCE_AUDIT_TEST_OK`
- `G2_SUKHBOLD2016_CANDIDATE_TEST_OK`
- `G2_SUKHBOLD_CHANNEL_PROJECTION_TEST_OK`
- clean aggregate audit exit code `0`;
- full `run_g2_preflight.sh` exit code `0`, ending in the intentional
  `G2_PREFLIGHT_BLOCKED` state;
- current projection: 26 ordinary review records plus 19 high-mass review
  records, zero canonical rows, zero runtime deposition.

The Opus final audit itself predates these last four edits; AGY then reviewed
the post-disposition tree and returned an engineering **PASS**. The scientific
**BLOCK** is unchanged and correct: no source package has project physics
approval, and the 40--120 M☉ physical fate/yield/remnant/lifetime/energy/
momentum prescription remains absent.
