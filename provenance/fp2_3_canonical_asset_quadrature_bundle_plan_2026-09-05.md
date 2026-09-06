# F-P2.3 canonical asset synchronization and source-closure quadrature control — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Repository: `kjhan0606/LagRamses`
Parent: F-P2.2 closure-record integrity and explicit canonical coverage
Approval: driver-approved on 2026-09-05 after the F-P2.2 Opus 5 conditional
pass audit

Status: implementation/evidence complete; the first bundled Claude Opus 5
audit returned CONDITIONAL PASS, its local findings were remediated, and the
follow-up disposition returned PASS. No runtime RAMSES coupling or physical
SED selection is included.

## Objective

Close the F-P2.2 audit findings that are local to canonical evidence and
source-closure numerics:

1. make every declared canonical asset checksum actionable;
2. make pilot and explicit artifact tests enforce the same meaningful
   provenance/freshness rules;
3. make the explicit-mode expectations independent of the ledger under test;
4. declare and bind the SED interpolation/quadrature convention, with a
   refinement/convergence guard; and
5. label the explicit SED as a synthetic non-physical wiring fixture and
   validate its sidecar contracts.

## Work packages

### Q1 — canonical asset synchronization

- Complete. Every one of the ten entries declared by
  `agn_nine_group_external_assets.json` is checked for exact ID, path, size,
  and current SHA-256; the stale explicit metadata digest was refreshed; and
  explicit rebind/transport reproduction commands were added.

### Q2 — artifact and attestation integrity

- Complete. Pilot and explicit artifact tests perform the same live-file hash
  sweep and use a repository-root-relative `simulation/snrt` attestation;
  the validator retains the attestation as provenance rather than treating its
  own recorded value as an independent pass criterion.

### Q3 — independent source closure

- Complete. Explicit source closure uses piecewise-linear interpolation of the
  declared energy-fraction density `f_E`, derives `q_E=f_E/(E*EV_ERG)`, and
  compares 2048- and 4096-subdivision logarithmic union grids at a declared
  `5e-6` relative tolerance. The convention is identity-bound and is shared
  by explicit photon, Verner, and source-weighted dust closures.

### Q4 — explicit metadata contracts

- Complete. Explicit support and rates are reconstructed from the SED; dust
  duplicate energy-fraction ledgers and source status enums are checked; and
  the canonical explicit ledger is marked as a synthetic, non-physical
  wiring fixture.

## Acceptance gates

- all declared asset entries have current path, size, and SHA-256 matches;
- pilot and explicit artifact tests both verify current provenance hashes and
  the scoped working-tree attestation;
- explicit expected group support and rates are independently reconstructed
  from the source SED;
- base/refined source quadrature is converged at the declared tolerance for
  photon moments, Verner closure, and source-weighted dust closure;
- negative status, duplicate-fraction, and payload/manifest tests fail closed;
- canonical pilot and explicit validators pass; and
- focused source/dust/RT tests, `py_compile`, and `git diff --check` pass.

## Exclusions

No physical AGN/stellar SED adoption, escape/obscuration calibration, mixed
STAR+AGN dust admission, `[40,120] M_sun` yield resolution, scattering/IR or
grain-temperature physics, live RAMSES coupling, production convergence, or
publication claim is part of this bundle.

The first disposition, remediation, and PASS follow-up are recorded in the
implementation evidence and
`provenance/claude_opus5_fp2_3_canonical_asset_quadrature_bundle_followup_audit_2026-09-05.md`.
No F-P2.3 acceptance gate remains blocking; the listed non-blocking items are
carried into later physical-source/publication work.
