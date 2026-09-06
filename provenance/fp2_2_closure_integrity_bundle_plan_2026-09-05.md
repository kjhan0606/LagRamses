# F-P2.2 closure-record integrity and explicit-path canonical coverage — 2026-09-05

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Parent: F-P2.1 source-closure verification; pre-approved by the driver after
the Claude Opus 5 conditional-pass audit on 2026-09-05.
Status: implementation/evidence complete; Claude Opus 5 end-of-bundle audit recorded as CONDITIONAL PASS; F-P2.3 approval pending.

## Objective

Close the remaining engineering weaknesses identified by the F-P2.1 Opus
audit without promoting an astrophysical source or widening into mixed
STAR+AGN admission:

1. bind the complete code dependency set used by the explicit AGN photon and
   source-weighted dust closures;
2. bind the serialized source-bound dust payload itself, not only its input
   files;
3. constrain dust status and propagate validated dust provenance into P4/P5
   outputs; and
4. give the repaired explicit nine-group path canonical artifact coverage and
   a truthful working-tree attestation.

## Work packages

### I1 — closure dependency manifests

- Record exact role/path/SHA-256 manifests for the AGN ledger and v2 dust
  closure, including `sed.py`, `primordial.py` where applicable, the loader,
  builder, and shared integrity helper.
- Re-hash every expected role at the runtime boundary and reject missing,
  substituted, extra, or stale entries.

### I2 — payload and status integrity

- Add a canonical JSON self-hash to source-bound photon/dust metadata.
- Reject modified closure arrays even when their input-file hashes remain
  unchanged.
- Restrict dust status to the versioned `reference_control` and
  `candidate_source_sed_matched` vocabulary with schema-consistent values.
- Record metadata, payload, source-table, and builder hashes in P4/P5 output
  attributes.

### I3 — explicit canonical nine-group control

- Add a deliberately synthetic, fully tabulated explicit SED control and
  regenerate its ledger, static input, and short transport artifact.
- Extend the canonical validator/artifact test with an explicit-source branch;
  retain the existing Sazonov pilot branch unchanged as a reference control.
- Record a current git-head plus clean/dirty working-tree status digest so a
  dirty canonical artifact is explicit rather than mislabeled.

## Exclusions

No mixed STAR+AGN aggregate dust admission, astrophysical SED/obscuration
selection, `[40,120] M_sun` yield decision, scattering/IR/grain physics, live
RAMSES coupling, production run, or publication claim is part of F-P2.2.

## Acceptance gates

- Source-bound AGN and v2 dust metadata carry exact dependency manifests and
  canonical payload hashes; all are revalidated before P4/P5 output creation.
- Array, manifest, status, and missing-file negatives fail closed.
- v1 null-identity reference controls continue to pass, while invalid status
  labels are rejected.
- Both pilot and explicit canonical nine-group validators/artifact checks pass,
  and the explicit artifact is marked synthetic candidate control rather than
  an approved physical SED.
- Focused RT/source/dust tests, Python compilation, and `git diff --check` pass.

After implementation, one bundled Claude Opus 5 read-only audit assessed
F-P2.2. The 2026-09-05 audit is recorded in
`provenance/claude_opus5_fp2_2_closure_integrity_bundle_end_audit_2026-09-05.md`
with verdict `CONDITIONAL PASS`. I1/I2 are closed; I3 remains conditional on
canonical asset synchronization, non-vacuous explicit artifact freshness
checks, and symmetric working-tree attestation. The recommended F-P2.3 bundle
is recorded there and awaits driver approval; no F-P2.3 implementation has
started.
