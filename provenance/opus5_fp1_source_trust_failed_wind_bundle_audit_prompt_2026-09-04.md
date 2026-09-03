# Claude Opus 5 bundle-end audit prompt — F-P1 source trust and LC18 cross-check

Conduct an independent, read-only implementation and physics audit of the
completed bundle in `/gpfs/kjhan/LRD_JWST` using Claude Opus 5 (`--model opus`).

Do not edit files, write generated artifacts, download data, contact authors,
launch RAMSES or other jobs, or assume that the dirty worktree is committed
source. Inspect files, hashes, git diffs, and narrowly scoped read-only tests
where useful. This is the bundle-end audit; there are no per-step audits.

## Final project purpose

The lagRamses project must become a production-ready and publication-ready
high-level hydrodynamics stack centered on radiative transfer, stellar/AGN
feedback, and dust. This bundle is specifically about stellar-feedback source
provenance, population/fate evidence, and an unresolved LC18 failed-wind
anomaly. Review whether the work advances that purpose without overstating a
review-only result as physical or runtime approval.

## Audit scope

Review the delta from base commit `296dd0c` and independently assess:

1. The code-owned source identity/rights lock profile and validator registry:
   exact candidate/release/DOI/Zenodo/license binding, exact nonempty inventory,
   byte and digest checks, composite fingerprint, scalar/date/path validation,
   rejection of symlinks and non-regular files, fail-closed exception behavior,
   and resistance to candidate substitution or a coherent rewrite of manifest,
   contract, and local bytes.
2. Fate sidecar and physical-package admission wiring: exact pinned
   repository-relative artifact paths, hash binding, registry-only validators,
   publication false while blocked, and whether any mutable or uncommitted
   evidence can open production.
3. LC18/ CDS cross-check science and implementation: parser reuse, one-to-one
   108-model join, 52 successful versus 56 failed rows, positive summary wind
   versus zero failed-release Wind tables, CDS endpoint comparison, duplicate
   phase-row policy, 3--8 unique phase nodes, monotonic age/mass checks, 96
   table7 records plus 12 explicit nulls, and no invented energy, structure,
   remnant, or wind composition values.
4. Whether the adversarial tests and JSON evidence establish the asserted
   properties rather than only reproducing a self-consistent fixture.

The expected result is fail-closed: physical node count 0 and false
`production_ready`, `publication_ready`, `canonical_conversion_allowed`, and
`runtime_deposition_allowed`. Treat this as intentional unless a bypass exists.
The four existing physical blockers must remain visible and unchanged.

## Required report

Provide a self-contained verdict of PASS, CONDITIONAL PASS, or FAIL. For every
finding give severity, exact file/symbol and line where possible, evidence,
reproduction route, and disposition. Check the numerical claims independently,
distinguish release-internal residuals from cross-source discrepancies, list
claims not reproducible from local staged bytes, and identify risks that belong
in the next implementation bundle. Conclude whether the driver may now draft
the next bundle plan.

Do not infer author confirmation, catalogue licensing, age-resolved winds,
per-node injected energies, canonical momentum/deposition, physical nodes, or
production/publication approval from a schema or synthetic test alone.
