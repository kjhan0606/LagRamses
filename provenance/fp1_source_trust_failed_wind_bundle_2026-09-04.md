# F-P1 source-trust and failed-wind cross-check bundle

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Status: implementation complete; bundle-end AGY/Claude Opus 5 audit pending

## Scope and purpose

This bundle implements the approved F-P1 follow-up for the production-ready,
publication-ready lagRamses high-level hydrodynamics stack. It is limited to
stellar-feedback source provenance, population/fate evidence, and the
Boccioli--Roberti LC18 failed-wind anomaly. It does not reopen completed
RAMSES topology, header, HDF5, restart, or generic distributed-runtime work.

This is an evidence and admission-boundary bundle, not a physical-source
approval. No canonical physical row, runtime deposition path, production
package, or publication package was enabled.

## Part A: trust-root repair

`fp1.source_identity_and_rights.v1` now uses a code-owned lock profile for the
candidate identity, release root, article/data DOI, Zenodo record and file,
CC-BY-4.0 metadata, exact five-file inventory, byte counts, per-file hashes,
published MD5 values, version-record hash, and composite fingerprint
`3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`.

The validator and registry now fail closed on malformed input, runner
exceptions, candidate/root substitutions, empty or extra inventories,
duplicate entries, non-regular files, symlinks, unsafe paths, wrong scalar
types, null identities, invalid dates, mutable terms-file changes, and
self-consistent manifest/contract/byte rewrites. The fate sidecar is pinned to
the exact repository-relative artifact set and cannot claim publication while
blocked. The genuine staged package remains unchanged by the adversarial
fixtures and passes the five source-rights requirements without blockers.

## Part B': LC18 failed-wind cross-check

The new review-only audit reads only manifest-pinned local bytes and reuses the
existing Boccioli--Roberti and Limongi--Chieffi parsers. Its evidence records:

- an exact one-to-one join of 108/108 summary and CDS models;
- 52 successful release rows with positive summary wind and nonzero release
  `Wind` tables;
- all 56 failed rows with positive summary wind but exactly zero release
  `Wind` tables, retained as unresolved;
- CDS endpoint wind positive for 53 failed rows and zero/quantized for 3;
- summary-to-CDS differences above the nominal 0.005 M_sun table5 half-bin
  for all 108 rows, with maximum absolute difference 1.5476 M_sun;
- summary-to-release-Wind residual maximum 0.007183005956193367 M_sun for
  the 52-row internal control;
- 845 unique phase nodes after 19 exact duplicate rows are collapsed, with
  3--8 unique phase points per model and no age or non-increasing-mass
  violations;
- 96 table7 structure records and 12 explicit null records, without inferring
  missing structure, binding-energy, or explosion-energy values.

The result deliberately reports `lc18_readme_consistency_pass=false`,
`failed_wind_anomaly_resolved=false`, and
`cross_source_difference_silently_reconciled=false`. The inquiry packet is
generated but has not been sent. CDS catalogue redistribution remains
review-only because no explicit catalogue licence was identified.

## Verification

The following checks passed after the final implementation patch:

- focused LC18 cross-check test;
- Python compilation of all changed tools and tests;
- `run_fp1_population_fate_contract.sh`;
- `run_g2_preflight.sh`, including all subordinate tests, terminating at the
  expected `G2_PREFLIGHT_BLOCKED` fail-closed state;
- JSON syntax validation and `git diff --check`.

## Current admission state

The candidate `boccioli_roberti2026_lc18` remains unqualified. The unchanged
hard blockers are:

1. failed-model Wind summary/table anomaly requires author confirmation or a
   corrected release;
2. age-resolved wind is missing;
3. per-node injected-energy mapping is missing;
4. canonical momentum and deposition are missing.

Physical node count is zero. `production_ready`, `publication_ready`,
`canonical_conversion_allowed`, and `runtime_deposition_allowed` are all
false.

## Next audit and planning boundary

This bundle is now complete at the implementation/evidence level. At bundle
end it is submitted independently to AGY (`gemini-3.8-flash-high`) and Claude
Opus 5. Their findings will be independently reproduced and triaged. The
driver will then write the next implementation-bundle plan, including any
accepted findings. Fable will evaluate that next plan for final-purpose fit,
scientific/technical justification, and feasibility before implementation
starts. No per-step audit is scheduled.
