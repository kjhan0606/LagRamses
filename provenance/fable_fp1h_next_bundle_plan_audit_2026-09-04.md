# Fable pre-implementation audit: next F-P1H bundle

Date: 2026-09-04
Model: Claude `fable`
Mode: read-only plan and repository review

## Verdict

**APPROVE WITH CHANGES.** Part A, the source trust-root repair, is approved.
The proposed coordinate/population validator is rejected for this bundle and
must be replaced by a source-byte-derived failed-wind cross-check and an
unsent author-inquiry packet.

## Reason for changing Part B

Fable verified that the LC18 mass, metallicity, rotation, and branch inventory
already exists in the candidate-grid contract/audit; the Boccioli--Roberti
summary table is already byte-parsed; and the fate resolver already implements
branch-specific half-open cells, a closed final edge, out-of-hull blocking,
and no nearest-node fallback. With no approved production metallicity domain,
no rotation/binary population decision, and zero physical nodes, a new
coordinate validator would only restate an existing deterministic block.

The locally staged Limongi--Chieffi 2018 CDS `table5.dat` instead provides
lifetime and total mass at eight phase ends for all 108 LC18 models.
`M_initial - M_total(PSN)` supplies an independent cumulative wind-mass check,
including the 56 models marked failed by Boccioli--Roberti. `table7.dat`
provides Fe-core mass, compactness, binding energy, and envelope H/He for 96
models. Cross-joining these source bytes therefore advances lifetime, wind,
and pre-SN-structure evidence and makes the failed-wind inquiry row- and
checksum-specific without resolving it by assumption.

## Required Part A changes

1. Use `passed = all(requirements) and not blockers` in both the rights
   validator and registry.
2. Add a code-owned lock profile for exact candidate/release identity, DOI and
   record identity, license, the non-empty five-file inventory, integer byte
   counts, SHA256 values, composite fingerprint, Zenodo-record SHA256, and
   Zenodo-published file MD5 values.
3. Derive redistribution from the pinned license record and pinned attribution;
   retain the mutable terms file as reported evidence only.
4. Rename `immutable_local_source_mirror` to
   `hash_locked_local_source_mirror` throughout the gate contract and tests.
5. Reject all symlinks/non-regular files, empty/extra/duplicate inventories,
   wrong scalar types, null identity fields, and invalid calendar dates.
6. Convert runner and registry exceptions into controlled fail-closed errors.
7. Reject blocked fate sidecars with publication enabled and reject absolute
   or escaping artifact paths.
8. Exercise all AGY/Codex bypasses in temporary fixtures and prove the staged
   package fingerprint is unchanged by the tests.

## Approved Part B' deliverable

- Add one read-only cross-check tool, one focused test, one generated JSON,
  one provenance inquiry packet, and one F-P1 runner entry.
- Join exactly 108 Boccioli--Roberti LC18 summary rows to exactly 108 CDS rows
  by mass, metallicity label, and rotation; report zero unmatched rows.
- Use 52 successful models as release-internal controls and report their
  independent CDS differences without assuming agreement.
- For all 56 failed models, record zero Wind-table sum, positive summary wind,
  CDS cumulative wind, and `resolution: unresolved`.
- Require strictly increasing cumulative phase ages and non-increasing total
  mass across the 3--8 unique source phases; preserve exact duplicate evidence
  and list violations without correction.
- Join 96 table7 structure rows and represent the other 12 as explicit nulls.
- Emit zero canonical rows, retain all admission blockers and false approvals,
  make no external-tree changes or downloads, and do not send the inquiry.

## Conditions that remain blocked

Physical nodes remain empty; canonical conversion, runtime deposition,
production, and publication remain false; the Fortran source-node consumer
remains false; compiled approval identity remains empty; the production
metallicity domain and rotation/binary population remain unselected; no
coordinate validator is registered; and the Boccioli--Roberti hard blockers
remain unchanged.

## Purpose-fit judgment

Part A is necessary publication provenance. Part B' is direct feedback-physics
progress because it turns the dominant failed-wind source anomaly into a
reproducible model-by-model comparison while adding the first source-derived
lifetime, cumulative-wind, and pre-SN-structure evidence for the leading
multi-metallicity, multi-rotation, redistributable candidate.

## Independent correction before implementation

The implementation-side reproduction accepted Fable's 108-row join and 52/56
outcome counts but rejected two numerical assumptions in the audit response.
All 52 successful-model summary winds differ from CDS
`M_initial - M_total(PSN)` by more than the nominal 0.005 M_sun half-bin, with a
maximum difference of 0.5842 M_sun. Three failed models also have a CDS wind
rounded to zero rather than positive. In addition, exact duplicate phase rows
leave 3--8 unique phases per model, not eight for every model. The approved
physical cross-check remains valuable, but its acceptance criterion is the
faithful measurement and exposure of these discrepancies, not a predeclared
agreement. No source value will be corrected or inferred.
