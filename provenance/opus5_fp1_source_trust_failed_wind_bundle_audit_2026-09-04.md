# Claude Opus 5 bundle-end audit — F-P1 source trust and LC18 failed-wind cross-check

Target: `/gpfs/kjhan/LRD_JWST`, delta from base `296dd0c` (which is HEAD; the
entire bundle is an uncommitted worktree).
Date: 2026-09-04. Mode: read-only. No files edited, no jobs launched, no
downloads, no author contact.

## Verdict: CONDITIONAL PASS

The bundle is genuinely fail-closed. The expected end state is reproduced and
enforced: physical node count 0; `production_ready`, `publication_ready`,
`canonical_conversion_allowed`, `runtime_deposition_allowed` all false; the four
`boccioli_roberti2026_lc18` hard blockers present verbatim and unchanged. I found
**no bypass** — not by candidate substitution, not by a coherent rewrite of
manifest + contract + local bytes, not through the fate sidecar, and not through
the physical-package contract. The trust-root repair does answer the AGY and
Codex `gpt-5.6-sol` findings that failed the previous bundle. The LC18 work is
honest science: it reuses existing parsers, invents nothing, and refuses to
reconcile a real cross-source discrepancy.

The conditions are structural, not corrections of a wrong result. They are listed
under Dispositions.

## Verification limitation, stated up front

**This session had no shell.** `Bash` is not available to me or to my
subagents (`No such tool available: Bash`). I therefore could **not**:

- run any test, tool, or `run_fp1_population_fate_contract.sh`;
- recompute any SHA-256, so every pinned digest is read, not verified;
- run `git diff 296dd0c`, so I cannot attribute any line to this bundle versus
  the base.

Base identity was established by reading `.git/logs/HEAD` (line 16 ends at
`296dd0c…` — so 296dd0c is HEAD and the delta is exactly the dirty worktree).
Tracked-vs-untracked was established by scanning `.git/index`. Everything below
is verified by reading source, configs, generated artifacts and provenance
documents and checking them against each other. Where a claim rests on
documentary rather than cryptographic evidence, I say so.

Partial compensation: several LC18 numbers were re-derived **by hand from the raw
staged bytes** (see §3), which is stronger than reading the generated JSON. And
`simulation/snrt/{tools,tests}/__pycache__/*fp1_lc18_failed_wind_crosscheck*.pyc`
exist, so the new tool and test have in fact been executed under Python 3.13 —
weak evidence, but real, that the run reported in the bundle document happened.

---

## 1. Source identity / rights lock profile and validator registry

`simulation/snrt/tools/validate_fp1_source_identity_rights.py`

**Binding is exact and code-owned.** `LOCKED_CANDIDATE_PROFILES` (lines 38–115)
pins candidate id, source candidate id, release root, article/manifest/attribution
citations, article DOI `10.1051/0004-6361/202557714`, data DOI
`10.5281/zenodo.19503168`, Zenodo record id `19503168` (int, `type(...) is int` at
:385 and :495, so `True` is rejected), record filename and SHA-256, title,
ordered creators, `cc-by-4.0`, licence URL, `retrieved_date` `2026-09-01`,
composite `3370571245be954b1330d0b91bae585ffed47b3a1c2d10ffa11fc4ef7b57065b`, and
a 5-file inventory with exact byte counts (1,646,297 total) and per-file SHA-256 /
published MD5.

**Inventory is exact and non-empty.** `_manifest_inventory` (:257–279) rejects
empty/malformed lists, unsafe or duplicate paths, and requires
`set(by_path) == expected_paths`. `_contract_inventory` (:282–292) requires
`set(records) == set(profile["files"])`. Extra, missing, empty and duplicate all
block.

**Byte and digest checks.** :408–467 compares declared bytes/sha256/md5 in both
the manifest and the contract against the lock, requires `type(...) is int` for
byte counts, and then hashes the file on disk and compares observed bytes, SHA-256
and MD5 to the lock. `md5: None` for the Zenodo record is load-bearing twice: it
excludes that file from the published-file set (:533–535) and makes any declared
`md5` on it a failure (:433–434).

**Composite fingerprint.** `_candidate_fingerprint` (:225–236) is SHA-256 over a
NUL/newline-framed payload seeded with the source id, records sorted by UTF-8 path
bytes, built **only from observed bytes** (:465–471), then compared to the code
literal (:472–474). This is the anti-substitution anchor and it works.

**Symlinks and non-regular files.** `_confined_regular_file` (:166–202) rejects
absolute paths, empty paths and `..`, `lstat`s the root and every path component
rejecting `S_ISLNK`, requires the final component `S_ISREG` and intermediates
`S_ISDIR`, and re-checks containment with `resolve(strict=True).relative_to(...)`.
`os.path.islink` is not used; `lstat` + `stat.S_IS*` is the correct primitive.

**Scalars and dates.** Strict `type(x) is int` (:385, :412–416, :495, :543)
rejects bool and float. `_calendar_date` (:216–222) round-trips through
`date.fromisoformat`, rejecting `2026-9-1` and `2026-99-99`.

**Fail-closed.** The public entry point converts *every* exception into a blocked
report with all requirements false (:622–636). The only swallow is
`_non_authoritative_terms_report` (:295–321), which is deliberate: it is computed
at :569, **after** `passed` is decided at :568, lands only in
`artifacts["source_use_terms_report"]`, and carries
`"authoritative_for_verdict": false`.

**Candidate substitution and coherent rewrite: both blocked.** An unregistered
candidate returns `candidate_has_no_code_registered_rights_profile` (:334–338).
A rewrite of README bytes plus matching manifest and contract records fails at
`source_file_not_lock_pinned:README`, because the expected SHA-256 lives in code,
not in the rewritten JSON. Confirmed by test `_coherent_rewrite`
(`tests/fp1_source_identity_rights.py:154–169`).

**Registry.** `fp1_gate_validator_registry.py` registers via a module-level dict
literal populated by a static import (:10–16, :26–33). `run_registered_validator`
(:64–136) takes three keyword-only **strings** and no callable; the `isinstance`
guard precedes the dict lookup so an unhashable id raises the domain error. It
converts runner exceptions to `GateValidatorRegistryError` (:81–86) and then
**re-derives the verdict** rather than trusting it: exact 12-key set equality,
identity echo, requirement key-set + `type(v) is bool`, non-empty string blockers,
and `passed = all(requirements.values()) and not blockers` cross-checked against
the reported `passed`/`status` (:126–135). Blockers are verdict-bearing; a runner
cannot force a pass while a blocker is present.

## 2. Fate sidecar and physical-package admission wiring

Paths are **exact repository-relative literals**, whitelisted in code and matched
by string equality, then hash-checked against a freshly streamed digest:
`audit_fp1_fate_admission.py:41–49, 91–103` (7 artifacts) and
`audit_fp1_physical_package_admission.py:36–42, 132–143` (5 artifacts). Artifact
*sets* must be exact (`:219–220`, `:122–123`). No globs, directories, prefixes or
absolute paths. Both tools reject escapes via `relative_to`.

**Validators are registry-only.** Contract evidence may declare exactly
`{"validator_id"}` and nothing else (`:163–166`), the id must be in the
contract-approved set (`:167–171`), and the contract's `approved_validator_ids`
must equal `registered_validator_ids()` from code (`:232–245`). There is no
`eval`, `exec`, `importlib`, `getattr` or `subprocess` anywhere in the F-P1
admission path. Cross-gate reuse is blocked by the registry's
`record["gate_id"] != gate_id` check.

**Publication false while blocked is enforced.**
`audit_fp1_fate_admission.py:299–304` requires `publication_ready` to be boolean
and, if not `production_ready`, to be `False`. Tested at
`tests/fp1_fate_admission.py:80–82`.

**Can mutable or uncommitted evidence open production? No.**
`qualified = not missing and not blockers` (`:303–304`) where `missing` is
`REQUIRED_GATES - passed` and `passed` is populated **only** by executed registry
validators. Nine gates are required; one validator is registered. `missing` is
non-empty by construction, so deleting every `hard_blockers` entry from the
editable contract still cannot qualify a candidate. That is the load-bearing
fail-closed property, and it holds.

**Current values, as staged.**
`config/fp1_physical_package_admission_contract_v1.json:199` inventory `[]`;
`:200–207` selection all `null`; `:208–214` all five approval flags `false`.
`data/fp1_physical_package_admission_audit.json`: `physical_node_count: 0`,
`production_ready/publication_ready/canonical_conversion_allowed/runtime_deposition_allowed`
all `false`, `status: "blocked_no_qualified_physical_package"`.
`data/fp1_fate_admission_audit.json`: `production_ready: false`,
`status: "blocked_review_only"`.

**The four blockers, verbatim** (contract :168–173, byte-identical in both audit
artifacts):
```
failed_model_wind_summary_table_anomaly_requires_author_or_corrected_release
age_resolved_wind_missing
per_node_injected_energy_mapping_missing
canonical_momentum_and_deposition_missing
```
That they are *unchanged* rests on the bundle document, not on a diff I could run.

## 3. LC18 / CDS cross-check

`tools/audit_fp1_lc18_failed_wind_crosscheck.py`

**Parser reuse is genuine at the byte→record layer.** It imports
`_summary_rows`, `_yield_table` and `audit_boccioli_roberti2026_candidate` from
`audit_g2_boccioli_roberti2026_candidate`, and `adapt_candidate(LIMONGI_ID)` from
`adapt_g2_candidate_sources` (:16–26). I read `_summary_rows` (:91–129 of the
pre-existing tool): `exploded` comes from the literal `Yes`/`No` column,
`wind_mass_msun` from the `M_wind` column. No new low-level text parsing is
written. (One layer above that, the aggregation *is* re-implemented — see F-17.)

**No network.** `adapt_g2_candidate_sources.py` contains no `requests`, `urllib`,
`http`, `socket` or `ftp`. The "CDS endpoint" is the locally staged
`table5.dat` / `table7.dat`, not a CDS service call. The CDS/VizieR URLs in
`acquisition_manifest_v1.json` are recorded provenance and are never dereferenced.
Staged bytes are fingerprint-verified (size **and** SHA-256) against the manifest
before parsing, `adapt_g2_candidate_sources.py:167–171`.

**One-to-one 108-model join is enforced, not asserted.** Coordinate key is
`(rotation, [Fe/H], initial mass)`. Duplicates raise on the release side (:281–284)
and the table7 side (:225–228); the CDS grouping must yield exactly 108 (:113–116);
and both unmatched sets must be empty with `len(release) == 108` (:395–398).

**Counts, verified in `data/fp1_lc18_failed_wind_crosscheck.json`:**

| Claim | Staged value | Where |
|---|---|---|
| join one-to-one, 108/108/108 | `summary_row_count 108`, `cds_model_count 108`, `joined_model_count 108`, `one_to_one true` | json :38–42 |
| 52 successful | `model_count 52`, `summary_wind_positive_count 52`, `release_wind_table_nonzero_count 52` | json :11694–11697 |
| 56 failed, positive summary wind, exactly-zero Wind table | `model_count 56`, `summary_wind_positive_count 56`, `release_wind_table_exact_zero_count 56`, `unresolved_count 56` | json :28–33 |
| CDS endpoint wind for failed rows | positive 53, zero 3 | json :28–29 |
| all 108 exceed the half-bin | `above_nominal_half_bin_count 108`, `nominal_cds_total_mass_half_bin_msun 0.005`, `max abs 1.5476` | json :19–25 |
| release-internal control residual | `0.007183005956193367` | json :11693 |
| successful-row CDS residual max | `0.5842000000000014` (matches the plan's pre-implementation 0.5842) | json :11692 |
| 845 unique phase nodes, 19 duplicates collapsed, 3–8 per model | `unique_phase_row_count 845`, `collapsed_duplicate_row_count 19`, min 3, max 8 | json :48, :241–250 |
| no age / mass violations, no missing PSN | all three counts 0 | json :244–248 |
| 96 table7 records + 12 explicit nulls | `available_model_count 96`, `explicit_null_model_count 12` | json :254–256 |
| review-only outputs | `canonical_rows_emitted 0`, `physical_nodes_emitted 0`, `lc18_readme_consistency_pass false` | json :17, :252, :325 |

**Independent corroboration from raw staged bytes — not just the generated JSON:**

- The CDS `ReadMe:44–46` file index gives `table5.dat` **864** records,
  `table6.dat` 108, `table7.dat` **96**. So 864 = 845 unique + 19 collapsed
  duplicates, and 96 of 108 models have table7 ⇒ 12 explicit nulls. The counts
  are attested by the source's own index, independently of the tool.
- `table5.dat` is `130 bytes/record × 864 + 864 newlines = 113,184` bytes —
  exactly the manifest-pinned size. Consistent.
- Worked model `(rot 0, [Fe/H] −3, 80 M_sun)`: `table5.dat:705–712` durations
  `2.55e4, 3.33e6, 2.69e5, 4.72e1, 6.10e-2, 2.65e-2, 2.27e-4, 2.29e-4` sum to
  **3624547.287956**, and total mass runs `80.00 → 79.90` ⇒ CDS wind
  `0.0999999…`. The JSON reports `cds_terminal_age_yr 3624547.287956001` and
  `cds_terminal_cumulative_wind_mass_msun 0.09999999999999432`. Exact.
- Structure passthrough: `table7.dat:82` is
  `0 -3 80 24.20 25.17 1.81 36.37 0.52 0.000 0.000 0.000 0.000 0.000 SNIIP`;
  the JSON carries `hydrogen_envelope_mass_msun 24.2`,
  `helium_envelope_mass_msun 25.17`, `iron_core_mass_msun 1.81`,
  `binding_energy_above_iron_core_1e44_j 36.37`, `compactness_xi_2p5 0.52`,
  `expected_supernova_type "SNIIP"`. Byte-faithful, nothing synthesized.
- The 19 collapsed duplicates are all `PSN` rows; coordinate `(300, −3, 120)`
  repeats 6× as six byte-identical lines at `table5.dat:859–864`.
- The `cds_terminal_wind_zero_count: 3` is exact float equality on
  `initial_mass − table5 PSN total_mass == 0.0` (:496–497) — i.e. models whose
  printed mass loss rounds to `0.00`, not a substituted zero.
- The inquiry packet's failed-row table
  (`…inquiry_packet_2026-09-04.md:84–139`) has exactly 56 rows and exactly 12
  `null` Fe-core entries, consistent with the JSON.
- The unit semantics are the CDS ones: `ReadMe:201` defines `Ebind` as
  `10+44J … Binding energy … Of the mass above the iron core`, matching the
  adapter's field name and the report's
  `binding_energy_is_not_injected_explosion_energy: true` firewall.

The only physics-flavoured constant in the new tool is the half-bin, and it is
read from `g2_limongi_phase_mass_history_contract_v1.json:18`
(`phase_endpoint_total_mass_precision_msun: 0.01`), which is consistent with
table5's `F5.2` two-decimal total-mass format. The remaining numeric literals are
`108`, `52`, `96`, `12` (cardinalities), `1024*1024` (I/O block), `1.0e-12`
(mass epsilon) and `0.5` (half of the contract precision).

**Computed versus asserted.** Every count above is *computed* from staged bytes
and then checked against an expectation, not read from a fixture. Of the
expectations: `108` is enforced three times (:113, :395) and also declared in the
BR26 contract; **`56` is data-derived**, taken from the pre-existing audit's
`lc18_failed_models_with_reported_wind_but_zero_wind_table_count`
(`audit_g2_boccioli_roberti2026_candidate.py:277–280`), which is a genuinely
independent computation over the same bytes; `52`, `96` and `12` are literals at
:523 and :531. Note `52` is over-determined rather than load-bearing — since
`len(successful) + len(failed) == 108` is separately enforced and `len(failed)`
must equal the data-derived value, the literal `52` cannot mask a drift. `3` and
`8` appear only in the test, not in the tool.

**Duplicate policy and monotonicity.** Duplicates are collapsed on the upstream
adapter's `phase_occurrence > 1` (:93–112) and recorded, not deleted; a
non-collapsed duplicate raises (:134–137). Age (:165–173) and non-increasing mass
with 1e-12 tolerance (:174–185) record violations rather than modifying values.

**Nothing invented.** table7 records are passed through verbatim minus the
coordinate (:233–237); absent models become `None` via `structures.get(coordinate)`
(:434), i.e. explicit nulls. No energy, remnant, structure or wind-composition
value is synthesized. `binding_energy_is_not_injected_explosion_energy: true` and
the six `scientific_limitations` entries (:598–605) state the boundaries correctly.

**The tool actively re-checks the fail-closed state** (:329–375): the phase
contract must be fail-closed; both reviews must emit zero canonical rows; the four
admission blockers must match exactly; `physical_node_inventory` must be `[]`; and
all five approval flags must be `false`. Any drift raises and exits 2.

**README claim verified.** The packet's basis for the anomaly is accurate. The
staged `README` line 33 reads "`*_Wind.txt` contains the wind yields as reported in
the original papers/pre-SN files", and line 38 reads "For failed explosions, we
decided to set both the Post and PreSN yields to zero since the mass cut is not
defined" — Wind is not covered by that policy, which is exactly what makes 56
all-zero Wind tables anomalous. The BHSN exception list is at README lines 44–46.

## 4. Do the tests establish the properties, or just reproduce a fixture?

They establish them.

- `tests/fp1_source_identity_rights.py` runs 18 adversarial cases against real
  copies of the staged package (`_fixture`, :54–72 uses `copy2`/`copytree`, so no
  write reaches `external/`), each asserting a **named** blocker. It asserts the
  genuine package passes on the real default paths **before** (:314–319) and
  **after** (:356–359) the whole matrix, with an unchanged fingerprint and
  identical per-file reports. It also proves the terms file is non-authoritative
  by substituting it and still passing (:298–310).
- `tests/fp1_lc18_failed_wind_crosscheck.py` hashes `LC18.zip`, `table5.dat` and
  `table7.dat` before and after the run and asserts equality (:51–56), then
  asserts every count above. The values are recomputed from real staged bytes on
  each run, not read back from a fixture.
- `tests/fp1_physical_package_admission.py` mutates the real contract in 12 ways
  and asserts a named error each time.

With the coverage gaps in F-6 and F-4 below.

---

## Findings

### F-1 — HIGH (structural). The code-owned trust anchor is uncommitted.
`validate_fp1_source_identity_rights.py`, `fp1_gate_validator_registry.py` and
`tests/fp1_source_identity_rights.py` are **untracked** (`.git/index` scan; git
status shows `??`), and 296dd0c is HEAD. Meanwhile everything they authenticate
lives under `external/`, which `.gitignore:49–50` excludes. The comment at
`validate_fp1_source_identity_rights.py:35–37` — "This profile is the independent
trust anchor… Updating it is a reviewed code change" — is the design intent but is
not yet true: there is no committed baseline, so no committed artifact would
detect a joint rewrite of the lock profile and the mirrored bytes.
**Evidence:** `.gitignore:49–50`; `.git/index` has no `external/g2_candidates`
entries and no entries for the three new files.
**Reproduction:** grep `.git/index`; compare against `git status`.
**Disposition:** commit the bundle. This is the single highest-value action and it
converts the anchor claim from aspiration to fact. Not a defect in the code.

### F-2 — HIGH (forward risk). Fate-admission approval is not coupled to the physical-package verdict.
`audit_fp1_fate_admission.py:318` sets
`production_ready = sidecar_ready and fate_report["production_ready"]`.
`audit_physical_package_admission` is invoked at :246–250 purely for its
raise-or-return behaviour — `blocked_no_qualified_physical_package` is a *valid*
return — and its report is embedded at :334 with no assertion on its
`production_ready`, `publication_ready`, `unique_hard_blockers` or
`physical_node_count`. Same for `source_node_report` and
`terminal_deposition_report`.
**Failure scenario:** once the fate map's two unresolved intervals are closed,
the sidecar can set `production_ready: true, publication_ready: true` while
F-P1H-E still reports 14 hard blockers and `physical_node_count: 0`. The only
thing preventing it today is the unconditional epoch pin at :292–295 — which is
precisely the line someone must delete in order to promote, and deleting it
removes the coupling with it.
**Currently inert:** `expected_unresolved` has 2 entries, forcing `sidecar_ready`
false.
**Disposition:** add an explicit assertion that
`physical_package_report["production_ready"]` and
`["publication_ready"]` are true before the fate sidecar may claim either.

### F-3 — MEDIUM. Evidence digests in the physical-package contract are self-certified.
`fp1_physical_package_admission_contract_v1.json:15–36` declares the SHA-256 of
each of the five evidence artifacts **inside the same editable file being
audited**. Only the *paths* are code-owned (`EXPECTED_EVIDENCE_PATHS`, :36–42).
Editing `data/fp1_high_mass_seam_review.json` and updating the contract's declared
digest passes. This is the same circular-evidence pattern the previous audit
rejected for the source package — and the source-rights gate is the one place it
was fixed.
**Currently inert:** `high_mass_ready` / `candidate_grid_ready` are only consumed
in the selection branch, which is unreachable (F-4).
**Disposition:** before any promotion, move the evidence digests into code
alongside `LOCKED_CANDIDATE_PROFILES`, or state explicitly in the contract that
they are self-certified and not a trust anchor.

### F-4 — MEDIUM. The production-opening branch is unreachable and untested.
`audit_fp1_physical_package_admission.py:212–213` raises unless
`required_birth_metallicity_domain_selected is False`, while :368 requires it to
be `True` for a selection. Both cannot hold, so `publication_ready: true` is
unreachable without a code edit. That is strongly fail-closed — but it means every
guard on that path is dead code with no coverage: the qualification check
(:343–344), non-empty inventory (:345–346), package/mapping SHA validity
(:347–353), source-node fingerprint agreement (:354–361) and the upstream-gate
conjunction (:362–371). `tests/fp1_physical_package_admission.py:172–174` only
exercises the *partial*-selection error. When the epoch pin at :212 is lifted, an
entirely untested code path becomes the production gate.
**Disposition:** add a test that builds a fully populated `selection` in a
temporary contract and asserts each of those five guards fires by name.

### F-5 — MEDIUM. `unique_hard_blockers` is reported but is not a gate.
`:396–402` unions `hard_blockers` across **all** candidates and emits it at :423.
Nothing checks it. A future admitted state could legitimately report
`publication_ready: true` alongside a non-empty `unique_hard_blockers` contributed
by unselected candidates. Cosmetic today, misleading later.
**Disposition:** rename to `unique_hard_blockers_all_candidates`, or scope it to
the selected candidate.

### F-6 — MEDIUM (test coverage). Gaps in the adversarial matrix.
`tests/fp1_source_identity_rights.py:330–349` does not cover:
`candidate_root_symlink_forbidden` (root-level symlink);
`candidate_source_contract_symlink_forbidden` (only the manifest symlink is
tested, :242–246); `acquisition_manifest_candidate_identity_not_unique`
(duplicate candidate records, validator :357–360); and a runner returning
`passed=True` with blockers present (registry :128–130).
Additionally, `_substitute_record_and_license` (:232–239) mutates the Zenodo
record file, which changes its SHA-256 and therefore trips the per-file byte lock
first — so `machine_readable_license_not_verified` (validator :517–519) and
`zenodo_published_file_identity_not_verified` (:549–550) are **never reached in a
failing state**. Those two blockers are untested defense-in-depth.
**Disposition:** add the four missing cases; exercise the licence/published-file
blockers by pointing the validator at a lock-consistent fixture whose record file
differs only in a field the byte lock cannot see (or accept them as unreachable
and document that).

### F-7 — LOW. `_regular_file_without_symlink` checks only the final component.
Validator :153–163 `lstat`s the target but does not walk its parents, so a
symlinked ancestor directory for the manifest, source contract or terms file is
not detected — asymmetric with `_confined_regular_file`'s full component walk
(:183–192). Bounded in practice by the field-level lock, which authenticates
content regardless of how the path was reached.

### F-8 — LOW. Two residual fail-crash (rather than fail-block) edges.
(a) `_blocked_report` calls `_sha256(TOOL_PATH)` at validator :253 from *inside*
the catch-all handler at :630–636; an unreadable tool file escapes uncaught.
(b) `fp1_gate_validator_registry._sha256` (:36–41) has no `OSError` handling and
is called at :133 outside the `try` at :81; an unreadable tool file escapes as a
raw `OSError`, which `audit_fp1_physical_package_admission.py:178` does not catch
(it catches only `GateValidatorRegistryError`). Neither yields a pass.

### F-9 — LOW. `hashlib.new("md5")` without `usedforsecurity=False`.
Validator :129 / :439. The test uses the flag (`tests/…:160`), the tool does not.
On a FIPS-restricted interpreter the gate blocks for an environmental reason
reported as `validator_exception:ValueError`. Fail-closed, but the blocker name
would misdescribe the cause.

### F-10 — LOW. Uncontrolled `KeyError` in the physical-package report assembly.
`audit_fp1_physical_package_admission.py:425–428` indexes
`high_mass["source_node_completeness"][...]` and
`high_mass["cross_engine_wind_review"][...]` directly, unlike the guarded `.get`
chains at :274–281. A missing key raises `KeyError`, which `main()`'s
`except PhysicalPackageAdmissionError` (:447) does not catch. Non-zero exit either
way, but Part A item 2 asks for *controlled* results.

### F-11 — LOW. Hex-shape validation asymmetry.
`audit_fp1_fate_admission.py:97–98` validates hexadecimality;
`audit_fp1_physical_package_admission.py:136` checks only `len == 64`. Harmless
(the digest comparison at :140 fails anyway) but the two tools should match.

### F-12 — LOW. Confinement root is `simulation/snrt`, not the repository.
Both admission tools confine to `SNRT_ROOT`. Separately,
`audit_fp1_fate_admission.py:37–40` reads
`patch/lagRamses/stellar_enrichment_config.f90` — outside that root — by regex,
with no digest pin. Semantic content *is* checked (interval bounds to 1e-12,
`compiled_fate_map_sha256` and `compiled_fate_approval_id` must be `""`,
`snii_source_node_fate_consumer_available` must be `.false.`), and both files are
git-tracked, so this is identity-unpinned but content-constrained.

### F-13 — LOW (science, precision of claim). The age-monotonicity check is tautological.
`cumulative_age` is built by accumulating `phase_duration_yr`
(`audit_fp1_lc18_failed_wind_crosscheck.py:155`), so the "strictly increasing
cumulative age" test at :165–170 reduces to "every phase duration is finite and
> 0". That is a real check, but it is not an independent monotonicity result and
should be stated as duration positivity. The **mass** check (:174–185) is a
genuine monotonicity test against the previous node.

### F-14 — INFORMATIONAL (science, resolved on inspection). Duplicate-row identity is delegated but sound.
The tool collapses rows on the upstream adapter's `phase_occurrence` and requires
`properties["all_duplicate_rows_physically_identical"] is True` (:383–386) before
collapsing, with a belt-and-braces post-check that no duplicate phase survived
(:133–137). The delegated flag is not a bare assertion: the adapter computes it
from an exact numeric-tuple signature set
(`adapt_g2_candidate_sources.py:301, 340, 353–359`, `len(signatures) == 1`), and
its stated policy is "any duplicate with different physical values must fail
closed". I confirmed the 19 collapsed rows are genuinely byte-identical in
`table5.dat`. No action needed.

### F-17 — MEDIUM. The phase-history aggregation is a second implementation, not reuse.
Part B' item 1 of the approved plan asks to "reuse existing source parsers rather
than introducing a second interpretation of the table formats". The byte→record
parsers *are* reused (§3). But `_build_phase_histories`
(`audit_fp1_lc18_failed_wind_crosscheck.py:87–216`) re-implements, nearly
line-for-line, the aggregation already present in
`audit_g2_limongi_phase_mass_history.py:93–152`: same 108-model grouping, same
`occurrence > 1` collapse, same `cumulative_age += duration`, same
`initial_mass − total_mass` cumulative wind, same `phase_rank` ordering, same
`0.5 × phase_endpoint_total_mass_precision_msun` half-bin (`:191–193` there
vs `:457–459` here). That module is not imported.
**Mitigating:** the two implementations agree exactly — the pre-existing
`data/g2_limongi_phase_mass_history_audit.json` reports
`collapsed_extra_row_count 19`, `phase_row_count_after_exact_collapse 845`,
min 3 / max 8, `model_count 108`, identical to the new artifact. So this is
redundancy, and the agreement is itself corroborating evidence — not divergence.
**Disposition:** import the existing aggregation, or state in the artifact that
the duplication is a deliberate independent reproduction. Do not leave two
copies to drift.

### F-18 — MEDIUM. The cross-check tool is less fail-closed than its predecessor on monotonicity.
Age violations, mass violations and missing-PSN are **collected and reported
only** (`:171–185`, surfaced in `diagnostics` at `:209–214`); the tool never
raises on them. The older `audit_g2_limongi_phase_mass_history.py:153–154` does
fail closed on invalid phase-derived wind mass. Today the only thing asserting
zero is the test (`tests/fp1_lc18_failed_wind_crosscheck.py:125–127`), and the
runner does execute that test — so the *contract* is fail-closed even though the
*tool* is not.
**Failure scenario:** anyone invoking `audit_fp1_lc18_failed_wind_crosscheck.py`
directly (it has a `main()` and a `--json-out`) on drifted table5 bytes gets exit
0 and a JSON artifact reporting nonzero violation counts, which reads as a
successful audit.
**Disposition:** raise on nonzero age/mass/PSN violation counts, matching the
predecessor. The approved plan's wording ("List duplicates and violations rather
than modifying or filling values") justifies *listing* them, not exit 0.

### F-19 — LOW (reporting hygiene). Two payload booleans are literals, not computed.
`"one_to_one": True` (:568) and `"hard_blockers_unchanged": True` (:590) are
hardcoded in the returned dict. Both are *justified* — they are unreachable
unless the raises at :395–398 and :356–375 passed — but as written the JSON field
is not evidence of the property, it is a restatement of an assumption a reader
cannot check from the artifact alone. Compute them (or emit the comparison inputs)
so a downstream reader can verify without reading the tool.

### F-20 — LOW. `LC18.zip` is opened twice, the second time without CRC validation.
`audit_boccioli_roberti2026_candidate` (:317) opens the archive and validates it
through `_archive_identity` (CRC + identity). `_release_rows` then re-opens the
same path at :246 with a bare `ZipFile(...)`, bypassing that validation. Same file,
same process, so the practical risk is negligible; but the second read is not
identity-checked. Pass the already-validated handle, or re-run `_archive_identity`.

### F-15 — LOW (science). CDS rights status is read from a mutable file and not labelled as such.
`limongi_chieffi_cds_redistribution_status` and
`…_production_license_status` come from the editable, non-hash-pinned
`config/g2_source_use_terms_evidence_v1.json` via the adapter (:536–538, :554–559).
The verdict does not depend on it (`review_use_only: True` is a literal, zero
canonical rows are emitted, and the status correctly reads
`no_explicit_catalogue_license_identified`), but the source-rights validator sets
a good precedent by tagging such reads `authoritative_for_verdict: false`. The
cross-check artifact should do the same.

### F-16 — INFORMATIONAL (next bundle). The cross-source discrepancy is reported but not characterised.
The tool computes signed differences per row (:426–432) but summarises only
`max |diff|`. Two facts visible in the inquiry-packet table
(`…inquiry_packet_2026-09-04.md:84–139`) are not in the evidence artifact:
every one of the 56 failed-row differences is **positive** (summary wind > CDS
wind), and for large winds the residual is ≈0.3–0.8 % of the wind mass
(e.g. 79.9214 vs 79.50 at B/300/120). A sign-consistent, roughly multiplicative
residual is strong evidence of a **definitional endpoint difference** rather than
random disagreement — which is exactly what inquiry question 3 is asking. Adding
signed and relative statistics would sharpen the question materially and costs
almost nothing.

---

## Release-internal residuals versus cross-source discrepancies

The bundle keeps these correctly separated, and the audit confirms the separation
is real:

- **Release-internal** (BR26 summary `M_wind` vs the sum of its own `*_Wind.txt`
  elements): 52 successful models, max residual **0.0072 M_sun**. This is a
  round-trip consistency check inside one release and is small, as expected.
- **Cross-source** (BR26 summary `M_wind` vs LC18 CDS `M_initial − M_total(PSN)`):
  **all 108** models exceed the 0.005 M_sun nominal table5 half-bin, max
  **1.5476 M_sun**, including all 52 successful controls. This is a genuine
  disagreement between two published sources, not a tolerance failure.
- **The anomaly itself** (56 failed models with positive summary wind and an
  exactly zero `*_Wind.txt`) is a third, separate thing, and it is the one that
  requires the authors.

The tool sets `agreement_required_for_this_review: false` (:511) and
`cross_source_difference_silently_reconciled: false` (:577). Correct.

## Claims not reproducible from local staged bytes

- **Author confirmation / corrected release** — none exists; correctly not
  claimed, and the packet is marked "not sent".
- **CDS catalogue licence** — correctly declared
  `no_explicit_catalogue_license_identified`; the derived artifact is review-only.
- **Zenodo CC-BY-4.0, published MD5s, title, creators** — reproducible only from
  the *local mirrored copy* of `zenodo_record_19503168.json`. It is hash-locked in
  code, which is the right control, but it is a local mirror of a remote
  assertion: it cannot detect that the upstream Zenodo record has since changed.
  Honest as stated; worth naming as a residual assumption.
- **Age-resolved winds, per-node injected energies, canonical momentum /
  deposition, physical nodes** — none inferred anywhere. Blockers 2–4 remain.
- **The `C, Ne, O, Si` phase ordering** used to sort phase rows comes from
  `g2_limongi_phase_mass_history_contract_v1.json:8`
  (`["MS","H","He","C","Ne","O","Si","PSN"]`). The CDS `ReadMe:164–166` defines
  only `MS`, `H`, `He` and `PSN` explicitly. The intermediate burning-stage order
  is physically standard and almost certainly right, but it is a **project
  assertion**, not a source-attested fact — and the age-monotonicity result
  depends on it, since a wrong order would reorder the durations. Worth stating
  in the artifact.
- **The `A/B/C/D → [Fe/H] 0/−1/−2/−3` join key mapping** lives in
  `g2_boccioli_roberti2026_candidate_contract_v1.json:57`. It *is* backed by the
  staged `README:42` ("A, B, C, D correspond to solar, 1/10 solar, 1/100 solar
  and 1/1000 solar"), so this one is source-attested.
- **That the tests pass** — not verifiable in this session (no shell).
- **Diff attribution against 296dd0c** — not verifiable in this session.

## Risks for the next implementation bundle

1. `audit_fp1_source_node_contract.py:112` — `APPROVED_RIGHTS_STATUSES =
   {"approved", "verified", "permitted"}`, accepted from a JSON string per node
   (:755–759). Inert only because `physical_node_inventory` is `[]`. The moment a
   physical node is added, node rights become self-certified JSON. This should be
   bound to the executed rights validator, not to a status string.
2. The epoch pins that currently do the work — `:212–213` in the physical-package
   tool and `:292–297` in the fate tool — will be deleted to promote. Every
   deletion must be paired with the coupling from F-2 and the tests from F-4.
3. The deferred `fp1.coordinate_hull_and_population.v1` validator remains correctly
   unregistered; it needs an approved metallicity domain and a rotation/binary
   decision first, per the approved plan.
4. The remaining 8 of 9 required gates have no registered validator. That is what
   keeps the system closed today; it is also the whole remaining workload.
5. Two implementations of the LC18 phase aggregation now exist (F-17) and agree.
   Fold them together before either is touched, or the next edit silently forks
   the 845/19/3/8 numbers.
6. `audit_fp1_lc18_failed_wind_crosscheck.py` exits 0 on nonzero monotonicity
   violations (F-18). Fix before the artifact is cited in a paper draft, since a
   reader will reasonably take exit 0 as "no violations".

## Ordered actions

1. **Commit the bundle** (F-1). Everything else labelled "code-owned" is
   conditional on this.
2. Couple fate-admission approval to the physical-package verdict (F-2).
3. Make `audit_fp1_lc18_failed_wind_crosscheck.py` raise on monotonicity
   violations (F-18).
4. Close the adversarial-matrix gaps (F-6) and add the selection-branch test
   (F-4).
5. Move the physical-package evidence digests into code, or label them
   self-certified (F-3).
6. Add signed/relative cross-source residual statistics (F-16) — cheap, and it
   materially sharpens inquiry question 3.
7. Deduplicate the phase aggregation (F-17); tidy F-5, F-7 … F-15, F-19, F-20 as
   routine hygiene.

## Answer to the closing question

**Yes — the driver may draft the next bundle plan.** The bundle does what it
claimed, does not overstate a review-only result as physical or runtime approval,
and leaves the system verifiably closed. The conditions above are inputs to that
plan, not blockers on writing it. F-1 (commit the bundle) should be done first,
since every "code-owned" claim in the plan depends on it.
