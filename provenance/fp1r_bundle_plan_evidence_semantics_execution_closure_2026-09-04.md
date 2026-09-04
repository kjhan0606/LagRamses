# Next bundle plan: F-P1R evidence, semantics, and execution closure

Date: 2026-09-04
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Status: Grok `APPROVE WITH CHANGES`; M1–M7 amendments applied; implementation
authorized by bundle-start plan governance
Parent evidence HEAD: `db1bb66`
Parent implementation: `25bd05f`
Parent verification boundary: `5aeb6d3`
Parent audit comparison: `fp1_identity_publication_closure_audit_comparison_2026-09-04.md`

## Authorization and governance

The user approved proceeding after the F-P1 identity/publication closure
bundle. This document is only the next driver plan. No F-P1R code or checked-in
evidence may be changed until Grok audits this plan and returns an acceptable
plan decision. After each F-P1R implementation step, Claude Opus 5 performs
the active implementation audit. AGY is retired and will not be called.

## Purpose and scope

The final purpose remains a production-ready and publication-ready lagRamses
high-level hydrodynamics stack focused on radiative transfer, stellar/AGN
feedback, and dust. F-P1R strengthens the evidence and semantics that protect
future stellar-yield admission. It does not select a source, create physical
source nodes, resolve the `[40,120] M_sun` fate seam, approve CDS
redistribution, or enable runtime feedback.

The real fail-closed state must remain unchanged: zero physical nodes,
unresolved `[0.8,1.0]` and `[40,120] M_sun` seams, zero canonical rows from
unapproved sources, false production/publication/conversion/deposition flags,
and the LC18 failed-wind anomaly unresolved.

## Findings driving this bundle

The F-P1 audits and independent triage agreed that the current implementation
has no live identity, mapping, rights, or runtime bypass. They identified
evidence/semantics/workflow gaps that should be closed before any future
positive physical-package admission:

- F1: the converter's admitted/generated mapping guard is sound, but the test
  suite never reaches a fully admitted positive converter path;
- F2: exact parsed CDS zero at `0.01 M_sun` table precision is not the same as
  an attested physical zero wind and the report needs outcome-specific
  qualification;
- F3: the population/fate runner validates the fate-admission report before
  regenerating the high-mass seam report consumed by that admission;
- F4: the publication-rights predicate accepts a caller-supplied terms digest
  and detached source-rights record instead of owning the terms-file read/hash
  boundary internally.

## Work packages

### R1 (P0) — exercise the admitted converter path in isolation

1. Extend `simulation/snrt/tests/yield_converter.py` with an explicitly
   synthetic, fully admitted package fixture. The fixture may patch only the
   converter-module seams (`DEFAULT_SOURCE_NODE_CONTRACT`,
   `DEFAULT_PHYSICAL_PACKAGE_CONTRACT`, and the bound
   `audit_source_node_contract` / `audit_physical_package_admission`). Use
   only temporary files and in-memory admitted selection/mapping seams; restore
   every mutable module seam in a `finally` block. Do not add a production skip
   flag, CLI contract override, or relax the converter's existing contract
   checks at `convert_yield_rows_to_canonical.py:379-380` and `:458-517`.
2. Reach the normal converter's positive path and write the table, sidecar,
   and source-node mapping only inside a temporary directory. Assert that the
   sidecar and mapping contain the exact generated asset/mapping hashes.
3. Build the admitted mapping from the non-writing proposal/generated
   document. Test a matching positive case, mapping-content mutation with the
   declared hash recomputed, content mutation without recomputing the declared
   hash, and hash-only mutation. Each disagreement must fail before any output
   path exists, demonstrating execution of the pre-write guard rather than
   merely inspecting it. Keep the existing blocked-path tests.
4. Assert the real repository's `simulation/snrt/config`,
   `simulation/snrt/data`, staged-source hashes, and blocked selection state
   are unchanged by the fixture. After the temporary fixture is restored, the
   real repository conversion must still fail closed. Do not weaken the
   production converter's requirement for the repository F-P1 contracts.

### R2 (P0) — qualify LC18 parsed-zero semantics

1. Replace ambiguous `cds_terminal_wind_zero_count` names in the review report
   with explicit parsed exact-zero terminology, or provide a versioned,
   backward-compatible alias whose interpretation is unambiguous. “Quantized”
   may describe the source print precision but must not imply that the pipeline
   rounded the winds. Apply the decision consistently to successful, failed,
   and all-model scopes and update all consumers/tests.
2. Emit the Table 5 precision (`0.01 M_sun`) and half-bin (`0.005 M_sun`) next
   to the counts, plus an explicit `physical_zero_inferred: false` field.
   Exact equality to zero must mean “parsed endpoint difference is zero at
   source precision,” never “physical wind is proven zero.”
3. Preserve and expose the outcome split: 4 successful parsed exact-zero
   endpoints with positive BR26 control wind, and 3 failed parsed exact-zero
   endpoints inside the unresolved BR26 zero-wind release anomaly. The 3
   failed CDS zeros must not be called the BR26 anomaly itself. Preserve
   52/56, 48/4, 53/3, and 101/7 accounting, including 52 successful models
   with nonzero BR26 Wind tables and 56 failed models with exact-zero
   unresolved Wind tables, plus all existing fail-closed blockers.
4. Add regression checks that reject accidental rewording which implies
   physical-zero inference or silently reconciles the failed-release anomaly.
   Do not rewrite historical provenance audit reports; update only live
   JSON/tools/tests and the unsent inquiry packet where needed.

### R3 (P0) — make the population/fate runner same-run fresh

1. Reorder `simulation/snrt/tests/run_fp1_population_fate_contract.sh` so the
   high-mass seam test/audit regenerates
   `data/fp1_high_mass_seam_review.json` before both physical-package and fate
   admission consume it. Alternatively, rerun both dependent audits after the
   final regeneration; the chosen implementation must cover both consumers.
2. In that same shell invocation capture the post-regeneration SHA-256 and
   assert that it equals `evidence_artifacts.high_mass_review.sha256` in both
   admission JSONs. A stale checked-in report must not make a green runner
   appear current. Do not add a per-run nonce: a byte-identical regeneration
   that matches the code-owned lock is the intended production state.
3. Keep the runner's expected terminal state and all production flags
   unchanged (`G2_PREFLIGHT_BLOCKED`, review-only high-mass seam, zero physical
   nodes). Do not broaden this into generic checkpoint or AMR work.

### R4 (P1) — make the publication terms boundary code-owned

1. Refactor `simulation/snrt/tools/fp1_publication_rights.py` so the
   production evaluator receives a candidate id and locked terms path/profile,
   reads the terms bytes itself, computes SHA-256 internally, parses
   `sources[candidate_id]` from those same hashed bytes, and compares path/hash
   to the code-owned lock. A caller-supplied digest or detached
   `source_record` must not be production authority. Missing keys, malformed
   JSON, wrong paths, and mutated bytes fail closed.
2. Keep isolated synthetic lock profiles available only as explicit in-memory
   test parameters. Put synthetic rights fields inside the hashed terms bytes.
   Test mutated bytes, wrong paths, malformed terms, missing candidate records,
   and review-only/approval failures; ensure all fail closed.
3. Update the LC18 cross-check to use the new API while retaining its current
   empty approval, review-only classification, and false publication gate.
   The terms file must not be edited to claim a license or redistribution
   permission that it does not contain; the multi-source catalog must remain
   the source of the parsed candidate record.
4. Ensure future export callers cannot consume a mutable JSON label without
   `require_publication_allowed()` and a fresh authoritative gate result. The
   bundle must not add a production export caller. `require_publication_allowed`
   must re-evaluate or refuse a forged gate dictionary; a hand-built
   `publication_ready` label is not permission.

## Acceptance gates

- **R1:** an isolated fully admitted converter success, content mutation, and
  hash mutation are all tested; mismatches fail before writes; all seams are
  restored in `finally`; real repository evidence remains byte-invariant.
- **R2:** report field names and interpretations distinguish parsed/quantized
  zeros from physical zeros; Table 5 precision/half-bin and
  `physical_zero_inferred: false` are explicit; 52/56, 48/4, 53/3, and 101/7
  remain exact; successful and failed release behavior is not conflated.
- **R3:** one runner invocation regenerates high-mass evidence before both
  consumers and proves both downstream reports bind the post-regeneration
  hash.
- **R4:** the code-owned terms reader/hash and bytes-derived candidate record
  are the sole production authority; altered terms/path/record/approval inputs
  cannot open publication.
- Full focused tests, population/fate contract, G2 preflight with expected
  `G2_PREFLIGHT_BLOCKED`, deterministic JSON regeneration, a second
  byte-identical high-mass regeneration, config/data and staged-source hash
  invariance (with the exception that the live LC18 JSON may change for R2
  names), compilation, and `git diff --check` pass.
- No source selection, physical-node creation, author inquiry, CDS
  redistribution, RAMSES run, or runtime feedback activation occurs.

## Explicit deferrals

The 40–120 M_sun physical fate decision, remaining validators, physical node
population, age-resolved source data, energy/momentum/deposition realization,
exact CDS/VizieR rights resolution, absolute-path portability, Python 3.9
compatibility, and generic restart/AMR infrastructure remain later bundles.

## Mandatory amendments from the Grok plan audit

The following M1–M7 edits are required before implementation and are now part
of this plan:

- **M1 — R1 isolation:** patch only the converter-module seams named in R1;
  use temporary outputs and in-memory admitted mapping/selection; restore all
  patches in `finally`; never add a production override or relax contract
  predicates; prove the repository path still blocks.
- **M2 — R1 mapping cases:** derive the admitted mapping from the non-writing
  proposal/generated document and cover matching, recomputed-hash mutation,
  unrecomputed mutation, and hash-only mutation, all before-write failures.
- **M3 — R2 semantics/split:** use parsed exact-zero names at Table 5 precision,
  publish the `0.005 M_sun` half-bin and `physical_zero_inferred: false`, retain
  the 52/56, 48/4, 53/3, and 101/7 split, and do not label the 3 failed CDS
  zeros as the BR26 anomaly; historical audit reports remain unchanged.
- **M4 — R3 freshness:** cover both fate and physical-package consumers and
  assert the same-invocation post-regeneration SHA-256 in both admission JSONs;
  do not use a per-run nonce; keep the blocked/review-only state.
- **M5 — R4 byte authority:** hash/read the locked terms bytes internally,
  parse the candidate record from those bytes, reject caller digest/record
  authority and malformed/missing/mutated inputs, and keep LC18 review-only.
- **M6 — acceptance hashes:** staged source bytes remain unchanged; the second
  high-mass regeneration is byte-identical; LC18 JSON may change only for R2
  naming/semantics; do not require all 248 config/data hashes to remain
  unchanged across R2.
- **M7 — lineage/order:** record `db1bb66` as parent evidence and `25bd05f` as
  implementation parent; execute R1 first, sequence R2 and R4 because both
  touch the LC18 tool, and integrate R3 last in the runner.

## Stop and audit rule

This plan requires Grok's bundle-start decision before implementation. After
Grok approval and each completed R1–R4 step, Claude Opus 5 audits that step.
If an Opus result is conditional or negative, run the GPT-5.6 Sol adjudication
path, incorporate confirmed findings into the driver plan, and obtain Grok's
review before starting the affected next bundle. After the whole F-P1R bundle
passes its step audits, record the final comparison and wait for a new user
approval before any subsequent bundle.
