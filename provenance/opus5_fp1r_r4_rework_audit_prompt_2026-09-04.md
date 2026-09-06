# Claude Opus 5 R4 rework implementation-stage audit

Act as the independent scientific, algorithmic, wiring, and implementation
auditor for the completed M-R4-1 rework in
`/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work read-only. Do not edit
files, run shell commands or jobs, select physical sources, contact authors,
or redistribute data. AGY is retired and must not be called or treated as an
active approval authority.

Audited implementation commit: `ac73d52` (`Close R4 wrong-path publication
attestation`). Previous Opus audit:
`provenance/opus5_fp1r_r4_publication_terms_audit_2026-09-04.md`.
GPT-5.6-Sol adjudication:
`provenance/gpt56sol_fp1r_r4_publication_terms_adjudication_2026-09-04.md`.
Grok plan approval:
`provenance/grok_fp1r_r4_rework_plan_audit_2026-09-04.md`.
Read the current bundle plan and all live implementation/evidence files:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/audit_governance_amendment_2026-09-04.md`
- `simulation/snrt/tools/fp1_publication_rights.py`
- `simulation/snrt/tests/fp1_publication_rights.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `simulation/snrt/config/g2_source_use_terms_evidence_v1.json`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. R4 remains pre-admission evidence hardening;
the unresolved source rights, failed-wind anomaly, fate seams, zero physical
nodes, zero canonical rows, disabled runtime deposition, and review-only
publication state are intentional and must remain unchanged.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for a live authority bypass, incorrect same-bytes semantics, a
publication opening, or a defect that invalidates R4. Use `CONDITIONAL PASS`
only for a material but non-blocking evidence gap. Do not treat optional
OSError coverage, the attestation-specific mutation message branch,
non-authoritative display labels, or absolute `/gpfs` paths as mandatory
unless you find a concrete safety issue.

## Required rework closure

The previous conditional/REWORK finding was:

- if `terms_path` mismatches the code-owned path, the evaluator opens no file
  but incorrectly reports `publication_terms_record_parsed=true` and
  `record_source=candidate_record_parsed_from_locked_terms_bytes`;
- the fix must initialize an explicit not-read/path-mismatch sentinel before
  the path check, let a successful locked-path read clear it (or replace it
  with a specific I/O/JSON/candidate error), and make the wrong-path test
  assert `publication_terms_record_parsed=false`,
  `record_source=not_available_due_to_terms_error`, and the named
  `publication_terms_record_unavailable` blocker.

Verify directly that:

1. the sentinel is initialized only for the skipped unlocked-path read and is
   overwritten by the locked-path reader result;
2. malformed JSON and missing candidates still report their specific errors,
   while a valid locked catalog still reports parsed=true and the positive
   record source;
3. wrong-path, mutated-byte, review-only, approval, forged-dictionary, and
   post-evaluation terms mutation cases remain fail-closed in the synthetic
   test, and no repository terms file is changed;
4. the LC18 caller still passes no digest or detached rights record, the live
   catalog and generated evidence remain review-only with false publication,
   and no export caller or physical/runtime activation was added;
5. the active governance contains Opus 5, Grok, and GPT-5.6-Sol only (apart
   from historical records), with AGY explicitly retired.

The driver reports that focused publication, LC18, physical-package,
converter, population/fate, and G2 preflight tests passed, with final
`G2_PREFLIGHT_BLOCKED`, and that the tracked LC18 JSON stayed byte-identical.
Assess those claims from the current files and provenance. Separate optional
maintenance observations from gate failures.

End with the verdict, checks actually established, mandatory fixes (if any),
non-blocking findings, and a direct statement on whether F-P1R R3 may begin.
Only an unconditional `PASS` opens R3. AGY is retired and must not be called.
