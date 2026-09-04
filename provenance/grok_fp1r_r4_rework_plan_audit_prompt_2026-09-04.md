# Grok plan audit: F-P1R R4 evidence rework

Audit `/gpfs/kjhan/LRD_JWST` read-only. Do not edit files, run jobs, select
physical sources, contact authors, or redistribute data. This is a plan audit
for the narrow completion rework inside the already approved F-P1R bundle.
AGY is retired and must not be called or treated as an auditor.

Read:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/opus5_fp1r_r4_publication_terms_audit_2026-09-04.md`
- `provenance/gpt56sol_fp1r_r4_publication_terms_adjudication_2026-09-04.md`
- `simulation/snrt/tools/fp1_publication_rights.py`
- `simulation/snrt/tests/fp1_publication_rights.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. The proposed rework only corrects evidence
semantics at the publication-rights boundary. It must not select a physical
yield source, create physical nodes, enable canonical conversion/runtime
deposition, edit the repository terms catalog, or broaden into RAMSES,
checkpoint, AMR, or generic infrastructure work.

## Proposed rework

Opus 5 returned `CONDITIONAL PASS`, and GPT-5.6-Sol independently returned
`REWORK R4`, for one confirmed defect: when the caller passes a path that does
not match the code-owned lock, no terms bytes are opened, but the current
payload reports `publication_terms_record_parsed=true` and
`record_source=candidate_record_parsed_from_locked_terms_bytes`. Publication
is still blocked, but the machine-readable rights attestation is false.

The driver proposes to:

1. initialize an explicit not-read/path-mismatch `terms_error` sentinel before
   the path check; clear it only after the locked-path reader successfully
   supplies a parsed candidate record;
2. make the wrong-path result report a false parsed-record requirement, an
   unavailable record source, a named blocker, and denied publication;
3. extend the temporary synthetic wrong-path regression to assert those two
   evidence fields;
4. rerun focused publication/LC18/physical-package/converter/population-fate
   tests and deterministic G2 preflight, then obtain a clean Claude Opus 5
   implementation audit before R3.

Evaluate whether this is the minimum scientifically and technically justified
repair, whether the acceptance tests prove the intended semantics, whether
the terms bytes remain the sole rights authority, and whether any missing
condition would make the plan unsafe or infeasible. Treat the following as
non-blocking unless you find a reason otherwise: extra OSError test coverage,
the attestation-mutation branch coverage, detached non-authoritative LC18
display labels, and absolute `/gpfs` paths in evidence.

Return exactly one plan decision: `APPROVE`, `APPROVE WITH CHANGES`, or
`REJECT`. Include concise rationale, mandatory changes if any, scope/feasibility
assessment, and whether the driver may implement the rework. The active audit
chain is Claude Opus 5 for implementation, Grok for bundle-plan review, and
GPT-5.6-Sol for conditional/negative Opus adjudication; AGY has no role.
