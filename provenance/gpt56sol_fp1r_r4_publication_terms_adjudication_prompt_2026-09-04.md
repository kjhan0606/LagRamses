# GPT-5.6-Sol adjudication: F-P1R R4 publication-terms boundary

Act as the independent adjudicator for Claude Opus 5's conditional
implementation-stage audit of F-P1R step R4 in
`/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work read-only: do not edit
files, launch jobs, select physical sources, contact authors, or redistribute
data. You may run bounded non-writing checks if useful. Inspect the actual
checkout and implementation commit, not only the audit summary. AGY is
retired and must not be called or treated as an active auditor.

Implementation commit: `3e4a4f5` (`Retire AGY and harden publication terms
boundary`). Claude Opus 5 audit:
`provenance/opus5_fp1r_r4_publication_terms_audit_2026-09-04.md`.
Accepted bundle plan:
`provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`.
Current governance:
`provenance/audit_governance_amendment_2026-09-04.md`.

The final project goal is a production-ready and publication-ready lagRamses
high-level hydrodynamics stack focused on radiative transfer, stellar/AGN
feedback, and dust. R4 is evidence/publication hardening before physical
source admission. It must preserve zero physical nodes, unresolved
`[0.8,1.0]` and `[40,120] M_sun` seams, no canonical conversion/runtime
deposition, review-only LC18 output, and blocked publication.

Return exactly one adjudication status: `ACCEPT R4`, `REWORK R4`, or
`BLOCK`. Use `BLOCK` only for an unsound publication authority boundary or a
live fail-closed bypass. Use `REWORK R4` for a material evidence or semantic
defect that must be fixed before F-P1R R3. Do not convert optional hygiene
findings into a block.

## R4 contract

Review:

- `simulation/snrt/tools/fp1_publication_rights.py`
- `simulation/snrt/tests/fp1_publication_rights.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `simulation/snrt/config/g2_source_use_terms_evidence_v1.json`

The evaluator must read and hash the locked terms bytes itself and parse
`sources[candidate_id]` from those same bytes. Caller-supplied digest and
detached source-record authority must be absent. Path mismatch, malformed or
unreadable terms, missing candidate, and mutated bytes must fail closed. The
synthetic lock profile must be private/test-only, with rights inside hashed
fixture bytes. LC18 must retain empty approval, review-only status, and false
publication permission. `require_publication_allowed()` must refuse forged or
mutated gates and re-evaluate current bytes.

Opus identified one finding: on the wrong-path branch the file is not opened,
but `terms_error` remains `None`, causing the blocked report to claim
`publication_terms_record_parsed=true` and
`record_source=candidate_record_parsed_from_locked_terms_bytes` despite an
empty record. Opus classified this as a mandatory non-bypass evidence fix:
use a not-read/path-mismatch sentinel and assert the false/unavailable fields
in the wrong-path test. Opus also listed as non-blocking the lack of coverage
for the post-evaluation mutation branch and the locked-file `OSError` branch,
plus display-only detached labels in LC18.

Independently decide whether the Opus finding is material and required under
the R4 plan, whether it can be deferred without weakening publication
authority, and whether the current implementation has any additional
bypass. Check the driver's focused-test and G2-preflight claims only as
read-only evidence; do not modify or regenerate tracked artifacts.

End with the exact status, a concise rationale, classification of M-R4-1 and
the other findings, required changes (if any), and a direct statement on
whether R3 may begin. AGY is retired and is not part of this adjudication.
