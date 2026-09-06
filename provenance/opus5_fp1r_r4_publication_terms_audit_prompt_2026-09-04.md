# Claude Opus 5 R4 implementation-stage audit request

Act as the independent scientific, algorithmic, wiring, and implementation
auditor for completed F-P1R step R4 in `/gpfs/kjhan/LRD_JWST`
(`kjhan0606/LagRamses`). Work read-only. Do not edit files, run shell
commands or jobs, select physical sources, contact authors, or redistribute
data. Use Claude Opus 5's own judgment. AGY is retired and must not be called
or treated as an active approval authority.

Audited implementation commit: `3e4a4f5` (`Retire AGY and harden publication
terms boundary`). Parent R2 pass: `00a48ac`. Read the F-P1R plan, current
governance, implementation, live test, and regenerated evidence:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/audit_governance_amendment_2026-09-04.md`
- `simulation/snrt/tools/fp1_publication_rights.py`
- `simulation/snrt/tests/fp1_publication_rights.py`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `simulation/snrt/config/g2_source_use_terms_evidence_v1.json`
- `provenance/fp1_lc18_failed_wind_inquiry_packet_2026-09-04.md`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. R4 is a provenance/publication boundary
inside the pre-admission path, not physical source activation. The real state
must remain fail-closed: zero physical nodes, unresolved `[0.8,1.0]` and
`[40,120] M_sun` seams, zero canonical rows/runtime deposition, unresolved
failed-wind anomaly, review-only LC18 output, and blocked publication.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for a caller-controlled rights/digest/record bypass, incorrect
same-bytes semantics, a live publication opening, a forged-gate acceptance,
or a code/data wiring defect that invalidates R4. Use `CONDITIONAL PASS` for a
material but non-blocking evidence gap. Do not treat the unresolved source
rights, anomaly, or lack of physical activation as defects; those are
intentional gates.

## R4 contract

1. The production evaluator must receive a candidate id and terms path, use a
   code-owned lock profile, read the terms bytes itself, compute SHA-256
   internally, and parse `sources[candidate_id]` from those same bytes. A
   caller-supplied digest or detached `source_record` must not be accepted as
   production authority. Missing keys, malformed JSON, missing candidates,
   wrong paths, unreadable files, and mutated bytes must fail closed.
2. Synthetic lock profiles are permitted only through an explicit private
   in-memory test seam. Synthetic rights fields must live inside the hashed
   terms bytes. The test must cover a positive admitted fixture, review-only
   and approval failures, mutated bytes, wrong paths, malformed JSON, and a
   missing candidate record. Check that the positive path reaches the normal
   evaluator and that no fixture edits a repository terms file.
3. The LC18 caller must use the new API without passing a digest or detached
   rights record. The live terms catalog must not be edited to claim rights;
   LC18 must remain review-only with an empty approval and false publication
   gate. The multi-source terms catalog remains the source of the candidate
   record.
4. `require_publication_allowed()` must not accept a hand-built or shallow
   copied dictionary whose labels say `allowed`, `publication_ready`, and
   `review_use_only=false`. It must detect mutation of the evaluated gate and
   re-evaluate current terms bytes before allowing an export. No production
   export caller may be introduced by this step.
5. The active governance must contain Claude Opus 5 as the implementation
   auditor, Grok as the bundle-start plan auditor, GPT-5.6 Sol as the
   conditional/negative Opus adjudicator, and no AGY role. Historical AGY
   reports may remain provenance and must not be treated as current approval.

Inspect for hidden semantic or wiring errors: any surviving old API call,
rights fields sourced from an object other than the bytes that were hashed,
an allowed result that ignores a terms parse error, test-only override
leaking into production, a mutable gate that can be forged, terms-file edits,
or an accidental change to the fail-closed physics/admission state. Check
whether the test assertions actually cover each required failure mode and
whether the generated LC18 evidence reflects the new gate schema. Separate
optional style, portability, and stronger type-validation suggestions from
gate failures.

End with the verdict, mandatory fixes (if any), non-blocking findings, and a
direct statement on whether F-P1R R3 may begin. AGY is retired and must not be
called.
