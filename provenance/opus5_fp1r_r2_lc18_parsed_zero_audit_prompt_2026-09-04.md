# Claude Opus 5 R2 implementation-stage audit request

Act as the independent scientific and implementation auditor for completed
F-P1R step R2 in `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work
read-only. Do not edit files, run shell commands or jobs, select physical
sources, contact authors, or redistribute data. Use Claude Opus 5's own
scientific and code judgment.

Audited implementation commit: `00a48ac` (`Clarify LC18 parsed-zero wind
semantics`). Parent R1 pass: `1a13457`. Read the F-P1R plan and the live R2
implementation:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `simulation/snrt/tools/audit_fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/tests/fp1_lc18_failed_wind_crosscheck.py`
- `simulation/snrt/data/fp1_lc18_failed_wind_crosscheck.json`
- `provenance/fp1_lc18_failed_wind_inquiry_packet_2026-09-04.md`
- `simulation/snrt/config/g2_limongi_phase_mass_history_contract_v1.json`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. R2 is a semantics/evidence correction inside
the pre-admission boundary, not physical source activation. The real state
must remain fail-closed: zero physical nodes, unresolved `[0.8,1.0]` and
`[40,120] M_sun` seams, zero canonical rows/runtime deposition, unresolved
failed-wind anomaly, and blocked publication.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for a scientific misinterpretation, count/accounting regression,
live fail-closed bypass, or a code/data wiring defect that invalidates the
R2 result. Use `CONDITIONAL PASS` for a material but non-blocking evidence gap.
Do not treat the unresolved anomaly or lack of physical activation as defects;
those are intentional gates.

## R2 contract

The previous ambiguous fields `cds_terminal_wind_positive_count` and
`cds_terminal_wind_zero_count` must be replaced in live report/tool/tests by
explicit parsed semantics. Confirm that:

1. `cds_terminal_wind_parsed_positive_count` and
   `cds_terminal_wind_parsed_exact_zero_count` are computed from the parsed
   `M_initial - M_total(PSN)` value, not from an unrecorded pipeline rounding;
2. every successful, failed, and all-model scope exposes Table 5 endpoint
   precision `0.01 M_sun`, half-bin `0.005 M_sun`, and
   `physical_zero_inferred: false`, with an interpretation that does not claim
   physical zero wind;
3. outcome accounting remains exactly 52 successful/56 failed, successful
   48 parsed-positive/4 parsed exact-zero, failed 53 parsed-positive/3 parsed
   exact-zero, and all-model 101/7, with the outcome map `{successful: 4,
   failed: 3}`;
4. the three failed parsed exact-zero endpoints are kept inside, but are not
   conflated with or used to resolve, the 56-model BR26 zero-Wind release
   anomaly; successful controls remain separate;
5. the phase contract's precision is the source of the half-bin, malformed or
   non-positive precision fails closed, all existing blockers/flags remain
   false, and no source selection/runtime deposition/redistribution is added;
6. live JSON is deterministically regenerated from the current tool, the
   focused R2 test and physical-package test pass according to the driver's
   claim, and no historical audit report is rewritten. The unsent inquiry
   packet may carry the corrected wording.

Inspect for hidden semantic or wiring errors: stale consumers, old ambiguous
   keys in live code/data, accidentally calling a parsed zero a physical zero,
   changing the release anomaly counts, treating a half-bin as a tolerance,
   or weakening the fail-closed publication/admission state. Check whether the
   test assertions actually reject a future rewording/inference. Separate
   optional style/portability/extra-field improvements from gate failures.
End with the verdict, mandatory fixes (if any), non-blocking findings, and a
direct statement on whether R4 may begin. AGY is retired and must not be
called.
