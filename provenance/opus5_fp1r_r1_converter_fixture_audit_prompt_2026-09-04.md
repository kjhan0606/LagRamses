# Claude Opus 5 R1 implementation-stage audit request

Act as the independent scientific and implementation auditor for completed
F-P1R step R1 in `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`). Work
read-only. Do not edit files, run shell commands or jobs, select physical
sources, contact authors, or redistribute data. Use Claude Opus 5's own
judgment. The audited implementation commit is `a514fd5`; the accepted
bundle plan and Grok plan audit are in:

- `provenance/fp1r_bundle_plan_evidence_semantics_execution_closure_2026-09-04.md`
- `provenance/grok_fp1r_bundle_plan_audit_2026-09-04.md`

The final project purpose is a production-ready and publication-ready
lagRamses high-level hydrodynamics stack focused on radiative transfer,
stellar/AGN feedback, and dust. R1 is only evidence hardening before physical
yield/fate activation. It must preserve the real fail-closed state: zero
physical nodes, unresolved `[0.8,1.0]` and `[40,120] M_sun` fate seams, no
canonical conversion or runtime deposition, and a blocked LC18 publication
gate.

Return exactly one top-level verdict: `PASS`, `CONDITIONAL PASS`, or `BLOCK`.
Use `BLOCK` for a defect that makes the converter admission boundary or
fail-closed behavior unsound. Use `CONDITIONAL PASS` for material but
non-blocking evidence or hardening gaps. Do not treat the intentional blocked
production state as a defect.

## R1 contract to audit

Inspect `simulation/snrt/tests/yield_converter.py` and the converter and
mapping code it exercises. The fixture is required to:

1. patch only converter-module test seams:
   `DEFAULT_SOURCE_NODE_CONTRACT`, `DEFAULT_PHYSICAL_PACKAGE_CONTRACT`, and
   the bound `audit_source_node_contract` /
   `audit_physical_package_admission` functions;
2. use temporary files and in-memory admitted mapping/selection data, derive
   the admitted mapping from the non-writing proposal rows, and restore every
   mutable seam in a `finally` block;
3. reach the real converter positive write path and write only a temporary
   table, sidecar, and source-node mapping;
4. verify sidecar/mapping asset and mapping SHA-256 values;
5. test a matching case plus mapping-content mutation with a recomputed hash,
   mapping-content mutation without recomputing the declared hash, and
   hash-only mutation. Every mismatch must fail before any output path exists;
6. retain the existing blocked-path tests and prove that, after seam restore,
   the genuine repository conversion remains fail-closed;
7. compare all `simulation/snrt/config` and `simulation/snrt/data` file hashes
   before and after the fixture, with no source/config/data changes.

Assess the actual implementation, not only comments. Check whether the
synthetic admitted package is a test-only seam rather than a production
override, whether the mapping equality/hash checks in
`simulation/snrt/tools/convert_yield_rows_to_canonical.py` remain authoritative
and before-write, whether the mutation cases genuinely exercise distinct
failure modes, whether exception/finally behavior restores state on failure,
and whether the test itself has hidden false positives or a mismatch with the
accepted plan. Consider Python/runtime compatibility and test maintainability,
but classify those separately from admission soundness.

The driver independently observed `YIELD_CONVERTER_TEST_OK` from the R1 test,
and the real config/data hash snapshot was identical before/after. Treat that
as a claim to inspect against the code; do not blindly accept it. No R2/R3/R4
work has started. End with the verdict and a concise list of mandatory fixes,
non-blocking findings, and whether R2 may begin.
