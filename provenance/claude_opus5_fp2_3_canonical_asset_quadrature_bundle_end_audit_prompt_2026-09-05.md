# F-P2.3 bundled end audit — canonical assets and source-closure quadrature

You are the single bundled Claude Opus 5 end auditor for F-P2.3 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`.

This is a read-only algorithm/implementation/provenance audit. Do not edit
files, do not launch simulations or tests, and do not use shell or external
network tools. Inspect the repository using only the permitted read/search
tools. The final purpose of this project is production-ready, publication-ready
high-level radiation transport and stellar/AGN feedback/dust physics in the
LagRamses/RAMSES workflow. Python helpers in this bundle are offline source
conversion, reference-oracle, and evidence tools; they are not being claimed
as the final Fortran RAMSES runtime.

This is a follow-up to the same-day F-P2.3 audit, whose report is included in
the reading list. Verify specifically that its MAJOR-1/MAJOR-2 repairs are
real and that the listed medium/minor consistency repairs did not introduce a
new mismatch. Return a fresh decisive verdict for the bundle.

Audit the complete F-P2.3 bundle, not just individual lines. Read:

- `provenance/fp2_3_canonical_asset_quadrature_bundle_plan_2026-09-05.md`
- `provenance/fp2_3_canonical_asset_quadrature_bundle_implementation_evidence_2026-09-05.md`
- `provenance/claude_opus5_fp2_2_closure_integrity_bundle_end_audit_2026-09-05.md`
- `simulation/snrt/snrt_core/sed.py`
- `simulation/snrt/snrt_core/primordial.py`
- `simulation/snrt/snrt_core/dust.py`
- `simulation/snrt/snrt_core/provenance.py`
- `simulation/snrt/tools/p4_build_agn_photon_ledger.py`
- `simulation/snrt/tools/build_draine_dust_opacity.py`
- `simulation/snrt/tools/validate_agn_nine_group_ledger.py`
- `simulation/snrt/tests/agn_nine_group_artifact.py`
- `simulation/snrt/tests/agn_nine_group_explicit_artifact.py`
- `simulation/snrt/tests/source_sed_dust_closure.py`
- `simulation/snrt/data/agn_nine_group_external_assets.json`
- `simulation/snrt/data/agn_nine_group_validation.json`
- `simulation/snrt/data/agn_nine_group_explicit_validation.json`
- `simulation/snrt/data/p4_explicit_agn_photon_ledger.json`
- `simulation/snrt/P4_AGN_RATE_LEDGER.md`
- `simulation/snrt/P4_DUST_OPACITY.md`

Assess at least these questions:

1. Does the external-asset manifest validate every declared asset with exact
   ID, canonical path, byte size, and current SHA-256, including both pilot
   and explicit controls? Are reproduction commands complete and truthful?
2. Are pilot and explicit artifact tests genuinely non-vacuous and symmetric
   about current live hashes, HEAD, and the `simulation/snrt` working-tree
   attestation? Is the attestation path evaluated from the actual repository
   root, and is a dirty development tree represented honestly?
3. Is the explicit SED convention scientifically and numerically coherent:
   piecewise-linear `f_E`, `q_E=f_E/(E*EV_ERG)`, exact group boundaries,
   base/refined quadrature, closed failure behavior, and identity binding?
   Does the independent validator reconstruct explicit group moments and
   Verner closure without simply trusting the serialized values?
4. Is the same declared convention correctly applied to source-weighted dust,
   including duplicated energy-fraction checks and source/status provenance?
5. Do the artifact and focused-test results support an engineering PASS while
   avoiding an astrophysical claim? Are the synthetic explicit SED, zero-dust
   controls, unresolved physical SED/obscuration, and absent live Fortran
   coupling stated clearly enough for publication provenance?
6. Identify any remaining correctness, auditability, reproducibility, or scope
   issue that must block this bundle, and distinguish it from a minor follow-up
   or a later high-level RT/feedback/dust bundle.

Return a decisive verdict exactly as one of `PASS`, `CONDITIONAL PASS`, or
`FAIL`, followed by severity-ranked findings with file/line references where
possible, the evidence supporting each finding, and explicit disposition for
each F-P2.3 acceptance gate. Do not invent runtime or astrophysical approval
that the inspected evidence does not establish.
