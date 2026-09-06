# F-P2 SNIa/HESMA bundle — independent Fable audit

- Model: actual `fable` through `claude -p --model fable`
- Date: 2026-09-03
- Scope: interval DTD, event ledger, HESMA review artifacts, approval sidecar, contract audit, and directly relevant tools/tests
- Mode: read-only; no files modified, no model selected, no runtime activated, no commit/push

## Verdict

**CONDITIONAL PASS as review-only evidence. Production activation is correctly BLOCKED.**

The review-only state is genuinely enforced at three independent layers, nothing physical is selected, and the recomputed checksums match. Conditions concern evidence wording, reproducibility, and physical blockers that the contract does not yet name.

## Verified

- Native and production SNIa Fortran mirrors are byte-identical; the relevant config mirrors differ only by the HDF5 schema constant and a comment.
- Sidecar hashes, the manifest package hash, and embedded tool/artifact hashes match current files.
- The DTD kernel has the normalized power-law interval integral, exact logarithmic branch, support clipping, transactional invalid-input rejection, and exact telescoping over a 2001-step split.
- The fail-closed chain for `use_snia=.true.` remains intact: config support check false → initialization error → unsupported driver path.
- Screening labels have no consumer outside tests. Comparison, estimator, packet, and contract audit require selection fields to remain null.
- n100 profile diagnostics are about 1.4532 M_sun and 1.51e51 erg; radial momentum budget is about 2.76e42 g cm/s, but that scalar is not present in the bundle.

## Findings

### High — production blockers not yet named by the contract

1. **WD remnant-reservoir debit is not represented.** `stellar_population_ledger.f90` computes living mass as initial − returned − remnant, while SNIa channel 4 has no remnant ownership. A SNIa ejecta event physically consumes the pre-existing WD reservoir; without an explicit debit or equivalent rule, wiring the event ledger could charge ejecta against living stars. This requires a contract and closure-test design before promotion.

2. **Momentum semantics are incomplete.** The contract has a signed vector deposited into gas momentum. A spherical event has approximately zero net vector by symmetry; a subgrid model needs either an explicit zero-vector convention plus a scalar radial momentum budget, or a documented deposition derivation from energy and ejecta mass. The current bundle has neither scalar radial momentum nor a finalized injection rule.

3. **There is no defined promotion path.** The admission audit accepts only `blocked_review_only`, so a future approval would require changing the audit logic. The adapter also normalizes a warning-bearing model without rejecting it. Define a v2 approval schema and validator, including warning handling and commit binding, before changing the gate.

### Medium — evidence and reproducibility

1. **The shell runner omits the negative-path Python tests.** It runs the audit tools but not `tests/fp2_snia_dtd_contract.py` and `tests/fp2_snia_event_source_admission.py`, so its “complete runner” claim currently covers the positive path only.

2. **`n300c` is a source-data anomaly, not merely a profile estimator warning.** Both profile estimators give the same approximately 641% mass discrepancy; the profile silicon mass alone exceeds the integrated ejecta mass. The source audit status/wording should quarantine or explicitly classify this model as anomalous.

3. **Artifact hashes depend on checkout paths and are not commit-bound.** Several hashed JSON payloads embed absolute paths, and the F-P2 bundle is currently untracked. Use path-free canonical payloads or repository-relative paths and bind the bundle to a commit before publication/production use.

4. **The “not linked into the production binary” statement is inaccurate as written.** The SNIa modules are listed in `stellar_enrichment_sources.mk` and are compiled into the production source set, even though runtime activation remains blocked. Correct the evidence wording or regenerate linked-build evidence.

5. **Promotion requirements omit normalization/IMF conversion, expectation-versus-Poisson realization, SNIa thermal coupling, and metallicity dependence.** These must be explicit decisions in the approval contract.

6. **The 5% screen is hardcoded and knife-edged.** All clean models show a common 3–5% profile excess; n1600c at 5.11% is flagged while n1600 at 4.68% is not. Make the threshold a contract field and report the common bias; retain it as diagnostic only.

### Low

- `kernel.alpha` is populated while other physics fields are null and is not checked as a fixture-only value.
- The DTD power-law branch has a cancellation risk for exponents just outside the 1e−12 logarithmic threshold.
- NaN/Inf, zero-width, α<−1, and additional ledger boundary cases need tests.
- `selected_model_id: n100` in a review fixture reads like a selection; rename to `review_fixture_model_id`.
- Source audit uses `sum` while adapter/comparison use `math.fsum`.
- The decay horizon may be cheaply closed by documenting the fully-decayed convention and sensitivity.

## Required before production approval

1. Design/test WD remnant debit, momentum convention, and a v2 approval schema.
2. Add negative-path tests to the runner; quarantine/classify n300c; remove absolute paths from hashed payloads; bind artifacts to a commit; correct the linking claim; and extend promotion prerequisites.
3. Then approve the progenitor/population model, DTD bounds and normalization/IMF conversion, decay mapping, returned mass/remnant, energy scalar, momentum realization, and named approval.

## Review-only assessment

Correctly enforced. Residual risk is primarily in evidence wording/reproducibility and in the unresolved physical contract, not in an accidental runtime activation.
