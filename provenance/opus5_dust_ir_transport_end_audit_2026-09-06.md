# Opus 5 DUST-3 bundle-end audit

Model: claude-opus-5. Session: 241e3d3b-932f-41e9-bba9-aa79d6509676.
Read-only audit; no edits/jobs/web/subagents. This record summarizes the
returned report and separately records the driver's disposition.

## Auditor verdict: CONDITIONAL PASS

No correctness defect found in algorithm, units, conservation or binding.
All five plan items and eight Fable conditions are implemented and exercised.
The audit certifies bounded bookkeeping evidence, not scientific accuracy.

Validated strengths:
- Shared total-power interpolation and segment slopes preserve weak CMB
  excess without subtractive cancellation. The temperature-dependent
  synthetic band fraction detects the old f(T)*delta_P approximation.
- Emission and opacity share the same reference-mixture density scaling;
  Planck band mean has the correct units and weighting.
- Transport uses energy density; photon counts are source diagnostics, so
  changes in temperature do not reprice radiation already in flight.
- Six-face outgoing flux matches telescoping upwind loss. The random
  anisotropic source-free test independently verifies boundary conservation.
- Every nonlinear source iterate uses the same old field. Reabsorption is
  recycled and finite iteration failure is rejected. Global residual and
  local fixed-point residual are related, rather than independent proofs;
  the independent boundary test makes the combined evidence sufficient.
- Actual static/P5/opacity/thermal/raw-table hashes are enforced. A valid
  but differently hashed static input is rejected before output creation.
- Work is proportionate: one module, runner, focused test, P5 hash attribute,
  documentation; no new gate framework or opacity schema.

## Findings and dispositions

1. Spectral-complement free escape omits nonzero far-IR absorption and
   underestimates trapping. **Adopted:** document the direction of the bias
   in the guide and output semantics; resolve before science/native IR.
2. Fixed-reference Planck absorption plus temperature-dependent emissivity
   does not ensure frequency-resolved detailed balance. **Retained promotion
   prerequisite:** temperature/frequency-consistent opacity before science or
   native IR. The measured 20/50 K spread is not a universal error bound.
3. Physical runs are optically thin; the approximately 2.89x reprocessing
   ratio merely propagates the opacity ratio. **Adopted:** evidence says so,
   and the existing physical runner test prints stationarity and cell tau.
   Nonlinear iteration evidence comes from the synthetic tau_box=0.8 cube.
4. Nonconvergence and thermal overflow had indistinguishable errors.
   **Implemented:** expose the reduced thermal-invalid flag; runner and
   tests distinguish thermal input/range invalidity from nonconvergence.
5. Missing negative tests for group upper edge >1 eV, CMB coverage and
   invalid angular weights. **Implemented:** three tests in the same suite.
6. Directly constructed closure could have a flat total-power segment.
   **Implemented:** strict finite positive increasing power validation plus
   a targeted flat-segment rejection test.
7. Diagnostic temperature rounds to T_CMB at very small excess, though
   emitted rates stay positive; inactive temperature is a zero sentinel.
   **Adopted:** guide/output explicitly prohibit reconstructing weak-excess
   source SEDs from that rounded temperature.
8. Output dataset units were missing. **Implemented:** energy, temperature
   and emitted photon datasets now carry explicit units attributes.

Driver judgment: the auditor's statement that a damped contraction converges
must be understood asymptotically, not as a guarantee within 128 iterations.
The implementation still rejects a result that exceeds the finite limit.
No extra audit round is needed for these bounded guard/test/label repairs;
rerun the two existing suites and record the results in bundle evidence.

No analytic radiative-transfer solution was added in this bundle; magnitude
accuracy inherits from the existing transport kernel. The dt/mesh/angular
comparisons remain bounded sensitivities, not spatial-order or production
accuracy claims. Native/live momentum and hydro coupling remain unapproved.
