# Fable DUST-5 plan review

Model claude-fable-5-1; session 134e1f0e-df6f-4f36-86f4-97030866b936;
duration 193 s. Read-only review; no jobs, edits or subagents.
This is a summary and driver disposition, not a verbatim transcript.

Verdict: **CONDITIONAL APPROVE**. Aligned with DUST-4's native handoff and
the final native RT/feedback/dust objective. Separate FP64 energy field is
correct: primary native state has nine FP32 photon groups, fixed direction
count, different angular normalization and atom budgets. Reusing it would
mix units, groups and semantics. One module, smoke driver and opt-in existing
test extension are proportionate; no new instrumentation framework.

Conditions adopted:

1. Native admission checks one-normalized weights, unit directions, positive
   quadrature/sigma, exact CMB node, strictly increasing T and total power,
   reciprocal neighbors, CFL<=1+1e-12; distinct table/state/config/shape/CFL/
   thermal-range/nonconvergence errors.
2. Match bath-relative weak source, log-T linear-in-power interpolation,
   thin-cell quadratic response below 1e-4, old-field vacuum escape,
   zero-primary old-inventory scale and half-relaxed bounded iteration.
   Both solves use 1e-9 stop tolerance; if stop-iterate differences alone
   violate 1e-8 comparison, rerun BOTH at 1e-12 rather than loosen the bound.
3. Force failure and prove bit-identical caller energy, temperature, photon
   diagnostics and scalar diagnostics. Zero/weak-source checks remain.
4. Add the native option explicitly; default existing test remains compiler-
   free. One plain temporary numeric fixture; no JSON/evidence framework.
   Module depends only on real64/IEEE, two-source smoke with GNU and Intel.
5. Only object-list/dependency changes in the dirty Makefile. Intel compile
   with production preprocessing/optimization in a private directory; do
   not claim full link/runtime qualification.

Deferred unchanged: abundance/depletion/grain evolution and physical source
approval, persistent live spectral layout/resolution, force/gas deposition,
AMR/MPI/restart/accelerated IR. No global source-to-binary or infrastructure
gate, no full production link, no additional Fable review before Opus end.
