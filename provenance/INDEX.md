# High-level RT/feedback/dust provenance index

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Index status: active engineering records; historical records remain in place
until an explicitly approved archive operation.
Record statuses: `active` = current contract/evidence, `superseded` = retained
but not a current gate, `long_term_record` = future G5/G6/science work.

This index is deliberately bundle-level. Individual native smokes and Python
validators are not separate project gates; they are referenced by the bundle
evidence that runs them. Audit prompts are request artifacts, not independent
scientific results. Scratch copies under `~/.claude/plans` are not project
evidence.

## Active bundle map

| Bundle | Status / purpose | Plan | Implementation evidence | Plan audit | End audit / disposition |
|---|---|---|---|---|---|
| F-P1.2 stellar feedback transaction | active, engineering contract | [plan](fp1_2_stellar_feedback_transaction_bundle_plan_2026-09-04.md) | [evidence](fp1_2_stellar_feedback_transaction_bundle_implementation_evidence_2026-09-04.md) | — | [Opus final](opus5_fp1_2_stellar_feedback_transaction_bundle_final_audit_2026-09-04.md) |
| F-P1.5 AGN effective efficiency | active, ledger convention | [plan](fp1_5_agn_effective_efficiency_convention_bundle_plan_2026-09-04.md) | [evidence](fp1_5_agn_effective_efficiency_convention_bundle_implementation_evidence_2026-09-04.md) | [Opus plan](opus5_fp1_5_agn_effective_efficiency_convention_plan_audit_2026-09-04.md) | [Opus final](opus5_fp1_5_agn_effective_efficiency_convention_bundle_final_audit_2026-09-05.md) |
| AGN native coupling | implemented; bounded closure PASS; live SNRT+AGN serial fresh-start only; no simultaneous legacy+SNRT | [approved plan](agn_native_coupling_bundle_plan_2026-09-05.md) | [native tests and caller compiles](agn_native_coupling_bundle_evidence_2026-09-05.md) | [Fable plan](fable_agn_native_coupling_plan_audit_2026-09-05.md) | [Opus end](opus5_agn_native_coupling_end_audit_2026-09-05.md), [Opus B1–B3 closure](opus5_agn_native_coupling_closure_audit_2026-09-05.md) |
| F-P1H-E validator admission | active, high-mass admission boundary | [plan](fp1h_e_validator_admission_bundle_plan_2026-09-04.md) | [closure](fp1h_e_validator_admission_bundle_closure_2026-09-04.md) | — | [Opus end](opus5_fp1h_e_bundle_end_audit_2026-09-04.md) |
| F-P2 source/SED/dust closure | active, source-to-dust boundary | [plan](fp2_source_sed_dust_closure_bundle_plan_2026-09-05.md) | [evidence](fp2_source_sed_dust_closure_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_source_sed_dust_closure_bundle_end_audit_2026-09-05.md) |
| F-P2.1 source closure verification | active, native source identity | [plan](fp2_1_source_closure_verification_bundle_plan_2026-09-05.md) | [evidence](fp2_1_source_closure_verification_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_1_source_closure_verification_bundle_end_audit_2026-09-05.md) |
| F-P2.2 closure integrity | active, ledger closure | [plan](fp2_2_closure_integrity_bundle_plan_2026-09-05.md) | [evidence](fp2_2_closure_integrity_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_2_closure_integrity_bundle_end_audit_2026-09-05.md) |
| F-P2.3 canonical quadrature | active, native angular/source quadrature | [plan](fp2_3_canonical_asset_quadrature_bundle_plan_2026-09-05.md) | [evidence](fp2_3_canonical_asset_quadrature_bundle_implementation_evidence_2026-09-05.md) | — | [Opus follow-up](claude_opus5_fp2_3_canonical_asset_quadrature_bundle_followup_audit_2026-09-05.md) |
| F-P2.4 native nine-group spectral contract | active, RT group contract | [plan](fp2_4_native_nine_group_spectral_contract_bundle_plan_2026-09-05.md) | [evidence](fp2_4_native_nine_group_spectral_contract_bundle_implementation_evidence_2026-09-05.md) | — | [Opus closure](claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md) |
| F-P2.5 native H/He thermochemistry | active, chemistry receiver | [plan](fp2_5_native_hhe_thermochemistry_bundle_plan_2026-09-05.md) | [evidence](fp2_5_native_hhe_thermochemistry_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_end_audit_2026-09-05.md) |
| F-P2.6 native RT transaction | active, transport/chemistry transaction | [plan](fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md) | [evidence](fp2_6_native_rt_chemistry_transaction_bundle_implementation_evidence_2026-09-05.md) | [Fable plan](fable_fp2_6_native_rt_chemistry_transaction_bundle_plan_audit_2026-09-05.md) | [Fable closure](fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md), conditional pass |
| F-P2.7 gate consolidation / initialized-RAMSES | active, current engineering bundle | [plan](fp2_7_gate_consolidation_initialized_ramses_bundle_plan_2026-09-05.md) | pending | [Fable plan](fable_fp2_7_gate_consolidation_initialized_ramses_bundle_plan_audit_2026-09-05.md) | pending; implementation in progress |

## Cross-cutting decisions and non-approvals

- [production/publication readiness plan](production_publication_readiness_plan.md)
  remains the governing high-level RT/stellar/AGN feedback/dust roadmap.
- [audit cadence amendment](audit_cadence_amendment_2026-09-05.md) governs one
  plan audit plus one end audit per bundle, with fallback only when needed.
- [Fable operational audit](fable_operational_instrumentation_gate_efficiency_audit_2026-09-05.md)
  classified the workspace `OVERINSTRUMENTED`; F-P2.7 is the approved
  consolidation response.
- Physical SED/yield/DTD/PISN approval, live production feedback, dust science
  approval, HDF5 restart, distributed AMR scaling, and publication convergence
  are not granted by any F-P2 record above.

## Classification rules

1. A bundle evidence file is the authoritative entry point for its native
   smokes, build/link, hashes, and conclusions.
2. A focused runner may be used for debugging, but its output is not a new
   gate or audit event.
3. Historical audit reports remain immutable records; later dispositions are
   appended to the bundle evidence/index rather than creating repeated audits.
4. Large HDF5, JAX, virtual-environment, build, and generated compiler
   products remain outside this index's active evidence unless a bundle
   explicitly pins a manifest/hash for them.
