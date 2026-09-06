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
| AGN event-level entrainment | implemented; Opus PASS within scope; no live/production qualification | [approved amended plan](agn_entrainment_bundle_plan_2026-09-05.md) | [native event checks and caller compiles](agn_entrainment_bundle_evidence_2026-09-05.md) | [Fable conditional approve and driver disposition](fable_agn_entrainment_bundle_plan_audit_2026-09-05.md) | [Opus end PASS](opus5_agn_entrainment_bundle_end_audit_2026-09-05.md) |
| AGN accepted fuel / overlapping deposition | implemented; bounded Opus closure PASS; no simultaneous RT/mechanical, live or restart qualification | [approved amended plan](agn_accepted_fuel_overlap_bundle_plan_2026-09-05.md) | [native tests, caller compiles and repair dispositions](agn_accepted_fuel_overlap_bundle_evidence_2026-09-06.md) | [Fable conditional approve](fable_agn_accepted_fuel_overlap_plan_audit_2026-09-05.md) | [Opus end](opus5_agn_accepted_fuel_overlap_end_audit_2026-09-06.md), [focused closure PASS](opus5_agn_accepted_fuel_overlap_closure_audit_2026-09-06.md) |
| SNRT AGN source-to-RT transaction atomicity | bounded ordering repair implemented; native rollback/commit and CPU/DUST_LIVE caller compiles pass; Opus conditional pass; driver-faithful failure coverage and physical/live qualification remain | [parent accepted-fuel plan](agn_accepted_fuel_overlap_bundle_plan_2026-09-05.md) | [implementation evidence](snrt_agn_rt_atomicity_bundle_implementation_evidence_2026-09-07.md) | — | [Opus end audit](opus5_snrt_agn_rt_atomicity_bundle_end_audit_2026-09-07.md); conditions retained |
| AGN accretion-powered partition reference | implemented; Opus CONDITIONAL PASS repairs applied; consolidated native/production build gate PASS; comparison-only, no restart/live qualification | [approved plan](agn_partition_reference_bundle_plan_2026-09-06.md) | [implementation evidence](agn_partition_reference_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_agn_partition_reference_plan_audit_2026-09-06.md) | [Opus bundle-end CONDITIONAL PASS and repair record](opus5_agn_partition_reference_bundle_end_audit_2026-09-06.md) |
| Dust scattering/state channel | implemented; static JAX/P4/P5 conditional candidate; Opus final CONDITIONAL PASS; v1/v2 absorption controls retained; no IR/live/native promotion | [plan](dust_scattering_state_bundle_plan_2026-09-06.md) | [implementation evidence](dust_scattering_state_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_scattering_state_bundle_plan_audit_2026-09-06.md) | [Opus initial](opus5_dust_scattering_state_bundle_end_audit_2026-09-06.md), [re-audit](opus5_dust_scattering_state_bundle_reaudit_2026-09-06.md), [final](opus5_dust_scattering_state_bundle_final_audit_2026-09-06.md) |
| Dust grain thermal / IR candidate | implemented; static JAX/P5 conditional candidate; one-pass recorded-not-transported; no recursive/native/live promotion | [plan](dust_ir_thermal_bundle_plan_2026-09-06.md) | [implementation evidence](dust_ir_thermal_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_ir_thermal_bundle_plan_audit_2026-09-06.md) | [Opus bundle-end conditional pass](opus5_dust_ir_thermal_bundle_end_audit_2026-09-06.md) |
| Dust excess IR transport | implemented static gray study; frozen P5 primary heating; local IR reabsorption; no native/live promotion | [plan](dust_ir_transport_bundle_plan_2026-09-06.md) | [evidence](dust_ir_transport_bundle_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_ir_transport_plan_audit_2026-09-06.md) | [Opus conditional pass and bounded repair disposition](opus5_dust_ir_transport_end_audit_2026-09-06.md) |
| Dust spectral IR transport | implemented static spectral candidate; full raw-domain transport and consistent node opacity; no native/live promotion | [plan](dust_spectral_ir_bundle_plan_2026-09-06.md) | [evidence](dust_spectral_ir_bundle_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_spectral_ir_plan_audit_2026-09-06.md) | [Opus conditional pass; all bounded repairs verified](opus5_dust_spectral_ir_end_audit_2026-09-06.md) |
| Dust native spectral IR operator | implemented FP64 Fortran candidate; GNU/Intel differential and rollback pass; no live driver activation | [plan](dust_native_ir_bundle_plan_2026-09-06.md) | [evidence](dust_native_ir_bundle_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_native_ir_plan_audit_2026-09-06.md) | [Opus conditional pass; two documentation conditions applied](opus5_dust_native_ir_end_audit_2026-09-06.md) |
| Dust primary-photon coupling | implemented native FP64 reference boundary; H/He saturation, returned photons, heating ledger and zero-dust regression pass; Opus conditional pass with bounded repairs applied; no CUDA/live activation | [plan](dust_primary_absorption_bundle_plan_2026-09-06.md) | [evidence](dust_primary_absorption_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_primary_absorption_plan_audit_2026-09-06.md) | [Opus end conditional pass](opus5_dust_primary_absorption_bundle_end_audit_2026-09-06.md) |
| Dust fourth-species CUDA boundary | implemented native FP32 H/He+dust ABI; legacy zero-dust bitwise regression and A10 GNU/Intel execution pass; Fable conditional pass with coverage/documentation conditions applied; no live RAMSES wiring | [plan](dust_cuda_species_dust_bundle_plan_2026-09-06.md) | [evidence](dust_cuda_species_dust_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_cuda_species_dust_bundle_plan_audit_2026-09-06.md) | [Fable conditional pass](fable_dust_cuda_species_dust_bundle_end_audit_2026-09-06.md) |
| Dust CUDA→FP64 transaction handoff | implemented native prepared ledger handoff; zero-dust scaffold only; production link, A10 CUDA, FP64 receiver, rollback/MPI and negative bundle gate pass; no nonzero live dust activation | [plan](dust_transaction_ledger_bundle_plan_2026-09-06.md) | [implementation evidence](dust_transaction_ledger_bundle_implementation_evidence_2026-09-06.md) | [Fable conditional approve](fable_dust_transaction_ledger_bundle_plan_audit_2026-09-06.md) | self end audit in evidence; nonzero dust/thermal/momentum remain deferred |
| Dust source-bound mapping / thermal receiver | implemented native FP64 opacity-binding check, explicit metallicity×dust-to-metal mapping, optical-depth construction, staged persistent thermal receiver and rollback smoke; live driver remains ZERO_SCAFFOLD | [plan](dust_mapping_receiver_bundle_plan_2026-09-06.md) | [implementation evidence](dust_mapping_receiver_bundle_implementation_evidence_2026-09-06.md) | — | [Codex end audit](dust_mapping_receiver_bundle_self_end_audit_2026-09-06.md); live dust state/physical grain closure remain deferred |
| Dust native contract admission | implemented bounded native opacity/thermal namelist admission, provenance/status checks, candidate-only runtime gate, invalid-load reset, GNU/Intel smoke and production-link integration; live driver remains ZERO_SCAFFOLD | [plan](dust_native_contract_admission_bundle_plan_2026-09-06.md) | [implementation evidence](dust_native_contract_admission_bundle_implementation_evidence_2026-09-06.md) | — | [Codex end audit](dust_native_contract_admission_bundle_self_end_audit_2026-09-06.md); physical dust approval/live state remain deferred |
| Cross-subsystem wiring/backend selection | implemented bounded wiring revalidation; dead `gpu_sink` autotune removed; SNRT CUDA-required, stellar/AGN CPU/OpenMP explicit, dust hybrid explicit; no automatic backend switching | [plan](wiring_backend_selection_bundle_plan_2026-09-06.md) | [implementation evidence](wiring_backend_selection_bundle_implementation_evidence_2026-09-06.md) | [Fable unavailable record](fable_wiring_backend_selection_bundle_plan_audit_2026-09-06.md) | [Codex end audit](wiring_backend_selection_bundle_self_end_audit_2026-09-06.md) |
| F-P1H-E validator admission | active, high-mass admission boundary | [plan](fp1h_e_validator_admission_bundle_plan_2026-09-04.md) | [closure](fp1h_e_validator_admission_bundle_closure_2026-09-04.md) | — | [Opus end](opus5_fp1h_e_bundle_end_audit_2026-09-04.md) |
| F-P2 source/SED/dust closure | active, source-to-dust boundary | [plan](fp2_source_sed_dust_closure_bundle_plan_2026-09-05.md) | [evidence](fp2_source_sed_dust_closure_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_source_sed_dust_closure_bundle_end_audit_2026-09-05.md) |
| F-P2.1 source closure verification | active, native source identity | [plan](fp2_1_source_closure_verification_bundle_plan_2026-09-05.md) | [evidence](fp2_1_source_closure_verification_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_1_source_closure_verification_bundle_end_audit_2026-09-05.md) |
| F-P2.2 closure integrity | active, ledger closure | [plan](fp2_2_closure_integrity_bundle_plan_2026-09-05.md) | [evidence](fp2_2_closure_integrity_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_2_closure_integrity_bundle_end_audit_2026-09-05.md) |
| F-P2.3 canonical quadrature | active, native angular/source quadrature | [plan](fp2_3_canonical_asset_quadrature_bundle_plan_2026-09-05.md) | [evidence](fp2_3_canonical_asset_quadrature_bundle_implementation_evidence_2026-09-05.md) | — | [Opus follow-up](claude_opus5_fp2_3_canonical_asset_quadrature_bundle_followup_audit_2026-09-05.md) |
| F-P2.4 native nine-group spectral contract | active, RT group contract | [plan](fp2_4_native_nine_group_spectral_contract_bundle_plan_2026-09-05.md) | [evidence](fp2_4_native_nine_group_spectral_contract_bundle_implementation_evidence_2026-09-05.md) | — | [Opus closure](claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_closure_audit_2026-09-05.md) |
| F-P2.5 native H/He thermochemistry | active, chemistry receiver | [plan](fp2_5_native_hhe_thermochemistry_bundle_plan_2026-09-05.md) | [evidence](fp2_5_native_hhe_thermochemistry_bundle_implementation_evidence_2026-09-05.md) | — | [Opus end](claude_opus5_fp2_5_native_hhe_thermochemistry_bundle_end_audit_2026-09-05.md) |
| F-P2.6 native RT transaction | active, transport/chemistry transaction | [plan](fp2_6_native_rt_chemistry_transaction_bundle_plan_2026-09-05.md) | [evidence](fp2_6_native_rt_chemistry_transaction_bundle_implementation_evidence_2026-09-05.md) | [Fable plan](fable_fp2_6_native_rt_chemistry_transaction_bundle_plan_audit_2026-09-05.md) | [Fable closure](fable_fp2_6_native_rt_chemistry_transaction_bundle_closure_audit_2026-09-05.md), conditional pass |
| F-P2.7 gate consolidation / initialized-RAMSES | active, current engineering bundle | [plan](fp2_7_gate_consolidation_initialized_ramses_bundle_plan_2026-09-05.md) | pending | [Fable plan](fable_fp2_7_gate_consolidation_initialized_ramses_bundle_plan_audit_2026-09-05.md) | pending; implementation in progress |

## Next bundle queue

| Bundle | Status / purpose | Plan | Plan audit | Authorization |
|---|---|---|---|---|
| F-P1H-F physical-source admission / high-mass seam | approved reduced scope; KL16/CK22 investigation parked by operator; physical admission remains unapproved, not an engineering-wide blocker | [plan and wrap-up](fp1h_f_physical_source_admission_bundle_plan_2026-09-05.md) | [Fable](fable_fp1h_f_physical_source_admission_bundle_plan_audit_2026-09-05.md), historical conditional approve | operator wrap-up 2026-09-05; do not automatically reopen AGB source research |

## Cross-cutting decisions and non-approvals

- KL16/CK22 AGB source discrepancies are parked, not resolved or approved.
  Evidence and readers remain in the [existing bundle record](fp1h_f_physical_source_admission_bundle_plan_2026-09-05.md).
  No author email was sent; contact is cancelled. The source-reader checks
  require explicit `--include-parked-agb`; standard selection/admission
  protections remain active. Continue source-independent native feedback work.
- [LC18 failed-wind inquiry and source review](fp1_lc18_failed_wind_inquiry_packet_2026-09-04.md)
  updated 2026-09-05: printed format is not physical precision; original-wind
  comparison and mixed-source mass deficit measured; no wind replacement,
  author contact, new gate or production activation. Supersedes the old
  parsed-zero precision interpretation, not the historical audit record.
- [production/publication readiness plan](production_publication_readiness_plan.md)
  remains the governing high-level RT/stellar/AGN feedback/dust roadmap.
- [audit cadence amendment](audit_cadence_amendment_2026-09-05.md) governs one
  plan audit plus one end audit per bundle, with fallback only when needed.
- [Fable operational audit](fable_operational_instrumentation_gate_efficiency_audit_2026-09-05.md)
  classified the workspace `OVERINSTRUMENTED`; F-P2.7 is the approved
  consolidation response.
- SNIa HESMA N100 and Maoz DTD retain their prior physical approval
  `FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ`; see
  [approval sidecar](../simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json).
  Runtime activation remains gated. Wind/AGB/SNII source admission, physical
  SED/PISN/dust choices, live production feedback, HDF5 restart, distributed
  AMR scaling and publication convergence require their own evidence.

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
