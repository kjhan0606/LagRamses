# High-level RT/feedback/dust provenance index

Project: `/gpfs/kjhan/chimes-cvode-gpu-20261001` (`kjhan0606/LagRamses`)
Index status: active engineering records; historical records remain in place
until an explicitly approved archive operation.
Record statuses: `active` = current contract/evidence, `superseded` = retained
but not a current gate, `long_term_record` = future G5/G6/science work.

Auditor policy (2026-11-05): internal Astra replaces every pending and future
Fable audit assignment, including plan/design reviews, scientific/technical
consultations, bundle-end/closeout audits, and fallback reviews. Completed
Fable links below remain historical records, including unavailable or
conditional outcomes; their verdicts are not rewritten and they do not assign
new work. Claude Opus 5 retains its separately assigned code
implementation/structure audits. See [project instructions](../AGENTS.md) and the
[audit cadence amendment](audit_cadence_amendment_2026-09-05.md).

This index is deliberately bundle-level. Individual native smokes and Python
validators are not separate project gates; they are referenced by the bundle
evidence that runs them. Audit prompts are request artifacts, not independent
scientific results. Scratch copies under `~/.claude/plans` are not project
evidence.

## Active bundle map

[M5 wiring and bounded coupling evidence](m5_l8_interpolation_boundary_evidence_2026-09-16.md):
source inspection records the M_N material-tile adapter connection. CPU
integration evidence for the current commit is summarized in
[real-source integration progress](real_source_integration_progress_2026-09-07.md):
AGN-only source attribution and effective stellar-population profiles both
completed with M5, live dust and CHIMES under reference-control inputs. The
effective-population run exercised M5 order 5 on levels 2 and 3; the AGN-only
run exercised the isolated AGN point source. These are bounded single-rank
CPU executions, not production qualification. An H100 CUDA-hybrid AGN-only
M5/live-dust/CHIMES reference integration also passed on source
`b1fc63c`; it verifies actual device RHS calls, not production assets, full L8,
cosmological, restart, multi-rank or performance qualification. An A40
CUDA-hybrid effective-population AMR-interface run on the same source also
passed: stellar births and a nonzero M5 ledger were observed across levels 2
and 3, with live dust and hybrid CHIMES. Its mixed M5 ledger does not
independently attribute photons to AGN; all source assets remain reference
controls. Details, limitations and hashes are in the same progress record. The
[M1--M5 dispatch evidence](mn_order_dispatch_efficiency_2026-09-27.md) adds a
source-active two-step coupled M5 run and bounded backend comparisons; none of
these establishes full L8 or cosmological production qualification. The earlier
`snrt_m5_full_connection_2026-09-15.md` entry referenced below is absent from
this checkout and must not be treated as evidence. Existing physical-asset
admission blockers remain in force.

[CHIMES dark mode-3 charge projection](chimes_dark_mode3_projection_2026-11-05.md):
bounded H-/C- roundoff repair with strict charge/element closure, native and
captured-input production-bridge evidence, plus a current-source CUDA-hybrid
M5/live-dust/CHIMES candidate run. Syntax job 451493 reached step 4 on both
continuous and resumed paths but failed the unchanged 626-dataset restart
comparison: 519 datasets bitwise-identical, with seven out-of-tolerance values
confined to five CHIMES species and two persistent SNRT He III fractions (the
payload header, not M5 angular moments). Its 22/31 CVODE warning counts are
not mapped to individual cells or retries. Raw HDF5 dumps were checksummed
and removed. The candidate remains nonproduction; the chemistry/ionization
restart divergence and warning attribution remain unresolved. Same-allocation
A/B job 451539 captured the prior cell and the exact comparator outlier in turn
on one A10; all 46 normalized input-manifest entries matched. Both strict
comparisons passed with zero tolerance violations. The prior-cell capture had
local CHIMES API differences within the integrated tolerance; the outlier-cell
capture was bitwise identical at steps 2 and 3. Per Astra's decision rule,
pass/pass does not implicate the capture target as the cause of job 451493's
failure, which remains unreproduced and unexplained. Raw HDF5 outputs were
checksummed then removed. Details and hashes are recorded in the bundle report.
Context replay job 451130 and integrated A40 replay job 451146 did not start;
Slurm rejected both with `QOSMaxGRESPerUser`.
The 2026-11-06 source-level fiber/retry review found no concrete thread-local
state leak: cooling is consumed inline after continuation resume, adapter
state is restored per fiber, the CHIMES rate cache is invalidated at switches,
and whole-cell CPU retry reinitializes from original inputs. Job 451493 remains
unexplained; no speculative patch or extra runtime matrix is justified. The
hybrid candidate remains nonproduction and CPU remains the configured default.

[Periodic sink antipode tie convention](sink_periodic_antipode_tie_bundle_2026-11-05.md):
the one-ULP negative-image rule is implemented and passed GNU/Intel native
fixtures plus CPU MPI2/OMP4 integration. The fixed-order sink-stat reduction
passes a native permutation test and CPU/CUDA source-object compilation.
Matched CPU replays 1271793 and 1271805 both passed the unchanged strict
626-dataset continuous/restart comparison with zero tolerance violations
(623 datasets bitwise-identical; three roundoff-level sink differences).
Together with the matched control failure, this supports sink-stat
accumulation as the cause of the tested CPU restart mismatch. The earlier
1271786 attempt stopped at provenance parsing before evolution. CUDA and
production admission remain unqualified; no production promotion is claimed.

[Event-weighted AGN SED assessment](agn_physical_sed_bundle_plan_2026-10-31.md):
Astra's plan review permits bounded offline work but keeps runtime admission
closed until an AGN-on source-state domain and out-of-group energy treatment
are justified. On 2026-11-06, the operator explicitly chose to keep production
AGN SED work blocked until a target AGN-on run is designated. The P4 Run0
snapshot remains diagnostic only. The CHIMES restart-context 2×2 replay and
integrated A40 replay were rejected at Slurm admission as jobs 451130 and
451146 (`FAILED`, `ExitCode=2:0`, `QOSMaxGRESPerUser`, 4 s and 29 s elapsed);
neither produced scientific evidence. Their 2026-11-06 terminal accounting
supersedes the earlier pending estimate; they were not resubmitted after the
same-allocation A/B and source review provided no new testable hypothesis. The
later integrated restart run is recorded separately under the
CHIMES dark mode-3 evidence. Current source inspection confirms that `mkrun.py`
rejects the non-admissible 25-group profile, while RAMSES permits only the
explicitly opted-in, collectively validated diagnostic fixture; diagnostic
permission does not change the production runtime-admission flag. CPU-only
grammar job 1271987 passed the existing native spectral-contract smoke for the
candidate: runtime admission stayed false, diagnostic permission stayed closed
without the full fixture gate, and the incompatible legacy grey H/He path was
rejected. It did not run RAMSES or promote the profile; details/hashes are in
the linked plan and `.snrt-admission-gate-20261106/`.
The later CPU-only MPI2 source-receipt integration, job 1272009, passed the
rank-local source-class/aggregate injection check and combined stellar plus
thermal-AGN attribution for one bounded M5 transaction. It produced no raw
dumps and leaves production AGN SED admission blocked; see the
[`M5 source attribution review`](m5_stellar_agn_source_attribution_2026-10-31.md).
Astra's 2026-11-06 next-bundle triage answered Q-GOAL: the receipt work advances
the effective-model objective, but current evidence identifies no further
source-independent code gap to bundle; Q-LEAN: do not add another runtime
matrix or instrumentation. Keep the strict CHIMES restart limitation and
AGN/Fe physical blockers unchanged. This is a sequencing disposition, not
project completion; resume the AGN SED bundle only after a target AGN-on run is
designated and its source domain is justified.

Historical current-source A40 integration smoke (2026-10-21,
Slurm job430823): a fresh binary built from source HEAD
`f459f6dc1febdaf1eb5c9ae25c97924be553e188` plus the recorded dirty working
tree (`DUST_LIVE=1`, `CHIMES=1`, `USE_CUDA=1`) completed four coupled order-5
M_N steps. The effective run exercised dust-live material/IR and
CHIMES photo/dark chemistry paths on 56 coarse and 64 fine leaf cells, with
CUDA material batches on one A40; `Run completed`, zero raw snapshots and the
`M5_COUPLED_AMR_INTERFACE_PASS` marker are recorded under
`/gpfs/kjhan/chimes-cvode-gpu-20261001/.m5-effective-amr-430823/` and the
binary/build hashes under `.m5-current-source-active-dust-chimes-build-430823/`.
This was a one-rank, four-step reference-control test
(CHIMES CPU backend, no restart comparison), not cosmological or production
qualification; physical-asset admission blockers are unchanged.

The current source does contain a compile-time `DUST_LIVE=1` RAMSES field,
HDF5 identity/restore contract, and M_N material-coupling path. The queued
`run_snrt_bundle_gate.sbatch` does not enable that profile (`DUST_LIVE` defaults
to 0), so its result is native/default-profile evidence only, not live-dust
runtime qualification.

[128^3 M5 runtime diagnosis](snrt_m5_cost_diagnosis_2026-09-14.md):
separate short kernel/closure cost measurement; original shadow run unchanged,
not a solver modification or a new physical gate.

[Density-selective ionization / shielding pilot](snrt_self_shielding_2026-09-13.md):
native Gaussian/lognormal density and luminosity comparisons; fixed temperature,
not a cooling-time experiment or angular-converged reference.

[CUDA/Tensor Core RT comparison and angular artifacts](snrt_cuda_and_angular_2026-09-13.md):
standalone S_N/M_N GPU backends, explicit mixed precision, and quadrature-rotation
diagnosis; not coupled RAMSES accelerator admission or angular-converged physics.

[Angular-method ionization comparisons](snrt_ionization_comparison_2026-09-13.md):
native fixed-grid evolving H tests and 128^3 reference execution; distinct from
coupled RAMSES production integration.

[M5 solver safeguards](snrt_m5_solver_2026-09-13.md):
covariance-metric step constraint with unchanged physical moments; native
near-boundary and directed-beam regressions, not coupled production approval.

[Selectable tensor M_N implementation](snrt_tensor_mn_implementation_2026-09-13.md):
9/16/25/36-moment native closure, conservative predictor/corrector and compact
state tests; historical standalone stage. Its blanket live-startup rejection
has been superseded by the full live connection above.

[SNRT runtime performance and hybrid measurements](snrt_performance_2026-09-12.md):
CPU/GPU IR improvements, integrated cost breakdown, and Tensor Core limitations;
stream-count tuning stopped per operator instruction.

[CHIMES photo/M5 high-thread tail](chimes_photo_omp_tail_2026-09-28.md):
bounded two-tile OpenMP preparation queue, 16/31-core A/B timing and physical
parity; speculative CVODE RHS skip rejected after integrated measurement.
The 2026-10-30 AGN-only CPU reference integration adds a three-step source and
ledger pass plus a measured cooling/coupling-dominated runtime; see
[CHIMES parallel efficiency](chimes_parallel_efficiency_2026-09-27.md). This
remains NONPRODUCTION reference-control evidence, not a whole-run speedup or
global energy-conservation qualification.

[CHIMES MPI/OpenMP layout comparison](mpi_omp_decomposition_2026-09-28.md):
fixed-total-core 1x28/2x14/4x7 comparison on 512- and 4096-cell controls;
4x7 was fastest on both. At 4096 cells, active-step times were
273.12/203.74/170.00 s. Per-cell photo-chemistry cost remains; production
scaling is not established by these single-run controls.

[CHIMES dark audit synthesis and P0/P1 implementation](chimes_dark_bundle_2026-09-12.md):
deterministic species-rate initialization, disabled explicit-candidate removal,
zero-photon-only work reduction, native parity and one-step integration evidence.

[CHIMES solver cost and optimized SUNDIALS](chimes_solver_cost_2026-09-12.md):
unoptimized numerical-library build identified; separate precise-FP64 release
build, native cost attribution and same-node integrated comparison.

[Exact-temperature dark rate reuse](chimes_rate_cache_2026-09-12.md):
same-cell temperature-only coefficient reuse; composition-dependent rates
and all illuminated chemistry remain recomputed.

[Three-step repeated-state verification](chimes_multistep_2026-09-12.md):
job542054 completed3 steps/96 CMB commits/3 IR commits without rejection;
initial-state validation only, no source-active galaxy claim.

Sequential effective-model implementation:
[v5/empirical-SNIa connection and integrated execution](effective_population_execution_2026-09-11.md).
Step1 native/live/restart passes; step2 configuration selected; step3 short
coupled execution/restart passes. The recorded long-macrostep CHIMES failure
is now [reproduced and repaired](chimes_long_interval_repair_2026-09-11.md):
unchanged~68.65Myr intervals complete the four-step coupled run, with the
original conservation thresholds. It was not deferred as microscopic physics.
[COLIBRE-type calibration and resolution convergence](galaxy_calibration_science_plan.md)
is a separate science campaign, not an implementation completion gate.
The operator now defers pc-scale zoom to long-term work and proceeds with the
small-box pilot. The previous25cMpc23--29-run design and its pilot matrix are
superseded. Report actual gas/star/DM mass and spatial resolution.
Local zoom results are not volume-complete galaxy-population calibration.
[New small-box hydro IC](galaxy_new_ic_2026-09-11.md):12.5cMpc,z99,128^3,
separate gas/DM, CAMB+2LPT/LLA; file checks pass. RAMSES reader/evolution and
halo-selected refined zoom is deferred. CPU costs include setup failures.
The scoped cosmological CMB/material repair and bounded MPI test passed;
the [actual128^3 coupled startup pilot](galaxy_coupled_pilot_2026-09-12.md)
is submitted as job540015. This is not galaxy calibration or source-active
performance qualification; its resource/startup result is tracked separately.
No science simulation has been launched.
[Approved initial Z=1e-10](galaxy_pilot_preenrichment_2026-09-11.md): new explicit
GRAFIC element fields prepared; actual yield/SED Z brackets verified. Native
reader and full cosmological coupled-model admission remain pending.

Latest completion standard (operator-approved2026-09-11):
[galaxy-scale effective physical model](general_physics_completion_plan_2026-09-11.md).
Declared literature-based approximations are allowed; conservation and native
coupling remain required. [Direct microscopic calculations](microscopic_physics_long_term_backlog.md)
are long-term research, not additional completion gates. This scope approval
does not itself constitute implementation or production qualification.

Approved remaining bundles1--5 (2026-09-11):
[execution and physical-data boundaries](remaining_physics_bundle_execution_2026-09-11.md).
Kind7 condensation/shock/sublimation/relative connection is
[complete within its bounded evaluated scope](kind7_dust_process_implementation_2026-09-11.md);
native tests, coadvected/relative live execution and bitwise restart pass.
Completed raw outputs were removed, with compact evidence retained. Do not infer
completion of PAH hard-source survival or a full matched stellar population.

Approved comparison-model execution (2026-09-11):
[PAH catalytic H2 complete within its stated domain](pah_catalytic_comparison_implementation_2026-09-11.md),
including shared H donor, thermal energy receipt and MPI2/OMP2 restart parity
for7,565 datasets. Four evaluated raw directories removed; inputs, binary,
hashes and compact results retained.
[Wind bi-stability velocity comparison](parsec_bistability_comparison_2026-09-11.md)
has physical-grid and native source/SED binding evidence; it does not add LBV
mass ejections. [Selective four-plan review and driver disposition](remaining_comparison_fable_2026-09-11.md)
records the comparison designs. [Mixed/truncated low-mass source](parsec_mixed_lowmass_evidence_2026-09-11.md),
[frozen-contact effective SNIa](snia_frozen_he_hybrid_plan_2026-09-11.md) and
[Al26/Fe60 transparent-decay material](stellar_radioactive_implementation_2026-09-11.md)
now have bounded actual MPI/live/restart evidence, with evaluated raw cleaned.
[Fe photon limits](fe_photon_comparisons_plan_2026-09-11.md) also pass bounded
MPI2/OMP2 live comparison. [PAH single-photon atomization](pah_single_photon_atomization_plan_2026-09-11.md)
passes actual weak/strong-source comparison, carbon return and7591-dataset
restart matching. These complete the approved named-comparison scope, not a
universal broadband source or self-consistent full stellar population.

Kind7 actual-source multilevel AMR repair (2026-09-11):
[completed empty-leaf and MPI ownership repair](snrt_multilevel_amr_fix_2026-09-11.md).
Consistent empty receiver batches are successful no-ops without skipping
collectives. Face classification uses current-grid ownership, not prospective
child-oct ownership. MPI2/OMP2 fully covered and mixed/derefining cases pass;
element/charge errors below1e-15, global mass2.147e-16, maximum IR4.544e-10.
This closes the recorded multilevel failure below within the tested scope.

Actual stellar/AGN sources with kind7 evolving grains (2026-09-11):
[completed bounded connection and driver evaluation](snrt_actual_source_transition_2026-09-11.md#final-driver-disposition).
Checked177-field AGN material map, narrow existing-sink/NENER0 admission and
both generators updated. Bondi empty-level/canonical-migration and coarse
donor-search defects repaired; trace-metal photo accuracy and negative-tail
retry fixed without relaxing budgets. MPI2/OMP2 actual source/thermal-to-jet
integration and restart pass: all hydro/chemistry/dust/RT/gravity/AGN queues
bitwise equal (441/449 datasets; eight near-zero velocity/momentum differences).
Native258 PASS; frontend48 PASS/1 skip. This record's uniform reference scope
is unchanged; its multilevel failure is resolved by the bounded repair above.
Excluded models are not declared production-ready.

SNRT high-temperature / evolving C-silicate grains (2026-09-11):
[completed connection and Fable disposition](snrt_transition_dust_evolution_2026-09-11.md).
Kind7 optionally admits existing growth/sputtering/coagulation/shattering,
with total-H collision density fixed across CHIMES mass paths. Current-bin
opacity/material/chemistry coupling passes dense cold/hot MPI2/OMP2 evolution
and two exact394-dataset restarts. Native254 PASS; exhausted-photon ODE and
bounded FP32-tail storage repaired. Existing fixed defaults and explicit
Fe/PAH/condensation/shock/drift/sublimation exclusions remain.

SNRT energy-aware high-temperature continuation (2026-09-11):
[approved implementation and live evaluation](snrt_hot_transition_plan_2026-09-11.md#approved-implementation-and-live-evaluation).
Explicit kind7 rapid-dissociation/atomic-NEQ comparison on 10--1e9 K;
ATcT binding costs, accepted-solution thermal root and ABI6. Native 251
checks pass; mixed-temperature and irradiated MPI2/OMP2 restarts each match
402 datasets bitwise. The subsequent mass-evolution connection is above;
fixed-grain results remain valid as that explicit comparison. First-order
topology assumptions remain; this is not unlimited molecular kinetics.

SNRT cold molecular / competing D03 spectra (2026-09-11):
[implementation and live completion](snrt_chimes_molecular_coupling_implementation_2026-09-11.md#live-completion).
Actual multi-cell scattering/photo/dark/material wiring, explicit kind6
data/restart identity and frontend opt-in completed. MPI2/OMP2 dark and
irradiated restarts match 402 datasets bitwise; neutral four-step run passes.
Bounded 10--10^4.98 K fixed-grain comparison, not general all-temperature
molecular/PAH/Fe/dynamic-dust admission. This supersedes the native-only
status of this connection, not other models' declared physical limits.

SNRT H/He intragroup spectrum (2026-09-10): [native implementation/evidence](snrt_band_implementation_2026-09-10.md).
Opt-in 64-node maximum-entropy N/E closure, evolving Verner absorption and
actual energy into chemistry; MPI2/restart exact. Default fixed retained;
not dust/CHIMES/GPU spectral qualification or completion of all medium physics.

SNRT angular refinement (2026-09-10): [live implementation/evidence](snrt_angular_implementation_2026-09-10.md).
80/320/720 choices, unchanged 80-ray default, refined MPI2 hydro/dust/feedback
restart and CUDA transport tests. The bounded H/He spectral extension is above;
full dust/network spectral coupling remains open.

Medium-term physics (2026-09-10): [scope, first Fe implementation and status](medium_physics_implementation_2026-09-10.md).
Eight groups approved for implementation; the first feature is optional
metallic-Fe seed growth plus Fe-specific thermal erosion, not completion of
the entire medium-term program. [Sources](medium_fe_sources_2026-09-10.md).
Continuation resolves the cold PAH/Fe 5 K bracket with native material/Planck
physics and four-step MPI2/OMP2 integration plus exact restart comparison.
Fixed-H PAH charge/photoelectron/IR and CHIMES gas+solid neutrality now connect
through 256 hydro carriers and exact MPI2/OMP2 restart, with native energy,
charge and element evidence. The optional H0--13 PAH comparison now adds
finite-H attachment/loss, binding energy and state-dependent mass through
3584 hydro carriers. Native/live conservation and exact MPI2/OMP2 restart
pass; the enlarged H-state physical identity uses checked HDF5 dataset I/O.
Carbon-skeleton destruction, H2 chemistry, higher charge states, general
hard-source survival and the other numbered extensions remain incomplete.
Group 6 now has a [native collision-aware size comparison and decision](dust_collision_comparison_2026-09-10.md):
coagulation/fragmentation, mass-ledger boundaries and 16/32/64-bin results.
Keep the live two-size default; no claim of multibin equivalence or a new
live multibin transport model. PAH carbon-fragmentation [source gaps](pah_fragmentation_sources_2026-09-10.md)
remain unresolved; no artificial atomization or inert-C22 closure was added.

Gas ideal MHD (2026-09-10): [implementation and native evidence](gas_mhd_implementation_2026-09-10.md).
CT/HLLD gas evolution, mixed AMR and same-size HDF5 restart verified;
dust-field and CR-field interactions explicitly deferred to long-term work.

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
| F-P2.7 gate consolidation / initialized-RAMSES | implemented: D2/D3 PASS; D4 PASS for recorded startup binaries; CPU-vs-CUDA binary-identity caveat retained | [plan](fp2_7_gate_consolidation_initialized_ramses_bundle_plan_2026-09-05.md) | [implementation evidence](fp2_7_gate_consolidation_initialized_ramses_bundle_implementation_evidence_2026-09-05.md) | [Fable plan](fable_fp2_7_gate_consolidation_initialized_ramses_bundle_plan_audit_2026-09-05.md) | [Codex D4 audit](codex_fp2_7_initialized_ramses_bundle_end_audit_2026-09-06.md) covers initialized runtime only; no physical-source or production approval |

## Next bundle queue

| Bundle | Status / purpose | Plan | Plan audit | Authorization |
|---|---|---|---|---|
| F-P1H-F physical-source admission / high-mass seam | approved reduced scope; KL16/CK22 investigation parked by operator; physical admission remains unapproved, not an engineering-wide blocker | [plan and wrap-up](fp1h_f_physical_source_admission_bundle_plan_2026-09-05.md) | [Fable](fable_fp1h_f_physical_source_admission_bundle_plan_audit_2026-09-05.md), historical conditional approve | operator wrap-up 2026-09-05; do not automatically reopen AGB source research |

## Cross-cutting decisions and non-approvals

- SNIa scope decision (2026-09-07): operator approved an explicitly effective
  SSP/remaining-stellar-mass accounting mode with unchanged empirical DTD and
  N100 source values; **implemented with bounded native/live/restart checks**. Preserve
  strict WD-supply mode and its existing rejection. Microscopic binary
  evolution/IFMR/retention work is medium-term, not the next implementation
  prerequisite. See the final operator-decision section of the
  [existing completion plan](production_completion_bundles_2026-09-07.md).
- KL16 auxiliary mass discrepancy: operator selected detailed yield-header .774 Msun over .744 (2026-09-07); printed article data take precedence where available. Other KL16/CK22 issues remain parked, not a full physical-package approval.
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
- [audit cadence amendment](audit_cadence_amendment_2026-09-05.md) governs
  selective Astra review: request a plan audit for a new physical model,
  substantive conservation/coupling change, unresolved scientific question,
  or explicit operator request. The driver performs end evaluation; there is
  no routine end audit or automatic fallback. Every requested plan audit asks
  **Q-GOAL** first, then **Q-LEAN**; this adds no audit or gate.
- [Fe edge quadrature comparison](fe_edge_quadrature_comparison_2026-11-06.json)
  verifies the current 128-node opacity averages against an edge-split
  2048-node reference for the checked group N/E means; maximum relative
  difference is 1.116%. It does not validate photoelectron heating or admit
  high-energy Fe at runtime.
- [Astra review of the Fe hard-photon response plan](astra_fe_hard_photon_response_plan_audit_2026-11-06.md)
  defers the proposed sub-edge extension: AGN25 groups 8–9 cover 25.827462% of
  `L_grid`, while occupied groups 10–25 cover 33.930488%. Keep hard-photon Fe
  fail-closed until an effective charged-grain electron/energy-partition model
  covers the intended source domain.
- Historical [Fable operational audit](fable_operational_instrumentation_gate_efficiency_audit_2026-09-05.md)
  classified the workspace `OVERINSTRUMENTED`; F-P2.7 is the approved
  consolidation response.
- SNIa HESMA N100 and Maoz DTD retain their prior physical approval
  `FP2-SNIA-PHYSICAL-2026-09-03-N100-MAOZ`; see
  [approval sidecar](../simulation/snrt/config/fp2_snia_event_source_approval_sidecar_v1.json).
  Runtime activation remains gated. Wind/AGB/SNII source admission, physical
  SED/PISN/dust choices, live production feedback, HDF5 restart, distributed
  AMR scaling and publication convergence require their own evidence.

## Classification rules

Latest M5 wiring evidence: [current-source wiring reconciliation and L8 boundary](m5_l8_interpolation_boundary_evidence_2026-09-16.md), supplemented by the
[2026-09-27 dispatch/coupled-run record](mn_order_dispatch_efficiency_2026-09-27.md)
and the current-source CPU, H100 CUDA-hybrid AGN-only, and A40 CUDA-hybrid
effective-population AMR source/dust/CHIMES executions described above. These
establish bounded source coupling, not full L8 or cosmological production
readiness. Full hybrid restart/multi-rank and
production-input qualification remain pending. The missing
`snrt_m5_full_connection_2026-09-15.md` artifact is not a valid current
evidence target. Earlier intermediate reports remain historical.

1. A bundle evidence file is the authoritative entry point for its native
   smokes, build/link, hashes, and conclusions.
2. A focused runner may be used for debugging, but its output is not a new
   gate or audit event.
3. Historical audit reports remain immutable records; later dispositions are
   appended to the bundle evidence/index rather than creating repeated audits.
4. Large HDF5, JAX, virtual-environment, build, and generated compiler
   products remain outside this index's active evidence unless a bundle
   explicitly pins a manifest/hash for them.
