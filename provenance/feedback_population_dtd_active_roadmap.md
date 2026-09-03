# Active population, DTD, and terminal-channel roadmap

Status: mandatory feedback work, 2026-09-02.

This file prevents the P0.4 fail-closed guards from being mistaken for a scope
reduction.  The guards only prohibit incomplete physics from entering a
production run.  They do not close or defer the following work to the generic
long-term validation backlog.

## Gate F-P1: stellar population and fate contract

Implementation status (2026-09-02): the first Claude Opus 5 audit returned
BLOCK with D1--D5.  Their engineering remediations are now complete as one
review bundle: continuous Chabrier normalization, production residual-ejecta
semantics, mandatory timestep ledger execution, strict production table
audit, and transactional legacy switches.  The consolidated re-audit returned
an engineering conditional pass for a code-mass tolerance floor and three
residual-inconsistent linked adapters; both conditions are now independently
fixed and await the next substantial bundled audit.  Scientific model/source
selection and the 40--120 Msun fate gap remain blocked and are not claimed as
a PASS.  See
`fp1_population_fate_contract_2026-09-02.md`.

- Select the single-star and/or binary population-synthesis basis and pin its
  version, license, checksum, IMF convention, initial-mass normalization,
  binary fraction and binary-parameter distributions.
- Define metallicity and age domains, stellar lifetimes, failed explosions,
  fallback, compact remnants, and mutually exclusive ownership of wind, AGB,
  SNII, SNIa, PPISN, and PISN returns.
- Distinguish per-star, per-event, instantaneous-rate, cumulative, and
  population-integrated source units.  Reject a second IMF convolution of an
  already population-integrated table.
- Require analytic normalization tests, channel-boundary tests, independent
  JAX/Fortran differentials, and a Claude Opus 5 gate audit.

## 중기 해결 계획: 40–120 M☉ terminal fate seam

상태 (2026-09-02): **중기 과제로 등록, production gate 차단 유지**.
이 구간은 하나의 ZAMS 질량범위로 direct collapse를 선언하지 않는다. 현재
문헌과 Fable의 `CONDITIONAL PASS`에 따라 source-node/최종 핵 구조 기반의
model-dependent resolver를 개발한다. 상세 근거와 이번 적용 기록은
[`fp1_mass40_120_literature_dossier_2026-09-02.md`](fp1_mass40_120_literature_dossier_2026-09-02.md)와
[`fp1_mass40_120_application_record_2026-09-02.md`](fp1_mass40_120_application_record_2026-09-02.md)에
있다.

1. **후보 source 선정:** Sukhbold W18/N20을 태양금속도 engine 비교용으로
   유지하고, Limongi–Chieffi/Ugolini 계열을 저금속도 production 후보로
   평가한다. 모든 후보는 source version/checksum, 금속도 정의, 회전,
   질량손실, explosion engine, lifetime, license를 함께 고정한다.
2. **fate resolver 구현:** `(source, M_ZAMS, Z, rotation/binary state,
   engine branch, pre-SN structure, PPISN criterion)` 키를 사용하고,
   source hull 밖 조회·nearest-node 대체·질량/금속도/회전 보간을 거부한다.
   실패붕괴, 봉투 방출을 동반한 direct collapse, fallback, PPISN, 완전
   붕괴·PISN을 서로 다른 outcome으로 표현한다.
3. **feedback 배선:** pre-terminal wind, terminal ejecta, remnant를 분리하고
   누적 age history와 terminal step을 중복 없이 ledger에 연결한다. PISN은
   완전파괴이므로 remnant owner가 없고, PPISN은 별도 pulse history를
   요구한다.
4. **runtime/검증:** 기존 선형 mass interpolation이 폭발·실패 노드를
   섞지 않도록 piecewise-constant source-cell 경로를 구현하고, 경계 포함
   규칙, lifetime timestep 해상도, unresolved mass bucket, restart/AMR/MPI
   closure를 검증한다. 현재 source-cell 선택·IMF-weighted diagnostic
   bucket·4개 계약 SHA256 admission gate와 native population 검증까지
   구현되었지만, production-linked build evidence는 아직 별도 차단이다.
5. **승인 후 전환:** age-resolved wind, decay inventory, energy/momentum/
   deposition, terminal remnant closure와 licensed source package가 모두
   충족될 때만 canonical rows를 생성한다. map checksum, sidecar, approval
   id가 일치하지 않으면 production을 열지 않는다.

**중기 종료 조건:** 40–120 M☉ source hull 전체에 대해 각 node의 fate,
wind/ejecta/remnant/energy/momentum/age/decay가 독립 재현되고, F-P1 및
F-P3 감사가 승인되며, 현재 `review_only_unresolved` 정책이 승인된 map
식별자로 교체될 때까지 이 과제는 열린 상태다.

## Gate F-P2: SNIa DTD and binary event model

Implementation status (2026-09-02): the interval-integrated mathematical
kernel, expected-event ledger, DTD/event-yield literature candidate matrices,
fail-closed contract, and source-normalized event-yield converter are
implemented and tested; the physical event model remains blocked.
The contract audit now checks the candidate matrix and dossier references in
addition to the source mirrors. See
`fp2_snia_dtd_contract_2026-09-02.md` and
`simulation/snrt/data/fp2_snia_dtd_contract_audit.json`.

The Keegans review package now has a source-semantic completeness audit in
`simulation/snrt/data/fp2_snia_keegans_format_audit.json`.  H, He, C, and N are
recorded as absent isotope rows, not explicit zero yields; the conversion
policy forbids inferring zeros.  This closes an audit ambiguity but does not
close the physical event-source gate: isotope aggregation, ejecta/energy/
momentum closure, and population weighting remain unapproved.

A second review-only source path is now staged for the HESMA
`yysd4-xap92`/Seitenzahl model record.  Its audit covers the 15 model archives,
384-column isotope profiles, and all 11 integrated project elements.  It also
reports profile-mass consistency warnings for `n300c` and `n1600c`; these are
not corrected silently.  HESMA remains an unselected candidate until a model,
decay convention, energy/momentum contract, and population weighting are
approved.  The review-only HESMA source adapter now requires an explicit model
id and produces only a checksum-bound intermediate document; it cannot emit a
converter input or runtime asset.  This advances source normalization without
silently selecting `n100` or inferring event energetics.
The checksum/admission sidecar and its audit are now also wired into the F-P2
contract; a checksum-consistent review extraction still remains blocked until
the physical event contract is approved.
An all-15-model HESMA comparison matrix is also generated and contract-audited;
it records no automatic ranking or population mixture, so model selection
remains an explicit physics decision.
The shell-boundary sensitivity audit also shows that the `n300c` mass
discrepancy persists under both documented edge policies; no silent profile
correction is permitted.

The 2026-09-03 blocker-reconciliation bundle hardens the review boundary:
HESMA models with any unresolved physical warning are rejected by the source
adapter and admission sidecar, `n300c` remains explicitly quarantined,
promotion fields are mirrored through one canonical 19-field list, all
review-audit paths are repository-relative, and the contract runner derives
its project root from its own location.  The negative admission cases and the
full F-P2 runner pass locally.  These changes do not approve a source or
enable SNIa at runtime.

The DTD kernel now evaluates the power-law primitive in a cancellation-safe
form near `alpha = -1`, with a logarithmic-series reference test at
`alpha = -0.99999999`; the native and production-mirror tests pass.

The next implementation bundle is the physical SNIa source contract: an
explicit WD-reservoir debit/closure policy, signed event-frame momentum and
its cell-deposition convention, and a versioned source/approval commit
binding.  Until those fields are populated
and independently reviewed, F-P2 remains blocked.

- Select and cite the DTD family, minimum/maximum delay, normalization per unit
  initial SSP mass, progenitor/binary assumptions, metallicity dependence,
  ejecta model, event energy, and stochastic or expectation-value realization.
- Implement interval-integrated events over `[age_old, age_new]`; never sample
  only the endpoint rate.  The cumulative event count must telescope under
  timestep subdivision and across restart.
- Test the analytic DTD integral, a single event, variable timesteps, zero-event
  intervals, restart/retry idempotence, and mass/species/energy closure.
- Keep SNIa disabled until the implementation and Claude Opus 5 audit pass.

## Gate F-P3: PISN/PPISN population eligibility

- Select a stellar-evolution/fate source and map metallicity, population tag,
  helium/core mass, mass loss, PPISN, direct collapse, and terminal yields.
- Do not treat a universal 140--260 Msun ZAMS interval as sufficient and do not
  let the generic yield interpolator create PISN events.
- Test both eligible and ineligible populations and all fate boundaries.  The
  reviewed science configuration may explicitly disable PISN, but the decision
  and its evidence are required.
- Require an independent Claude Opus 5 gate audit.

## Gate F-P4: integrated channel realization

- Combine approved winds, AGB, SNII, SNIa, PPISN/PISN decisions, remnants,
  energy, momentum, and delayed-cooling ownership without overlap.
- Demonstrate timestep/refinement/restart/MPI invariance and source-to-cell
  closure in the production RAMSES binary.
- Compare against the stopped transitional-feedback baseline only after the new
  model's own ledgers close.

No gate above can be closed by a synthetic fixture, a parser-only test, or a
silent default.  Each gate ends with a distinct Claude Opus 5 physics and code
audit before work advances.
