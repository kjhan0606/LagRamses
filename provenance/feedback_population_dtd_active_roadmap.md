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

### 2026-09-03 AGY/Opus 통합 실행 게이트

고질량 seam은 더 이상 하나의 포괄적인 중기 항목으로만 추적하지 않고,
다음 순서의 독립 승격 게이트로 수행한다. 이 순서는 두 감사의 공통
필수자료를 기존 F-P1 계획에 합치고, Opus가 발견한 runtime/schema 배선
위험을 source 승인보다 먼저 제거한다.

1. **F-P1H-A — build-bound admission identity:** production binary가
   namelist에 적힌 임의의 64자 문자열을 승인 digest로 받아들이지 않도록,
   fate-map SHA256과 approval id를 컴파일된 승인 identity와 일치시킨다.
   Review build의 컴파일 identity는 비어 있어야 하며 항상 fail closed다.
   **상태: 구현 및 단위/F-P1 계약 테스트 통과.**
2. **F-P1H-B — lossless source-node schema:** canonical 32-field payload와
   별도로 resolver의 source/Z/rotation/binary/engine/lifetime/criterion 축을
   보존하는 immutable node sidecar를 정의한다. 축을 줄일 경우 freeze 또는
   population marginalization의 분포·가중치·승인 id를 기록하고, 누락과
   물리적 zero를 구분한다. Failed/direct-collapse 모델은 누락 row가 아니라
   명시적인 wind-only/zero-terminal/remnant record여야 한다.
   **상태: 84개 필수 필드, 12개 resolver 축, zero/null 및 축 축약 규칙을
   계약·변환기·asset audit에 구현하고 회귀검사 통과. 실제 physical node는
   의도적으로 0개다.**
3. **F-P1H-C — channel and deposition ownership:** SNII terminal candidate
   domain이 40 M☉에서 잘려 60--120 M☉ 폭발 node를 버리지 않도록 fate
   resolver와 channel window의 역할을 분리한다. Wind, terminal ejecta,
   remnant, PPISN pulse, PISN disruption의 owner를 배타적으로 정하고,
   scalar radial momentum은 별도 deposition contract 없이 canonical vector에
   넣지 않는다.
   **상태: channel 3을 `[8,120]` candidate domain으로 확장하고 에너지·운동량·
   deposition 소유권 및 exactly-once 요구사항을 계약화했다. 실제 source-node
   fate/deposition consumer는 아직 없으므로 driver가 40 M☉ 초과 경로를
   거부하며 runtime deposition은 false다. 물리 node 기반 exactly-once 실행
   검증은 F-P1H-E 이후 승격 조건이다.**
4. **F-P1H-D — coverage and closure:** `[40,120]` mass cells의 gap/overlap,
   모든 resolver 축의 source hull과 ±epsilon 경계, branch별 wind consistency,
   decay epoch/duplicate isotope, source 정밀도별 mass closure, age telescoping,
   exactly-once terminal energy, restart 및 population closure를 자동화한다.
   **상태: 내부 review gate 완료. Branch별 실제 candidate node를 120 M☉까지
   보존하고 flattened-union 보간을 금지했다. W18/N20 18개 outcome, source
   정밀도 질량수지, wind branch 차이, radioactive epoch, K-40 duplicate,
   failed-node completeness를 자동 검사한다. Physical node가 없으므로 age/
   energy/restart closure의 실제 자료 실행은 F-P1H-E 이후 승격 조건이다.**
5. **F-P1H-E — physical package admission:** 교정되고 재배포 가능한
   multi-Z/multi-rotation package를 승인한다. Fate/pre-SN structure,
   lifetime/age-resolved wind, terminal/remnant/fallback, decay, energy kind와
   injected-energy mapping, momentum/deposition, PPISN/PISN 자료가 없는 node는
   승인하지 않는다. Boccioli--Roberti failed-model Wind anomaly는 저자 확인
   또는 수정 release 전까지 차단한다; Sukhbold W18/N20은 validation branch다.
   **상태: 9개 필수 gate와 4개 후보 qualification matrix를 가진 기계식
   admission contract를 구현하고 최상위 sidecar에 checksum-bound artifact로
   연결했다. 현재 후보 4개 모두 미승인, physical node 0개, production/
   publication/runtime deposition은 차단 상태다. 해시된 임의 validator 파일의
   자기승인을 막는 code registry와 Boccioli--Roberti 2026 source identity/rights
   validator의 초안을 구현했다. 이전 독립 감사에서 확인된 candidate substitution,
   self-consistent package rewrite, mutable rights evidence, sidecar
   path/publication invariant 위험은 2026-09-04 admission-closure 묶음에서
   code-owned evidence lock, 순수 coupling/selection predicate, controlled-error
   및 adversarial fixture로 보강했다. LC18 failed-wind cross-check에는 phase
   invariant diagnostics, differential G2 check, signed/relative residual과
   non-authoritative CDS rights 표시를 추가했다. 구현 commit은
   `033799a2d2ea8618877596122f02a2007d8d64bb`이다. 집중/F-P1/G2 증거와 248개
   config/data fixture 불변성 검사는 통과했고, AGY는 PASS, Opus 5는
   CONDITIONAL PASS를 반환했지만, Opus 5가 selection hash의
   검증 fingerprint 결속과 CDS-derived publication 경계를 추가 gating item으로
   지적했다. 이 두 항목과 성공 control의 zero-CDS 통계는 다음 묶음 후보로
   기록하며, 현재 묶음은 physical source 승인이나 runtime 승인이 아니다.
   독립감사와 재현/triage는 완료되었고, 다음 구현 묶음은 사용자 명시
   승인이 있을 때까지 시작하지 않는다.
   나머지 8개 executable validator, 실제 교정 자료, physical node, runtime
   consumer가 남았으며 production 승인은 계속 닫혀 있다.**
6. **F-P1H-F — promotion and bundled audit:** 물리 node를 채운 뒤에만
   converter, map, source package, sidecar, compiled identity를 다시 hash하고
   unresolved bucket이 정확히 0인지 검증한다. 전체 F-P1/F-P3 묶음 테스트와
   독립 감사를 통과하기 전에는 canonical conversion과 runtime deposition을
   열지 않는다.

상세 공통 요구와 모델별 추가사항은
`fp1_high_mass_required_data_comparison_2026-09-03.md`를 기준으로 한다.

## Gate F-P2: SNIa DTD and binary event model

Implementation status (2026-09-03): the interval-integrated mathematical
kernel, expected-event ledger, DTD/event-yield literature candidate matrices,
fail-closed contract, and source-normalized event-yield converter are
implemented and tested.  The Maoz field DTD plus HESMA yysd4-xap92/n100
physical baseline is now approved, while runtime activation remains gated on
the AMR/MPI caller.
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

The review-only SNIa physical contract now has a guarded population-ledger
path: an explicitly supplied WD reservoir is debited transactionally, mass
closure is recomputed, and source-frame, isotropic-zero, or radial momentum
conventions are validated.  A radial budget requires an explicit unit cell
direction.  The new
`stellar_snia_cell_deposition` adapter converts a validated event budget into
one-cell mass, momentum, event-energy, bulk-kinetic-energy, and total-energy
density increments.  Its only admitted thermal policy is explicitly approved
all-to-total-energy; fractional thermalization is rejected until a declared
non-thermal receiver exists.  The adapter is native/production hash-matched
and unit-tested.  The production `stellar_ramses_bridge` now wraps it with
explicit RAMSES code-unit conversion and normalized multi-cell weighting;
pre-write validation keeps failed SNIa policies transactional.  The bridge
is now callable from the runtime-facing AMR leaf-cell path.  The runtime uses
the particle's located AMR leaf cell as a one-cell NGP target and passes the
local RAMSES owner rank to the row-major `unew` bridge.  SNIa activation is
still disabled by the independent production gate, so no SNIa event is sent
to AMR cells in a production run yet.

The follow-up implementation bundle (2026-09-03) adds the actual
RAMSES-facing `unew(local_cell,variable)` adapter: it validates local target
indices, MPI owner rank, target uniqueness, and weighted multi-cell geometry,
deposits through variable-major scratch storage, and scatters only after the
scratch transaction succeeds.  The production bridge and the actual
`stellar_ramses_runtime`/`feedback.kjhan3` objects now compile together under
the `/gpfs` Makefile; the runtime now selects the particle's AMR leaf cell
and calls the adapter behind the still-closed production gate.  The
source-identity manifest now includes all SNIa objects actually linked into
the production binary, and the historical Fable reproduction test records
F3/F4/F7/F8 as resolved in the current tree rather than replaying them as
active defects.

The runtime-caller bundle is now implemented and source/build-linked.  Its
normal retry/restart invariant, weighted multi-cell bridge conservation, and
runtime-disabled negative path are now explicit in automated evidence and pass
against the linked `/gpfs` binary.  The production audit distinguishes
`physical_baseline_ready` from actual `production_ready`; the latter remains
false while activation is disabled.  The declared staging commit is an
ancestor of the current HEAD, but the current dirty worktree is explicitly
reported as not source-bound until committed.  A hard crash exactly-once
guarantee is not claimed because RAMSES `unew` writes and the existing
`indtab` checkpoint commit are not one atomic transaction; a persistent
pending-event journal is a separate production-hardening item.  Full
net-yield and metallicity sensitivity remain separate publication extensions.

The next F-P1 staging bundle records the Huscher et al. 2025 AGB release as a
review-only candidate for the `[0.8, 1.0)` M☉ lifetime seam.  The endpoint
models and source fingerprints are present, but the release is
lifetime-integrated rather than an age-resolved per-star history.  The new
`fp1_low_mass_seam_review.json` therefore keeps the seam unresolved and emits
no canonical source row.  This is evidence hardening, not a physical source
selection or a runtime activation.

- Select and cite the DTD family, minimum/maximum delay, normalization per unit
  initial SSP mass, progenitor/binary assumptions, metallicity dependence,
  ejecta model, event energy, and stochastic or expectation-value realization.
- Implement interval-integrated events over `[age_old, age_new]`; never sample
  only the endpoint rate.  The cumulative event count must telescope under
  timestep subdivision and across restart.
- Test the analytic DTD integral, a single event, variable timesteps, zero-event
  intervals, restart/retry idempotence, and mass/species/energy closure.
- Keep SNIa disabled until the caller qualification evidence and Claude Opus 5
  audit pass.

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
silent default.  Audits are no longer run after every individual step.  Each
implementation bundle ends with independent AGY (`gemini-3.8-flash-high`) and
Claude Opus 5 audits.  Their findings are reproduced and triaged before the
driver writes the next bundle plan; Fable then evaluates that plan for final-
purpose alignment, scientific/technical justification, and feasibility before
the next bundle starts.
