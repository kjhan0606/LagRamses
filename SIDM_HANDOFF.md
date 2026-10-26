# SIDM 프로젝트 업무인수인계서

작성 시각: 2026-08-24 KST

## 2026-10-21 최신 lagRamses 기준 소스 전환 및 쌍성 수치 게이트

사용자 지시에 따라 GitHub `kjhan0606/LagRamses`의 최신 기본 브랜치 `main`
커밋 `9bec08cd5e3fa570d881a86da067b8fcfbb6e37d`를 별도 worktree
`/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/latest_main`
에 받았다. 현재 작업 브랜치는 `codex/sidm-cloud-ledger`이며 기존 `sidm`
브랜치(01b6ead)와 그 미커밋 변경은 별도 유지한다. 앞으로 1번 연구의 코드
기준은 이 최신 main이고, 과거 브랜치의 검증된 안전성 수정은 필요한 것만
개별적으로 이식·재검증한다. 양쪽을 자동으로 병합했다고 주장하지 않는다.

LagEunha에서 최신 코드의 비-HDF5 Intel MPI/FFTW 초도 빌드를 완료하고 sink
2개의 짧은 중력 진화를 수동 실행했다. 이 초도 smoke run에서는 소스 수정·
커밋·푸시를 하지 않았다. 첫
HDF5=1, SNRT=0 빌드는 링크 누락으로 실패했으며 별도 기록했다. 첫 적응격자
시험의 심한 격자 토글은 `m_refine`을 한 단계에만 준 입력 오류였다. 단계별
입력을 교정한 뒤에는 토글이 사라졌으나 쌍성 궤도량과 질량중심 오차가 아직
크다. 균일 64^3 및 128^3 대조군은 각각 약 3초와 10초 벽시계로 두 덤프를
만들고 두 BH 질량과 궤도 운동을 확인했다. 128^3 덤프는 약 125 MB씩이다.
코드 메모리 로그는 64^3에서 1.7 GB, 128^3에서 368 MB로 단조적이지 않아
별도 peak-RSS 계측이 필요하다. RAMSES `econs` 진단과 sink-cloud 중력의
수렴은 아직 검증되지 않아 실제 SIDM+쌍성 연구의
생산 실행을 승인하지 않는다. 전체 입력·빌드·로그·수치 비교와 출력 시계의
함정은
`/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/LATEST_MAIN_GATE.md`
에 기록했다. 현재 `/goal`은 사용자 확정 목표로 active이며 1번에 집중한다.

## 2026-10-21 phase-resolved PM-cloud gate 재감사

목표 1의 현 단계는 live binary+SIDM halo 계산이 아니라, 그 전에 sink-cloud
알고리즘의 수치적 hardening/energy work를 분리하는 것이다. 위의 최신-main
전환 항목은 초기 checkout 시점의 기록이다. 현재 기준 checkout은
`/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/latest_main`,
HEAD `9bec08cd5e3fa570d881a86da067b8fcfbb6e37d`, branch
`codex/sidm-cloud-ledger`다. 별도 build worktree에는 reset-boundary ledger를
기록하는 sink-source 변경이 있으며, 그 binary는 준비된 2-step runner에서만
사용하도록 고정되어 있다. 기존 `sidm` worktree와 자동 병합하지 않는다.

보존된 오프라인 결과
`pm_cloud_reset_energy_uniform128_v2.json`은 SIDM-off, halo 없는 uniform
128^3 두-cloud 상태에서 가상 same-time reset map을 적용한다. 진단상
`delta_K_reset=-1.9539839348e-5`, force-squared PM 진단 변화
`delta_U_F=+8.0448242693e-5`, 합은 `+6.0908403345e-5`다. 절댓값 Kepler
proxy 대비 2.98은 단위·규모 비교일 뿐 orbital-energy 주입률이 아니다.
이 결과는 실제 kick/drift/다음 reset을 함께 닫은 시간스텝 예산이 아니며
쌍성 hardening 또는 SIDM 물리 결과로 해석하지 않는다.

Astra의 독립 재감사 결론은 **계측 설계 조건부 승인, 현재 실행파일은 full-step
검증에 불충분**이다. 기존 writer는 `reset_pre / canonical_average /
reset_post` 경계 inventory만 기록하고 force-source epoch를 계측하지 않는다.
기존 초기·최종 dump도 실제 내부 gather, `synchro_fine` half-kick,
`move_fine`의 sink-cloud velocity overwrite 및 후속 reset을 독립 관측하지
않아 이를 실행된 closed-step ledger로 대체할 수 없다. 따라서 기존
`run_cloudledger_uniform128_nstep2_20261021/submit.sbatch`는 energy gate로
제출하지 않는다.

다음 기술 단계는 reset inventory를 유지하고 별도 opt-in operation ledger를
추가해 force-ready/source epoch, sync kick, 실제 projection 및 move kick,
drift, reset endpoint를 phase·level/subcycle별로 기록하는 것이다. 고정 질량
조건에서 endpoint `delta K`가 sync/projection/move-kick/reset 항의 합과
일치하는지, drift가 kinetic energy를 바꾸지 않는지 검증한다. source-matched
PM force-squared 진단은 계속 Hamiltonian이 아닌 별도 수치 진단으로 표기한다.
Astra가 권고한 다음 pilot은 새 binary로 fresh uniform-128^3, `nstepmax=2`
조건의 1-rank 및 2-rank 실행을 각각 해 비교하는 것이다. 초기 seed→cloud와
terminal partial interval은 정상 closed interval 판정에서 제외한다.

준비된 현 boundary-only bundle의 read-only 사전검사는 통과했다: `bash -n`,
`checksums.sha256` 전체, live-gate namelist contract. 이번 재감사에서는 추가
소스 수정·빌드·실행·Slurm 제출은 하지 않았다. 특히 기존 runner는 2-rank
한 조건뿐이므로 operation ledger 구현 뒤 1/2-rank 두 조건에 맞춰 다시
감사·준비해야 한다. 프로젝트의
최신-main `AGENTS.md`는 코드 구현 역할을 Claude Opus 5로 지정하고 있다.
현재 runner의 제출에도 별도 명시 승인이 필요하다. 이 gate는 PM-cloud 수치
결합만 검사하며 SIDM 열수송/경화 가설을 검정하지 않는다.

## 2026-10-21 finite-crossing 및 cloud-ledger 단위 회귀 확인

승인된 finite-time crossing 진단 코드와 테스트 파일은 대상 검증 디렉터리에서
이미 Astra 검토본과 일치했다. 중복 수정 없이 단일 프로세스로 테스트를 실행해
`14 passed in 1.12s`를 확인했다. 코드 SHA256은
`92239db1668fb712415a5d3a38b8d5aba7d6ff4a8c8864f5e40e3ce135de07e2`, 테스트는
`ab676a42bcb0ec4cc4cb42d51757b63ec15a160c58b24063bd786a799d4d3039`다. 기존
v3 대규모 진단 JSON은 `432d12934a4d073f30a378ecd4756bc15acc3edb61a97d9aa2b31fcb5eebf52b`
로 전후 동일하다. 이는 finite-time crossing 진단의 회귀 확인이며 RAMSES
실행이나 새 대규모 진단이 아니다.

최신-main `tests/sidm/test_cloud_ledger.py`도 단일 프로세스 unittest로 재확인해
`19 tests, 2.343s, OK`였다. 이는 reset-boundary cloud inventory 검증기의 fixture,
rank/event/namelist 계약 및 실패 검출을 확인한다. 소스가 `FORCE_SOURCE_EPOCH
not_instrumented`로 기록하는 것처럼 실제 gather, half-kick, cloud projection,
drift를 잇는 full-step energy closure를 증명하지 않는다. 이번 두 테스트 실행 중
소스 수정, 빌드, 시뮬레이션, Slurm 제출은 없었다.

## 2026-10-21 Astra operation-ledger 설계 재검토

Astra는 opt-in, operation-resolved ledger를 **설계 조건부 승인**했다. 이는
현재 PM-cloud 쌍성 적분의 수치 선행조건이며 SIDM 열수송 가설을 직접 시험하지
않는다. 최신-main call order는 `create_sink → rho_fine → force_fine →
synchro_fine → move_fine`이다. `move_fine`에서 sink 속도는 기존 입자 속도가
아니라 `vsink + 0.5*dtnew*ff`로 투영되므로 이를 별도 kinetic-energy 항으로
기록해야 한다. Force-source epoch, level/subcycle, 실제 gather 가속도, 입력·출력
속도, 실제 `dteff/dtnew`, drift 위치를 연산 순번과 함께 남긴다. 재생성 경계를
넘는 입자 ID로 rank-local slot을 오인하거나 replicated sink 상태를 중복 합산하지
않는다.

최소 닫힘식은 고정 질량 구간에서
`ΔK_endpoints = ΔK_sync + ΔK_projection + ΔK_move_kick + ΔK_reset`이며 drift의
kinetic-energy 변화는 0이어야 한다. 첫 1-rank pilot의 모든 rank/event 완결성,
실제 gather의 독립 reference 일치, 정규화 잔차의 roundoff 규모(예비 목표
`1e-10`)를 확인한 뒤 동일 입력 2-rank에 간다. 진단 off/on은 물리 상태가 같아야
한다. nstepmax=2는 plumbing 검증일 뿐 장기 궤도 수렴이 아니다. 다음 matched-time
mesh/Δt 비교의 orbit/work 오차가 사전 1–2% 예산을 만족하지 못하거나 reset/
projection 효과가 해상도에 따른 바닥값으로 남으면 target-scale SIDM science run을
보류한다. PM의 `U_F ∝ -Σg²`는 Hamiltonian이 아니므로 그 변화의 잔차만으로 pass/fail
하지 않는다. 구현 및 새 run 승인은 별도 역할·실행 절차를 따른다.

## 2026-10-21 Astra 3개월 실행 가능성 및 규모 재검토

Astra의 Q-GOAL/Q-LEAN 검토 결론은 **다음 작은 pilot은 타당하지만 생산 규모는
아직 추천 불가**다. 현재 64³/128³ 벽시계와 dump 크기는 SIDM 없는 두-cloud
실행이고, 128³는 `ngridmax`의 99.9%를 차지한다. RAMSES의 보고 메모리도
64³ 1.7 GB, 128³ 367.6 MB, corrected adaptive 약 6.1 GB로 비일관적이며,
adaptive Kepler proxy 변화 약 −3.9%의 원인도 아직 분리되지 않았다. 실패한
Slurm pilot은 MPI 초기화에서 멈춰 RAMSES halo 실행시간·메모리·물리량을 측정하지
못했다. 따라서 이 자료에서 백만 입자 SIDM halo의 비용, 필요한 N, 완료시간을
외삽하지 않는다.

Astra가 제시한 3개월 분기는 다음 순서다: (1) operation ledger를 구현·독립 검증한
뒤 uniform 128³, SIDM-off, `nstepmax=2`를 1-rank, 이어 같은 입력 2-rank로 실행;
(2) matched-time 64³/128³ 및 Δt orbit/work 비교, 필요할 때만 256³를 추가하고
목표 cloud/separation 비율의 moving-AMR를 검사; (3) 준비된 halo의 collisionless
안정성과 SIDM rate·열응답·peak 자원을 낮은 규모에서 측정; (4) 이 gate와 정밀도가
통과할 때만 binary/monopole × SIDM-off/active A/active B의 6조건을 짧은 구간부터
실행한다. `2^18→2^20`은 측정용 탐색 계단 예시일 뿐 생산 추천 N이 아니다.
`a_h/3`는 잠정 판정 구간이며 `a_h/10`은 3개월 go/no-go의 필수조건으로 두지 않는다.

확대는 matched orbit/work 오차가 예산 약 1–2% 이내이고, collisionless drift,
SIDM sampled-pair/radial shell energy, seed 변동, 경계 및 heat-response 보정이
분리될 때만 한다. 열응답 시간은 binary-free perturbation에서 독립 정의하며
`P_bin≈L_χ` 또는 동일 장부 잔차에서 정의하지 않는다. 잠정 20% 효과를 주장하려면
해당 대비의 총 불확실도가 약 6.7% 이하여야 하고, 수치 편향에는 별도 예산을 둔다.
저정밀도나 미닫힘은 물리 가설의 반증이 아니라 feasibility no-go다. 실행비용은
노드/job peak memory와 rank RSS, wall time, grid/particle peak, output bytes로
직접 잰다. ledger plumbing pilot은 full dump 0개, halo 비용 pilot은 제한된
초기/최종 dump로 실제 크기를 잰다.

이번 세션의 기존 2-rank boundary-inventory bundle은 shell syntax, manifest 내
9개 checksum, effective namelist contract 사전검사를 통과했다. 그러나 이는 실행
준비성일 뿐 operation-level gate가 아니며 디렉터리에 실행 산출물도 없다. 코드
구현, MPI probe/Slurm 실행, 물리 비교는 수행하지 않았다.

## 2026-10-21 사용자 확정 연구 목표와 현재 우선순위

사용자는 Astra의 상위 연구 제안 1·2·3번을 모두 연구 계획에 포함하고,
현재는 1번에 주력하도록 지시했다. 새 시뮬레이션 수행도 허용되어 있다.
다음 문단은 `/goal`에 넣을 목표 문안이다. 기존 JCAP 방법론 논문의 완료
기준을 충족했다고 간주하거나 현재 원고의 주제를 자동으로 교체하지 않는다.

> SIDM에서 독창적이고 검증 가능한 물리 결과를 얻어 높은 수준의 물리학
> 저널 논문으로 발전시킨다. 연구 계획에는 (1) 살아 있는 초대질량 블랙홀
> 쌍성과 SIDM 헤일로의 시간 의존 열수송 병목, (2) 내부상태 반응 지연이
> 만드는 동적 열역학과 중력계 반응, (3) 회전하는 SIDM 붕괴핵의 각운동량
> 장벽을 모두 포함한다. 현재 주력은 1번이다. 기존 연구가 제시한 유한
> 열저장고의 에너지 예산을 충족하는 계에서, 유한한 열수송률과 살아 있는
> 쌍성·헤일로의 상호작용이 쌍성 궤도 경화의 지속 또는 정체를 결정하는지
> 판정한다. 대조군과 보존량·수송량 진단을 갖춘 작은 시험부터 시작하고,
> 실제 측정한 실행시간·저장량과 수치 수렴에 따라 생산 계산 규모를 정한다.
> 첫 3개월에는 1번의 구현 가능성과 핵심 물리 가설을 검증 가능한 근거로
> 평가한다. 2번은 소규모 반응계 시험, 3번은 회전 평형계와 각운동량 수송의
> 설계를 후속 단계로 유지한다. 기존 JCAP 방법론 논문은 검증 범위에 맞춰
> 계속 완성한다.

1번의 첫 판정은 선행연구의 누적 energy-budget 부족과 별개인 수송률 효과를
검출할 수 있어야 한다. 동일한 초기 에너지 예산에서 finite halo와 외곽
thermostat를 비교하고, CDM 쌍성·SIDM 단일 블랙홀·정적인 배경을 대조한다.
궤도 경화율, 중심 spike의 진화, 궤도-헤일로 에너지 교환 및 바깥쪽 열유량을
동시에 기록한다. 블랙홀 궤도 적분, 평형 초기조건, 중력·산란·시간보폭,
Brownian motion, 경계 조건을 검증하고 영향반경의 실입자수로 해상도를
판정한다. 탐색 규모는 제안서의 수치일 뿐 확정된 생산 사양이 아니다.
첫 달의 작은 검증과 계측 결과로 3개월 범위와 go/no-go 기준을 갱신한다.
2번과 3번의 상세 대조군·단계별 규모·중단 기준은 아래 Astra 제안서에 있다.

현재 `/goal` 시스템에는 위 1·2·3번과 기존 JCAP 방법론 논문을 함께 담은
목표가 active로 설정되어 있다. 아래 초기 문안과 실제 `/goal`은 동일한
연구 범위를 지향하되 표현이 완전히 같지는 않을 수 있다.

## 2026-10-21 Astra 물리 종합저널 관점의 독립 연구 제안

사용자가 현재 원고/진행방향에 얽매이지 않는 새 SIDM 물리연구를 Astra에게
제안받도록 요청했고, 새 시뮬레이션의 수행도 명시적으로 허용했다. 새
gpt-6-astra agent는 기존 원고·이전 심사평을 읽지 않고 최신 primary 문헌을
조사했다. 결과는 아래 지정 경로의 research_proposals.txt에 보존했다.
/gpfs/kjhan/Hydro/Sidm/validation/astra_physics_editor_20261021.FLymij
request.txt가 의뢰문이다. 실제 저널 편집자의 결정이나 게재 보장이 아니다.

후보 8개 중 우선순위는 (1) 살아 있는 SMBH 쌍성의 시간 의존 열수송 병목,
(2) 내부상태 반응지연의 동적 열용량/벌크점성/중력계 진화,
(3) 회전하는 SIDM 붕괴핵의 각운동량 장벽이다. 나머지는 angular closure의
한계, halo 이력의 최소 상태변수, 물리적 분열질량, 초기/후기 비탄성 evolution의
일관성, bar resonance response다. 후보별 선행연구·차별성·반증조건과 상위
3개 pilot/해상도/stop-go/3개월 및 장기경로를 담았다. 측정하지 않은 실행시간은
발명하지 않았고 저장량 예시는 입자 record 크기만의 산술로 한정했다.

1번에서 유한 열저장고의 energy-budget 조건 자체는 기존 PRL 연구에 있다.
그 조건을 이미 만족하는 baseline에서 시간 의존 수송과 live backreaction이
추가로 hardening을 막는지 판정하는 것이 새 질문이다. 2번은 detailed balance와
총 에너지 보존을 전제로 하며 수동 반응에서 자발 성장모드를 가정하지 않는다.
3번은 Newtonian 수축 endpoint를 BH formation으로 주장하지 않는다.

이 제안서 작성 과정에서는 해당 물리주제의 시뮬레이션을 시작하지 않았다.
이후 사용자가 상단의 1·2·3번을 연구 계획으로 채택하고 1번을 현 주력으로
확정했다. 새 시뮬레이션 허용은 유효하다. 기존 methods 논문의 검증 과제는
미완료다. 이 절 작성 당시에는 `/goal`의 기존 항목이 paused였지만, 현재는
상단 목표가 active이며 별도 두-sink 수치 게이트의 진행 상황은 이 문서의
최신 항목에 기록했다.

## 2026-10-21 실제 산란 연산자의 균질 이완 검증

Goal continuation에서 현재 worktree를 다시 확인했다. 코드 HEAD는 01b6ead이고
paper HEAD는 d31a077이며 기존 미커밋 변경은 보존했다. 현재 squeue에 SIDM
계산은 없다. 기존 L17 output_00005는 CDM a=0.573638632210679,
SIDM3 a=0.250336146250867로 여전히 epoch가 다르다. 이를 저적색편이 matched
science result로 취급하지 않았다. 이전 goal turn의 독립 심사는 progress였다.

새 tests/sidm/homogeneous_relaxation.f90 및 run_homogeneous_relaxation.py는
production sub_sidm_scatter 전체를 무변경 추출하고 실제 probability/angular
helpers와 amr/random.f90를 링크한다. 균질 periodic leaf와 free streaming은
통제된 environment가 제공한다. Maxwell-molecule형 sigma*v=constant 조건의
독립 moment 방정식과 quadrupole/energy 이완을 비교한다. 그 조건은 수치 검증용
모델이며 실제 dark-matter 모델의 추천이 아니다. 1 km/s floor의 moment 보정
상한을 별도 기록한다. 실제 AMR tree, gravity, MPI 또는 halo 검증은 아니다.

LagEunha foreground 수동, nice19, 각 실험 1 thread로 실행했다. 17 조건 ×
8 production RNG substreams의 matrix 결과는
/gpfs/kjhan/Hydro/Sidm/validation/homogeneous_matrix_20261021b/verification.json
이며 SHA256 a566725d7d9cfee9a8a27cfb792741886bbe94d18b8f42777e1b48ce2e32767c다.
같은 root의 interpretation.txt와 PDF/PNG에 derivation, scope 및 결과를 기록했다.
8192 입자로 8/64/512 cells, dt=0.04/0.02/0.01을 대조했고 2048/32768 입자,
grid translation, Rutherford epsilon=0.1/0.01의 equal-total/equal-viscosity,
f_diss=0.2도 측정했다. 모든 ensemble mean이 초기 stress 대비 절대 0.02의
사전 기준 안에 들었다. 최대 discrepancy 0.01650765, Pmax 0.0925,
momentum residual 1.38e-16, kinetic+removed energy ledger 1.29e-14 이하다.
Floor의 moment 보정 상한은 최대 7.8e-7이다. 평균의 통계 오차를 함께 보존했고
이 결과로 일반 AMR convergence 또는 시간 적분 차수를 인증하지 않는다.

homogeneous_lowN_128_20261021b는 mean occupancy=4 조건을 미리 정한 128개
substream ensemble로 확인했다. 최대 mean discrepancy 0.00375310,
최대 standard error 0.00262339다. 첫 8 streams는 matrix와 겹치는 nested
sample이며 독립 재현이라고 부르지 않는다. 동일한 tolerance를 유지했다.
homogeneous_negative_20261021b는 별도 build에서 Rutherford를 isotropic
sampler로 의도적으로 바꿨다. conservation은 유지되지만 이완 discrepancy가
0.48123으로 기준을 실패하여 동역학 오류를 검출했다. 생산 소스는 수정하지 않았다.

초기 a-series는 새 driver의 seed 설정 오류 때문에 INVALIDATED로 표시했다.
큰 integer를 base-4096 digit에 넣어 일부 run의 정수 overflow 및 128-run 시도
timeout을 일으켰다. 정식 rans(256,0,seed_bank) 부분수열과 digit 범위 검사를
사용하도록 driver를 수정하고 b-series를 새 경로에 전부 재실행했다. 256은
serial seed-bank의 dimension이며 MPI ranks가 아니다. 기존 RNG 소스는 그대로다.
실패 자료는 삭제하지 않았다. 원고·기존 production binary·출력·Overleaf는
변경하지 않았다. 실제 AMR/중력 halo/인증된 matched low-z zoom은 계속 필요하다.

## 2026-10-15 Astra 천문학 저널 독립 심사 완료

사용자가 천문학 탑 저널 referee 관점의 평가를 요청했다. 새 gpt-6-astra agent가
이전 심사평을 읽지 않고 현재 원고 전체를 독립 평가했다. 원고와 PDF를
`/gpfs/kjhan/Hydro/Sidm/validation/astra_astro_referee_20261015.B7MO95`에 고정했고,
의뢰문은 request.txt, 한국어 상세 심사평은 referee_report.txt에 보존했다.
심사평 SHA256은 22933aac4a1d64d16b8a60128aa1cd8950c93a9a6242246233298edc5bfba2ed다.
원고/PDF SHA는 아래 운전자 polishing 완료 기록과 동일하다.

공식 저널 결정이 아닌 AI referee-style 판정은 Nature Astronomy급 현 상태
비추천, MNRAS/ApJ/A&A 전문 방법론 논문 substantial major revision이다.
보고서에 major 6개, minor 10개, 원고 위치, primary 문헌/공식 정책 링크,
외부 manifest와 원고 증거의 구분, 필수/선택 검증 및 12주 연구 경로를 포함했다.
핵심 보완은 AMR/continuum collision 수렴, 중력과 결합한 고립 halo benchmark,
중심 기여의 차별화, 유지할 dissipative 기능의 dynamical accuracy,
timestep/full restart 및 공개 재현 패키지다. ADM은 열역학 closure 재설계 또는
experimental 범위로 분리하도록 권고했다. 이미 수정한 coefficient/guard를
미해결 코드 오류로 취급하지 않았다. 새 대형 zoom은 필수조건으로 요구하지 않는다.

이번 요청은 평가이며 연구 방향 변경의 실행 승인이 아니다. 심사평을 근거로
원고나 코드를 수정하지 않았고, 새 시뮬레이션·commit·push도 수행하지 않았다.

## 2026-10-14 사용자 지시에 따른 운전자 문장 교정 완료

사용자가 "Fable로 하지말고 자네가 수행해"라고 지시하여 Codex가 추가·수정된
논문 문장의 polishing을 직접 수행했다. 아래의 Fable 대기 기록은 이전 시점의
기록이며 이번 지시와 완료 결과로 대체된다. Fable에는 의뢰하지 않았다.
`/home/kjhan/WRITING.md`에 따라 초록, 검증표, iSIDM/ADM 설명과 결론의
문장 구조와 용어를 다듬었다. 적분 결과 및 재결합 reference의 차이가 상대량임을
명시했다. 과학적 의미 변경은 없으며 기존의 미검증 물리와 적용 제한을 보존했다.

증거 경로는 `/gpfs/kjhan/Hydro/Sidm/validation/paper_polish_20261014.7OXDHv`다.
작업 전 원고/PDF, 작업 후 원고, polishing.diff, checks.json, 컴파일 로그와
PDF 표·본문 확인 이미지를 보관했다. 자동 비교에서 수식·inline math·수치의
출현 횟수·인용·label·참조·그림 경로가 보존되었고 수정 행의 금지 대명사는 0개다.
ADM trajectory 표/driver 및 zoom profile JSON의 SHA도 작업 전과 동일하다.
Lageunha에서 수동 latexmk 빌드를 완료했고 paper/main.pdf를 갱신했다.
최종 PDF는 21쪽이며 undefined reference/citation, overfull box, oversized float는
없다. 기존 hyperref/font 경고와 runtime 표의 underfull 3개는 그대로다.

최종 main.tex SHA256은
9f6f00d50da662a708cf914a0caee2b2a603c01e641df14f2c46b8fa1e63f4cc,
main.pdf SHA256은
3f07bdaedc0b760ea07f811b5e3e28c9d1c92eb6b85220a453ed9d6a1a77ebd0다.
과학 검증의 추가 완료나 투고 준비 완료를 뜻하지 않는다. 코드와 시뮬레이션은
변경하지 않았고 commit/push 및 Overleaf push는 수행하지 않았다.

## 2026-10-14 승인된 P0 결과의 논문 갱신

사용자가 논문 갱신을 명시적으로 지시했다. `/home/kjhan/WRITING.md`를 읽고
기존 paper repo의 미커밋 변경을 보존하여 사실·수치와 검증 범위를 반영했다.
현재 paper HEAD는 d31a0779c9db71e943f8d42933251c686d4c8e36이며 commit/push,
GitHub push 및 Overleaf push는 하지 않았다. 아래 P0 절의 paper 미변경 상태는
이전 작업 시점의 기록이다. 이번 갱신은 코드나 production을 추가 변경하지 않았다.

원고·PDF·기존 ADM 표·driver의 작업 전 사본은
`/gpfs/kjhan/Hydro/Sidm/validation/paper_p0_update_20261014.RmT7ga/before`에
보관했다. 새 원고/PDF/표는 같은 root의 after에 보관했고 `verification.json`과
`source_change_final.diff`에 변경 범위와 근거를 기록했다.

`paper/main.tex`의 초록, 구현 설명, 검증표, 결론 및 runtime 표에 동일 질량
iSIDM 제한, all-ground restart 보존, HPM의 hydro=F/star=F/sink=F 조건,
HPM force의 synchronization 전 ordering을 반영했다. ADM 제동복사 질량
스케일링과 재결합의 kinetic-energy coefficient/cgs 식을 명시했다.
냉각 적분기를 기존의 ten-iteration 설명에서 실제 80-iteration bisection과
1e-10 bracket tolerance로 교정했다. ADM의 nuclei-only energy와 Saha closure,
나머지 원자 냉각항·chemistry·PdV·shock 및 full restart의 미검증 범위를 남겼다.

ADM toy trajectory 기준은 960.77869 K, 2000-step implicit 상대 차이는
1.1796655e-7이다. 원고는 각각 960.78 K와 1.2e-7로 표기한다.
`paper/adm_timestep_convergence.txt`는 GNU/Intel P0 원시 출력 양쪽과 byte-for-byte
일치한다. Single-step 차이는 약 0.024%다. 수치 적분의 일치를 ADM 전체 물리
검증으로 확대하지 않는다. 계수 비교의 978개 point와 tail 비교의 5개 y^2
point를 구분했다. 검증표는 page overflow를 피하도록 두 표로 나눴다.

원래 zoom JSON의 20-500 physical h^-1 kpc, 양쪽 count>=50인 공통 17 bins를
재계산했다. Median absolute density differences는 SIDM1 0.81267777%, SIDM3
1.73201950%이며 dispersion은 0.11152862%, 0.59785958%다. Signed median을
오차 크기로 읽지 않도록 교정했다. Astra 요약의 첫 absolute-density 숫자를
그대로 복사하지 않고 원래 배열에서 계산했다. 원문 심사평은 수정하지 않았다.
Legacy IC/gravity provenance 미인증을 본문·표·caption·결론에 명시했다.
Drag cap의 50% 조건도 relative-speed cap이 아니라 각 component kick임을
실제 식에 맞춰 교정했다. 기존 profile figures와 simulation payload는 재계산하지 않았다.

Lageunha foreground 수동 latexmk 빌드가 통과했다. 새 PDF는 21쪽이며 undefined
reference/citation, overfull box 및 oversized float는 없다. 기존 style의 hyperref
경고, small-cap font fallback 및 runtime 표의 underfull 3개는 nonfatal이다.
물리적 PDF page 8/12/17의 검증표·냉각식·temperature paragraph를 이미지로
확인했다. `paper/main.pdf`도 새 PDF로 갱신했다. 원고 SHA는
7b7fda1775458b6e7c5b0032f41cd1855109736173d19330ad12360a6c026255,
PDF SHA는 6485fad0047a3e87e452a83c586cc0780b4d97ed5f92b1430cfedc9a54489ca9다.

현재 호출 가능한 Fable 도구/경로가 없어 Fable polishing은 미완료다.
운전자가 Fable을 가장하거나 임의의 다른 model을 대체하지 않았다.
`fable_polishing_request_final.json`에 변경 문장과 보존 조건을 준비했지만
외부로 발송하지 않았다. 추가 문장의 WRITING.md 금지 대명사 자동검사와
git diff --check는 통과했다. 투고/푸시 또는 polishing 완료로 표시하지 않는다.

## 2026-10-14 Astra 심사 후 승인된 P0 교정·소규모 검증

사용자의 승인으로 우선 코드 교정과 bounded validation을 수행했다. 원문 Astra
심사평과 기존 production/IC/checkpoint, 원래 bin executable은 보존했다.
현재 작업은 GitHub commit/push 또는 Overleaf push를 하지 않았다. paper/main.tex도
동결 심사본과 같은 SHA aa18c7ef3c748ead0d365c5890d9dab2c949e93fc44ac646191be715955eed3d다.

증거 root는 `/gpfs/kjhan/Hydro/Sidm/validation/astra_p0_20261014.FmTymC`이며
`verification_manifest.json`에 범위·source SHA·미완료 항목을 기록했다.

- 양쪽 dark_cooling_mod의 제동복사 mass scaling을 m_e'^(-3/2)로 교정하고
  덮어써지던 임시 대입을 제거했다. 재결합은 Rosenberg–Fan(2017) 식(13)의
  g_bf=1 전체 bound-level 합/적분의 kinetic-energy loss coefficient로 교체했다.
  81-node 고정 log table과 양 끝의 해석적 점근식을 사용하므로 mutable 초기화나
  매 입자마다 quadrature가 없다. m_p' >> m_e', optically thin, nonrelativistic
  적용 영역이며 이 변경이 Saha/EOS 또는 나머지 원자항을 인증하지 않는다.
- lagRamses/cuRamses read_params에서 HPM+hydro를 금지했다. 기존 star/sink
  금지를 유지한다. dark-only pressure force 분리는 구현하지 않았으므로 HPM은
  hydro=F, star=F, sink=F의 제한적 prototype이다.
- iSIDM 첫 호출에서 전체 rank의 eligible DM macro mass를 검사한다. 양의 유한
  동일 질량(64 machine-epsilon 상대 허용오차)만 허용하며, 상태/RNG 변경 및
  OpenMP region 전에 명시적으로 중단한다. 서로 다른 weights의 상태 전환을
  보존적으로 구현한 것이 아니라 미지원 입력을 거부한 것이다. Elastic SIDM
  zoom의 multi-mass 사용을 제한하지 않는다.
- nrestart>0이면 excited fraction 요청과 무관하게 checkpoint 상태를 보존한다.
  all-ground restart를 새 IC로 오인해 다시 여기시키지 않는다. 새 IC의
  pre-populated excited state를 보존하는 기존 경로도 유지했다.

모든 실행은 Lageunha foreground 수동으로 수행했다. GNU13.2/ifx2025.3에서
각 patch당 978개 coefficient point(SM+40keV+2MeV 전자질량, 두 coupling,
온도 전이 구간/점근영역)를 독립 cgs reference와 비교했다. 제동복사 최대
상대 차이는 0.03110%(유지한 SM 정규화의 반올림), 재결합은 0.38283%다.
Reference의 bound-level tail 64/256 비교 차이는 9.84e-11이었다.
`rates_gnu/verification.json`, `rates_ifx/verification.json`에 결과와 SHA가 있다.
기존 adiabatic/heating/empty-cell 단위시험도 새 source와 링크하여 통과했다.
보존한 toy one-zone closure에서 200000-step explicit reference는
960.77869 K, 2000-step implicit 상대 차이는 1.1796655e-7이다.
기존 원고의 681.05492 K는 교정 전 냉각함수 결과이므로 현 source의 검증값으로
재사용하지 않는다. 새 결과는 `rates_*/corrected_timestep.txt`에만 저장했다.

`safety_ifx_final/verification.json`은 실제 production routine/block을 compile한
33개 MPI(1/2/4 ranks) guard/init test, 양쪽 identical HPM guard의 8개 consumer
조합, 동일 질량 산란의 3개 momentum/kinetic+internal-energy case를 기록한다.
전체 checkpoint-reader 또는 AMR evolution의 end-to-end 인증은 아니다.
별도 clean source에서 CPU HDF5+FFTW 전체 build가 94.752s에 통과했다.
`build_serial/SOURCE_SHA256SUMS`, `BINARY_SHA256SUMS` 확인이 통과했다.
새 binary SHA는 3678b014b0ebfb4898950d1e2a9e62ebf14f28c2a0818450d57cc7b0d30fd6c6다.
Estimator regression, angular integrals 7개, HPM wiring도 통과했다.
실제 binary의 HPM+hydro namelist는 MPI1/2에서 init 전에 거부됐다.
`binary_hpm_reject_final/effective.nml`은 noutput=1, aout=1.1, tout=1e100,
foutput=fbackup=1e9이고 예상/실제 dump 모두 0이다. 진화 또는 IC 읽기는 안 했다.

실패한 초기 시도는 별도 경로에 보존했다. 기존 Makefile의 불완전한 module
dependency로 fresh -j2 build가 race해 -j1로 검증했다; Makefile은 수정하지 않았다.
Runtime에는 `/home/kjhan/local/hdf5/lib`가 Conda보다 먼저 와야 한다.
VPATH-selected cuRamses/update_time.f90의 clean_stop은 zero exit를 반환하므로
launch gate는 exit code만 믿지 않고 ERROR text와 진화/출력 여부를 검사해야 한다.

아직 ADM science production을 허용하는 결과가 아니다. 충돌 이온화·여기의
항별 계수 감사, Compton 정량 reference, Saha/CIE·열용량·총에너지/PdV/shock
closure·분자 bookkeeping, iSIDM 상세평형, full restart integration, AMR/halo
수렴을 남긴다. 다음 우선 작업은 나머지 원자 냉각항과 에너지 ledger 검증이며,
논문의 ADM 수치/그림 갱신과 Fable polishing은 아직 수행하지 않았다.

## 2026-10-14 Astra 독립 사전심사 의뢰

사용자가 Astra에 물리학 최상위 저널 심사자 기준의 검토를 명시적으로 요청했다.
`gpt-6-astra` 독립 reviewer task `/root/astra_referee`에 원고 전체, PDF, 그림,
최신 검증 자료와 source 경로를 제공했다. 실제 저널의 공식 심사·판정이 아니라
AI 사전심사이며 PRL의 범용 중요성과 PRD/JCAP의 전문적 출판 가능성을 구분한다.
심사는 완료됐다. 원고·코드·production·Overleaf를 변경하지 않았다.

고정 입력과 의뢰 원문은
`/gpfs/kjhan/Hydro/Sidm/validation/astra_referee_20261014.rqHwq4`에 보존했다.
`INPUT_SHA256SUMS`, `request.txt`, `review_manifest.json`으로 원고 버전과
심사 범위를 기록했다. main.tex SHA는
aa18c7ef3c748ead0d365c5890d9dab2c949e93fc44ac646191be715955eed3d다.
원문 심사평의 지정 보관 위치는 같은 root의 `referee_report.txt`다.
평가 요청은 독창성·물리 모델·수치 검증·주장별 근거·major/minor issues,
명시적 출판 판정 및 3개월 내 우선 연구/시뮬레이션과 통과 기준을 포함한다.
이전 Fable 심사평 원문을 복구하거나 재현한 결과로 표시하지 않는다.

최종 심사평은 336행/53,645 bytes이며 SHA256은
d0819dd0d8914c5d878eb094d0e81d8ac5bc28d5fbf015d8424fc70a9fa498c5다.
PRL reject, PRD/JCAP은 현 원고 reject 후 대폭 재구성한 재투고 권고다.
주요 의견12개·세부 의견10개와 측정 가능한 통과 기준/12주 계획을 포함한다.
우선 지적은 ADM 원자 냉각계수/질량 스케일링, 허용 hydro 설정의 HPM gas-force
경로, unequal macro-mass iSIDM 총에너지 보존, all-ground restart와 상세평형,
AMR transport/time convergence, legacy science provenance 및 signed-median
표현이다. 실제 과거 결과의 영향 크기까지 측정했다고 주장하지 않는다.
핵심 권고는 검증된 elastic/phenomenological-angular SIDM 중심으로 범위를 좁히고
independent isolated-halo 및 AMR benchmark를 완성하는 것이다. 저적색편이
그림만 추가하면 PRL/JCAP 기준을 통과한다는 평가가 아니다. 보고서의 범위
축소·코드 수정·새 계산 제안은 승인/시행 결과가 아니며 본 의뢰로 실행하지 않았다.
원문 보고서는 수정하지 않고 보존했다. 선택 소스8개를 audit_sources로 복사하고
SOURCE_SHA256SUMS로 source evidence를 보존했다.

## 2026-10-14 coarse transfer-kernel OpenMP 버그 교정 및 pilot 종료

이 절은 아래의 전체 OMP4 pipeline 미통과 기록을 대체한다. LagMUSIC
`paper-Ib-total-velocity` branch에 `c65da8fd23b147f8e29ed790c1acffcca83cd9fa`
(`Serialize overlapping coarse transfer-kernel writes`)를 commit/push했고,
GitHub 원격 SHA와 로컬 HEAD 일치를 확인했다. 수정은
`src/convolution_kernel.cc`의 OLD_KERNEL_SAMPLING coarse averaging loop와
신규 `tests/test_kernel_thread_replay.py`만 포함한다. 기존 미커밋 writer는
변경하거나 commit하지 않았다. master/upstream으로 merge/push하지 않았다.

원인은 40^3 fine kernel을 16^3 coarse periodic cache로 averaging할 때
ACC_RC의 periodic mapping이 같은 출력 cell로 겹치는 것이다. 병렬 ix loop의
대입/누적 ordering에 race가 있었다. 유효 kernel cells만 비교하면 기존
base-L4 kernel의 OMP1/4 차이는 relative L2 6.447e-5였고 상위 kernel 및
white noise는 roundoff 수준이었다. FFTW padding 2 cells는 미초기화일 수
있으므로 물리 kernel 검사에 넣지 않는다. OLD coarse loop는 upstream에서
이어진 코드로, 로컬 git blame에서 Oliver Hahn의 2019 commit을 확인했다.
수정은 fine extent <=2*coarse extent인 모든 축에서만 parallelize하며, 그 외에는
기존 직렬 overwrite/add 순서를 유지한다. 입력을 바꾸거나 기준 오차를 완화하지
않았다. 이전/수정 binary의 OMP1 noise/kernel/24 output fields 차이는 0이었다.

LagEunha foreground 수동 시험만 수행했다(nice19, 최대4 total threads,
build 최대2 processes). 검증 root는
`/gpfs/kjhan/Hydro/Sidm/lm_gate/kernel_fix_1013.KDkEPW`다. clean detached
f7be6f8 source + 위 2 files로 full rebuild했으며, 공통 외부 white-noise anchor를
지정해 MPI/OMP (1,1), (1,2), (1,4), (1,4) repeat, (2,2), (4,1)을 비교했다.
wrapped L4–6 24 fields의 최대 relative L2는 2.722e-8, 정상 L5–7은 1.968e-6로
기존 2e-5 gate를 통과했다. 유효 kernel의 차이는 각각 <=3.06e-15/3.95e-15다.
수정 전 f7be6f8 binary로 같은 외부-input 시험을 수행한 negative control은
88개의 kernel/field 비교에서 실패했고 input anchor는 보존됐다. 기존 RNG
회귀시험도 모두 통과했다. 증거는 final_wrapped_replay/final_normal_replay의
thread_replay.json, negative_replay/thread_replay.json, rng_replay/regression.json,
build/raw logs 및 verification.json이다. production/target/extended-reader
인증은 아니며, 기존 출력·다른 프로젝트 계산은 보존했다.

SIDM timestep pilot 421295는 COMPLETED 0:0, elapsed 01:21:27이다. 20 coarse
steps 중 warmup2 제외 18 steps는 4332.54s(240.697s/step), a=0.250382823에서
0.250736434로 진행했다. FFT direct solve21, fine solves2568, 최대 final
relative residual9.997e-5, convergence gate=true, warning/nonconvergence0이다.
운영 출력 gate도 통과했고 accepted_as_science_trajectory=false다. 이 국소
속도를 고정한 외삽은 pilot 끝에서 z=2까지11.71 compute days, z=1까지35.35일,
z=0까지106.25일이다. queue/AMR 성장/초기 단계·여러 physics/resolution variants는
포함하지 않으며, 90일 생산 완료를 보장하지 않는다. 따라서 신규 certified IC와
reader gate 뒤 우선 matched z=1 이하 결과를 얻는 축소 행렬을 검토해야 한다.
원시 pilot 결과는 validation/l17_timestep_pilot/421295/pilot_summary.json이다.

## 2026-10-13 외부 white-noise phase-anchor 수동 검증

사용자 지적에 따라 MPI rank별 PRNG 재생성은 판정 대상에서 제외하고, 하나의
저장 white-noise 실현값을 읽는 경로를 직접 점검했다. LagEunha에서만 `nice 19`,
최대 4 MPI ranks의 작은 L5 IC를 생성했다. production source, 기존 IC, 대형
production 실행, 제출 및 push는 변경하지 않았다. 입력 anchor는
`/gpfs/kjhan/Hydro/Sidm/lm_gate/manual_1013b/z5r1/wnoise_0004.bin`
(16^3 float64, SHA-256
`d7c582abdc1c529b2b032ae40fd8028c0790264fce03ed5134feeba8b8400760`)다.

- `manual_1013c_white_noise`의 1/4-rank run은 모두 `seed[4]`에 위 파일을
  지정했고, log의 “read from file” 및 compact-file read marker를 확인했다.
  각 run 뒤의 `wnoise_0004.bin` SHA-256도 rank 사이에 동일했다. 저장본이
  원본 SHA와 다른 것은 load 뒤 작업 field를 cache로 다시 쓴 결과이므로
  원본 불일치로 해석하지 않는다. 그러나 rank당 32 OpenMP threads 설정에서
  density/Poisson 뒤 IC field 차이는 최대 relative L2 `1.416022e-3`였다.
- `manual_1013d_white_noise_omp1`은 같은 binary·anchor·1/4 ranks에서
  `OMP_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`로만 바꿨다. phase-anchor 및
  rank-cache 일치는 모두 통과했고, 16 IC fields의 최대 relative L2는
  `3.1555207e-8`로 기존 `2e-5` 기준을 통과했다. 따라서 external white-noise
  read와 MPI 경로는 이 작은 fixture에서 통과하며, 앞선 차이는 rank별 난수
  생성이 아니라 높은 OpenMP thread 설정의 후속 병렬 수치 경로에서 발생한다.

증거는 각 root의 `manifest.json`, `music.log`, `white_noise_replay.json`에
있다. 재현 도구는 `tests/sidm/lagmusic_white_noise_replay.py`이며 기본 실행은
`OMP_NUM_THREADS=1`을 명시한다. 단위검사 3개를 통과했다. 이 결과는 작은 L5
fixture의 phase-anchor/MPI 재현성 검증일 뿐이며, production IC·target halo·과학
결과 또는 production release를 인증하지 않는다.

## 2026-10-13 실행·중력 수렴 재감사

### 11:22 KST LagMUSIC 난수 버그 수정·커밋·푸시

사용자가 명백한 코드 버그의 수정·commit/push를 승인했다. 독립 LagMUSIC
저장소의 현재 branch `paper-Ib-total-velocity`에
`f7be6f87268548f9723d1c9e51277542b824978b`를 커밋하고 GitHub로 정상
fast-forward push했다. `ls-remote`로 원격 branch SHA와 로컬 HEAD 일치를
확인했다. master/upstream Bitbucket에는 push/merge하지 않았다.

수정은 `src/random.cc`와 신규 regression tests 3개만 포함한다.
periodic cube 중복 fill을 제거해 OpenMP allocation/fill race를 막고,
FFT 확장 크기가 PRNG period를 넘는 unpack은 직렬 last-write 순서를 유지한다.
F.2b auto slab RNG는 이러한 alias geometry를 eligibility check에서 제외해
rank0 refinement로 안전하게 전환한다. 정상 크기의 SPMD는 유지하며,
사용자가 부적격 slab 경로를 강제하면 noise 출력 전에 명시적으로 거부한다.
기존 phase를 바꾸는 새 난수 알고리즘을 도입하거나 tolerance를 완화하지 않았다.

LagEunha에서 clean detached worktree(35a4eb9 + 정확한 수정/시험 files)를
Intel MPI로 build하고 foreground 수동 시험했다. build 최대2 processes,
MPI 최대4 ranks × OMP1, 내부 RNG boost 최대4 threads다. root는
`/gpfs/kjhan/Hydro/Sidm/lm_gate/rng_fix_1013.htc1Fm`이며 `replay_c`가
최종 성공 시험이다. wrapped L4–6 MPI 및 repeat/forced legacy path의
24 fields 최대 relative L2는 1.276e-8, normal L5–7은 4.136e-8로
기존 2e-5 기준을 통과했다. pre-fix 1-thread noise는 세 level 모두 bitwise
유지했다. 직접 C++ unpack test는 7개 geometry(축별 alias/다중 wrap/경계
crossing/exact period)에 대해 독립 cell-order reference 및 OMP1/4,
4회 반복을 통과했다. 같은 시험에 pre-fix random object를 연결한 negative
control은 `OMP periodic writeback changed`, exit1로 실제 버그를 검출했다.

별도 전체 pipeline OMP4 trial은 noise 차이 ≤5.14e-16인데 density 차이가
base에서 6.053e-5로 strict 2e-5를 넘었다. 난수 외 convolution/solver
pipeline의 원인은 아직 격리하지 않았다. 이 실패는 `regression.log`와
`replay/wrapped_omp4`에 보존했고 tolerance를 면제하거나 전체 OMP4 검증
통과로 표시하지 않는다. 직접 RNG unit test와 MPI4×OMP1 full-field replay를
구분한다. 초기 `replay_b`는 실제 explicit rejection을 caught-fatal 문자열로만
검사한 test assertion 때문에 실패했으며, nonzero exit도 인정하는 올바른
검사와 core-dump suppression을 적용한 새 `replay_c`에서 통과했다.

기존 미커밋 grafic writer는 수정·stage·commit하지 않았다(SHA
efb108e7f94df120e3d1f4554def21ed6fee68fdada51c6b04143d28ac02c936 유지).
따라서 clean HEAD writer의 44-byte GRAFIC을 포함한 standalone RNG replay이며,
확장48-byte SIDM reader/full physics gate를 새로 통과했다는 뜻은 아니다.
RAMSES 초기 dump/MG 로그 정책 수정, 생산 계산, Overleaf 작업은 하지 않았다.
`verification.json`, `replay_c/regression.json`, 모든 positive/negative raw logs
및 source/binary hashes를 위 root의 `SHA256SUMS`로 검증했다. production
release=false 및 기존 출력 보존을 유지한다.

### LagMUSIC 수동 small-gate 실패 원인 격리

사용자의 원인 파악 지시에 따라 LagEunha에서만 foreground, nice 19,
MPI 최대 4 ranks로 작은 비교시험을 수행했다. production 수정·제출·push는
하지 않았다. 진단 roots는 `/gpfs/kjhan/Hydro/Sidm/lm_gate/cause_1013a`와
`/gpfs/kjhan/Hydro/Sidm/lm_gate/cause_1013b`다. 원래 `manual_1013a/b` 및
기존 production 데이터는 변경하지 않았다. 신규 진단 도구는
`tests/sidm/lagmusic_cause_probe.py`이며 관련 unit tests 총 6개를 통과했다.

1. **MPI 차이는 작은 시험의 periodic alias/writeback 경로 차이**다.
   L5 PRNG periodic resolution은 32인데 확장 FFT 작업영역은 40³이다.
   serial `random.cc:1100`은 `unpack_subvolume_from_buffer()`로 40³ 결과를
   periodic 32³ cube에 되쓴 다음 다시 export한다. 동일 periodic cell에
   여러 값이 들어가 마지막 쓰기 값이 남는다. F.2b는 extended FFT slab을
   직접 export하므로 이 되쓰기와 일치하지 않는다. averaging을 제거한
   L4–5 진단에서 base noise는 bitwise 동일, L5 raw noise relative L2는
   0.1409122111였다. MPI raw field에 serial periodic-last-write mapping을
   적용하면 차이가 1.982e-16, max abs 9.992e-16으로 사라졌다. 이로써
   FFT 반올림이나 grafic writer/header 문제가 아니라 export semantics가
   원인임을 격리했다. F.3 correct_avg는 fine 차이를 base로 전파한다.
   같은 binary에서 `slab_rng_refinement=no`와 `rng_boost_threads=no`를
   함께 쓰면 모든 noise level이 bitwise 같고 24 output fields의 최대
   relative L2는 1.276e-8로 기존 2e-5 tolerance를 통과한다. 반면
   boost만 끄면 기존 density 14.68% 차이가 그대로 남는다.
2. **자동 OMP boost는 별도의 중복 writeback 문제를 노출**한다.
   F.2b를 끄고 boost만 유지하면 output 최대 relative L2가 0.1314276이다.
   `unpack_subvolume_from_buffer()`의 OpenMP `collapse(2)`는 wrap된 cube
   block들이 같은 destination cell을 쓰는 것을 막지 않는다. 40>32
   geometry에서 serial과 multi-thread overwrite ordering이 다르며
   동시 쓰기 race가 존재한다. 같은 4-thread 조건 반복의 noise는 이번에
   bitwise 같았으므로 run-to-run 비결정성이 관측됐다고 주장하지 않는다.
   반복시험 config의 출력 경로가 첫 `legacy_ref_boost/ic`를 가리켜 해당
   **진단용** IC가 다시 쓰였다. 첫 noise/log/JSON과 반복 noise/log는 각각
   보존했고 원래 gate IC/production에는 영향이 없다. 첫 output JSON은
   최초 실행의 측정값이며 이 재쓰기 caveat를 포함해 해석한다.
   L5–7(base 5) non-wrapping 비교는 RNG 작업 크기 32³/56³/88³,
   각 level periodic resolution 32/64/128로 wrap alias가 없다. auto MPI 및
   auto boost를 유지해도 noise 차이 ≤4.45e-16, 24 output fields 최대
   relative L2 9.367e-9로 통과했다. 따라서 현재 작은 fixture의 geometry가
   취약 경로를 자극한 것이며 일반적인 모든 MPI IC가 틀렸다는 결론은 아니다.
3. **step-zero dump는 주기 출력 guard 누락**이다. compiled SIDM
   `amr_step.jaehyun.f90:304`의 `mod(nstep_coarse,foutput)==0`에서
   nstep=0이면 foutput=1e9여도 참이다. 실제 output info의 coarse step=0,
   aexp=0.00990099로 확인했다. aout=1.1/tout=1e100은 미래지만 periodic
   조건만으로 dump가 발생한다. `nstepmax=0`도 initialization 호출을 막지
   않는다. 이번에는 고치지 않았으며 기존 zero-dump gate 실패는 유지한다.
4. **L5 residual 초과는 cap 부족이 아니라 acceptance 정책 차이**다.
   실제 VPATH 우선 source는 `patch/lagRamses/multigrid_fine_commons.f90`
   (cuRamses 사본 아님)이며 relative epsilon 외 absolute floor도 허용한다.
   별도 copied source에 read-only MG_CAUSE print만 추가하여 기존 objects와
   relink했고 상태·iteration/acceptance를 변경하지 않았다. L5 두 번째 solve
   iter=2에서 err=1.3691891333e-4, epsilon=1e-4,
   res_norm2=3.8965203087e-21, i_res_norm2=2.0785000339e-13,
   absolute floor=1e-20였다. relative=False, absolute=True, cap=False로
   종료했다. cap80/abort failure가 아니다. relative-only audit는 이를
   거부하지만 solver는 absolute 조건으로 수용한다. 이를 force accuracy의
   독립 인증 또는 gate 합격으로 확대 해석하지 않는다. diagnostic reader는
   예상 초기 dump 1개(예산 상한 2 MiB)를 명시하고 실행했고 보존했다.
   원래 dump sizes를 재측정하거나 삭제하지 않았다.

증거 JSON은 `cause_1013a/rank_path_report.json`,
`cause_1013b/geometry_report.json`, 절대 잔차 raw log는
`cause_1013a/reader_l5/reader.log`다. production release=false를 유지한다.
후속 권고는 non-wrapping small fixture(예: base5/high7)로 전체 gate 재시험,
SIDM의 step-zero periodic-only output guard 검증, absolute MG 수렴값 및
stop reason의 구조화 로그/감사 정책 일치다. 사용자 요청은 원인 파악이므로
이 production fix들을 이번 작업에 자동으로 적용하지 않았다.

### 10:31 KST LagMUSIC 소규모 검증 승인 및 제출

후속 직접 지시로 소규모 검증을 Lageunha 수동 실행으로 전환했다. 자동 제출
`421344`는 PENDING 상태에서 취소됐고 elapsed 0을 확인했다. 입력 bundle은
변경하지 않았다. 수동 host는 `LagEunha`, 실행 root는
`/gpfs/kjhan/Hydro/Sidm/lm_gate/manual_1013a`, 고정 manual runner와 로그는
`/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_manual_20261013.MoX1Zz`다.
수동 경로는 host guard를 추가해 syntax에서 실행을 거부하고, 가짜 Slurm ID를
사용하지 않는다. 4 ranks × 1 thread, compile 최대 2 processes,
nice 19 및 전체 timeout 2700초로 기존 계산을 보존한다. GPFS 여유 약 279 TB,
새 full dump 예상 0개/0 bytes, small IC 상한 128 MiB 정책을 유지한다.
두 clean build는 성공했다. 첫 시도 `manual_1013a`의 generator는 작은 level-4
base에서 padding 8/overlap 4를 적용한 subgrid가 half-box보다 커서 중단됐다.
실패 로그와 binary를 보존했다. small test에만 padding 2/overlap 2,
nexpand 0,1,1을 사용해 새 root `manual_1013b`에서 재시도했다. 실제 target의
production padding/mask 설정을 바꾸지 않았다. binary hash 검증 후 동일 host의
새 root에 두 바이너리를 복사해 재사용했고 빌드를 다시 하지 않았다. replay
tooling은 `/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_retry_20261013.sSyLrU`다.
아래의 421344 queue 상태는 수동 전환 이전 기록이다.

수동 small gate의 최종 판정은 **실패/production 미해제**다. 검증된 부분과
실패한 부분을 분리해 다음 결과를 보존했다.

- 확장 48-byte header, Omega_b 0.04, 모든 검사 field의 header geometry,
  유한한 payload, binary refinement mask 및 writer namelist metadata는 통과했다.
- 기본 auto MPI RNG 경로의 1/4 rank replay는 24 fields 중 21개에서 실패했다.
  세 refmap은 같지만 density/velocity/displacement는 다르다. level-4 density
  relative L2는 0.1468160916, correlation은 0.9892141289다. tolerance를
  완화하거나 같은 seed만으로 동일 realization이라고 주장하지 않았다. log는
  4-rank에서 F.2b refinement RNG 및 correct_avg_slab, 1-rank에서 기존 RNG
  경로의 실행을 보여준다. 차이를 유발한 구체적 routine은 아직 격리하지 않았고
  LagMUSIC production 코드를 수정하지 않았다.
- level-5/6 1-rank IC를 4-rank reader로 읽은 실제 particle mass sum은 각각
  정확히 1이고 nonfinite mass/position/velocity는 0이다. counts는 각각
  4,971 / 13,497이며 모든 예상 mass tier가 존재하고 합산 count/mass가 같다.
  header/namelist Omega_b check와 실제 FFT direct solve는 두 경우 모두 통과했다.
- level-6 strict MG gate는 max relative error 5.958e-5로 통과했다. level-5는
  한 solve의 relative error 1.369e-4가 epsilon 1e-4 이상이어서 보수적 gate가
  실패했다. warning 또는 cap-80 초과는 없었다. 로그의 absolute-density floor
  조건을 확인하지 않았으므로 이를 코드의 명시적 nonconvergence abort 또는
  확정된 물리 미수렴으로 확대 해석하지 않는다. 앞선 commentary의 두 경우
  모두 MG 통과라는 표현은 전체 strict log 감사 뒤 정정했다.
- zero-dump gate는 실패했다. 두 cold-start가 `output_00001`을 한 개씩 생성했다.
  실제 compiled `patch/lagRamses/amr_step.jaehyun.f90` line 304는
  `mod(nstep_coarse,foutput)==0`을 step-zero guard 없이 검사한다. 큰 foutput도
  coarse step 0에서는 출력을 막지 못했다. snapshot은 삭제하지 않고 보존했으며
  새 두 output의 한 번 측정한 allocated sizes는 2,804 / 1,284 KiB,
  합계 약 3.992 MiB다. dummy future output clocks만으로 suppression을 인증하지
  않는다. 기존 recovery/pilot은 restart gate이며 이번 cold-start failure와 구분한다.

최종 JSON/log 및 검사 코드 snapshot은
`/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_report_20261013.IVusXd`다.
`ic_gate.json`은 header_payload_pass=true, rank_replay_pass=false,
`reader_gate.json`은 mass_header_checks_pass=true, reader_physics_pass=false,
no_dump_gate_pass=false, new_full_dumps=2, production_release=false다.
diagnostic report 생성의 exit 0은 실패들을 기록했다는 뜻이며 numerical gate의
합격을 의미하지 않는다. 네 unit-test method를 다시 통과하고 results.sha256으로
보고서와 두 raw reader log를 검증했다. earlier diagnostic snapshot
`lagmusic_diagnostic_20261013.C3Yq5u`의 level-5 strict gate failure traceback도
보존했다. 최고 IC level 5/6 간 base density correlation 0.9594778578은 진단값이며
같은 phase/target certification이 아니다.

LagMUSIC binary SHA256은
`9e0f493442a65c466cc054055fd7ef332e1bba56e6c92b3d56d177626c0016b9`,
diagnostic reader는
`de83f8c497d98eebe18318b82c25e88ec1f4d53528adb3f2d34e3af518faa5e5`다.
수동 계산은 모두 종료됐다. 대형 IC 또는 생산 실행은 제출하지 않았고 원래
LagMUSIC 미커밋 writer, production source, 기존 데이터 및 다른 프로젝트 작업을
변경하지 않았다. 병렬 RNG 재현성 원인 분리, SIDM cold-start periodic-output
guard의 검증, relative/absolute MG acceptance 정책 확인이 후속 과제다.

사용자는 LagMUSIC 대안을 검토한 뒤 기존 로컬 변경을 보존한 소규모 IC 생성 및
SIDM reader 연동 검증을 명시적으로 승인했다. 아래의 MUSIC2 writer 승인 대기는
해당 경로의 과거 상태이며, 이번 승인을 MUSIC2 저장소 수정 승인으로 해석하지
않는다. 전체 JCAP 목표는 여전히 미완료다. 대형 IC generation, production 실행,
원격 commit/push, Overleaf 변경은 이번 범위에 포함하지 않는다.

LagMUSIC 독립 프로젝트는 `/home/kjhan/BACKUP/LagMUSIC/music`, GitHub 원격은
`https://github.com/kjhan0606/lagMUSIC.git`, source commit은
`35a4eb978ae0e84673d86a0e4722421118277700`이다. 기존 미커밋
`src/plugins/output_grafic2.cc` 변경은 48-byte Omega_b header, 생성 namelist의
cosmological Omega_b 및 slice buffer 처리 교정을 포함한다. 해당 변경을 그대로
분리 bundle에 복사했으며 원래 저장소를 수정하지 않았다.

Slurm `421344` (`sidm_lagmusic_gate`)를 제출했다. 10:31 조회에서
`PENDING (QOSMaxGRESPerUser)`이며 계산 결과는 아직 없다. 4 MPI ranks ×
2 threads, GPU 1개, memory 32 GiB, time limit 45분이다. 고정 bundle은
`/gpfs/kjhan/Hydro/Sidm/validation/lagmusic_submission_20261013.UFaMyF`,
실행 root는 `/gpfs/kjhan/Hydro/Sidm/lm_gate/421344`이다. bundle의 코드,
source override, writer diff 및 unit-test 결과를 SHA256SUMS로 묶고 읽기 전용으로
보존했다. 기존 `421284`, `421295` bundle이나 legacy checkpoint는 변경하지 않았다.

검증 대상은 box 128 h^-1 Mpc, zstart 100, Omega_m 0.3111, Omega_b 0.04의
작은 직육면체 zoom이다. IC base level 4, 최고 level 5/6이며 실제 target halo나
production L10/L11 IC가 아니다. Eisenstein total-matter analytic fit과 기존
density-derived 2LPT velocity 경로를 명시해 CAMB 전환을 분리했다. level-6
IC를 1/4 MPI rank로 생성해 payload replay를 비교하고 level-5 base density와의
상관은 별도 진단값으로 보존한다. 같은 seed만으로 기존 MUSIC2 IC의 위상 또는
다른 highest-level hierarchy의 위상 일치를 인증하지 않는다. binary refinement
mask 검사는 실제 target containment/low-resolution contamination 검사를 대신하지
않는다. 후자는 대형 target-specific IC를 승인받기 전에 별도 gate가 필요하다.

새 isolated reader source는 `eadfc9f`와 기존 reader/timestep 수정 일곱 파일을
사용한다. `init_part` 끝의 read-only MPI mass/finite/tier diagnostic은 고정
테스트 snapshot에만 추가했다. 사용자 worktree의 init_part에는 diagnostic을
추가하지 않았다. 로그의 실제 입자 질량 합 1, 유한한 position/velocity/mass,
각 IC level mass tier/count, extended-header 일치, 실제 FFT dispatch 및 strict
MG convergence를 요구한다. 두 reader cold-start는 `nstepmax=0`인 초기 단계
gate이며 evolution benchmark나 science trajectory 인증이 아니다.

effective namelist는 각 case의 `run.nml`이고 80-character initfile 한계를
넘지 않는 경로만 허용한다. `noutput=1`, `aout=1.1` (z=-0.090909...),
`tout=1e100`, `foutput=fbackup=1e9`로 새 full dump 예상 0개/0 bytes다.
예상 small-IC storage의 보수적 상한은 128 MiB이고 제출 전 GPFS 여유 공간은
약 279 TB였다. runner는 예상하지 않은 `output_*` directory를 거부한다.
실행 전 manifest 검증과 실행 후 effective-input SHA 비교를 유지한다.

`test_lagmusic_small_gate.py`의 네 unit-test method가 grammar에서 통과했다.
계산 실행·payload 결과·full reader build 성공을 의미하지 않는다. 10:26 이후
조회한 기존 CDM recovery `421284_0`는 COMPLETED 0:0이며 cap 40/80 JSON은
각각 511 fine solves, max error 9.347e-5, FFT direct solves 3회,
warning/threshold failure 0을 기록했다. SIDM recovery `421284_1`는 10:31에
RUNNING이고 timestep pilot `421295`는 dependency 대기 상태였다. CDM 단독
통과로 paired recovery 또는 pilot gate를 해제했다고 주장하지 않는다.

최신 제출 상태: 아래의 초기 제출 `421268`과 임시 hold 제출 `421283`은
시작 전에 취소했다. 수정본 `421284`를 제출하고 hold를 해제했으며,
후속 조회 상태는 `PENDING (QOSMaxGRESPerUser)`였다. 기존 생산 런을
재시작하거나 다른 프로젝트 작업을 중단하지 않았다. 수정 및 검증 근거는
아래 「재제출 및 timestep 진단」 절에 기록한다.

2026-10-13 09:47 KST 후속 scheduler 확인에서도 `421284`는
`QOSMaxGRESPerUser`, `421295`는 `afterok:421284_*` 미충족으로 대기했다.
`sacctmgr show qos normal`의 `MaxTRESPU`는 `gres/gpu=8`이다. 09:48 KST에
`sacct AllocTRES`로 확인한 실제 할당은 FDM 작업 1개, LUMINA 작업 2개,
parallel generation array의
다섯 task 각 1개로 합계 8개였다. CPU-only 데이터 작업의 `gres/gpu=0`을
GPU 점유로 계산하지 않았다. 이 관측은 사용자별 한도 대기를 설명하며
SIDM binary/runtime 실패나 예약된 시작 시각을 의미하지 않는다. 다른
프로젝트 작업을 취소하거나 queue priority, QOS, GPU 요청을 바꾸지 않았다.
MUSIC2 writer 교정은 별도 승인 요청 뒤 명시적인 사용자 승인을 기다린다.
goal continuation 자체를 해당 승인으로 해석하지 않는다.

09:50 KST 확인에서도 두 Slurm handle은 여전히 PENDING이고 종료 상태가
아니었다. MUSIC2 writer의 별도 저장소 변경 승인 미수신이 세 연속 goal turn에
걸쳐 지속됐다. reader block/전체 build, IC 해상도/단위, angular normalization 및
현 scheduler 상태의 안전한 검증은 마쳤다. 다음 IC 준비 단계는 writer 변경의
명시적 승인이 필요하므로 전체 목표는 미완료인 승인 대기 상태로 넘긴다.
GPU queue 대기만을 실패나 blocker로 해석하지 않는다. 이미 제출한 `421284`,
`421295`는 취소하거나 수정하지 않았으며 scheduler dependency와 runtime
gate를 유지한다. 승인 요청 범위는 MUSIC2의 Omega_b 확장 writer와 소규모
검증이고 대형 IC generation, production launch, 외부 push는 포함하지 않는다.

이 절은 아래의 8월 초기 상태보다 우선하는 현황 기록이다. Slurm 조회에서
기존 SIDM 생산 런은 실행·대기 중이지 않았으며, L17 연속 실행 네 개는 모두
정상 job-control 정지와 `Run completed`, 종료 코드 `0:0`을 남겼다.

| 런 | 최신 출력 | scale factor | redshift | 마지막 job / 종료 |
|---|---|---:|---:|---|
| L17 CDM | `output_00005` | 0.573638632 | 0.743258 | `324571` / 2026-08-28 |
| L17 SIDM3 | `output_00005` | 0.250336146 | 2.994629 | `324572` / 2026-08-28 |
| L18 SIDM1 | `output_00012` | 0.344282459 | 1.904592 | `323699` / 2026-08-25 |

L17의 마지막 두 출력에는 AMR·gravity·particle 파일이 각각 24개씩 있고 빈
파일은 없었다. 다만 legacy binary가 생성한 출력에는 `COMPLETE` marker가
없다. 파일 개수와 정상 종료는 full payload 무결성이나 물리 수렴의 증명이
아니므로, 재시작 reader 검증과 solver 검증을 별도로 거쳐야 한다.

기존 로그 여섯 개를 `epsilon=1e-4`, iteration cap 10으로 재감사했다.

| 로그 | warning 수 | cap 도달·허용오차 초과 수 | 최대 final relative error | FFTW enabled banner |
|---|---:|---:|---:|---|
| z4 CDM `323683` | 17 | 1 | 4.206e-3 | 없음, 요청 무시 |
| z4 SIDM3 `323771` | 17 | 1 | 4.204e-3 | 있음 |
| L17 CDM `323844` | 1 | 1 | 1.706e-4 | 없음, 요청 무시 |
| L17 CDM `324571` | 0 | 0 | 9.998e-5 | 없음, 요청 무시 |
| L17 SIDM3 `323845` | 1 | 1 | 2.043e-4 | 있음 |
| L17 SIDM3 `324572` | 1 | 1 | 1.106e-4 | 있음 |

L17의 세 실패는 모두 첫 fine step 전 AMR level 11 solve에서 발생했다.
CDM과 SIDM3의 기존 binary는 서로 다른 7월 source revision에서 빌드됐으며,
CDM에는 FFTW support가 없었다. 기존 z4 분석 runner는 MG summary를
기록할 뿐, 실패를 이유로 정지하지 않았다. 따라서 기존 profile 비교가 수치를
생성했더라도 science convergence 인증은 보류한다.

`tests/sidm/analyze_mg_convergence_log.py`에 strict gate를 추가했다.
`--require-converged`는 warning, 비유한 residual, 허용오차 이상의 relative
residual, solve 기록 누락을 거부한다. absolute-density floor는 로그만으로
판정할 수 없으므로 모호한 high relative residual도 보수적으로 거부한다.
`--require-fftw`는 enabled banner와 실제 `FFT direct solve DONE` 기록을 모두
요구한다. z4 분석 runner는 profile 분석 전에 두 검증을 실행한다. 회귀 시험
`test_mg_convergence_gate.py` 10개와 shell syntax check가 통과했다.

짧은 operator/restart 검증 array `421268`을 제출했다. 기존 z4 checkpoint를
입력으로 현재 source의 첫 Poisson reconstruction과 iteration headroom을
측정한다. 이 검증만으로 science seed의 타당성을 인증하지 않는다.
각 task는 24 MPI ranks × 2 threads, 200 GiB, GPU 1개, 최대 2시간을 요청하며,
array concurrency는 1이다. 제출 직후 상태는 `PENDING (QOSMaxGRESPerUser)`였다.

- 고정 source revision: `eadfc9fdac3e10cac3bf3415c05ab67fb061ae96`
- 고정 runner/analyzer: `/gpfs/kjhan/Hydro/Sidm/validation/l17_mg_submission_20261013.aU9BeU`
- 실행 위치: `/gpfs/kjhan/Hydro/Sidm/validation/l17_mg_recovery/421268/{cdm,sidm3}/cap{40,80}`
- effective namelist: 위 각 directory의 `run.nml` (job 시작 시 생성)
- build: compute node의 source archive에서 `HDF5=1 USE_FFTW=1`, 두 모델에 같은 source/options
- solver: `epsilon=1e-4`, `maxiter_fine=40,80`, `abort_on_mg_nonconvergence=T`, predictor restart
- stopping: seed coarse step부터 2 step. CDM `nstepmax=103`, SIDM3 `nstepmax=222`
- output: `noutput=1`, dummy `aout=1.1` (z=-0.090909), `tout=1e100`, `foutput=fbackup=1e9`
- 신규 full dumps: 0 (seed link 이외의 output이 생기면 gate failure)
- seed 크기: CDM 6.571 GiB, SIDM3 6.364 GiB. 제출 전 GPFS free space는 `df -h` 기준 279T

runtime 통과는 아직 확인하지 않았다. 다음은 array의 build/runtime summary,
실제 MG iteration 수, FFTW solve, 종료 조건, 출력 억제를 확인하는 것이다.
이후 측정한 iteration headroom과 source를 고정하고, 기존 seed의 solver 차이가
profile에 준 영향을 검증해 공통 epoch의 CDM/SIDM 재계산을 설계한다. 기존
checkpoint의 restart 성공만으로 과거 trajectory의 물리 수렴을 증명하지 않는다.

### FFTW Poisson 실행파일 근거 추가 검사

동일 날짜의 후속 audit에서 legacy SIDM3의 `FFTW enabled` banner를 실제
FFTW dispatch 증거로 볼 수 없음을 확인했다. 보존된 L17 SIDM3 실행파일의
regular symbol table에는 일반 FFTW plan/execute 함수가 있지만 RAMSES의
`fftw_poisson_solve_uniform_` entry point는 없었다. 실제 로그에도
`FFT direct solve DONE`이 없고 base level 8의 MG iteration 기록이 있었다.
따라서 `ldd`에 동적 FFTW library가 없다는 사실만으로 FFTW library 미포함을
판정하면 안 된다. 정적 library 포함과 Poisson operator 포함은 별개의 조건이다.

`tests/sidm/audit_poisson_binary.py`는 binary를 실행하지 않고 `nm -g`,
`nm -D`로 Poisson entry point와 FFTW plan/execute symbol을 함께 검사한다.
정적 및 동적 linking을 지원하고, library만 존재하거나 feature variable만
존재하는 binary는 거부한다. stripped binary 등 근거가 부족한 경우도 통과로
취급하지 않는다. 이 검사는 linked-code 증거이며 mesh 조건, 실제 dispatch,
solver 정확도 또는 과거 trajectory를 인증하지 않는다.

실제 binary 세 개를 검사했다. CDM과 SIDM3의 legacy binary는 해당 gate를
통과하지 못했다. 두 SHA256은 저장된 `binary_324571.sha256` 및
`binary_324572.sha256`과 일치했다. 현재 worktree의 `bin/ramses_final3d`는
linked-code gate를 통과했다. 현재 binary를 production으로 실행하지 않았으며
source/options 일치나 runtime 검증 완료를 주장하지 않는다. 구형 binary의
symbol/banner 불일치가 발생한 정확한 build 경위는 확인하지 못했다.

- binary audit 결과: `/gpfs/kjhan/Hydro/Sidm/validation/poisson_binary_20261013.V212pt/`
- legacy CDM SHA256: `dfdb23e406da06472e7c076d5264963233714f69440def1bfb215335167099eb`
- legacy SIDM3 SHA256: `fabd3a888c386e342e6c1b81aad6b0e4d3f98fc3bed80aa25c84d5218e6cbd1e`
- 현재 local binary SHA256: `b15abca01adc4c5017f720fa03a261d979fc8614fe2d94be9f2d677e0942e981`

binary gate의 synthetic 회귀 시험 8개, 기존 log gate 시험 10개 및 수정된
runner의 shell syntax check가 통과했다. 로컬 recovery runner의 차기 제출
version에는 build 직후 binary gate를 넣었다. 이미 제출한 `421268`의 고정
snapshot은 수정하지 않았다. 기존 제출본은 source-only archive에서 새로
빌드하며, 실제 FFT solve 기록과 MG 수렴을 요구하는 runtime gate를 유지한다.
새 archive에는 tracked `.o`/`.mod`/`ramses_final` binary가 없음을 확인했다.

### 재제출 및 timestep 진단

빌드 검증에서 두 실행환경 문제를 실제로 재현했다. 현재 Lmod의 Intel compiler
module은 `intel/tbb`, `intel/compiler-rt`, `intel/umf`를 선행 조건으로 요구한다.
`module purge; module load intel/compiler/2025.3.0 intel/mpi/2021.17`은
`syntax`와 `grammar` 모두에서 실패했다. 수정 runner는 다음 module을 명시한다.

`intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2
intel/compiler/2025.3.0 intel/mpi/2021.17`

비어 있는 source archive의 `make -j2`도 실패했다. Makefile에 완전한 Fortran
module dependency가 없어 `amr_commons`를 `amr_parameters.mod` 생성 전에
컴파일했기 때문이다. 검증 runner는 `make -j1`으로 변경했다. 수정 환경에서
별도 source archive의 전체 `HDF5=1 USE_FFTW=1` build가 성공했고 binary gate도
통과했다. 이 archive에는 아래 timestep 진단의 로컬 source patch 네 개도 포함된다.
기존 worktree binary와 사용자 소유 module 파일은 덮어쓰거나 제거하지 않았다.

`421268`은 `2026-10-13 08:47:36 KST`, 중간 hold 제출 `421283`은
`08:52:14 KST`에 취소됐으며 둘 다 elapsed `00:00:00`이었다. 최종 `421284`는
다음 고정 snapshot을 사용한다.

- runner/analyzers: `/gpfs/kjhan/Hydro/Sidm/validation/l17_mg_submission_20261013c.y29dUI`
- runner SHA256: `698cb3ef284ef7f5ead816e24416b6dc4d79c43c9ec29e2c98fee608d38bf35f`
- source revision: 기존과 같은 `eadfc9fdac3e10cac3bf3415c05ab67fb061ae96`
- 실행 root: `/gpfs/kjhan/Hydro/Sidm/validation/l17_mg_recovery/421284/{cdm,sidm3}/cap{40,80}`
- effective namelist: 위 각 directory의 `run.nml`, task 시작 후 생성
- 24 ranks × 2 threads, 200 GiB, GPU 1개, task당 최대 2시간, array concurrency 1
- dummy `aout=1.1`, `tout=1e100`, `foutput=fbackup=1e9`, 신규 full dumps 0
- build 직후 binary gate와 실행 후 실제 FFT solve 및 strict MG gate를 모두 요구

최종 제출은 순수한 고정 commit에서 중력/restart gate를 측정한다. 로컬 timestep
진단 source patch는 `421284`에 포함하지 않는다. 기존 고정 snapshot 세 개는
덮어쓰지 않고 보존했다. 구형 제출의 실패를 새 제출의 runtime 실패로 해석하지 않는다.

SIDM timestep 계측은 기존 `sidm_estimator_diagnostics=.true.`에서만 켜지며
기본값은 계속 false다. `newdt_fine`은 비-SIDM 제약 후의 timestep, 이전 timestep,
이전 `Pmax`, courant, SIDM 후보와 적용 후 timestep을 `SIDM_DT`로 기록한다.
기존 `MIN(dt, sidm_courant/sidm_Pmax*dtold)` 식은 그대로 유지했다.
`amr_step`은 뒤이어 적용되는 parent/subcycle bound와 선택 결과를
`SIDM_DT_PARENT`로 기록한다. 두 patch의 같은 경로를 함께 수정했다.

`analyze_timestep_diagnostics.py`는 두 phase의 coarse/fine/level/epoch 일치,
후보식, minimum, 유한한 값 및 pair completeness를 검증한다. legacy 로그처럼
해당 기록이 없는 경우 원인을 역추정하지 않고 거부한다. summary의 call 수는
독립 물리 step 수가 아니다. child recursion이 부모 dt를 나중에 조정할 수 있고,
비-SIDM 제약은 하나로 묶으므로 이 기록만으로 최종 coarse dt 또는 wall-time
speedup을 증명하지 않는다. 진단은 짧은 pilot에서 사용하며 장기 로그에는 기본
off를 유지한다.

`test_timestep_diagnostics.py`는 실제 두 source의 표시된 block을 추출해 실행했다.
GNU Fortran 13.2.0과 Intel ifx 2025.3.0에서 각각 256개 scalar/array 통제 사례의
결과가 기존 식과 bitwise 동일했고 32개의 paired telemetry record가 검증됐다.
log parser의 회귀 시험 8개도 통과했다. 이는 full cosmological/MPI smoke를
대신하지 않는다. 전체 clean build 외에 격리된 directory에서 fresh module을
사용해 두 patch의 수정된 Fortran unit 네 개를 ifx `-O3`로 컴파일했다.
project root의 오래된 `amr_parameters.mod`는 별도 compile의 module 검색을
오염시켰으므로 보존한 채 작업 directory를 격리해 검증했다.

근거 root는 `/gpfs/kjhan/Hydro/Sidm/validation/sidm_dt_20261013.5tDoHn`이다.
`gnu_current/`와 `ifx_blocks/`에 테스트 JSON 및 로그, `instrumentation.diff`에
source 변경, `build_serial.log`에 전체 build, `clean_build_binary_gate.json`에
새 실행파일 근거를 저장했다. 실패했던 병렬 build log도 삭제하지 않았다.
새 binary SHA256은 `c1a404451cafc2db402c5640e82458690227c18eb35c9ba536156811bf36893b`다.
진단을 포함한 새 binary로 cosmological runtime pilot을 실행하지는 않았다.

### 조건부 성능 pilot 제출 `421295`

중력 gate 이후에 진행 속도와 timestep 선택을 측정할 short diagnostic evolution을
제출했다. `421295`의 scheduler dependency는 `afterok:421284_*`이며, 초기 상태는
`PENDING (Dependency)`였다. 기존 생산 런을 재개하는 job이 아니다.
실행 전에 CDM/SIDM3의 cap 40/80 네 결과를 다시 검사한다. 저장된 green summary만
읽지 않고 각 원시 로그의 strict MG 및 실제 FFT solve도 재검증한다.
cap 하나라도 실패하거나 source revision이 다르면 pilot을 실행하지 않는다.

- frozen bundle: `/gpfs/kjhan/Hydro/Sidm/validation/l17_dt_submission_20261013.IZDpzN`
- 실행 root: `/gpfs/kjhan/Hydro/Sidm/validation/l17_timestep_pilot/421295`
- effective namelist: 위 root의 `run.nml`, task 시작 시 bundle의 `run.nml`에서 복사
- seed: 기존 L17 SIDM3 `output_00005`, `a=0.250336146250867`, coarse 838, 24 ranks
- stopping: coarse 858, 20 coarse steps, 최대 2시간
- 자원: 24 MPI ranks × 2 threads, 200 GiB, GPU 1개
- output: `noutput=1`, dummy `aout=1.1` (z=-0.090909), `tout=1e100`, `foutput=fbackup=1e9`
- 새 full dumps: 0. seed는 6.680 GiB의 원본에 대한 link이며 새로 복제하지 않음
- 제출 시 GPFS free space: `df -h` 기준 279T
- `epsilon=1e-4`, cap 80, nonconvergence abort, predictor restart, constant cross section 3,
  `sidm_courant=0.1`, isotropic angular kernel, opt-in estimator/timestep diagnostics

bundle은 앞 절의 clean-build binary `c1a404...893b`를 복사한 것으로 새로 production
binary를 덮어쓰지 않았다. 빌드 근거, 네 source unit, source diff, GNU/ifx 통제 시험,
namelist와 분석 code를 SHA256 manifest로 묶고 파일을 read-only로 보존했다.
`validate_pilot_bundle.py`는 실행 전에 hash, seed metadata, output/physics 계약과
operator gate를 검사한다. 현재 recovery artifact가 없을 때 admission이 거부됨을
확인했다. Scheduler의 afterok 조건과 artifact 검증을 모두 요구한다.

`analyze_timestep_diagnostics.py`의 pilot 모드는 initial adjustment 2 coarse steps를
제외한다. coarse-boundary의 high-precision base-level `aexp`와 해당 completed step의
wall-time 기록을 맞춰 국소 속도를 계산한다. 네 자리 Fine-step epoch를 정확한
끝점으로 쓰지 않는다. final base record가 없으면 정확한 마지막 base epoch까지만
측정한다. pair completeness, stop step, 연속한 epoch, positive/finite timing 및
strict gravity evidence가 모두 필요하다. 조건부 실행/분석 회귀 시험 13개와 기존
parser 시험 8개 및 두 production block의 256개 bitwise 통제 사례가 통과했다.

성공 시 `pilot_summary.json`은 timestep 제약별 call 통계, warmup을 제외한
seconds/coarse step 및 seconds/delta-a, 동일 국소 속도 가정의 z=0 운영 외삽을
기록한다. 짧은 구간의 외삽은 steady-rate 보장이나 90일 완료 약속이 아니다.
legacy seed의 과거 trajectory가 인증되지 않았으므로 결과는 halo science profile,
CDM matched-epoch 상대 speedup 또는 force/mass 수렴의 증거로 사용하지 않는다.
이 시점에는 pilot의 실제 runtime과 결과를 아직 확보하지 못했다.

## 2026-10-13 angular-kernel 검증 확장

`grammar`에서 production helper `patch/cuRamses/sidm_angular.f90`를
GNU Fortran 13.2.0으로 직접 컴파일해 검증했다. Rutherford regulator는
`epsilon=0.001,0.01,0.1,1,10`으로 확장했으며, 각 regulator에서 50,000개
표본을 검사했다. CDF의 KS 통계량은 각각 `0.00243,0.00593,0.00381,0.00504,
0.00412`이며 공통 기준 `1.95/sqrt(50000)=0.00872`보다 작았다.

회전 기하 검증은 좌표축 여섯 방향과 `abs(v_hat_z)=0.9` 분기 경계의 양쪽
네 방향에서 극각 다섯 개와 방위각 네 개를 사용한다. 200개 조합에서 unit
direction의 norm과 incoming direction에 대한 cosine이 `2e-13` 이내로
일치했다. 기존의 oblique direction에 대한 표본별 회전 검사도 유지했다.

다음 세 angular weight의 평균을 검사했다. `mu=cos(theta)`, `b=2*epsilon`,
`L=log(1+2/b)`일 때 normalized angular density는
`p(mu)=b*(b+2)/(2*(1-mu+b)^2)`이다.

- 운동량 weight `1-mu`: 평균 `b*(b+2)*L/2-b`
- 점성 weight `1-mu^2`: 평균 `b*(b+2)*((b+1)*L-2)`
- identical-particle weight `1-abs(mu)`: 평균 `b*(b+2)*log(1+1/(b*(b+2)))/2`

세 sampled mean은 모든 regulator에서 해석값과 5 standard errors 이내로
일치했다. `test_angular_transport_integrals.py`는 별도 adaptive quadrature로
정규화와 세 평균을 재계산하며, 20개 적분 비교가 absolute difference `3e-10`
이내로 통과했다. `1-abs(mu)`에는 추가 factor of two를 적용하지 않았다.
다른 문헌의 수송 단면적과 비교할 때 weight 정의와 factor를 함께 확인한다.
기본 `epsilon=0.01`에서 total cross section을 각 weighted cross section으로
나눈 값은 운동량 `13.6564570`, 점성 `9.1424487`, identical-particle weight
`15.2391817`이다.

측정 JSON과 PDF/PNG는
`/gpfs/kjhan/Hydro/Sidm/validation/angular_kernel_20261013.8KwC79/`에 보존했다.
JSON에는 production helper, Fortran driver, Python verifier의 SHA256과
compiler version, 각 측정의 표준오차를 기록했다. Figure의 선은 해석값이며,
점과 오차 막대는 표본 평균과 5 standard errors이다. 우측 panel은 KS 통계량과
통과 기준을 비교한다.

검증 범위는 Rutherford inverse-CDF helper와 lab-frame rotation이다.
Standalone driver의 uniform input은 Fortran `random_number`로 생성하며,
production operator의 thread-local `ranf`나 velocity-dependent regulator
선택까지 인증하지 않는다. 이 검증은 cosmological angular halo response 또는
mass/force resolution convergence를 대신하지 않는다.

원고 근거 대조에서 `main.tex`의 angular table과 본문은 세 regulator와 운동량
weight만 설명하고 있음을 확인했다. 원고의 Yukawa angular 설명에도 코드의
`max(v_rel_kms^2,1)` floor를 명시할 필요가 있다. 원고는 기존 로컬 수정이
있어 문장을 추가로 바꾸지 않았다. 추가·수정 문장의 polishing은 사용자가
지정한 Fable 경로로 진행한다. 현재 세션에는 호출 가능한 Fable 도구가 없다.

### 등방 모델과 angular transport 정규화 대조

2026-10-13 후속 계산에서 `sigma_total/sigma_weighted`와 동일 weighted cross
section을 가진 등방 모델에 대한 `sigma_total,R/sigma_total,iso`를 구분했다.
등방 density는 `p_iso(mu)=1/2`이고 세 weight의 등방 평균은 각각
`1`, `2/3`, `1/2`다. 같은 weight를 쓰는 두 모델의 total cross section 비율은
`<weight>_iso/<weight>_R`이다. 따라서 앞서 기록한 점성 inverse moment
`9.1424487`은 등방 모델 matching multiplier가 아니다.

| Rutherford epsilon | 운동량 matching | 점성 matching | identical-particle matching |
|---|---:|---:|---:|
| 0.001 | 84.5213723 | 33.8238134 | 45.2081047 |
| 0.01 | 13.6564570 | 6.0949658 | 7.6195909 |
| 0.1 | 3.0530906 | 1.7267189 | 1.9169044 |
| 1 | 1.2943497 | 1.0489894 | 1.0612734 |
| 10 | 1.0328023 | 1.0009087 | 1.0011359 |

표의 숫자는 dimensionless `sigma_total,R/sigma_total,iso`다. 예를 들어
등방 `sigma_total/m=3 cm^2/g`, Rutherford `epsilon=0.01`에 대해 동일
점성 수송 단면적을 얻으려면 Rutherford total은 `18.2848974 cm^2/g`다.
동일 ordinary momentum weight는 `40.9693711 cm^2/g`, 동일
identical-particle weight는 `22.8587726 cm^2/g`를 요구한다. 동일 total
`3 cm^2/g` 비교는 동일 수송 강도 비교가 아니다.

문헌의 계수를 원 논문에서 확인했다. Robertson, Massey & Eke (2017)의
[angular scattering 논문](https://arxiv.org/html/1612.03906), §2.2는
modified weight `1-|mu|`를 사용하며 등방 limit은 total의 절반이다.
Yang & Yu (2022),
[Gravothermal evolution](https://arxiv.org/html/2205.03392), Eq. (5)는
점성 weight에 `3/2`를 곱해 등방 cross section을 total과 같게 정의한다.
Arido, Fischer & Garny (2025),
[hybrid scattering 논문](https://arxiv.org/html/2410.07175v2), §2.2는
modified weight에 `2`, 점성 weight에 `3/2`를 사용한다. 같은 계수를 양쪽
모델에 곱하면 위 matching 비율은 유지된다. JSON에는 unscaled moment와
등방 limit을 1로 만든 moment를 구분해 기록했다.

점성 matching은 준정상 isolated-halo 비교의 근거가 있지만 임의의 observable에
대한 equivalence 증명이 아니다. 위 2025 논문의 §6.3은 merger의 DM-galaxy
offset을 modified transfer나 viscosity 하나로 일반적으로 matching할 수
없음을 보인다. Ordinary `1-mu` weight를 identical-particle halo response의
유일한 비교 기준으로 채택하지 않는다. angular helper 검증을 새로운
cosmological angular-response 결과나 최초 구현으로 주장하지 않는다.

코드의 확률과 timestep 비용도 구분한다. 입자 위치·질량·속도·cell volume·dt를
고정하면 pair probability는 total cross section에 선형이다. SIDM 제한이
지배하고 같은 이전 상태를 사용하는 국소 timestep 후보는 total 증가에 반비례한다.
따라서 기본 epsilon의 점성 matching은 국소 확률을 6.095배 높일 수 있다.
AMR/중력 timestep과 evolution이 달라지므로 wall time이 6.095배라는 측정이나
예측은 아니다. 기존 110.59일의 비인증 trajectory 외삽에 표의 비율을 곱해
새 production 완료일을 제시하지 않는다. pilot `421295`는 등방 SIDM3이며
matched Rutherford throughput을 측정하지 않는다.

`test_sidm_angular.py`의 JSON과 verification figure에 matching 비율을
추가했다. figure의 세 번째 panel은 위 total 비율이며 첫 두 panel은 원래의
표본 평균 및 CDF 검증이다. analytic moment 계산은 `b>=20`에서 24항의
적분 series를 사용해 큰 `b`에서 점성식의 cancellation을 피한다. 생산
Fortran sampler나 산란 법칙은 바꾸지 않았다. quadrature 비교는 기존 다섯
regulator에 `b=19.999,20.001,200,20000,2e8`을 추가했다. 등방 평균,
equal-transport 단면적, 문헌 정규화, isotropic limit, invalid regulator 및
아래 velocity-law 차이의 회귀 시험 7개가 `grammar`에서 통과했다. 기존
25만 angle 표본과 rotation 200 cases도 다시 통과했다.

최신 artifact는
`/gpfs/kjhan/Hydro/Sidm/validation/angular_matching_20261013.yy2D2v/`의
`angular_validation.json`, PDF/PNG, `integral_tests.log`, `angular_test.log`다.
이전 artifact 두 개는 덮어쓰지 않았다. 원고 및 Overleaf는 변경하지 않았다.

### Yukawa 명칭의 모델 범위 확인

같은 audit에서 source의 `sidm_type='yukawa'`는 total normalization을
`sigma0/[1+(v/v0)^2]^2`로 계산하고 pair probability에 넣음을 재확인했다.
`sidm_angular='yukawa'`는 normalized angular shape에
`b=2*v0^2/max(v_kms^2,1)`을 사용한다. 두 옵션은 독립적인 선택이다.
Robertson et al. angular 논문의 Eqs. (12)–(13)에 제시된 Rutherford Born
differential law를 각도로 적분하면 total은
`sigma0/[1+(v/v0)^2]`다. 동일한 sigma0와 v0를 두면 코드의 total은 Born
total의 `1/[1+(v/v0)^2]`다. `v=v0`에서는 절반이고 `v=10*v0`에서는
`1/101`이다. 따라서 현재 옵션 조합을 그 Born differential model의 full
normalization으로 인증하지 않는다. identical-particle Møller amplitude까지
계산한 물리 모델로도 인증하지 않는다.

원고 Eq. `sidm_sigma`는 코드의 제곱 velocity law를 이미 명시한다. 법칙을
조용히 바꾸지 않는다. 현재 옵션을 phenomenological parameterization으로
한정할지 별도 Born-consistent 옵션을 도입할지는 후속 모델 선택이다. 원고의
물리 해석과 참고문헌 범위를 검토한 뒤 결정해야 한다. JSON의 velocity table은
`v0=100 km/s`에서 `v=0.1,1,10,100,1000 km/s`를 대조하고 1 km/s angular
floor를 기록한다. source SHA256과 확인한 formula contract를 보존했다.
Rutherford Born total을 독립 quadrature로 확인했으나 실제 velocity-dependent
operator의 runtime 검증으로 취급하지 않는다. velocity-dependent angular
shape에 고정 epsilon의 단일 matching multiplier를 적용하지 않는다.

다음 angular 물리 시험은 동일 seed의 등방 기준, 같은 total의 Rutherford,
같은 점성 수송의 Rutherford를 구분해 설계한다. 작은 controlled 시험에서
response와 timestep 수렴을 확인한 뒤 cosmological 비교의 실행 예산을 판단한다.
`epsilon=0.1`의 상대적으로 작은 matching 비용은 초기 operator 시험에 쓸 수
있지만 `epsilon=0.01`의 검증을 대체하지 않는다. 새 런, output schedule 변경,
새 microphysics 옵션이나 hybrid algorithm 구현을 이번 audit에서 승인하지 않았다.

## 2026-10-13 L17 시간·저장 용량 예산

`tests/sidm/audit_l17_runtime_budget.py`로 CDM/SIDM3의 완료 구간 네 개를
재계산했다. 시작 checkpoint의 coarse step, 종료 checkpoint의 coarse step과
epoch, 실제 로그의 rank 수, 저장된 effective namelist의 restart 번호가 서로
일치해야 계산한다. 종료 및 elapsed 기록 누락도 거부한다. 설정된 24 ranks ×
2 threads는 보존된 runner에서 확인한다. synthetic 회귀 시험 16개가 통과했다.
이 audit는 계산 실행이나 science gate 통과를 승인하지 않는다.

| 모델 / job | 측정 scale factor 구간 | elapsed | 현재 endpoint에서 z=0까지 동일 속도 추정 |
|---|---|---:|---:|
| CDM `323844` | 0.201033 → 0.413541 | 22.16 h | 2.55 compute days |
| CDM `324571` | 0.413541 → 0.573639 | 22.03 h | 2.44 compute days |
| SIDM3 `323845` | 0.200094 → 0.244123 | 21.96 h | 15.71 compute days |
| SIDM3 `324572` | 0.244123 → 0.250336 | 22.00 h | 110.59 compute days |

SIDM3의 단위 `delta_a`당 비용은 두 구간 사이에 7.10배 증가했다. 이 수치는
동일 epoch의 CDM/SIDM 상대 성능 측정이 아니다. 마지막 SIDM3 속도를 고정한
단순 외삽에서도 90일을 넘으며, 가상 allocation availability가 50%이면
calendar time은 221.18일이다. z=4의 공통 seed부터 같은 속도로 다시 계산하는
가상 시나리오는 약 118일의 compute time을 요구한다. 이러한 숫자는 신뢰구간도
완료 보장도 아니다. solver/source 변경, AMR 증가, timestep 변화 및 실제 queue
대기는 별도 benchmark가 필요하다. 기존 trajectory는 중력 gate를 통과하지 않았다.

SIDM3의 마지막 checkpoint에서 `z=2,1,0.5,0`까지 추정 compute time은 각각
`12.24,36.83,61.42,110.59`일이다. `z=2,1,0.5,0`의 네 공통 epoch를 비교용
출력 후보로 기록했으며, 기존 output schedule을 바꾸거나 새 런을 제출하지 않았다.
현재 legacy checkpoint를 이어서 만든 결과를 인증된 science 결과로 쓸 수는 없다.

최신 dump를 한 directory 깊이에서 한 번씩 측정했다. CDM은 6.685 GiB,
SIDM3는 6.680 GiB이며 각각 77 files, AMR/grav/part 각 24 files, 빈 파일 0개다.
이는 payload 무결성 검증이 아니다. 공통 z=4 seed에서 각 모델의 네 epoch를
저장할 때 dump 크기를 현재의 1/2/4배로 가정한 추가 용량은
`53.46/106.92/213.84 GiB`다. 크기 증가율은 가정이며 상한이 아니다. IC, build,
로그, 비상 checkpoint, 수렴 실험 및 재계산 용량은 별도로 더해야 한다.

JSON 근거는
`/gpfs/kjhan/Hydro/Sidm/validation/l17_budget_20261013.Mf8pip/runtime_storage_budget.json`이다.
입력 로그와 audit code의 SHA256, 각 구간의 solver gate, epoch별 추정과 용량
시나리오를 보존했다. snapshot size는 하루 안에 반복 측정하지 않는다.

다음 실행 판단은 대체 제출 `421284`의 operator gate 이후에 한다. 같은 source와
solver로 공통 epoch의 짧은 runtime benchmark를 수행해 timestep과 wall time을
측정하고, 90일 예산을 다시 평가한다. 두 coarse step의 recovery test는 장기
throughput benchmark가 아니다. L17/L18 cap 변경은 force/AMR 시험이며 동일
입자 질량의 시험으로 mass-resolution 수렴을 대체하지 않는다.

## 2026-10-13 zoom IC 질량·격자 해상도 감사

`tests/sidm/audit_zoom_resolution.py`는 기존 IC의 첫 Fortran header만 읽고
기존 입자 QC 및 snapshot 단위와 대조한다. IC payload, 모든 particle 파일,
noise cache를 새로 읽거나 디렉터리를 재귀 순회하지 않는다. IC level 8–11의
density, 세 velocity, 세 displacement, refinement mask, particle-variable
field 36개의 header 총 `1872 bytes`를 읽었다. 각 level에서 모든 field의
geometry/cosmology가 같고, level 간 cell width가 1/2씩 감소했다.

| IC level | 실제 header dimensions | box의 equivalent uniform nside |
|---|---|---:|
| 8 | 256 × 256 × 256 | 256 |
| 9 | 100 × 120 × 94 | 512 |
| 10 | 168 × 208 × 152 | 1024 |
| 11 | 304 × 382 × 272 | 2048 |

`music_zoom_l11.used.conf`와 선택한 Lagrangian point file의 SHA256은
보존된 provenance manifest의 두 항목과 일치했다. density/velocity의 실제
transfer function 및 IC generator binary/source revision은 이 검사로
인증하지 않는다. 공통 box는 nominal `128 h^-1 Mpc`, IC 시작은 `z=100`,
`baryons=no`, transfer plugin은 `eisenstein`, random seeds는 level 8–11에서
`20260401,20260402,20260403,20260404`다. 같은 seed 설정만으로 새 IC의
phase matching을 증명하지 않는다.

세 런의 보존된 effective namelist `run_324571.nml`, `run_324572.nml`,
`run_323699.nml`에서 같은 네 IC path와 `levelmin=8`을 확인했다.
AMR cap은 L17 CDM/SIDM3에서 17, L18 SIDM1에서 18이다. 세 최신 snapshot의
`unit_d*unit_l^3`가 상대오차 `1e-12` 이내로 같았으며 cosmology도 IC header와
일치했다. 즉 기존 L17/L18 비교는 입자 질량을 바꾸지 않는다. AMR cap 변화는
힘 격자와 leaf-cell SIDM interaction locality를 함께 바꿀 수 있으므로 순수한
입자 질량 수렴 시험으로 제시하지 않는다. nominal comoving cell width는 cap
16/17/18에서 각각 약 `1.953/0.977/0.488 h^-1 kpc`다. cell width는 실제
occupied finest cells 또는 calibrated gravitational softening 측정이 아니다.

cached `output_00003_particle_qc.json`은 legacy CDM snapshot에서 입자 파일
32개, 전체 `31,981,734` particles, 최고 해상도 `13,699,720` particles를
기록한다. mass tiers는 IC level 11/10/9/8의 `2^(-3*level)`과 정확히 일치하고
합산 code mass는 1이다. 이 QC를 새로 계산하지 않았다. 별도 L17/SIDM 출력의
payload 검사나 low-z science acceptance로 취급하지 않는다.

snapshot cgs units와 `M_sun=1.98847e33 g`으로 환산한 질량은 다음과 같다.
`h=H0/100`을 적용한 h^-1 단위 숫자를 별도 열에 기록했다. 과거 catalogue의
근사 cosmological-constant 환산과 섞지 않는다.

| 최고 IC level 후보 | particle mass / M_sun | particle mass / (h^-1 M_sun) | 실현 상태 |
|---|---:|---:|---|
| 9 | 1.98463e9 | 1.34280e9 | 새 mask/IC 필요 |
| 10 | 2.48078e8 | 1.67850e8 | 새 mask/IC 필요 |
| 11 | 3.10098e7 | 2.09812e7 | 기존 IC 존재, science 인증 미완료 |

현재 base tier 입자는 `16,679,773`개로 전체의 52.15%다. 최고 tier 42.84%만
eight-to-one grouping하고 나머지를 유지하는 count 시나리오는 총
`19,994,479` particles, 즉 기존 대비 `1.5995`배 감소다. 실제 regenerated
mask나 contamination 검증 결과가 아니다. 기존 base 입자만 그대로 남긴
절대적인 count 하한도 1668만 개이므로, base tier 고정 조건에서 count 감소의
이론적 ceiling은 약 `1.9174`배다. 이 ceiling은 wall time에 대한 상한이 아니다.
AMR 수준·timestep·gravity 비용이 바뀌면 runtime 감소는 입자 count와 다르다.

기존 `mass_cut_refine=2.3283064365e-10`은 level-11 mass의 약 두 배다.
생산 `rho_fine.kjhan.f90`는 더 무거운 DM을 입자 수 기반 refinement 계산에서
제외한다. 이 cutoff를 재사용하면 level-10/9 후보의 최고 해상도 입자도 제외된다.
cutoff를 finest와 next mass 사이에 두는 예시는 level 10에서
`1.8626451492e-9`, level 9에서 `1.4901161194e-8`이다. 새 effective namelist의
cutoff, occupancy 및 실제 AMR layout를 함께 검증해야 한다. 관련
[RAMSES zoom 문서](https://ramses-fork-nbrucy.readthedocs.io/en/developper_guide/wiki/Music.html)는
목표 Lagrangian region으로 MUSIC IC를 재생성하고 최고 IC level에 맞춰
`mass_cut_refine`을 바꾸는 절차를 설명한다. 최고 IC directory만 빼거나
RAMSES `levelmax`만 낮춰 mass convergence를 얻었다고 주장하지 않는다.

### legacy baryon header와 total-matter DMO 계약

36개 IC header는 모두 legacy `44-byte` 형식이고 `omega_b` field가 없다.
현재 reader는 legacy fallback 뒤 namelist 값을 사용하며, 위 DMO snapshot은
`omega_b=0`이다. MUSIC cosmology의 `Omega_b=0.04`는 density transfer
계산에 관한 값이고, DMO particle mass가 total matter를 대표하는 조건과 구분해야
한다. 교정 전 `init_part.f90`의 GRAFIC particle mass 식은 hydro guard 없이
`2^(-3*level)*(1-omega_b/omega_m)`다. extended header의 일치 조건을 맞추려고
DMO namelist의 omega_b만 0.04로 올리면, 나머지가 같은 code path에서 particle
질량이 `0.87142398`배가 된다. total-matter one-fluid DMO에서 약 12.86%의
질량을 제거할 수 있는 변경이므로 채택하지 않는다. 기존 omega_b=0 런에
12.86% mass deficit이 있었다는 주장은 아니다.

새 IC를 만들기 전에 cosmological baryon metadata와 gas/particle mass split을
구분하는 reader/writer contract, total-matter density 및 velocity transfer
선택을 확인해야 한다. 기존 IC의 header를 임의로 바꾸거나 새 IC generation을
시작하지 않았다. legacy restart를 쓰는 현재 recovery/pilot의 diagnostic 범위는
유지하며, restart 성공으로 이 신규 IC 계약을 인증하지 않는다.

### 후속 DMO reader 교정과 writer 범위 확인

같은 날 후속 작업에서 `patch/{lagRamses,cuRamses}/init_part.f90`의 GRAFIC
mass initialization을 교정했다. one-fluid DMO는 `2^(-3*level)`의 total-matter
particle mass를 사용하고, hydro에서만 기존 `1-omega_b/omega_m`를 곱한다.
hydro mass 식과 legacy DMO `omega_b=0` 결과는 bitwise 유지했다. cosmological
baryon metadata를 가진 새 DMO IC를 쓰면서 total matter를 중복 차감하지 않도록
한다. ASCII/Gadget/restart mass records를 재작성하지 않는다. 기존 checkpoint와
제출된 고정 recovery/pilot bundle을 수정하지 않았다.

`patch/lagRamses/init_time.kjhan.f90`는 extended read의 성공 여부로 header
presence를 판정한다. 기존의 `omega_b_header>0` 판정과 달리, 확장 header에
명시된 0도 namelist와 비교한다. 확장 header 및 namelist의 baryon 값이
비유한 값, 음수 또는 Omega_m보다 크면 거부한다. legacy 44-byte fallback은
호환성을 위해 유지하며 absent-header warning을 내므로 신규 science IC의
확장 header 요구를 대신하지 않는다. 기존 hydro default 0.047 정책은 이
수정에서 바꾸지 않았다.

`tests/sidm/test_dmo_grafic_mass.py`는 표시한 실제 생산 Fortran block을
추출해 컴파일한다. 각 patch에서 hydro on/off, 네 Omega_m, 네 baryon fraction,
네 IC levels의 128 cases를 검사한다. hydro의 이전 식 및 legacy DMO 결과는
bitwise 같고, DMO의 full-box mass는 baryon metadata와 무관하게 1이다.
별도 driver는 실제 unformatted 44/48-byte native header를 쓰고 실제 read,
fallback 및 metadata gate를 실행한다. zero-header mismatch, positive mismatch,
NaN, 음수, 상한 초과와 tolerance 경계 조건을 포함한다. GNU 13.2와 Intel
ifx 2025.3 각각에서 256 mass cases 및 header 14 cases가 통과했다.
검증 범위는 compiled blocks이며 cosmological cold-start, MPI collective,
모든 GRAFIC payload와 generator/reader end-to-end 검증은 아직 아니다.

검증 artifact는
`/gpfs/kjhan/Hydro/Sidm/validation/dmo_reader_20261013.q28kvO/`다. source SHA,
test SHA, compiler flags/version을 `gnu_blocks.json`, `ifx_blocks.json`에
보존했다. source archive는 baseline `eadfc9f`에 기존 timestep 진단 네 파일과
reader 교정 세 파일만 복사했다. 별도 archive의 전체 HDF5+FFTW serial build가
성공하고 FFTW linked-code gate도 통과했다. binary SHA256은
`fc3eed35aafad343f2df40415601efdd878a7fdcb83ad3931de0e79bfae7c8ae`다.
실제 FFT dispatch/MG 수렴과 cold-start를 아직 검증하지 않았다. 사용자
worktree의 binary/module 파일을 덮어쓰지 않았다.

기존 IC 생성기에 지정된 `/gpfs/kjhan/Hydro/MUSIC/MUSIC2`는 독립 저장소다.
원격은 `https://github.com/cosmo-sims/MUSIC2.git`, 조회 revision은
`0aa769debbe786bf2d562469abba3826b711abc4`다. main/transfer_eisenstein/
output_grafic2의 세 inspected source는 HEAD와 diff가 없다. 현재 source의
DMO density 선택은 `delta_matter`이고 2LPT velocity 선택은 `theta_matter`다.
Eisenstein plugin은 `tf_distinct=false`, `tf_withvel=false`이며
`(Omega_b/Omega_m)*T_b + (1-Omega_b/Omega_m)*T_c`의 analytic total fit을
반환한다. 이 backend는 별도 Boltzmann velocity transfer를 제공하지 않으므로
analytic fit과 성장 법칙 범위로 해석한다. 현재 경로가 과거 생성 binary와
일치하거나 과거 IC payload의 transfer를 인증했다고 주장하지 않는다.

- 현재 main source SHA: `2a3dc446e18745f5d70383b9ea3a9f8966c9f27402e694b54028d65cc6641c12`
- 현재 grafic2 writer SHA: `cbbc95f9f70dd787b204433d52aa59dbc285f0d9f0c09069fea2c2b4602f83c0`
- 현재 Eisenstein source SHA: `de1be6fe9d3c49e0840e1dd4df998f08bcb4aabb2d150a02aaf347080cbc8a8f`
- 현재 MUSIC binary SHA: `11a0d7a54b52d0ea227345917a72f01e8e26b011b61a655b8c06187cb22da75b`

현재 writer의 header struct는 44 bytes이고 cosmological Omega_b를 쓰지 않는다.
새 science IC에는 Omega_b를 포함하는 48-byte writer와 corrected reader의
small end-to-end gate가 필요하다. 독립 MUSIC2 저장소의 writer 변경은 승인
대상으로 분리한다. 해당 저장소를 수정하거나 새 IC를 생성하지 않았다.
승인 범위는 writer 교정과 소규모 검증이며 대형 IC generation/production
launch, 원격 commit/push 및 Overleaf 변경을 포함하지 않는다.

### 수렴 실험 크기에 대한 후속 판단

입자 질량 수렴의 우선 후보는 같은 target/phase를 가진 최고 IC level 10과 11의
비교다. level 9를 유일한 저적색편이 과학 결과로 대체하지 않는다. level 10의
약 8배 무거운 입자 mass는 host-scale benchmark 후보지만, core 내부의
입자 수·resolved radius·contamination·timestep 및 AMR layout 검증을 대신하지
않는다. 기존 halo에서 단순 eight-to-one scaling으로 남는 입자 수를 새 halo의
측정값으로 제시하지 않는다.

[Mace et al.의 SIDM 수렴 시험](https://arxiv.org/html/2402.01604), §V는
isolated core-collapse 시험에서 `N<1e5`의 realization noise와 일부 `N=1e5`
case의 미수렴을 보고한다. 이 숫자를 모든 cosmological core 형성이나 RAMSES
host/subhalo의 보편적인 합격 기준으로 사용하지 않는다. 저질량 subhalo의
collapse 주장은 실제로 더 높은 mass resolution과 timestep 시험이 필요하다.

다음 실행 판단은 recovery gate와 SIDM timestep pilot 뒤에 한다. 신규 IC 계약과
phase/mask 검증을 먼저 완료하고, level-10 후보의 짧은 same-epoch benchmark로
비용을 측정한다. level 10/11 mass 비교와 AMR cap 비교를 구분하고, science
profile 결과가 필요한 epoch를 유지한 채 90일 예산을 평가한다. 아직 3개월 내
완료를 보장하거나 production matrix를 제출하지 않았다.

header format/endian, bounded read, corruption, cached tier/count normalization,
count-scenario 범위, unit conversion, saved namelist mismatch 검사의 회귀 시험
7개가 `grammar`에서 통과했다. 실제 audit의 첫 시도는 float32 IC 값과 출력
문자열의 exact equality에서 정지했다. 문자열 반올림만 허용하는 상대오차
`2e-13` 비교로 고쳤고 실제 감사가 통과했다. 초기 로그는 보존했다. 최신
artifact는 `/gpfs/kjhan/Hydro/Sidm/validation/zoom_resolution_20261013.8b8GTo/`의
`resolution_audit.json`, `unit_tests.log`, `audit.log`다. JSON은 header-only
SHA256과 전체 payload 인증을 구분하고, 새 science seed/production 승인을
false로 유지한다. 원고와 Overleaf는 변경하지 않았다.

## 1. 이관 결과와 저장소 구조

SIDM 관련 업무의 로컬 기준 위치를 다음과 같이 통합했다.

- 프로젝트 루트: `/home/kjhan/BACKUP/lagRamses-SIDM`
- 코드: 프로젝트 루트 자체가 LagRamses의 `sidm` Git worktree이다.
- 코드 기준점: `main`의 `0373ef7`에서 `sidm` 브랜치를 만들었다.
- 논문: `paper/`에 독립 Overleaf 저장소를 유지한 채 Git submodule로 등록했다.
- 논문 원격: `https://git@git.overleaf.com/6a64afa4af99ca2536a1f4c6`
- 논문 기준 커밋: `b40372c`, 2026-08-17 공저자 Overleaf 수정까지 포함한다.
- 기존 논문 경로 `/home/kjhan/paper_sidm_overleaf`는 새 `paper/`를 가리키는
  호환 심볼릭 링크로 남겼다.
- Lageunha 원시 데이터: `/gpfs/kjhan/Hydro/Sidm`
- 로컬 실행 메타데이터: `runs/lageunha/snapshot-2026-08-24`

원시 시뮬레이션 출력과 초기조건은 대용량이며 재시작의 기준 자료이다. 해당
자료는 GPFS에 그대로 보존했다. 로컬 스냅샷에는 namelist, 실행 및 분석 스크립트,
JSON 진단, 종료 코드, job-control 상태, 그리고 각 출력의 `info_*.txt`만 포함했다.

## 2. 서버 운영 원칙

1. 모든 작업을 시작하기 전에 `hostname`으로 현재 서버를 확인한다.
2. `syntax`에는 CPU 계산을 올리지 않는다. 컴파일과 검증은 `grammar`에서 한다.
3. Lageunha에서는 실행 중인 다른 작업과 CPU 및 NUMA 배치를 먼저 확인한다.
4. 원격 런을 수정하기 전에 PID, 작업 디렉터리, 최신 로그 시각, 체크포인트 구성,
   `jobcontrol.txt`, launcher return code를 함께 확인한다.
5. 로컬 메타데이터 스냅샷은 읽기용 기록이다. 재시작은 반드시 GPFS의 원본
   출력에서 수행한다.

2026-08-24 확인 당시 Lageunha에서 실행 중인 64-rank RAMSES 계산은 SIDM이
아니라 VoidSim의 `compact726` 런이었다. 활성 SIDM 프로세스는 없었다.

## 3. 구현 범위와 검증 현황

주요 SIDM 연산자는 `patch/cuRamses/sidm_scatter.f90`에 있다. 관련 매개변수,
호출 순서, 시간 간격 제한, 입자 상태, 냉각, 그리고 load balancing은 다음 파일에
분산되어 있다.

- `patch/cuRamses/amr_parameters.jaehyun.f90`
- `patch/cuRamses/read_params.jaehyun.f90`
- `patch/cuRamses/amr_step.jaehyun.f90`
- `patch/cuRamses/newdt_fine.kjhan.f90`
- `patch/cuRamses/pm_commons.f90`
- `patch/cuRamses/dark_cooling_mod.f90`
- `patch/cuRamses/dark_cooling_fine.f90`
- 대응하는 `patch/lagRamses/` 파일
- `patch/lagRamses/load_balance.kjhan.f90`와 `patch/cuRamses/bisection.f90`

현재 구현은 constant, velocity-dependent, anisotropic, inelastic, 그리고
dissipative self-interaction을 지원한다. 입자와 baryonic gas 사이의 drag 및
atomic dark matter 냉각도 포함한다. 핵심 구현 커밋은 다음과 같다.

- `4f77365`: SIDM 단위 변환, 입자 수 clamp, iSIDM channel, drag energetics 수정
- `3268da9`: iSIDM 미시적 운동학과 논문 검증 절 추가
- 이후 main의 load-balance 및 process-time 수정도 현재 worktree에 포함된다.

논문에 기록된 검증 결과는 다음과 같다.

- 과거 원고에 기록된 Monte Carlo event ratio는 `a=0.3`에서 `1.07 +/- 0.06`,
  `a=0.6`에서 `1.006 +/- 0.044`였다. 현재 paper bundle에는 원시 count,
  epoch별 expectation, 표본 집계 방식과 uncertainty 정의가 없어 검증할 수
  없으므로, 2026-10-22에 원고와 그림에서 이 수치 주장을 보류했다. 원자료를
  복구하면 다시 검증할 수 있다.
- elastic event의 운동량과 운동에너지는 machine precision에서 보존된다.
- `delta=10 keV`, `m_chi=1 GeV`인 two-state test는 up-scattering threshold와
  de-excitation energy를 재현한다.
- dark-matter--gas drag smoke test는 결합 운동량과 에너지를 보존하며 음의 gas
  internal energy를 만들지 않는다.
- ADM implicit cooling은 2,000-step 결과가 200,000-step explicit reference와
  `1.3e-5`의 상대 차이로 수렴한다.
- 신규 ADM 실행은 `adm_T_init` (기본 `1 K`)에서 particle internal energy를
  명시적으로 초기화한다. 이 값은 과학 실행에서 반드시 의도한 thermal initial
  condition으로 설정해야 한다. HDF5 checkpoint는 이를 보존하지만 legacy binary
  particle backup은 `edp`를 저장하지 않으므로 ADM binary restart는 명시적으로
  거부된다.
- `tests/adm/test_adm_leaf_density.py`는 leaf 내 macro-particle 질량 합과 proper
  volume의 density 변환을 검사하고,
  `tests/adm/test_adm_initial_temperature.py`는 신규 초기화와 restart guard를
  검사한다. `tests/adm/run_adm_amr_leaf_temperature_smoke.sbatch`는 1, 2, 4, 8
  particle/leaf에서 총 leaf mass를 고정한 runtime smoke다.
  2026-08-26의 job `324433`은 HDF5 출력에서 각각 4,096, 8,192, 16,384, 32,768
  macro-particle과 동일한 총질량 `4.096e-3`을 확인했고, fully refined cooling
  step 뒤의 질량가중 평균 온도는 모두 `941.7222186--941.7222188 K`였다.
- `tests/adm/run_adm_hdf5_restart_smoke.sbatch`는 연속 실행과 HDF5 restart를
  particle ID별로 비교한다.  2026-08-26의 job `324439`은 16,384 particle에서
  `dark_energy_int`와 `dark_h2_frac`가 bitwise 동일함을 확인했다. 위치와 속도는
  restart 후 force reconstruction의 순서 차이만 남았고 최대 상대 차이는
  `6.499e-9`였다.
- `tests/adm/adm_thermal_history.py`는 dark radiation과의 kinetic decoupling
  뒤에 $T_D\propto a^{-2}$라고 두고 `adm_T_init` 후보를 계산한다. 이 도구는
  decoupling redshift `z_kd`를 미시 파라미터만으로 추정하지 않는다. 그 값에는
  별도의 dark recombination 및 kinetic-coupling 계산이 필요하다. Production zoom의
  $z_{\rm init}=100$, $\xi=0.5$에서 job `324452`는
  `z_kd=100,300,1000,3000,10000`을 각각 `137.64,46.18,13.89,4.63,1.39 K`로
  변환했다. 다섯 AMR runtime case는 모두 완료했고 두 step 뒤의 온도 변화는
  $8\times10^{-16}$ 이하의 round-off 수준이었다. 현재 $1\,\mathrm{K}$ floor는
  이 가정에서 `z_kd=1.39e4`에 해당한다.

직접 관련된 테스트와 입력은 `bin/test_darkcool.f90`,
`tests/load_balance/preflight_sidm_work.nml`,
`tests/load_balance/preflight_sidm_ab.nml`, 그리고
`tests/load_balance/test_domain_leaf_cost.f90`에 있다. 논문 저장소에는
`test_adm_timestep.f90`와 figure 재생성 스크립트가 있다.

## 4. Paper-I 줌 시뮬레이션 현황

세 런은 동일한 zoom initial condition을 사용하고 refinement를 `levelmax=18`까지
허용한다. Collisionless CDM 기준 런은 이미 `z=0`에 도달했다.

| 런 | 물리 모형 | 최신 출력 | scale factor | redshift | MPI | 상태 |
|---|---|---:|---:|---:|---:|---|
| CDM | collisionless LambdaCDM | `output_00009` | 1.002018 | 약 0 | 32 | 완료 |
| SIDM1 | `sigma/m=1 cm^2/g` | `output_00010` | 0.344185 | 1.905 | 32 | 정상 정지, `rc=0` |
| SIDM3 | `sigma/m=3 cm^2/g` | `output_00009` | 0.263555 | 2.794 | 24 | 정상 정지, `rc=0` |
| SIDM0.3 | `sigma/m=0.3 cm^2/g` | 없음 | 없음 | 없음 | 예정 | 시작되지 않음 |

SIDM1과 SIDM3는 2026-08-02 23:07 KST에 `jobcontrol` 요청을 받아 coarse-step
경계에서 체크포인트를 쓰고 종료했다. 두 launcher return code는 모두 0이며 최신
로그에 `Run completed`가 있다. 최신 출력의 파일 수는 SIDM1이 101개이고 SIDM3이
77개이다. 구성은 각 MPI rank의 AMR, gravity, particle 파일과 공통 metadata로
일관된다. 별도의 non-empty error log나 치명적 오류 기록은 발견되지 않았다.

마지막 12개 coarse step의 평균 시간은 SIDM1이 약 367.0초이고 SIDM3이 약
366.3초였다. SIDM1에는 434.5초의 load-balance 인접 step이 하나 포함되어 있다.

### 재시작 전 반드시 고칠 항목

live namelist의 `nrestart`가 최신 체크포인트보다 뒤처져 있다.

- SIDM1의 현재 값은 `nrestart=7`이다. 재시작 값은 10이어야 한다.
- SIDM3의 현재 값은 `nrestart=8`이다. 재시작 값은 9여야 한다.
- 두 `jobcontrol.txt`에는 현재 `0 -1`이 남아 있다. 실행을 확정하기 전에는
  지우지 않는다.

원격 경로는 다음과 같다.

- SIDM1: `/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm1`
- SIDM3: `/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_run_sidm3`
- 공통 zoom IC: `/gpfs/kjhan/Hydro/Sidm/Agn/Zoom0/zoom_l11/ic_zoom_l11`

기존 실행 배치는 SIDM1이 NUMA1에서 32 MPI ranks x 2 OpenMP threads이고 SIDM3가
NUMA0에서 24 MPI ranks x 2 OpenMP threads이다. SIDM3 executable은 SIDM1의
파일과 연결되어 있었으므로 새 binary를 배치할 때 shared target을 한 번만 바꾸고
기존 binary를 먼저 보존해야 한다.

## 5. Run0 시뮬레이션 현황

Run0은 `128 h^-1 Mpc`, `1024^3` dark-matter particles와 초기 gas elements를 쓰는
hydrodynamic cosmological suite이다. 모든 Run0 모델은 중간 적색편이에서 정지해
있으며 `z=0` 결과가 아니다.

| 런 | 최신 출력 | scale factor | redshift |
|---|---:|---:|---:|
| `run_cdm` | 16 | 0.208365 | 3.799 |
| `run_dsidm` | 17 | 0.210738 | 3.745 |
| `run_idm` | 17 | 0.207156 | 3.827 |
| `run_isidm` | 19 | 0.212623 | 3.703 |
| `run_isidm_endo` | 20 | 0.237594 | 3.209 |
| `run_sidm` | 19 | 0.236877 | 3.222 |
| `run_sidm5` | 16 | 0.203728 | 3.908 |
| `run_sidm10` | 17 | 0.211781 | 3.722 |

Run0 자료는 구현 검증과 초기 matched-halo 진단에 사용됐지만 현재 Paper I의
저적색편이 core 측정을 완성하지 않는다.

## 6. 논문 현황

논문 주 저장소는 `paper/`이다. 2026-08-17 Overleaf의 최신 공저자 커밋은 제목의
`cuRAMSES AMR code`를 `Cosmological AMR code`로 바꿨으며 이관 전에 fast-forward로
반영했다. 논문 저장소는 이관 시점에 clean 상태이다.

현재 원고는 다음 결과를 포함한다.

- SIDM 및 ADM 연산자와 runtime parameter 설명
- Monte Carlo rate, conservation, iSIDM threshold, drag, ADM cooling 검증
- collisionless LambdaCDM zoom의 `z=0` density map
- 공통 초기 epoch 부근의 임시 matched SIDM profile

원고의 `z` 약 4 profile은 large-scale matching과 초기 수치 진단용이며 물리적인
SIDM core 측정이 아니다. 다음 주 그림은 아직 완성되지 않았다.

1. 공통 `z=0`에서 CDM, SIDM1, SIDM3의 중심 projected-density map
2. 같은 epoch의 density, SIDM/CDM ratio, one-dimensional velocity dispersion,
   그리고 가능하면 circular velocity profile
3. `z=4,3,2,1,0`에 걸친 core와 cumulative scattering의 진화
4. force 및 mass resolution에 대한 core-size convergence
5. 부록의 `sidm_courant` 또는 `P_max` timestep convergence

대규모 Mpc density map은 CDM과 거의 구별되지 않으므로 SIDM science figure로
추가하지 않는다. SIDM 비교 map은 resolved halo core에 집중한다.

## 7. 물리적 한계

ADM 부분은 아직 full cosmological validation을 거치지 않았다. 현재 모형은
macro-particle에 internal energy를 부여하고 leaf cell 안의 eligible
macro-particle 질량 합으로 density를 추정하며 Saha ionisation equilibrium을
사용한다. Dark pressure, photo-heating, resolved dark shock, 그리고 fully
coupled chemical network는 포함하지 않는다.
따라서 fragmentation과 dark-disc vertical structure는 검증 범위 밖이다.
`adm_T_init`의 과학적 선택과 cosmological thermal history의 일관성도 아직
별도 검증이 필요하다. Job `324452`가 탐색한 physically labelled initial
temperatures는 atomic cooling을 전혀 유발하지 않았다. Source audit에서 `edp`는
초기화, checkpoint와 particle migration, 그리고 `dark_cooling_fine`에만 연결된다.
Dark pressure, adiabatic compression, shock heating, 그리고 Compton heating이
없으므로 이 particle-only implementation에서 그 초기온도 후보들은 gravitational
zoom 결과를 바꾸지 않는다. 따라서 현 구현으로 ADM cosmological zoom을 science
production으로 제출하지 않는다. 의미 있는 다음 단계는 dark kinetic thermal history와
pressure 또는 energy-coupling prescription을 설계하고 별도 검증하는 것이다.
2026-08-26에 `adm_adiabatic`과 `adm_T_floor`를 추가했다. Cosmological ADM run은
기본적으로 step 시작과 끝의 Friedmann scale factor 사이에서
$u_D\propto a^{-2}$를 적용한다. `adm_T_floor`는 기본 `1 K`를 보존하지만 더 이른
kinetic decoupling history를 위해 양의 값으로 낮출 수 있다. Compton term은 이제
signed net rate로 계산하고 implicit solve가 heating branch도 bracket한다. 이 변경은
dark pressure force나 $P\,dV$ work를 아직 제공하지 않는다. 그 둘은 AMR 경계를
넘는 보존적 dark-fluid solver 또는 명시적으로 근사한 HPM closure를 선택한 뒤에만
추가한다.

2026-08-26에 그 HPM closure의 첫 구현을 추가했다. `adm_hpm=.true.`이면 leaf cell에
$P_D=(\gamma_D-1)\rho_D e_D$를 NGP로 deposit하고, 중앙차분
$-\nabla P_D/\rho_D$를 기존 particle-mesh force에 더한다. pressure working field는
새 대형 배열을 만들지 않고 이미 예약된 `rho_star`를 재사용한다. 이 초기 경로는
모든 collisionless macro-particle을 ADM으로 취급하므로 `poisson=T`, `star=F`,
`sink=F`에서만 허용된다. `adm_hpm_courant`는 $c_s^2=\gamma_D(\gamma_D-1)e_D$의
sound-crossing CFL을 추가한다. `tests/adm/run_adm_hpm_unit.sh`는 pressure law와
gradient sign을, `tests/adm/test_adm_hpm_wiring.py`는 force/KDK/CFL 연결을 검사한다.
독립 실행 smoke (`tests/adm/run_adm_hpm_gradient_smoke.sbatch`)는 8-cell periodic
pressure gradient에서 HPM on/off의 최대 $x$-velocity 차이
$2.586\times10^{-7}$와 총 momentum kick 차이 $7.024\times10^{-19}$를 확인했다.
이는 force-operator 검증일 뿐, production cosmological result는 아니다. 이 HPM은
Riemann flux, shock capture, compressional $P\,dV$ heating, 그리고 full chemical
network를 풀지 않는다. 따라서 dark-disc thickness나 fragmentation의 정량 예측에는
사용하지 않고, 우선 pressure-suppression 및 controlled convergence 연구에 한정한다.

## 8. 권장 작업 순서

1. Lageunha의 자원 가용성과 다른 작업의 NUMA 배치를 확인한다.
2. `grammar`에서 SIDM binary를 재현하고 smoke 및 load-balance gate를 통과시킨다.
3. GPFS 원본의 최신 checkpoint를 rank별로 다시 검사한다.
4. live namelist의 `nrestart`를 SIDM1은 10, SIDM3은 9로 수정한다.
5. 기존 executable을 보존하고 검증된 binary를 한 번만 배치한다.
6. 실행을 확정한 직후에만 `jobcontrol.txt`를 비우고 기존 launcher로 재시작한다.
7. PID, 최신 log timestamp, fine/coarse step 증가, memory, `Pmax`, 그리고 error
   log를 시작 직후 확인한다.
8. 두 런이 공통 output epoch에 도달할 때마다 profile 및 scattering diagnostic을
   생성한다.
9. `z=0` 후 Paper-I figure 3부터 6까지 교체하고 convergence suite를 수행한다.

이 순서는 상태 조회와 문서 작업에 관한 승인만으로 원격 namelist나 실행 상태를
바꾸지 않는다. 실제 재시작은 별도 실행 지시를 받은 뒤 수행한다.

## 9. Git 사용

코드 저장소와 논문 저장소는 서로 독립적이다.

```bash
git -C /home/kjhan/BACKUP/lagRamses-SIDM status -sb
git -C /home/kjhan/BACKUP/lagRamses-SIDM submodule status
git -C /home/kjhan/BACKUP/lagRamses-SIDM/paper status -sb
```

논문을 먼저 commit하고 Overleaf에 push한 뒤 상위 `sidm` 브랜치의 submodule
pointer를 갱신한다. 코드와 인수인계 문서는 LagRamses GitHub 저장소의 `sidm`
브랜치에서 관리한다. 이관 시점에는 로컬 브랜치만 만들었으며 원격 push는 하지
않았다.

## 10. 2026-10-22 finite-time crossing endpoint audit

Astra가 지적한 유한 시간창의 두 endpoint 결함을 승인된 범위에서 재확인했다.
검증 디렉터리
`/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/`의
`eddington_finite_time_crossing_diagnostic.py`에는 이미 입력 float 표현값에 대한
정확한 `Fraction` 기반 crossing count와, 실제 초기 inward-boundary crossing만
한 주기 뒤로 넘기고 그 밖의 초기 위상 모호성은 거부하는 delay helper가 있었다.
테스트 파일에는 인접 float 경계, 비유한 clock, 초기 경계, 양의 지연 회귀 테스트도
이미 추가되어 있어 이번 승인 실행에서 소스 파일은 수정하지 않았다.

단일 프로세스로 다음만 실행했고 **25개 테스트가 1.43초에 통과**했다:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m pytest -q -p no:cacheprovider \
  test_eddington_finite_time_crossing_diagnostic.py
```

검증 시 SHA-256은 진단 코드 `bbd44c24562644fce9e2e650d9a70c08cc94b0a169f206e63b2f9908b9618ac4`,
테스트 `05300a1e7e771046c7b455f345fd56783a715ef41652becbb0e59b9fcba121fa`였다.
보호 대상 v3 산출물 `eddington_finite_time_crossings_N1048576_seeds20261021_28_v3.json`
은 읽기만 했으며 SHA-256은
`432d12934a4d073f30a378ecd4756bc15acc3edb61a97d9aa2b31fcb5eebf52b`이다.
이번 작업에서는 해당 JSON을 수정하지 않았고 MPI/Slurm, 대규모 진단, 새 시뮬레이션,
commit 또는 push를 수행하지 않았다.

## 11. 2026-10-22 Objective 1 scientific and implementation handoff

현재 물리 질문은 기존 열저장고 에너지 예산을 만족하는 계에서조차 유한한 열수송
시간이 살아 있는 SMBH 쌍성의 경화를 제한하는지 여부다. 독창성 주장은 단순한 live-halo
backreaction, SIDM-altered dynamical friction, 또는 core stalling 감소로 두지 않는다.
ACD식 쌍성 에너지 주입과 실제 SIDM 산란을 연결하고, 수렴된 binary hardening 및
시간·반경 분해 에너지 재분배와 encounter supply를 통해 인과성을 보여야 한다.

선행연구 경계:

- ACD, *Phys. Rev. Lett.* **133**, 021401 (2024), DOI `10.1103/PhysRevLett.133.021401`는
  쌍성 에너지 예산, SIDM core reservoir, 빠른 재열화와 spike 보충을 다룬다.
- Sabarish et al. (2025), arXiv `2505.14779`는 단일 SMBH SIDM spike의 N-body,
  accretion 및 heat-conduction calibration을 다룬다.
- van den Bosch & Dattathri (2026), arXiv `2511.14912`는 live SIDM halo 속 단일
  massive perturber, core stalling 감소 및 gravothermal response와 겹친다.
- Tiruvaskar & Gordon (2025), arXiv `2506.18153`는 SMBH-binary population의 SIDM-spike
  drag 및 PTA fit을 다룬다.

따라서 이들과의 구별점은 binary-period-resolved heating과 SIDM 에너지 전달·수송,
실제 경화율의 수렴 및 대조군이어야 한다. ACD Fig. 3의 축은
`ΔE_SMBH orbit / U_SIDM core`이며 Eq. 7은 쌍성 궤도 에너지 손실을 정의한다. 다만
Eq. 8의 spike binding energy와 isothermal-core reservoir의 적분 범위·정규화 대응은
곡선 수준 재현 전까지 완전히 확정하지 않는다. 그때까지는 “ACD-inspired” test라
부르고, 이 참고문헌 해석 문제만으로 feasibility 판단을 지연하지 않는다.

현재 최소 hard interval 후보는 기하학적으로 `L >= 19`를 요구하지만, 실제 endpoint
encounter sample과 cell-local SIDM scattering 통계는 아직 없다. 정적 pericenter
확률은 유한시간 encounter flux가 아니다. 기존 PM force 비교도 고정 격자 force-field
범위에 한정되어 particle gather/kick, 반복 reset, AMR 또는 binary hardening 수렴을
검증하지 않는다.

다음 구현 단계는 source-level, operation-resolved kinetic ledger다. Claude Opus 5가
fresh-start, box-interior, uniform-grid, SIDM-off, fixed-mass two-cloud 조건에서 실제
sync half-kick, `vsink` projection, move half-kick, drift 및 averaging/recreation의
순서와 각 단계의 kinetic-energy delta를 기록한다. delta 합은 같은 phase의 endpoint
`ΔK`와 비교하되 이를 전체 energy-conservation proof라 부르지 않는다. 이어 MPI startup
probe를 gate 시작부에 합친 단일 bounded executable gate로 같은 입력의 1-rank/2-rank
결과만 비교한다. 별도 MPI probe, 구 inventory-only gate 또는 별도 SIDM preflight는
새 의무로 추가하지 않는다. 이 gate 통과 후 target-scale force/sampling/runtime·storage
feasibility를 측정한다. 코드 구조 감사·구현은 현 project `AGENTS.md`의 Claude Opus 5
역할, 계획·평가·검증은 Astra 역할을 따른다. 이 gate의 MPI/Slurm 실행은 아직 승인되지
않았으며, 별도 명시 승인 전에는 제출하지 않는다.

JCAP 방법론 원고는 검증 범위에 맞춰 `f_diss=0.2`의 단일 homogeneous relaxation
사례만 표기하고, `A = ρ(σ/mχ)v` 및 총 운동에너지 `E`를 정의했다. 측정 최대
ensemble-mean standard error는 `0.013259206`, floor-moment 보정 상한은
`7.79013e-7`이며, 이 보정값을 전체 velocity-distribution 오차로 해석하지 않는다.
노이즈가 있는 timestep sequence는 결정론적 수렴 차수를 입증하지 않는다. 논문은
새 임시 디렉터리에서 22페이지로 빌드했고 unresolved references/citations 및
overfull boxes가 없었으며 `git diff --check -- main.tex`도 통과했다. 이 변경은 기존
dirty paper submodule 안에 남아 있고 commit 또는 Overleaf push는 하지 않았다.

## 12. 2026-10-22 operation-resolved cloud-ledger gate 및 원고 증거 정리

§11의 “구현 단계/Slurm 승인 대기”는 작성 당시 상태로, 아래 완료 상태로
대체한다. 사용자는 구현 역할을 Codex로 이관하고 grammar `normal`의 bounded
CPU/MPI gate를 승인했다. Astra 설계 감사는 재구현·재제출을 불필요하다고 판단했다.

최신-main 기준은 `/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/latest_main`
(branch `codex/sidm-cloud-ledger`, HEAD `9bec08cd5e3fa570d881a86da067b8fcfbb6e37`)다.
완료 계약은 `run_cloudledger_kinetic_normal_20261022/RUN_CONTRACT.md`에 있고,
build `1206413`, 1-rank job `1206414`, 2-rank job `1206415`가 모두 `normal`
partition에서 exit `0:0`으로 완료됐다. 재제출은 하지 않았다. 현재 네 계측 소스는
build manifest의 SHA256과 모두 일치하고 n1/n2 각각의 `checksums.sha256` 전체가
검증됐다. 재실행하지 않은 보존 출력만 분석했다.

fresh-start, uniform L7, sink cloud당 2,109개 고정질량 입자, SIDM/hydro off,
OMP=1, `nstepmax=2` 조건에서 event 0·1·2와 두 closed interval이 확인됐다. 첫
sync의 `dt=0`과 `ΔK=0`은 의도한 fresh-start 경계다. 각 run은 RAMSES 4/6초와 ledger
4,324,867/4,327,562 bytes를 사용했고 full dump를 만들지 않았다. n1/n2의 저장된
kinetic/reset/operation 값 최대 차이는 `9.49e-20` code-energy이며 sink endpoint
최대 운동에너지로 정규화하면 `5.20e-15`다. 네 interval의 최대 절대 closure residual은
약 `3.27e-19` code-energy다. drift 항은 코드가 위치만 갱신하는 계약에 따라 `0`을
기록하므로 독립적인 drift-energy 계측으로 부르지 않는다.

이 gate는 source-level operation accounting의 plumbing·자기일관성 검증이다. 잘못된
force도 닫힐 수 있어 독립 force 정확도나 orbit convergence를 증명하지 않는다.
AMR/subcycling, 다중 thread, SIDM 산란/열수송, distributed FFT, hardening 또는
물리적 전체 에너지 보존은 검증되지 않았다. 다음으로는 보존 출력의 비교와 target-scale
force/encounter/runtime/storage feasibility를 설계하되, 이 gate만으로 production
실행을 시작하지 않는다.

보존된 strict validator의 실행 SHA는 manifest에 그대로 남아 있다. 그 validator의
마지막 안내 문구가 strict 모드에서도 boundary-only라고 잘못 출력하는 것을 확인해
현재 worktree의 문구만 조건부로 고쳤다. 기존 실행 파일·로그·manifest는 수정하지
않았다. current source의 boundary 23개 및 kinetic 11개 synthetic regression이
통과했다.

JCAP 원고에서는 원시 event count와 uncertainty 정의를 찾을 수 없는 과거 두 epoch
event-ratio 수치와 오차막대 주장을 보류하고, 검증 가능한 leaf-local estimator/MPI
시험으로 범위를 제한했다. 그림을 두 패널로 다시 생성했다. hom-relaxation nested 행
표현 및 README의 실행 경로·축약 데이터 한계도 교정했다. 논문은 임시 경로에서 24쪽으로
빌드됐고 원고 디렉터리의 `main.pdf`는 덮어쓰지 않았다. 코드·논문 commit/push 및
Overleaf 반영은 하지 않았다.

## 13. 2026-10-26 L10 force-probe 실행 제출 상태

준비된 L10 live-binary + active-SIDM force-probe bundle은
`/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/run_binary_forceprobe_l10_20261026/`
이다. 사용자는 이 bounded 검증을 Grammar의 Slurm `normal` 파티션으로 실행하도록
승인했다. 이번 제출 시도에서 현재 로그인 노드 `syntax`의 `sinfo -s`에는 `normal`
파티션이 없었고, Grammar alias `lageunha` 접속은 `Connection closed by 10.0.200.25
port 22`로 끝났다. Bundle의 `bash submit.sh`는 파티션 확인 단계에서 exit 2로
정상 거부됐다. 따라서 `sbatch`는 호출되지 않았으며 build 또는 RAMSES 작업은
제출·실행되지 않았다. 기존 산출물은 변경하지 않았다. `syntax`의 GPU 파티션에
대체 제출하지 않았다.

실행 계약상 clean build job 할당 한도는 2시간, 진화 job 한도는 1시간이다. 이는
대기시간이나 실측 소요시간이 아니다. 유사한 27-step L10 integration gate는
144초였으나 force-probe 계측 오버헤드는 아직 측정되지 않았으므로 이번 simulation의
실행시간은 미확정이다. 현재 남은 조치는 Grammar `normal` 로그인 접근이 복구된 뒤
해시·입력 검사를 거쳐 해당 bundle의 `bash submit.sh`를 실행하는 것이다.

## 14. 2026-10-26 JCAP 증거 패키지와 남은 원고 주장

paper submodule은 HEAD `d31a0779c9db71e943f8d42933251c686d4c8e36`에서 여러 원고·그림·
분석 파일이 수정된 상태다. 아래 검증 자료는
`paper/data/README.md`의 claim-to-evidence 표와 함께 uncommitted로 보존한다. 코드
저장소 HEAD는 `01b6ead49be96d97a6d3366dfffe4fde36726fe0`이며 각 bundle은 실행 당시의
dirty diff와 소스 스냅샷을 별도로 pin한다. 아직 immutable public release, commit 또는
Overleaf push는 만들지 않았다.

- `data/pair_probability_regression_20261026/`: production pair-probability 함수의
  compiled Fortran regression은 `N_c=2…129`의 128개 점유도를 double-precision epsilon
  이내로 통과했다. GNU Fortran 13.2, 명령, 결과, 소스 스냅샷과 SHA-256을 보존했다.
- `data/sidm_angular_validation_20261026/`: Rutherford angular helper는 다섯 조절값에서
  각 50,000 표본의 CDF와 세 수송 모멘트 검사, 두 기저 분기의 회전 200건을 통과했다.
  독립 적분 7건도 통과했다. Yukawa 부분은 source-formula audit와 독립 적분이며 compiled
  scattering run이나 Born 정규화 검증으로 해석하지 않는다.
- `data/profile_epoch_contract_20261026/`: matched-epoch 분석 계약 합성 시험 3건이
  통과했고 수정된 paper analyzer의 별도 repo provenance와 Python 의존성을 pin했다. 실제
  halo IC·중력 provenance나 profile 수렴을 인증하지 않는다.
- `data/adm_atomic_rates_20261026/`: GNU Fortran 13.2와 Intel Fortran 2025.3에서 두
  source variant 각각 978개 coefficient 비교가 통과했다. 최대 차이는 bremsstrahlung
  0.031097%, recombination 0.382823%다. reference-tail 변화는 `9.84e-11`, 2,000 대
  200,000 timestep 상대 차이는 `1.18e-7`이다. 이는 두 coefficient와 고정밀도 local
  integration 결과일 뿐 complete ADM/EOS/chemistry 검증은 아니다.
- `data/adm_hpm_unit_20261026/`: source 기반 closure unit과 wiring regression은 통과했다.
  이상기체 압력, 1D pressure-gradient 부호, uniform-pressure zero-force, sound speed 및
  source call order/field exchange/timestep hook/입력 guard 범위다. periodic AMR gradient와
  MPI executable-rejection smoke 결과는 없다.

Astra가 각 bundle과 metadata를 감사해 checksum 및 기술 범위 불일치가 없음을 확인했다.
Fable이 다듬은 HPM·molecular replacement text를 현재 `main.tex`의 표와 본문에 적용했다.
HPM 주장은 보존된 closure unit과 source-wiring 범위로 좁혔고, molecular fits는 독립
검증을 하지 않았다고 명시했다. 새 초안은 임시 경로
`/tmp/jcap_scope_rebuild_20261026_v3/main.pdf`에서 24페이지로 빌드됐다. 최신 build log에
overfull box와 unresolved reference/citation은 없다. 기존 `paper/main.pdf`는 덮어쓰지
않았다.

출판 전에 남은 주요 주장 간극은 inelastic threshold/33-case MPI와 routine restart
(`main.tex` §abstract, validation table, §iSIDM, conclusion), DM-baryon drag smoke
(`§abstract`, validation table, conclusion), 그리고 식별 가능한 immutable release다.
이 항목은 아직 source/output provenance가 없거나 release bundle에 고정되지 않았다.
현재 JCAP 원고는 여전히 release-ready가 아니다.

## 15. 2026-10-26 L10 force-probe gate 실행 결과

§13에 기록된 syntax 제출 실패는 Grammar 제출 성공으로 대체한다. 제출 전 frozen bundle은 source, input, bundle hash 검사를 모두 통과했다. Grammar Slurm `normal` 파티션에서 build job 1219347은 exit code 0으로 2:31 만에 완료했다. 성공한 build에 종속된 run job 1219348은 exit code 0으로 2:25 만에 완료했다. `srun` 통합 실행은 143초 걸렸고 monitor sample step 27에서 종료했다. 앞선 유사 gate 실행은 144초 걸렸다. 143초와 144초의 차이는 logging overhead를 통제해 측정한 결과가 아니다.

Build batch의 Slurm MaxRSS는 26,363,636 KiB, 약 25.1 GiB였다. RAMSES 단계의 MaxRSS 602,520 KiB는 MPI task 하나의 최대 RSS이며 allocation 전체 메모리 합계가 아니다. Unit 및 source-contract 테스트 46개가 모두 통과했다. 실행 결과는 force generation 28회, force-probe row 112개, central sample 112개, cloud sample 236,208개, force-level fallback row 0개였다. Force-probe row 112개는 두 sink, move와 sync phase 모두, coarse step 0부터 27까지를 포함했다. Paired orbit checkpoint 5개는 initial checkpoint 1개, interior sample step 3개, terminal checkpoint 1개로 구성됐다. `metrics.txt`는 interior sample 3개를 별도로 기록했다. Monitor terminal sample time은 2.000345225e-5였고 ledger closure는 2.070209278e-5까지 이어졌다.

확률 guard는 sampled pair 12,853,776개와 accepted scatter event 454개를 기록했다. 최대 preflight probability는 한도 0.1보다 낮은 0.05246104였고 최대 actual sampled probability는 0.03797227이었다. 28개 interval의 최대 bin residual은 1.8112132734479844e-14 code-energy, 최대 state gap은 0, 최대 global scatter delta는 1.35525271560688e-20 code-energy였다. `run.complete`, checksum manifest 2개, 필수 validator 전부, Slurm exit code 2개가 통과했다.

RAMSES log는 768,787 bytes, force-probe CSV는 113,172 bytes, 실행 파일 사본은 9,497,800 bytes였다. Run directory 총 저장량은 11,270,764 bytes로, IC symlink target과 별도 build tree를 제외한다. Force-probe CSV는 interval당 약 4 KiB였다. Full snapshot은 생성되지 않았다. 저장량은 짧은 진단 gate 비용이며 production run snapshot 비용을 측정하지 않는다.

Gate는 force-source 및 force-use epoch 계측, kick accounting, probability safety, scatter-to-ledger coverage를 검증했다. 실행 조건은 frozen L7-L10 fixture, 입자 1,048,576개, MPI rank 2개, thread 2개였다. 물리적 PM-force bias, cloud-reset error, 궤도 수렴, binary-halo 전체 에너지 보존, 유한 열수송 응답은 측정하지 않았다. Accounting residual과 validator 통과만으로 PM-force bias, cloud-reset error, orbital convergence, total binary-halo energy conservation 또는 finite heat-transport response의 상한을 정할 수 없다. 종료 시점까지의 진화 시간은 초기 궤도 주기의 약 2.8e-5였다.

Astra는 수정된 Grammar `normal` CPU contract와 검증된 matched inputs를 사용하는 bounded target-scale force, sampling, cost screen을 다음 단계로 권한다. Collisionless 조건부터 시작하고 matched short SIDM 구간을 추가한다. 기존 `a10`/GPU L17 script는 수정 없이 재사용할 수 없다. Production 규모를 정하기 전에 applied-force sensitivity, timestep distribution, encounter statistics, active scattering leaves, physical-time throughput, peak memory, diagnostic storage를 측정해야 한다.

## 16. 2026-10-26 hard-binary endpoint pilot 제안

목표는 단일 hard-binary endpoint에서의 공학적 실행 가능성을 확인하는 일이다. 물리적으로 준비된 binary나 열수송 결과를 제시하는 단계가 아니다. 비교 run에는 입자 1,048,576개의 periodic-spherical-mean cusp realization 하나를 사용한다. `M_halo=1`, `M_sink=0.0099`, `q=1`이며 sink 두 개의 질량은 각각 `0.00495`다. Matched run 사이에서 DM bytes와 boost를 동일하게 유지한다. IC는 hash와 row 검사를 통과했으나 주기 경계 평균 potential을 근사적으로 사용한다. PM-cloud regularization 아래에서 평형화된 IC는 아니다.

Candidate hardening radius는 `a_h=0.00111125`다. Endpoint는 `a=a_h/3=3.7042e-4`이고 점질량 주기는 `4.50e-4`, 상대 속도는 `5.17`이다. L7부터 L17까지 cloud radius는 `7.6294e-5`이며 `a/rcloud=4.86`이다. `a/rcloud=4.86`은 필요한 기하 조건만 보여주며 수렴을 입증하지 않는다. L17과 L18을 비교한다.

1. **Stage 1.** Halo와 SIDM이 없는 진공 조건에서 L17 및 L18 force calibration을 `0.1P=4.50e-5`까지 수행한다. 제안된 최대 timestep은 `P/200=2.25e-6`이다. Applied relative acceleration을 순간 Newtonian binary reference와 비교하고 tangential error, sink 자체 가속도와 중심질량 가속도, 해상도 민감도, 짧은 궤도 구간을 측정한다. 임시 admission 기준은 companion-force bias `≤2%`, L17에서 L18로 해상도를 바꿀 때 변화 `≤1%`, timestep 민감도 `<1%`다. 기준값은 아직 측정되지 않았으며 확립된 정확도 주장이 아니다.

2. **Stage 2.** 동일한 DM bytes와 binary sinks를 사용해 L17 live-halo run을 collisionless 조건과 `\sigma/m=1` 조건에서 각각 수행한다. 두 run은 동일한 `0.1P` 구간을 사용하며 step ceiling은 2048이다. Physical-time throughput, probability와 ledger, active mesh, timestep, radial occupancy, eligible-pair fraction, encounter identity와 crossing, scatter counts를 측정한다. 2048-step ceiling 도달은 완료로 인정하지 않는다.

3. **Stage 3.** 앞선 단계가 실행 가능하면 timestep ceiling을 절반으로 낮춘다. Vacuum test로 거동을 설명하지 못하면 live-metric spatial comparison을 추가한다.

CPU pilot의 제안 자원은 Grammar `normal` 파티션의 MPI 2개 × OpenMP 2개, 메모리 16 GiB, run당 2시간이며 GPU는 사용하지 않는다. 별도 clean build가 필요하면 task 1개, CPU 2개, 메모리 48 GiB, 2시간을 제안한다. `ngridmax=1.5M`, `npartmax=2M`, run당 1 GiB reserve도 제안값이며 아직 측정되지 않았다. Capacity와 reserve 적정성을 제출 전에 확인한다. Full snapshot은 쓰지 않고 output clock은 미래로 설정한다. Ledger closure 뒤 monitor stop에 도달한 실행만 성공으로 인정한다. Slurm timeout이나 step ceiling 도달은 미완료다.

기존 `run_l17_timestep_pilot.sbatch`는 `a10`/GPU, 24×2 자원, 인증되지 않은 cosmological `z=3` seed를 사용하므로 수정 없이 재사용할 수 없다. Bundle을 고정하기 전에 exact sink generator와 boost, source 상태와 L17 capability, 실제 cloud support와 sink ownership, force-reference convention, sampling diagnostics, capacity와 timestep semantics를 확인한다. 다음 pilot은 수정된 Grammar `normal` CPU contract와 검증된 matched inputs를 사용해야 한다. 계획은 아직 frozen runnable bundle로 고정되지 않았다.

통과 결과는 PM-consistent halo preparation, matched SIDM monopole, 문서화된 binary preparation, 최소 binary 대 monopole × collisionless 대 active-SIDM 2×2 control experiment으로 이어갈 근거가 된다. Production 계산 전 usable reservoir를 확인하고 random-energy flux를 advection 및 gravity work와 분리해야 한다. Hardening과 flux 정밀도도 충분해야 한다. 3개월 안에 20% 효과를 3σ로 판별할 계산이 불가능하면 method 또는 model을 전환한다. Computational no-go를 physical null로 부르지 않는다. L10 gate만으로 direct repetition이나 production 실행을 정당화할 수 없다.

## 17. 2026-10-26 JCAP 범위 수정과 v4 검토 번들

원고에서 inelastic suite의 MPI 통과 주장을 삭제했다. `test_referee_safety.py`는 1, 2, 4 ranks에서 guard 및 초기화 MPI 사례 33개, HPM consumer 사례 8개, equal-mass kinematics 사례 3개를 정의한다. 실행 출력은 보존되지 않았다. README도 실행 결과를 확인할 수 없는 범위로 수정했다.

Drag 관련 원고는 구현만 기술하며 수치 검증을 주장하지 않는다. Smoke run의 manifest, output, analysis, source provenance는 보존되지 않았다.

v4 검토 PDF는 `/tmp/jcap_scope_rebuild_20261026_v4/main.pdf`에 컴파일했다. PDF는 26쪽, 1,321,066 B다. 빌드 로그에는 overfull box와 미해결 참조 또는 인용이 없다. 비치명적인 hyperref 및 font-substitution 경고는 남았다. `paper/main.pdf`는 덮어쓰지 않았다.

`paper/data/jcap_scope_draft_20261026_v4/` 번들에는 원고 소스와 README, project-local TeX 입력 8개, log와 FLS, paper HEAD, 대상 dirty 상태와 patch, checksum 파일이 포함된다. 불변 release ID와 DOI는 없으며 commit이나 Overleaf push도 하지 않았다. 공개 release 준비는 끝나지 않았다.

## 18. 2026-10-26 Objective 1 Stage 1 force·timestep 결과

Grammar `normal`에서 수행한 진공 쌍성 보정 세 건이 모두 진단 게이트를 통과했다. L17 job 1219630은 149.82초, L18 job 1219701은 280.77초, L17 시간 간격 절반 job 1219884는 285.58초 걸렸다. 진단 산출물 크기는 각각 1,599,163 B, 2,918,711 B, 2,741,726 B였고 full snapshot은 생성하지 않았다. 보존된 L17/L18 비교 JSON과 시간 민감도 JSON은 각각 SHA-256 sidecar 검증을 통과했다. 두 분석 스크립트의 해시도 보고서 기록과 일치한다.

공통 source epoch 21개와 sink-time 표본 42개에서 순간 Newtonian 점질량 쌍성 가속도를 기준으로 한 중앙 PM 가속도 residual의 median/p95/max는 L17에서 0.153897/0.192002/0.243177%, L18의 보간 vector 비교에서 0.061454/0.157814/0.158548%였다. L18 residual은 42개 중 38개에서 더 작았다. 두 residual vector 차이는 Newtonian reference로 정규화한 크기 기준 최대 0.2752%였다. 검증 범위는 짧은 진공 궤도의 reference-force 비교이며 formal convergence order를 정하지 않는다. L18의 실제 시간 간격도 L17보다 짧았다. 따라서 궤도 차이를 순수한 공간 해상도 효과로 해석할 수 없다. Mesh, cloud regularization, self-force와 solver의 개별 기여도는 분리하지 못했다.

L17 시간 간격 ceiling을 P/200에서 P/400으로 낮춘 실행은 41개 구간에서 실제 median `dt`가 기준 실행의 0.5였다. 공통 시각에서 최대 차이는 separation/초기 separation 3.30e-5, energy/초기 절대값 9.56e-6, angular momentum/초기값 4.82e-6, phase 6.43e-5 rad였다. 기준 및 절반 간격 실행의 에너지 drift는 각각 약 −0.1339%, −0.1343%로 거의 같았다. 짧은 시험에서는 시간 간격을 절반으로 줄여도 누적 drift가 개선되지 않았다. 검증 범위는 temporal sensitivity 진단이며 PM Hamiltonian 보존, SIDM 수송, live-halo 경화 결과를 포함하지 않는다.

Astra는 추가 진공 matrix를 중단하고 L17 live-halo Stage 2로 전환하라고 권고했다. 동일한 초기 DM realization과 쌍성 조건으로 collisionless binary+halo를 먼저 실행한다. Collisionless 계산의 실행 가능성을 확인하고 OMP=1 난수 guard의 runtime 검증을 마친 뒤 matched active-SIDM 실행을 진행한다. 잠정 실행 조건은 초기 궤도 주기의 0.1, 시간 간격 상한 P/400, Grammar `normal`의 2 MPI × 1 OpenMP, 16 GiB, 2시간 ceiling, full snapshot 금지다. 이는 생산 계산 승인이 아니라 공학적 feasibility 단계다. 후보 IC는 목표 PM-cloud 쌍성 퍼텐셜에서 평형 상태임이 검증되지 않았다. 본 계산으로 publishable hardening을 주장하지 않는다. 생산 물리 결과를 얻으려면 PM-consistent halo 준비, binary insertion transient 평가, monopole 대조군이 필요하다.

3개월 feasibility를 위한 L10 선별 수치도 분리해 기록한다. Active-SIDM monopole runtime job 1218013은 RAMSES 4,942초 동안 code time `7.684410e-4`까지 1,000 coarse interval을 진행했고, accepted scatter 17,567건 및 sampled pair 459,588,146건을 기록했다. 실행은 확률 guard와 kinetic ledger의 sustained-runtime 검증에 한정되며 radial heat flux, 열수송, 쌍성 경화는 측정하지 않았다. 이 계측률을 유지한다고 가정해 code time 0.25까지 선형 외삽하면 18.61일이다. 별도의 L10 binary force-probe job 1219348은 143초에 code time `2.00034522534e-5`까지 진행했다. 이 fixture의 초기 sink separation은 0.05이고 총질량은 0.0099이므로 Kepler period는 약 0.706 code time이며, 실행 구간은 그 주기의 약 `2.8333e-5`다. 이는 Stage 1의 hard-binary fixture(P=`4.501984e-4`)와 다른 넓은 쌍성이다. L10 binary gate 처리율을 code time 0.25까지 선형 외삽한 20.69일은 넓은 쌍성 계산에 기반하며 hard-binary 처리율로 해석하지 않는다. L10 monopole 처리율을 target hard-binary의 0.1P (code time `4.501984e-5`) 및 0.25P (`1.125496e-4`)로 환산하면 각각 약 289.5초와 723.8초다. 이 역시 L10 code-time 처리율을 일정하다고 가정한 시나리오이지 L17 live-halo 예측은 아니다. L17 진공 실행은 각각 0.1P에서 149.82초(P/200)와 285.58초(P/400)였으나, DM 및 SIDM 계산비용은 포함하지 않는다. code time 0.25까지 네 대조군을 각각 계산하는 L10 단순 산술은 약 74–83일이며, 빌드·반복 seed·분석과 L17 추가 비용까지 고려하면 3개월 일정은 여유롭지 않다. L17 live-halo의 물리시간당 throughput과 peak memory를 직접 측정하기 전에는 생산 규모나 완료일을 확정하지 않는다.

현재 latest-main source worktree의 active-SIDM OMP=1 guard는 Astra 정적 검토를 통과했지만 runtime 검증은 하지 않았다. 전용 Grammar `normal` 검증 승인 요청이 대기 중이다. 과거 L10 hybrid SIDM 결과는 난수 stream 독립성이나 재현성의 증거로 사용하지 않는다. 보고서와 실행 자료는 `/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/` 아래의 force-probe, timestep, `run_sidm_runtime_1000steps_20261026/` bundles에 보존한다.

## 19. 2026-10-26 JCAP RNG/restart 범위 수정 및 v5 검토 번들

Astra의 현재 원고 감사 결론은 “범위가 제한된 JCAP 방법론 논문으로는 일관되지만 아직 제출 준비가 끝나지 않음”이다. RNG와 checkpoint/HDF5 문구를 소스 구현과 보존된 검증 범위에 맞게 좁혔다.

- `main.tex`는 per-thread seed vector가 rank-local state에서 jump offset으로 구성된다는 구현만 기술한다. MPI/OpenMP 배열 전반에서 독립·비반복 stream 또는 hybrid stochastic reproducibility가 검증됐다는 주장은 삭제했다. `paper/data/README.md`에도 해당 증거 간극을 기록했다.
- iSIDM 초기화 문구는 `nrestart>0`일 때 initializer가 state resampling을 건너뛴다는 source-level 동작만 말한다. 전체 checkpoint round trip, 특히 all-ground population 보존은 미검증으로 명시했다.
- HDF5 경로의 `dark_energy_int` write/read는 구현으로 표현하고 end-to-end restoration은 미검증으로 명시했다. 원고 내 추가 restart parameter 문구도 동일하게 정리했다.
- Astra가 수정 문구와 증거표를 재감사해 PASS했다. 추가적인 deterministic-stream, universal restart-preservation, HDF5 round-trip 실증 주장은 남지 않았다.

수정된 원고는 `/tmp`의 별도 output directory에 빌드했다. 26쪽, 1,321,845 bytes, SHA-256 `a6a736ed785cd45316e05783e4c03a2d67267f98b7b4e3753a3990245aa1c2bc`다. 최종 log에 미해결 citation/reference 및 overfull box는 없다. 기존 v4에도 있던 underfull table-cell 3건과 비치명적인 hyperref/font 경고는 남는다. `paper/main.pdf`는 덮어쓰지 않았다.

현재 source snapshot, README, project-local TeX inputs, build log/FLS, Git HEAD 및 dirty patch를 `paper/data/jcap_scope_draft_20261026_v5/`에 고정했고 20개 checksum 검사를 모두 통과했다. 원고 SHA-256은 `3e1665b9c076c1662d94c2a45e309ef9ad2c3c81c5db0e641aef6202e573961f`, 증거표 SHA-256은 `08db93bbc4d31add9d725042d8bdc3ae55fb8a5ae70d35c74802e585ba032b78`이다. v4는 보존했다. Paper 저장소 HEAD `d31a0779c9db71e943f8d42933251c686d4c8e36`의 미커밋 작업도 보존했고 commit, GitHub push, Overleaf push는 하지 않았다.

남은 제출 전 핵심은 untracked/dirty production source까지 포함하는 immutable software/data release를 만들고 manuscript의 Data Availability 식별자와 맞추는 것이다. Inelastic/MPI, drag, full restart, cosmological ADM 검증을 원고 주장으로 추가하지 않는 한 새로운 physics campaign은 필수 조건이 아니다. 이 JCAP 방법론 논문은 live-binary 열수송 가설을 검증하지 않는다.

## 20. 2026-10-26 Series 2 stop-fix gate와 열수송 연구 다음 단계

Grammar의 기존 Series 2 collisionless ledger 작업 1222155는 105초에 RAMSES를 종료했고, task 최대 RSS는 685,500 KiB였다. 두 완결 구간의 최대 상대 M/P/K 잔차는 각각 `3.0573454672e-14`와 `2.4893749232e-14`였다. 그러나 로그는 `step=2`의 세 번째 ledger 구간을 열고 대응 상태 기록 없이 끝났다. 기존 validator의 `PASS`는 두 구간의 회계 결과를 뜻하며 종료 수명주기까지 통과했다는 뜻이 아니다.

종료 수명주기 수정은 별도 worktree `/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/latest_main_series2_stopfix_20261026_v1/`에 보존했다. 기반 HEAD는 `f2b5aa101ed6426503c0ba0ccbd978ddd706c4e0`, 감사된 source delta SHA-256은 `166db4f16e64d7a7aa46acb65b2a052856f23772df9c0e1c907ee183c8a7c493`다. Astra의 정적 검토는 monitor off, `nstepmax=2`, 출력 요청이 없는 dump-free 2-rank gate 범위에서 PASS했다. 사용자가 compile-only job을 승인해 Grammar `normal`에 job 1222259를 1 node, 1 CPU, 32 GiB, 최대 1시간으로 제출했다. 제출 직전 HEAD와 source delta 및 16개 build-input checksum을 재검증했다. 당시 `bin/ramses_final3d`는 없었으며 `bin_ledger/` 아래 과거 후보 실행 파일은 있었지만 새 source delta를 반영한 검증 실행 파일로 취급하지 않았다. Slurm job은 2026-10-26 20:43 KST경 exit `0:0`으로 완료됐다. elapsed는 2분 39초, batch MaxRSS는 `26357692K`였다. Intel ifx 2025.3.0으로 생성된 `bin/ramses_final3d`의 크기는 9,626,640 B, SHA-256은 `429e06231c68dd6eac4b1e70862d600b18f3d5b373619c42b4ac034732179018`이다. `ldd`에서 누락 라이브러리가 없었다. Build 결과와 source/build/log artifact 24개 checksum을 `BUILD_RESULT.json`, `BUILD_ARTIFACTS.sha256`에 보존하고 검증했다. 이 job은 compile/link 전용이며 RAMSES, MPI, 테스트를 실행하지 않았다.

검토 가능한 run package에 빌드된 실행 파일을 `ramses_ledger_smoke`로 복사했다. 실행 파일 hash는 build binary와 같다. `RUN_MANIFEST.json`에는 현재 gate 설정을 고정했고, `SOURCE_PILOT_MANIFEST.json`은 과거 IC engineering-pilot 출처를 명시적으로 분리해 보존했다. 현재 effective namelist hash는 `367aac21088bc342fd220eaa9e41ad1ac45e1d3391d5c1a4853c2ed63dcc60a6`이다. Build-artifact, run-package, binary 및 input checksum 검사를 통과했다. Particle IC는 기존 생성 bundle의 절대 경로로 symlink되며 전체 파일 SHA-256은 `e2a1a15933e2009866708dd024426cc056b7473eb99bf79d217cb8bdf2047f32`다. 사전 실행 일반 파일과 sink IC 합계는 9,648,993 B이며 symlink된 particle IC는 제외한다. 설정은 2 MPI × 1 OpenMP, Grammar `normal`, 16 GiB, 최대 1시간, `nstepmax=2`, `sidm=.false.`다. `noutput=1` 및 미래 `aout=1.1`, `tout=delta_tout=delta_aout=1.0d100`, `foutput=fbackup=100000`을 고정했다. full dump 예상 수는 0이고 runner는 snapshot 생성과 총 산출물 100 MiB 초과를 실패 처리한다. 메타데이터 교정 뒤 package checksum을 다시 만들었으며 현재 `RUN_MANIFEST.json` SHA-256은 `7b049da091016d30e91df81ff90facfd0961d1444a9eba459f3cdaa5db03242f`, `package.sha256` SHA-256은 `fca94c32bca89ec30fdabbb9fb489f1d71434f7d17ea24a3b0d1602cd561ed5b`다. 새 run은 아직 제출하지 않았다. 앞선 16-GiB 승인은 이 수정본 gate의 별도 새 제출을 승인하지 않으므로 사용자 승인을 기다린다.

Stop-fix candidate의 `sidm_scatter.f90`는 수정된 `cell_index_from_part`를 산란 및 사전 점유 분류에서 호출하고 `sidm_probability_guard` 사전 rate gate와 `OMP_NUM_THREADS=1` fail-closed 검사를 포함한다. 감사 범위의 kinetic ledger와 같은 source tree에 존재한다. 준비된 stop-fix MPI gate의 설정은 `sidm=.false.`이므로 통합된 classifier·확률 보호·난수 제한의 실제 active-SIDM 운용을 검증하지 않는다. 이후 calibration은 이 source identity와 모든 기능의 결합을 별도 확인해야 한다.

Astra의 Q-GOAL 판정에 따르면 stop-fix gate는 필요한 회계·정상 종료 검증이지만 열수송 결과는 아니다. Gate 통과는 저장 M/P/K 회계의 완결성을 지지할 뿐 물리적 에너지 보존, SIDM 수송, 평형 또는 쌍성 경화를 입증하지 않는다. Gate의 목적은 다음 SIDM feasibility 계산을 해석할 수 있도록 회계 경계를 닫는 것이다.

Gate 이후의 권고 순서는 다음과 같다.

1. 동일한 frozen source, 입자 realization, live-binary IC, mesh와 cloud 조건에서 collisionless 및 active-SIDM 보정을 수행한다. 제안된 첫 active 구간은 2 MPI × 1 OpenMP, 20–50개 완결 interval, 최대 1시간, full snapshot 없음이다. 동일한 step 수는 SIDM에 의한 timestep 변화 때문에 같은 물리 시간을 뜻하지 않으므로 실제 code time과 event 수를 함께 비교한다. 진단 바이트에도 사전 상한을 둔다.
2. Calibration이 입자 공급과 진단 정밀도를 통과하면 동일 초기조건의 live binary 대 단일 monopole, collisionless 대 active SIDM의 2×2 비교를 물리 시간에 맞춰 단계적으로 수행한다. 단일 monopole 대조군은 쌍성 가열 없이 SIDM 완화와 초기 drift를 분리한다. 모든 대조군을 한 번에 제출하지 않는다.
3. 확대 전에 timestep 민감도, live-halo force/cloud 해상도, 독립 난수 stream 또는 realization의 산포를 측정한다. 수렴 판단은 보존 잔차만이 아니라 열수송과 hardening 관측량의 오차를 대상으로 한다.

첫 열수송 비교에는 반경별 random-motion 에너지 유량과 bulk advection의 분리, 중력 일, 수용된 산란의 구각 에너지 교환, 이동 경계·재중심화·입자 inventory 항, 실제 산란 적격쌍·수용 사건·고유 참여/경계 통과 ID, 쌍성 separation·eccentricity·binding-energy 변화·hardening rate가 필요하다. 사용 가능한 열저장고 에너지 용량은 요구되는 쌍성 에너지 전달량과 별도로 계산한다. 총 헤일로 kinetic energy만으로 저장고 충분성을 판정하지 않는다.

Astra의 이전 응답에 있던 `2a_h` 내부 interval당 `0.0107` pair-slot 추산은 보존된 계산으로 재현되지 않아 철회한다. Astra의 독립 초기 IC 기하 감사는 `/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/run_eddington_periodic_cusp_ic_writer_N1048576_seed20261026_20261026/ics/ic_part` 파일(SHA-256 `e2a1a15933e2009866708dd024426cc056b7473eb99bf79d217cb8bdf2047f32`)을 사용했다. 중심 기준 좌표로 `r=\sqrt{x^2+y^2+z^2}<2a_h`, `2a_h=0.0022225023917858823`를 선택하고 각 좌표에 `1.25`를 더한 뒤 `dx=2.5/2^{16}` 격자로 분류했다. 반경 내부 입자는 46개였고 서로 다른 L16 cell도 46개였다. 두 입자가 모두 반경 안에 있는 초기 무순서 동일-cell 쌍과 disjoint pairing slot은 0개였다. 반경 밖 입자와의 쌍은 집계에 넣지 않았다. 이 수치는 초기 IC의 기하 공급만 나타낸다. 실제 AMR leaf·입자 목록·속도 적격 조건을 처리하는 runtime pairing, 산란 적격률 또는 사건률을 대신하지 않는다. 따라서 짧은 active-SIDM calibration은 시간에 따라 변하는 실제 cell 점유도·적격쌍·수용 사건과 flux 오차를 측정해야 한다. 외곽 산란 검출만으로 내부 열수송이 분해된다고 단정하지 않는다. 생산 N을 늘리기 전에 local pairing 또는 유한 입자 공급이 목표 신호를 측정할 수 있는지 평가한다. 대규모 실행은 transport 통계가 목표 신호를 분해하고, 수치 오차가 신호보다 작으며, 측정된 물리시간당 throughput·peak memory·interval당 저장량이 필요한 비교군과 기간을 지원할 때만 제안한다.

3개월 산출물은 검증 자료에 근거한 진행 또는 연구 방향 전환 결정이다. 이 thermal-transport gate와 SIDM-on/off calibration은 연구방향 ②의 내부상태 반응 지연 및 ③의 회전 붕괴/각운동량 장벽을 검증하지 않는다. 세 연구 방향은 모두 계획에 유지하며 현재 우선순위는 ①이다. 기존 JCAP 방법론 논문의 재현성 릴리스와 누락된 figure 입력도 별도 미완료 항목으로 유지한다.

## 21. 2026-10-26 JCAP inelastic-safety source snapshot

단독 JCAP 자료 릴리스에서 빠졌던 `test_referee_safety.py`와 해당 하네스가 컴파일하는 production source 세 파일을 `paper/data/source_reproduction_20261026_v1/`에 고정했다. 코드 원본은 `sidm` branch HEAD `01b6ead49be96d97a6d3366dfffe4fde36726fe0`이며 선택된 네 파일의 working-tree diff SHA-256은 `6fbb24900769ee40383a28263fd10c39f49ad740da71a5d6ba6e0faad5fdad81`이다. Bundle의 `checksums.sha256` 전체 검증이 통과했고 manifest SHA-256은 `da6260a5b7eec969697de1bcc9465e255c596fa2436ffb618ad8d6189efb415f`이다. Paper repo HEAD `d31a0779c9db71e943f8d42933251c686d4c8e36`는 유지했다. Top-level `paper/data/README.md`와 release-candidate inventory를 갱신해 exact source path, dirty patch와 범위 한계를 기록했다. 이 source-only bundle은 테스트를 실행하지 않았고 inelastic suite 통과, restart 또는 gas coupling을 주장하지 않는다. Astra에 read-only source·provenance·범위 감사를 요청했다. License/접근 조건, 단일 clean release source, zoom-density reduced inputs, DOI는 아직 없으며 commit, push, Overleaf 반영, 공개 release를 하지 않았다.

## 22. 2026-10-26 Series 2 stop-fix gate 제출

사용자가 stop-fix gate의 MPI/Slurm 실행을 승인해, 사전 패키지·입력·실행파일·빌드 산출물 checksum을 재검증한 뒤 Grammar `normal`에 job `1222398`을 제출했다. 제출 시각은 2026-10-26 20:57:52이며 실행 디렉터리는 `/gpfs/kjhan/Hydro/Sidm/validation/smbh_two_body_20261021.LRhN1I/sidm_series2_collisionless_ledger_smoke_20261026_stopfix_v1/`이다. 동일 이름의 활성 job과 기존 `run.complete`, `ramses.log`, `ledger-validation.json`은 없었다.

첫 상태 확인에서 job은 `PENDING (Priority)`였다. 따라서 RAMSES/MPI 실행은 아직 시작되지 않았고 결과는 없다. 고정된 `RUN_MANIFEST.json`은 제출 전 package hash에 포함되므로 수정하지 않았다. 이 gate는 2 MPI × 1 OpenMP, 16 GiB, 최대 1시간의 2-step collisionless ledger/lifecycle 검사이며 SIDM 열수송이나 경화를 시험하지 않는다. 예상 full snapshot은 0개이고 실행 runner가 저장량 100 MiB 상한을 검사한다.

## 23. 2026-10-26 stop-fix gate 실행 권한 수정 및 재제출

첫 제출 job `1222398`은 1초 뒤 Slurm exit `13:0`으로 실패했다. Slurm stderr와 RAMSES redirect log 모두 `execve(): .../ramses_ledger_smoke: Permission denied`를 기록했다. RAMSES 본체는 시작되지 않았고 ledger 결과, snapshot, `run.complete`는 생성되지 않았다. 원인은 실행 파일 내용이 아니라 run package에 복사된 바이너리 권한이 `0644`였던 것이다. 실패 stdout/stderr/log와 원래 package는 보존했다.

Astra의 재시도 package 정적 감사는 PASS했다. 새 디렉터리 `sidm_series2_collisionless_ledger_smoke_20261026_stopfix_retry1/`에서 실행 파일을 `0755`로 설정했으며 바이너리 SHA-256은 동일한 `429e06231c68dd6eac4b1e70862d600b18f3d5b373619c42b4ac034732179018`이다. 입력, effective namelist, 2 MPI × 1 OpenMP, Grammar `normal`, 16 GiB, 최대 1시간, 무덤프 정책 및 100 MiB 상한은 기존 승인과 같다. 새 경로·작업명 `sidm-s2-fix-r1`에 맞춰 package checksums와 manifest 경로/준비 용량을 다시 고정했다. 빌드·입력·실행파일·package checksum이 모두 통과했고 준비 파일 크기는 9,649,130 B다. 재시도 디렉터리에는 실행 산출물이 없었다.

같은 승인 범위의 재시도 job `1222480`을 Grammar `normal`에 제출했다. 첫 상태 조회에서는 `PENDING (Priority)`였으므로 아직 RAMSES 결과는 없다. 계산이 시작되면 별도로 job 결과와 ledger validator를 확인한다.

## 24. Astra의 Objective 1 인과 검정 설계 보강

Astra의 재충전/수송 가설 검토는 “받아들일 총 에너지가 충분함”과 “필요한 속도로 열을 내보냄”을 별개로 측정하라고 권고했다. 쌍성의 사전 지정된 초기·최종 반장축에 대해 근사 Keplerian 고정질량 조건에서는 `ΔE_req=G m1 m2 (1/a_f−1/a_i)/2`를 계산하되, 비-Keplerian 보정과 다른 에너지 sink를 별도 장부에 둔다. 저장고 용량은 사전에 고정한 반경 영역과 경계조건 안에서 허용 가능한 헤일로 상태열의 `E_cap=E_halo,lim−E_halo,initial`로 정의한다. 중력에너지는 해당 periodic/PM 모델과 일관되게 포함하고, 용량의 보수적 하한이 요구량을 넘는지 판정한다. 전체 kinetic energy를 이용 가능 용량으로 간주하지 않는다.

수송 응답은 쌍성의 경화 결과와 독립적으로 측정한다. 준비된 SIDM monopole에서 cloud/상호작용 오염 영역 바깥에 작고 알려진 등방 에너지 pulse를 가하고, 고정 반경을 통과하는 누적 에너지와 `q_r=½ρ⟨|v−u|²(v_r−u_r)⟩`를 추적한다. bulk advection, pressure work, gravitational work, 필요시 nonlocal collision exchange를 구분하고, collisionless pulse control과 무자극 monopole로 비교한다. 응답시간은 예를 들어 사전 지정 에너지 분율이 지정 반경을 통과하는 시각으로 정의하며, ledger 잔차나 저장 에너지를 관측 flux로 나눈 값으로 정의하지 않는다.

가장 작은 인과 대조는 저장고 용량, 중력, insertion 이력, 국소 쌍성 상호작용이 같고 독립 측정된 수송 응답만 다른 두 live-binary 계산이다. 상호작용 영역 바깥의 transport-only 변경이 입자 위치와 쌍성 힘을 직접 바꾸지 않고 선언된 보존량을 지킨다면 인과 대조에 유용하지만, 그 조작 자체를 검증해야 하며 인공 제어로 표기한다. 전역적으로 `σ/m`만 바꾸거나 무한 thermal bath를 두는 비교는 relax/replenishment 또는 저장 용량까지 함께 바꿀 수 있어 단독으로 수송률 인과를 분리하지 못한다. 그러한 대조가 불가능하면 독립 pulse·단일 monopole·collisionless 대조에 근거한 결론은 “수송 병목과 일치”로 제한하고 고유한 인과 증명이라고 부르지 않는다.

지지 근거는 용량 하한의 충분성, 국소 가열 축적과 그보다 느린 독립 측정 수송, 빠른 수송에서 축적 감소 및 유의한 경화 회복, 그리고 수치오차·삽입 transient보다 큰 차이를 함께 요구한다. 반대로 측정 수송이 필요한 열을 충분히 운반하는데도 경화가 멈추거나, 수송을 크게 빠르게 해도 경화가 변하지 않고 궤도 공급이 제한적이면 해당 조건에서 병목 가설을 기각한다. 희박한 산란 쌍, 반복 상관, streaming을 heat로 오인, 변하는 중력 퍼텐셜, 각운동량의 인공 재분배, noisy flux는 결과를 inconclusive하게 만들 수 있다. 이 제안은 아직 실행되지 않았고 현재 stop-fix gate의 범위를 바꾸지 않는다.

## 25. 2026-10-26 Series 2 stop-fix gate 결과 및 Astra 감사

재시도 job `1222480`은 Grammar `normal`에서 Slurm/batch/MPI step 모두 `COMPLETED`, exit `0:0`으로 끝났다. 경과는 Slurm 98초, RAMSES 97초였고 MPI task 최대 RSS는 682,016 KiB였다. 최종 `run.complete`, 10개 output checksum 및 고정 입력/build checksum이 통과했다. 실행 파일은 승인된 SHA-256 `429e06231c68dd6eac4b1e70862d600b18f3d5b373619c42b4ac034732179018`이다.

Ledger validator는 `PASS`를 기록했다. 완결된 구간은 step 0과 1이며 각각 양의 `dt`는 `2.974e-6`, `2.947e-6` code time, 최대 상대 M/P/K 잔차는 `3.0573454672e-14`, `2.4893749232e-14`다. 종료 로그는 `SIDM_KLEDGER_TERMINAL step=2 interval_active=F`를 남겼고 step 2의 새 BEGIN은 없다. Full snapshot은 0개다. 이는 collisionless 저장 M/P/K 장부와 종료 수명주기 검증이지, SIDM 산란·열수송·평형·쌍성 경화 또는 restart/출력 활성화 경로의 검증이 아니다.

Astra의 독립 read-only 결과 감사도 PASS했다. `artifact-sizes.json`은 후속 완료 산출물 생성 전 9,992,008 B를 기록한다. 완료 후 직접 파일과 sink IC를 합산한 최종 크기는 9,993,453 B이며 차이 1,445 B는 뒤에 생성된 `artifact-sizes.json`(197 B), `output-checksums.sha256`(820 B), `run.complete`(228 B), Slurm stdout의 최종 checksum 확인(200 B) 합과 정확히 같다. 최종 크기 역시 100 MiB 상한 아래이고 symlink particle IC payload는 제외한다. 같은 retry job은 종료됐고 queue에 남아 있지 않다.

다음 연구 gate는 이 승인된 collisionless run과 동일한 source/IC/mesh에서 active SIDM을 켜 실제 산란·cell 분류·확률 guard·OMP=1 제한·ledger 결합을 확인하는 짧은 calibration이다. Astra에 20–50 완결 interval의 구체적 stop 조건, 계측 및 판정 기준을 read-only로 설계 의뢰했다. 이는 별도 active-SIDM 계산이며 이번 결과만으로 실행 승인이 되지는 않는다.

## 26. Active-SIDM ledger calibration의 코드상 선행조건

완료된 stop-fix build worktree의 namelist 기본값과 guard를 다시 확인했다. `sidm_probability_guard` 기본값은 `.false.`이며 `sidm=`.true.`와 `sidm_kinetic_ledger=`.true.`를 함께 켤 때 이 guard를 명시적으로 `.true.`로 설정하지 않으면 초기화 단계에서 중단된다. Guard는 rate preflight 뒤 실제 산란을 수행하고 `SIDM_RATE_ACTUAL` 카운터를 남긴다. `sidm_estimator_diagnostics=.true.`는 별도로 rate-preflight 및 leaf estimator 로그를 기록한다. Kinetic ledger는 `OMP_NUM_THREADS=1`을 강제하므로 다음 run도 MPI 2 × OpenMP 1이어야 한다. 이 source-level 조건은 다음 package의 namelist·runner acceptance에 반영해야 한다. 활성화 후 기록량은 아직 측정하지 않았으므로 별도 크기 상한을 유지한다.

## 27. 2026-10-26 active-SIDM 20-interval calibration: 확률 guard 중단

Astra가 수정된 실행 패키지를 정적으로 PASS 판정했고 사용자가 제출을 승인했다. package SHA-256은 4e5733b1becec855a5bd849540fad5f239d3f74dab23284e9ea31e74ee21e853, 실행파일 SHA-256은 429e06231c68dd6eac4b1e70862d600b18f3d5b373619c42b4ac034732179018이다. Grammar `normal`에 job 1222562를 2 MPI × 1 OpenMP, 8 GiB, 1시간 한도로 제출했다.

Job은 1:59 후 FAILED로 끝났다. RAMSES runner exit code는 137이고 MPI step은 CANCELLED+ (1:58, 0:9), batch exit은 9:0이다. sacct task MaxRSS는 684,152 KiB, batch MaxRSS는 7,300 KiB다. 이는 OOM이나 시간제한이 아니라 의도된 사전 산란 확률 guard의 MPI abort다.

Ledger step 0과 1은 각각 PASS였고 최대 상대 M/P/K 잔차는 1.8609928931e-14와 7.6131368912e-15다. Coarse step 2는 BEGIN을 기록했지만 STATUS와 terminal certificate 없이 중단됐다. Level 14의 preflight rate는 coarse step 1에서 0이었고 step 2에서 1.999140e-8 s^-1로 증가했다. 이때 dt_phys=2.947497e8 s, P_bound=5.892458236, P_limit=0.1이어서 해당 level의 산란 전에 MPI_ABORT가 발생했다. 따라서 20구간 gate는 실패이며 validator와 run.complete는 생성되지 않았다. 이는 앞선 정상 종료 lifecycle 실패와는 다른, 안전 guard가 만든 의도된 중단이다.

실행 전 bootstrap coarse step 0은 각 level의 dt=0, P_bound=0, 실제 Pmax=0으로 기록됐고 accepted event는 없었다. 완결된 step 1에서 7 events가 수락됐다. 미완결 step 2의 level 8/9에서 2 events가 더 수락됐으나 이들은 닫힌 ledger 구간의 결과로 세지 않는다. 총 sampled-pair 수는 step 1에서 315,980, 미완결 step 2에서 중단 전까지 315,991이었다. 기대 사건 합은 각각 6.55021과 6.54953이다.

Full snapshot은 없었다. 부분 산출물을 보존했고 실행파일과 sink IC를 포함한 직접 측정 저장량은 10,092,037 B로 100 MiB 아래이며 particle IC symlink payload는 제외한다. 이 job은 열수송, 평형 또는 쌍성 경화를 측정하지 않았다.

Astra의 독립 결과 감사는 이 gate를 FAIL로 확인했다. 원인은 이전 timestep 선택이 과거 rate bound에 의존하는 동안 새 leaf/입자 점유가 커져 rate가 급증한 것이다. Guard는 과도한 확률로 산란하기 전에 계산을 멈췄다. 다음 설계는 guard를 완화하지 않고 새로 활성화된 cell population에서도 적절한 timestep을 산란 및 해당 sweep의 write 전에 보장하는 방법을 다뤄야 한다. collisionless_courant_factor는 중력·입자 CFL을 낮추지만 절대 timestep 상한은 아니므로, 임의 값을 넣은 재실행을 보장된 해결로 간주하지 않는다. 새 timestep 설정이나 source 변경을 사용하는 계산은 별도 패키지와 승인을 거쳐야 한다.

Astra는 재실행에 앞서 pinned geometry classifier, rate-bound helper, timestep helper를 함께 사용하는 결정론적 new-pair/zero-to-positive-rate regression으로 이 전이를 재현하라고 권고했다. 이 회귀 검증은 아직 만들거나 실행하지 않았다. 장기 계산을 위한 안전한 처리는 과거 rate를 사용하는 timestep predictor가 새 점유를 어떻게 처리할지 또는 유효한 rollback/subcycling이 가능한지를 설계해야 하며, 이미 실행된 하위-level 연산을 버리거나 확률 기준을 완화해서는 안 된다.
