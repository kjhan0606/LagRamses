# Block grid-major 프로젝트 최종 결론

- 최종 상태: **CLOSED**
- 폐문일: 2026-08-22
- 대상 저장소/브랜치: `github.com/kjhan0606/LagRamses`, `main`
- 검증된 구현/실행 기준선: `c861f2a17e417e46834dc05fc790d4e5eb1face3`
- 승인 범위: 사용자가 확정한 기본 기능·물리 동등성·기능 smoke 범위

## 1. 최종 판정

block grid-major 개편은 승인 범위에서 완료되었다. `B=64` block grid-major
cell index, positive growth를 검증한 grid dynamic capacity, no-overflow 경로를
검증한 particle dynamic-capacity 코드, production load balance headroom 사전 확장,
canonical AMR 비교, CUDA block index 교정과 nGR 기능 경로를 `main`에 통합했다.
legacy output은 기본 형식으로 유지한다.

추가 대형 시뮬레이션이나 VoidSim 실행은 완료 조건에서 제외했다. 아래의 알려진
제한은 숨기거나 합격으로 재해석하지 않고 그대로 동결한다. 이 문서 작성 이후
추가 필수 Slurm 실행은 없으며, 새 요구가 생길 때만 별도 범위로 재개한다.

## 2. 완료 근거

| 항목 | 근거 | 최종 판정 |
|---|---|---|
| collective grid growth | job 321554, 1 node/32 MPI, tight run `90688 -> 113408`, hydro/part/grav와 physics log exact | PASS |
| production LB headroom | job 321942, growth가 첫 real LB 전에 발생, AMR/info 포함 runtime output exact | PASS |
| canonical AMR topology | job 321554와 321942의 6개 output pair | PASS |
| CUDA/nGR 기본 기능 | job 322420 (Slurm/batch gate FAILED); CPU/GPU child run은 각각 rc=0·Run completed, CUDA MG·nGR scalar·particle gather/deposit marker positive; strict comparator는 transient topology와 velocity relL2=`5.9469e-3 > 2e-3`로 FAIL | LIMITED FUNCTIONAL-PATH EVIDENCE / STRICT PARITY FAIL |
| HDF5 restart 기본 I/O | job 322427 (Slurm/batch gate FAILED); continuous/restart child run은 각각 rc=0·Run completed 및 restore marker 확인; one-step exact `h5diff` FAIL | OPERATIONAL RESTORE PATH OBSERVED / EXACT CONTINUITY FAIL |
| layout 물리 동등성 | job 322428, ABBA 두 pair의 초기 payload와 8-step physics trace exact, canonical topology exact | PASS |
| layout 성능 | job 322428, shared node/4 MPI/8 step | `NOISY_INCONCLUSIVE` |

Phase 4 base gate는 auto grid capacity가 실제로 성장한 뒤에도 동일한 물리 payload를
보존함을 확인했다. production-LB gate는 `lb_grid_headroom=0.85`에서 candidate가
첫 load balance 전에 control과 같은 capacity로 collective 성장하고, 이후 LB state와
전체 legacy output이 일치함을 확인했다.

CUDA job 322420의 배치 gate 자체는 통과하지 않았다. 다만 stale host residual norm이
device norm을 덮던 multigrid 오류를 수정한 뒤 두 child run의 완주와 네 CUDA 경로의
실제 실행을 좁은 기능 증거로 수용했다. strict CPU/GPU 수치 동등성은 승인하지 않는다.

## 3. 최종 성능 관측

job 322428은 같은 현재 source에서 pristine `B=64`와, 정확히 두 assignment만
`B=ngridmax`로 바꾼 legacy-equivalent control을 비교했다. 두 실행 모두 fixed
capacity, CPU-only, 1 node/4 MPI/OMP 1, 256^3 IC와 8 coarse step을 사용했다.

| 지표 | legacy 평균 (s) | block 평균 (s) | block/legacy |
|---|---:|---:|---:|
| program elapsed | 811.248 | 846.374 | 1.0433 |
| timer total | 816.719 | 853.836 | 1.0454 |
| external wall | 911.107 | 951.991 | 1.0449 |

관측점은 block layout이 약 4--5% 느린 방향이다. 그러나 shared node의 co-tenant가
legacy 구간에서 바뀌었고 paired 차이도 `+3.16%`, `+5.48%`로 갈렸다. 따라서 이
결과로 5% 이상 성능 회귀나 성능 동등성을 확정하지 않는다. 결론은
`correctness PASS / performance characterization complete / NOISY_INCONCLUSIVE`다.

원시 보고서는 다음 위치에 보존한다.

```text
/home/kjhan/BACKUP/lagramses_qa/perf_layout_jobs/322428/reports/performance.md
/home/kjhan/BACKUP/lagramses_qa/perf_layout_jobs/322428/reports/performance.json
/home/kjhan/BACKUP/lagramses_qa/perf_layout_jobs/322428/reports/canonical_output1.json
```

## 4. 폐문 시 알려진 제한

다음 항목은 이번 완료 범위에 포함하지 않으며, 필요하면 새 과제로 연다.

1. CUDA smoke의 strict CPU/GPU numerical-parity gate는 통과하지 않았다. 정확히
   `m_refine=8`인 fixture에서 CPU 직렬 합산과 GPU `atomicAdd` 순서 차이가 refinement
   tie를 다르게 선택해 transient topology와 velocity 차이를 만들었다. block index나
   수정된 MG residual 계산의 오류 증거로 판정하지 않았지만 strict parity를 승인하지도
   않는다.
2. HDF5 checkpoint의 기본 restore/I/O는 동작하지만 one-step restart output의
   zero-tolerance exact continuity는 통과하지 않았다. AMR topology와 particle ID/mass,
   phi/scalar는 강하게 일치했으나 L6 force와 velocity, 일부 history/capacity state가
   달랐다. production-correct exact restart를 요구하면 별도 blocker로 재개한다.
3. Particle dynamic-capacity 구현과 no-overflow 경로는 통합·검증했지만, 실제
   particle-capacity growth가 발생한 positive runtime event는 확인하지 못했다.
   따라서 particle-growth runtime 승인은 주장하지 않는다.
4. RT/ATON/radiation, optional FDM/MG bundle, repeated/refine-coarse growth,
   fixed-capacity negative semantics, multi-node와 different-ncpu restart는 전수 runtime
   검증하지 않았다.
5. 성능 결과는 shared-node 방향성 측정이다. 전용 노드 반복 benchmark를 수행하지
   않았으므로 확정적인 performance claim에 사용하지 않는다.

## 5. 운영 결론

- `main`이 이 프로젝트의 최종 통합 브랜치다.
- VoidSim과 다른 프로젝트 자산은 완료 근거에 포함하지 않는다.
- 추가 기본 테스트, 예약 작업, 필수 merge는 남아 있지 않다.
- 향후 변경은 이 프로젝트의 미완료 작업이 아니라 명시적으로 승인된 새 과제로
  취급한다.

이로써 block grid-major 프로젝트를 폐문한다.
