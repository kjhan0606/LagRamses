# RT / feedback / dust 비교 실행본 — 종료 인계

2026-09-07. 사용자 승인: 종료 범위를 다시 넓히지 않고, 이미 실행 확인된
선택 모형을 고정하여 인계한다. **이 범위의 구현은 종료한다.**
production 장시간 운전이나 publication 수준의 과학 검증 완료를 뜻하지 않는다.

## 고정한 실행 범위

- 프로젝트/작업공간: `kjhan0606/LagRamses`, `/gpfs/kjhan/LRD_JWST`.
- 구현 소스: `073daf62db0ea26106131a1c66474700564febbd` (GitHub main에 푸시됨).
  바이너리는 커밋 직전 빌드이므로 식별은 아래 SHA256을 사용한다.
- 단일 MPI rank, OpenMP 2 threads, 고정 level 3, 4 coarse steps.
  `SNRT_BACKEND=openmp`, reduced-c=0.01c. GPU 전환이나 다중 rank로 확대하지 않는다.
- HDF5=1, SNRT=1, DUST_LIVE=1, USE_CUDA=1, USE_FFTW=0, NVAR=30.
  GPU 실행은 필요 없지만 Intel MPI/OpenMP·HDF5·CUDA 공유 라이브러리는 필요하다.
  현재 작업공간에서 `ldd`의 누락 라이브러리가 없음을 확인했다.
- 피드백: Kroupa 0.08–120 Msun / binary SSP fraction 0.5, KL16 AGB,
  LC18 `wind_only_collapse`, effective-SSP SNIa. 기존 질량 회계와 DTD 유지.
- 복사: 별 입자의 초기 질량·나이·Z → 독립 BPASS 복사 인구 → native SNRT.
  AGN은 기존 BH의 실제 accretion → reference partition 경로. spin/MAD는 비활성.
- dust: Draine optics + DL01 bulk graphite/silicate 30/70, 단일 온도 U(T), live IR.
- 이 입력은 비우주론적 수치 비교용이다. `cooling=.false.`는 별도 gas cooling
  모듈 비활성 설정이며 SNRT H/He 열화학을 끄는 뜻이 아니다.

## 보존해야 할 실행 파일과 입력

아래 경로는 모두 `/gpfs/kjhan/LRD_JWST/` 기준이다. 숨김 디렉토리도 임시 폐기
대상이 아니라 **현재 인계의 로컬 의존성**이다. 옮기거나 삭제하지 않았다.
GitHub clone만으로 바이너리와 로컬 yield 입력까지 배포되는 것은 아니다.
531 MB BPASS 원본 HDF5는 오프라인 변환용이며 이 실행에서는 읽지 않는다.

| 파일 | SHA256 |
| --- | --- |
| `.bpass-native.v0ZwR6/ramses_bpass_native3d` | `5555c1ba1fb428aa643f88332ccafe8eac277c0e37eae3e0ce9c027c2a598dba` |
| `.bpass-native.v0ZwR6/live/run.nml` | `04be3ae220035c4ace435fc701d564da2a96c5b60ea3a9510d7a17bf9b9c896c` |
| `.bpass-native.v0ZwR6/live/ic_sink` | `8e59783a64e8baeba5d090b781a6c44501e62dc619926ceeefff6581afd4b53d` |
| `.agb-physical.4LAOTJ/snia-input/history.nml` | `894339fd9e1b1238b1a658cb5829b20c25cdb94b992c308e405012792091690d` |
| `.agb-physical.4LAOTJ/snia-input/yields.dat` | `de15e4e62b757acaac581b4bf536c733cc866d9db49f88a2d4393991cb7c87d2` |
| `simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml` | `3a167ebaf393bd2ac11217cb0a8b278e2e8da793a91c4b36e8a57e4464e4dff5` |
| `simulation/snrt/config/snrt_stellar_sed_bpass_independent_v2.nml` | `ebd5e097e0e9eda50ca029747cb12d89f0a1d89ac02dda665f1381fcccbab840` |
| `simulation/snrt/config/dust_dl01_bulk_030_reference_v4.nml` | `8072b09a592a2f13e9c3607f4d3d0a2c5d6e328072b3bfb8429eee55cced077e` |
| `simulation/snrt/config/snrt_group_contract_reference_control_v1.nml` | `5825dae4d55d5f4448880c0ebd5b9b727ecfa97cbcf573b609b6e9bfa6e6ca94` |
| `simulation/snrt/config/snrt_secondary_table_contract_v1.nml` | `558cbc6114ba31d45dd068fc19f131499fecf366e09cba80d6a7e648aa43854c` |
| `simulation/snrt/config/snrt_agn_driver_faithful_smoke_yields.dat` | `d88d62a619a8ed3c2e1f05395af86a85728a24d6d2ff3d880585a58414e8d69c` |

## 재현 실행 안내 — 이번 종료에서는 다시 실행하지 않음

기존 출력 디렉토리에서 재실행하지 않는다. 새 실행이 필요할 때만 아래처럼
새 디렉토리를 만든다. 복사하는 namelist는 placeholder가 이미 해결된 실제 입력이다.

```bash
set -euo pipefail
comparison_root=/gpfs/kjhan/LRD_JWST
comparison_run=$(mktemp -d "$comparison_root/.comparison-replay.XXXXXX")
cp -pn "$comparison_root/.bpass-native.v0ZwR6/live/run.nml" "$comparison_run/run.nml"
cp -pn "$comparison_root/.bpass-native.v0ZwR6/live/ic_sink" "$comparison_run/ic_sink"
cd "$comparison_run"
pwd
rg -n 'nstepmax|noutput|aout|tout|foutput|fbackup' run.nml
df -h .
```

실행 전 위 **새 절대 경로의 namelist**와 여유 공간을 다시 보고한다.
`noutput=1, aout=2, tout=1e30`은 종료 전 도달하지 않는다.
`foutput=2, fbackup=1000000, nstepmax=4`: 덤프 2개, 약 56 MB/개,
총 약 112 MB. 새로운 계산의 출력 정책으로 무심코 복제하지 않는다.

다음은 보고 후 별도로 실행할 명령이며, 설정 변경 없는 비교 재현에 한정한다.

```bash
env -u SNRT_DRIVER_TEST_SEED_SOURCE -u SNRT_RT_TX_DIAGNOSTIC_MODE \
 OMP_NUM_THREADS=2 I_MPI_FABRICS=shm SNRT_RT_ENABLE=1 SNRT_BACKEND=openmp \
 SNRT_AGN_MODEL=partition_reference_v1 SNRT_REDUCED_C=.01 SNRT_RT_LEVEL=3 \
 SNRT_ALLOW_REFERENCE_CONTROL=1 SNRT_P1_DIAGNOSTIC=0 \
 SNRT_GROUP_CONTRACT="$comparison_root/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml" \
 SNRT_SECONDARY_TABLE_CONTRACT="$comparison_root/simulation/snrt/config/snrt_secondary_table_contract_v1.nml" \
 SNRT_DUST_CONTRACT="$comparison_root/simulation/snrt/config/dust_dl01_bulk_030_reference_v4.nml" \
 SNRT_STELLAR_SED="$comparison_root/simulation/snrt/config/snrt_stellar_sed_bpass_independent_v2.nml" \
 PHASE0_YIELD_TABLE="$comparison_root/.agb-physical.4LAOTJ/snia-input/yields.dat" \
 PHASE0_SNIA_RUNTIME_CONTRACT="$comparison_root/simulation/snrt/config/fp2_snia_effective_ssp_runtime_v1.nml" \
 "$comparison_root/.bpass-native.v0ZwR6/ramses_bpass_native3d" run.nml > run.log 2>&1
```

재시작은 `output_00001`(coarse step 2)을 **다른 새 디렉토리에 복사**하고
`nrestart=1`로 설정하여 step 4까지 이어간다. 나머지 입력·환경은 동일하게 유지한다.
이 경우 새 덤프 1개를 예상한다. source-free checkpoint에 BPASS를 추가하거나
인구·출처·dust 물질을 바꾸는 것은 허용된 재시작이 아니다.

## 종료 판단에 사용한 기존 증거

새 빌드·계산·회귀검사·외부 감사를 추가하지 않았다. 기존 결과 파일과 기록,
현재 바이너리·입력 SHA256 및 동적 라이브러리 가용성을 읽기 전용으로 확인했다.

- `.bpass-native.v0ZwR6/stellar-tests.log`: 기존 v1 유지, BPASS 적분/분할/범위,
  9개 잘못된 입력 거부 통과.
- `live/run.log`: 4-step 완료, 실제 광원 활성화, 화학 실패 없음.
  최대 dust/IR 수지 오차 6.1550e-10. 기존 출력은 총 110145532 bytes.
- `restart/run.log`: 완료. 기존 비교에서 RT·dust는 정확히 동일;
  94개 hydro/RT 배열 중 91개 동일, 운동량 3개만 반올림 수준 차이.
- `restart-source-mismatch/run.log`: 출처만 변경해도 진화 전 HDF5 거부.

상세 숫자와 실패 이력은 [기존 진행 기록](real_source_integration_progress_2026-09-07.md)에 있다.
`Run completed`만으로 과학적 검증 통과를 주장하지 않는다.

## 알려진 제한 — 종료를 다시 미루는 새 필수조건이 아님

BPASS와 feedback은 서로 다른 인구이며, BPASS 0–1 Myr는 첫 스펙트럼 유지,
미측정 파장 꼬리는 0, escape fraction=1이다. 수송의 에너지·단면적은 공통 AGN
기준 grey 근사다. dust는 bulk 단일 온도이며 PAH 확률 가열/승화 모형이 아니다.
전체 feedback 입력의 공통 Z 범위는 .007–.01345이며, 실제 프로파일은 Z=.01이다.
큰 BH seed와 긴 초기 timestep은 수치 비교용이다. 이 짧은 결합 실행은 첫 terminal
AGB 방출까지 도달하지 않았다. 기존 비우주론적 SFRD 출력의 NaN은 남아 있으나,
기록된 진화 상태의 NaN_CHK 개수는 모두 0이었다.

Primary SNRT의 auto OpenMP/CUDA 선택은 이미 구현되어 있지만 **이번 고정 실행본은
OpenMP만 확인한 범위**다. mechanical feedback과 dust thermal/IR은 host 구현이다.
전체 조합의 MPI/AMR/GPU 장시간 자격이나 CUDA 없는 빌드를 주장하지 않는다.

## 별도 검증 계획 제안 — 자동 실행하지 않음

사용자가 후속 검증을 요청하면 다음 세 항목을 하나의 검증 작업으로 묶는다.
기존 검사를 그대로 반복하거나 단계별 새 감사 체계를 만들지 않는다.

1. 선택 모형의 적분 질량·원소·광자·에너지 회계 및 AGB 방출을 포함하는 시간 범위.
2. timestep·공간/각 해상도·reduced-c에 따른 주요 물리 관측량의 민감도/수렴.
3. 실제 사용할 규모의 MPI/AMR·백엔드·재시작 운영 검증. 목표 규모와 저장량은 사전 확정.

스펙트럼별 수송 계수/혼합 복사장 개선, 같은 인구의 stellar SED,
미시적 SNIa/BPS, 고질량 대안과 고급 dust 모형은 **별도 물리 개선 과제**다.
별도 승인 없이 구현을 시작하거나 이번 종료 조건으로 되돌리지 않는다.
