# CHIMES GPU+CPU 하이브리드 속도 업무인수인계서

작성일: 2026-10-06 (Asia/Seoul)

## 1. 결론

다크 CHIMES의 CVODE 적분에서 GPU+CPU 하이브리드는 8스레드 CPU보다 의미 있게 빠르지 않다. H200에서 같은 셀, 같은 공차로 잰 최고 기록은 512셀 **1.02배**다. 128셀은 0.54배, 256셀은 0.84배로 CPU보다 느리다. 생산 경로로 켜지 말 것. 공차를 느슨하게 하지 말 것. 물리 항을 추가하지 말 것.

이 루틴은 GPU가 힘을 쓸 자리가 거의 없다. 냉각 커널은 셀당 스레드 하나의 긴 `pow`/`exp` 사슬이라 H200의 SM을 채우지 못한다. RHS만 보면 GPU 스레드 하나가 CPU 스레드 하나보다 약 4배 빠르다(1.90 µs 대 7.62 µs). CPU 8스레드가 RHS를 동시에 하면 그 처리량이 GPU RHS 파이프라인과 비슷하다. 적분 스텝의 약 60%는 호스트 CVODE라 GPU가 건드리지 않는다. RHS를 공짜로 만들어도 상한은 약 1.14배다.

## 2. 저장소

| 항목 | 값 |
|---|---|
| 작업 저장소 | `/gpfs/kjhan/chimes-cvode-gpu-20261001` |
| origin | `git@github.com:kjhan0606/LagRamses.git` |
| branch | `agent/chimes-cvode-gpu-20261001` |
| 이 문서의 커밋 | 이 파일을 추가한 커밋. 부모는 `89be0e8` |
| 로컬 사본 | `/home/kjhan/BACKUP/LRD_JWST/CHIMES_GPU_HYBRID_HANDOVER_2026-10-06.md` |

`/home/kjhan/BACKUP/LRD_JWST`는 다른 Git 저장소다. 그 트리를 커밋하거나 푸시하지 말 것. `.cvode-gpu/`, `.chimes-real-ramses-128-build-*`, `.mg-rb-dynamic-build-*`는 빌드·로그·체크아웃이다. 커밋하지 말 것.

이전 인수인계 `CHIMES_GPU_CVODE_HANDOVER_2026-10-06.md`와 패치 `CHIMES_GPU_CVODE_FIXED_ROUTE_2026-10-06.patch`는 고정 worker/backend 진단용이다. 그 패치는 GPFS 소스에 적용하지 않았다. 이번 속도 결론과 섞지 말 것.

## 3. 현재 스케줄

기본은 **OpenMP 스레드 8개 + CUDA 스트림 8개**다. 스트림 1개 + 스레드 8개가 아니다.

`SNRT_CHIMES_COMPUTE_BACKEND=cuda_integrated`이고 `SNRT_CHIMES_CELL_SCHEDULER=level_queue`이면, 레벨 셀이 스레드당 16개 이상일 때 워커가 모두 GPU 브로커가 된다. 각자 `cuda_acquire_stream()`으로 스트림 하나를 임대해 자기 셀의 RHS 배치를 띄운다. 풀은 `cuda_pool_init(0, 8)`이다. 브로커는 버퍼 두 개를 가지지만 둘 다 그 스트림 하나에 순서대로 넣는다. 런치는 `cudaMemcpyAsync`와 커널, 이벤트까지이고, `cudaEventSynchronize`는 그 브로커의 섬유가 모두 멈췄을 때만 한다.

`SNRT_CHIMES_GPU_BROKERS`를 두면 그 수만큼만 브로커가 된다. `SNRT_CHIMES_GPU_WORKER=thread0`이면 브로커는 1개다. 셀이 `16 * workers`보다 적으면 브로커는 0이고 CPU만 쓴다. 32셀, 8스레드가 여기에 해당한다.

코드의 `shared_team`은 `false`다. 스트림 하나로 모으는 경로는 소스에 남아 있으나 켜지지 않는다. 작업 414163에서 느렸기 때문이다.

CVODE 제어와 157×157 dense LU는 호스트에 있다. GPU는 `N_spectra==0`인 다크 반응률과 다크 순냉각만 계산한다. 디바이스 LU를 타일 경로에 켜지 말 것.

## 4. 측정

공차는 생산 다크 셀과 같다. `relativeTolerance=1e-10`, `absoluteTolerance=1e-17`, `explicitTolerance=1e-10`, `scale_metal_tolerances=1`. dt는 `1e10` s. 비교는 8스레드 CPU 대 8스레드 하이브리드이고, 어느 쪽도 참값으로 두지 않는다.

드라이버는 `.cvode-gpu/hybrid-fixture/chimes_hybrid_integrate.cpp`다. 브리지(MPI/OpenSSL)에 링크하지 않는다. 냉각 훅은 RHS 픽스처와 같다. `hybrid_cooling_mode=1`, `hybrid_cooling_fn=hybrid_dark`, `snrt_chimes_molecular_coefficients=dummy_mol`. 함량은 브리지 보존 장부와 같지 않다. `chimes_network`는 그 장부를 검사하지 않는다.

H200 잡에는 `--partition=h200 --exclude=syn102 --gres=gpu:H200`이 필요하다. `gpu:1`은 H200 자격이 되지 않는다. 자격 GPU는 H200/H100/A100이다. A10 결과는 자격 측정이 아니다.

### 4.1 RHS만, 작업 414110, syn104, H200 NVL

로그: `.cvode-gpu/hybrid-fixture/chimes-hybrid-414110.log`

32셀 RHS 픽스처. 순냉각과 157종. 최대 상대편차 `3.184e-14`. 비트 일치는 아니다.

| 비교 | 배속 | GPU µs | CPU µs |
|---|---:|---:|---:|
| gpu1 / cpu1 | 4.00× | 1.904 | 7.615 |
| gpu4 / cpu4 | 3.39× | 0.570 | 1.936 |
| hyb8 / cpu8 | 2.00× | 0.524 | 1.047 |

이것은 RHS 마이크로벤치다. 적분 벽시계가 아니다.

### 4.2 적분, 브로커마다 스트림, 작업 414170

로그: `.cvode-gpu/hybrid-fixture/chimes-cvode-414170.log`

상태 0, retries 0, fallback 0, errors 0. 온도 최대 상대오차 `4.066e-15`, 함량 `1.227e-11`(종 147, 뜨거운 셀), 절대오차 `7.216e-16`.

| 셀 | CPU 벽시계 | 하이브리드 | 배속 | 평균 배치 | GPU RHS |
|---:|---:|---:|---:|---:|---:|
| 32 | 0.0986 s | 0.0971 s | 1.02× | 0 (CPU) | 0 |
| 128 | 0.3543 s | 0.6508 s | 0.54× | 6.83 | 161409 |
| 256 | 0.7072 s | 0.8421 s | 0.84× | 13.56 | 323154 |
| 512 | 1.4025 s | 1.3765 s | 1.02× | 27.10 | 646658 |

32셀 하이브리드는 CPU와 비트 단위로 같다. GPU를 쓰지 않기 때문이다. 128셀 이상은 GPU가 RHS를 전부 계산한다. 함량이 `1e-11`까지 다른 것은 그 증거다.

같은 정책의 앞선 작업 414142는 512셀에서 1.06배였다. 414170의 1.02배와 같은 자리다. 512셀에서 CPU와 같다는 결론을 바꾸지 않는다.

작업 414120은 옛 4브로커 스케줄이다. 32셀에서 배치가 6.8이고 브로커가 셀을 모두 가져가 CPU 워커가 쉬었다. 벽시계는 0.29배였다. 그 숫자를 현재 코드의 성능으로 쓰지 말 것.

### 4.3 한 스트림으로 모은 배치, 작업 414163

로그: `.cvode-gpu/hybrid-fixture/chimes-cvode-414163.log`

0번 스레드만 런치하고 나머지 7스레드는 호스트 CVODE를 이어 갔다. 배치는 약 56–64셀까지 찼다. 정확도는 4.2절과 같다.

| 셀 | 배속 | 평균 배치 |
|---:|---:|---:|
| 32 | 1.02× | 0 (CPU) |
| 128 | 0.31× | 56 |
| 256 | 0.36× | 62 |
| 512 | 0.34× | 64 |

8스트림 × 약 27셀은 GPU 위에 약 200셀이 있다. 1스트림 × 64셀은 64셀이다. 배치를 키워도 `<<<n,1>>>` 냉각 커널의 지연은 거의 줄지 않는다. 한 덩어리로 합치면 8스트림의 동시성이 빠진다. 이 정책을 다시 켜지 말 것.

### 4.4 멈춘 작업

- 414150: 공유 배치의 조건 변수와 `level_launch_busy`가 `cudaStreamSynchronize` 동안 겹쳐 128셀 하이브리드에서 멈췄다. Slurm 시간 초과. 그 대기 코드는 빼 두었다.
- 414168: 8브로커인데 런치와 대기가 `active_brokers<=1`에만 묶여 섬유가 RHS에서 멈췄다. WorkDir를 확인한 뒤 이 작업만 `scancel`했다. 그 조건은 제거했다.

## 5. 비동기로 다른 셀을 하는 구조

이미 들어가 있다. 섬유가 `submit_at`개 막히면 배치를 큐에 넣고, 같은 스레드가 아직 막히지 않은 섬유의 호스트 CVODE를 계속 돈다. 512셀이면 브로커 하나가 셀 64개를 미리 가져가고 `submit_at`은 32다. 뒤 32개의 호스트 일이 GPU 왕복보다 짧다. 8스레드가 모두 브로커라, GPU를 기다리는 동안 순수 CPU 셀을 집어 갈 워커도 없다. 그래서 겹치기가 맞아도 벽시계는 1.02배에 머문다.

## 6. 버스 데이터

테이블과 반응 목록은 스트림 초기화 때 한 번만 올라가고 디바이스에 남는다.

RHS 한 번, 셀 하나:

| 방향 | 내용 | 크기 |
|---|---|---:|
| CPU→GPU | 함량 157개 | 1256 B |
| CPU→GPU | 온도, 밀도, 원소, 기둥밀도, 플래그 | 약 230 B |
| GPU→CPU | 생성률·파괴률 157쌍 | 2512 B |
| GPU→CPU | 우주선 반응률 136개와 분자 계수 | 약 1160 B |

전체 출력 18936 B는 이미 컴팩트 3680 B로 줄여 두었다. 512셀 적분의 왕복은 약 3.1 GB이고 초당 약 2.3 GB다. H200 대역의 일부라, 이 크기를 줄여도 1.02배가 크게 움직이지 않는다.

셀 적분 동안 안 변하는 밀도·원소·기둥밀도·플래그는 셀마다 한 번만 보낼 수 있다. CPU→GPU의 약 15%다. 디바이스 순냉각을 호스트가 그대로 쓰면 우주선 반응률 귀환 약 1160 B를 생략할 수 있다. 함량과 생성·파괴률은 RHS마다 바뀌고 호스트 CVODE가 오차 판정에 쓴다. 이 바이트를 버스에서 빼려면 상태와 `ydot`을 디바이스에 두고 적분기도 거기 있어야 한다. 델타 압축은 공차 `1e-10`에서 수용·거부를 바꿀 수 있다. 하지 말 것.

## 7. CVODE 제어를 GPU로 옮기는 일

시작하지 않았다. SUNDIALS CVODE의 스텝 크기, 오차 판정, Newton 반복은 호스트 루프다. 플래그로 GPU로 가지 않는다. 옮기려면 연산자 분할 시간이 끝날 때까지 BDF를 디바이스에서 돌리고 마지막 함량만 내려받는 적분기가 필요하다.

제어 분기 자체가 빨라지는 구조가 아니다. 이득은 RHS마다 호스트로 돌아오지 않는 것이다. 커널이 셀당 스레드 하나이면 지연은 더 길어진다. 512셀로는 SM이 차지 않는다. 셀마다 수용 스텝 수와 Newton 횟수가 달라 워프는 가장 긴 셀을 기다린다. 큰 레벨에서 그 불균형을 나누어야 수 배를 바라볼 수 있고, 루프만 옮기면 지금과 같은 자리다. 공차를 유지한 채 CPU CVODE와 수용·거부를 맞추는 일은 이번 스케줄러와 별개다.

## 8. 하지 말 것

- 공차를 풀지 말 것. `explicitTolerance==0`이면 explicit Euler 관문이 사라진다. 생산 값은 `1e-10`이다.
- 128³를 이 하이브리드가 빨라졌다는 근거로 돌리지 말 것. 옛 1브로커 128³(작업 413581, 디바이스 냉각 이전, A100)는 수치 회귀는 통과했고 가속은 없었다.
- `scancel -u kjhan`과 `pkill -f`를 쓰지 말 것. 취소 전에 `sacct -j <id> -o JobID,JobName,WorkDir -X`로 WorkDir를 확인할 것. 이 세션이 제출하지 않은 작업은 남의 것이다. 작업 414071은 이 세션의 것이 아니다.
- 로그인 노드에서 무거운 적분을 돌리지 말 것. 로컬 `nvcc -c`는 A10 호스트에서 문법 확인용으로만 했다.
- 스케줄러 폴링은 60초보다 자주 하지 말 것.
- `peak_concurrent_batches`를 고치려 하지 말 것. 금속 합산 순서와 `pow`/`exp` 식을 바꾸지 말 것.
- 하이브리드를 생산 기본값으로 켜지 말 것. CPU가 기본이다.

## 9. 빌드 메모

모듈: `intel/tbb/2022.3 intel/compiler-rt/2025.3.0 intel/umf/1.0.2 intel/compiler/2025.3.0 intel/mpi/2021.17 cuda/13.0.2`.

H200 잡은 `sm_90+sm_80`, `--fmad=false`, `-ffp-contract=off`, `-DCHIMES_USE_DOUBLE_PRECISION`으로 `patch/lagRamses/snrt_chimes_rhs_cuda.cu`를 컴파일한다. 디바이스 `pow10`은 `powten`이다. `chimes_proto.h`를 포함하는 번역 단위는 같은 double 매크로가 있어야 한다.

잡 스크립트: `.cvode-gpu/hybrid-fixture/chimes_hybrid_fixture.sbatch`. 출력은 `chimes-cvode-%j.log`. 이 디렉터리는 커밋하지 않는다.
