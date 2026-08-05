# CUDA-stream 기반 동적 OMP+GPU 작업분배 — 적용 가능 루틴 분석

작성: cuRamses 소스 조사 (`patch/cuRamses/`) 기반. 대상 독자: 개발자.

이 문서는 (1) 현재 하이드로/중력 단계에서 쓰고 있는 "CUDA 스트림
획득 실패 시 CPU로 폴백"하는 실시간 동적 작업분배 메커니즘을 정리하고,
(2) 같은 패턴을 적용할 수 있는 다른 루틴들을 판정 기준과 함께 열거한다.

---

## 1. 현재 메커니즘

### 1.1 스트림 풀

`patch/cuRamses/cuda_stream_pool.{h,cu}` 는 OpenMP 스레드가 동시에
접근할 수 있는 스레드-안전 CUDA 스트림 풀을 제공한다.

- `cuda_pool_init(local_rank, n_streams)` — MPI local rank → GPU 매핑,
  `n_streams`개의 non-blocking 스트림 생성.
- `cuda_acquire_stream()` — 놀고 있는 스트림 슬롯을 원자적
  (`__sync_lock_test_and_set`)으로 하나 잡아 슬롯 번호를 돌려준다.
  모두 사용 중이면 **-1**을 돌려준다 (`cuda_stream_pool.cu:278`).
- `cuda_release_stream(slot)` / `cuda_stream_sync(slot)` — 반납/동기화.
- 각 슬롯은 자기 전용 device 버퍼(hydro I/O, 중간 배열, stencil,
  scatter-reduce 출력)와 pinned host staging 버퍼를 갖는다.

### 1.2 핵심 디스패치 패턴 (the pattern)

동적 분배의 본질은 다음 한 덩어리다 (예: `force_fine.kjhan.f90:249`
`force_gradient_hybrid`).

```fortran
!$omp parallel private(igrid, ngrid, stream_slot, gs)
  stream_slot = cuda_acquire_stream_c()          ! 스레드마다 스트림 획득 시도
  if (stream_slot >= 0) then                      ! 잡았으면 이 스레드는 GPU 워커
     gs => force_gstates(stream_slot)
     call force_gstate_ensure(gs, FORCE_SUPER_SIZE)
     gs%off = 0
  end if

  !$omp do schedule(dynamic)                       ! 배치를 동적 스케줄로 분배
  do igrid = 1, ncache, nvector
     ngrid = MIN(nvector, ncache - igrid + 1)
     if (stream_slot >= 0) then
        call force_gpu_gather_batch(gs, ilevel, icount, igrid, ngrid, stream_slot)  ! GPU 경로
     else
        call gradient_phi(ilevel, igrid, ngrid, icount)                              ! CPU 경로(기존 벡터 루틴)
     end if
  end do
  !$omp end do nowait

  if (stream_slot >= 0) then
     if (gs%off > 0) call force_gpu_flush_scatter(gs, stream_slot, ilevel)
     call cuda_release_stream_c(stream_slot)
  end if
!$omp end parallel
```

동작 원리:

1. `!$omp parallel` 진입 직후 **각 스레드가 스트림 하나를 잡으려 시도**한다.
   스트림이 `n_streams`개면 그 수만큼의 스레드가 GPU 워커가 되고,
   나머지 스레드는 `stream_slot=-1`로 CPU 워커가 된다.
2. `!$omp do schedule(dynamic)`이 배치(`nvector`개 grid 묶음)를 먼저 끝난
   스레드에 계속 넘긴다. GPU 워커는 배치를 빠르게 소화하므로 자연히 더 많은
   배치를 가져가고, CPU 워커는 더 적게 가져간다 → **실행시간 기준 자동 부하분산**.
3. GPU가 없거나(`cuda_pool` 미초기화) 메모리 부족으로 업로드가 실패하면
   모든 스레드가 -1을 받아 **그대로 순수 OMP-CPU 경로로 우아하게 강등**된다.

즉 "몇 %를 GPU에 줄지"를 정적으로 튜닝하지 않는다. 스트림 개수 = GPU 워커
스레드 수이고, dynamic 스케줄이 나머지를 실시간으로 맞춘다.

### 1.3 GPU-상주 메쉬 (gather 계열의 전제)

`cuda_mesh_upload(uold, f, son, ...)` 가 step당 한 번 메쉬를 GPU에 올려두면
(`cuda_stream_pool.cu:374`), gather 계열 커널은 이 상주 배열에서 stencil을
직접 읽는다. `godunov_fine` 진입 시 업로드하고
(`godunov_fine.kjhan.f90:1138`), 하이드로가 끝나면
`cuda_mesh_free_c()`로 반납한다 (`amr_step.jaehyun.f90:762`).

### 1.4 두 가지 GPU 전략 구분

| 전략 | 방식 | 사용처 |
|------|------|--------|
| **A. 스트림 동적분배** | 스레드가 스트림 획득 or CPU 폴백, `schedule(dynamic)` | 본 문서의 대상 (하이드로/중력) |
| **B. 레벨 통째 오프로드** | `use_mg_gpu` 토글로 레벨 전체를 GPU로 (CPU와 배치 공유 없음) | multigrid Gauss–Seidel/residual (`multigrid_fine_commons.f90:352`) |

MG(멀티그리드 Poisson)는 전략 B다. 레드-블랙 Gauss–Seidel이 grid간 의존성이
있어 배치 단위 CPU/GPU 혼합이 까다롭기 때문이다. 본 문서에서 "적용 가능
루틴"이라 할 때는 전략 A를 뜻한다.

---

## 2. 이미 전략 A가 적용된 루틴

| 루틴 | 파일:라인 | 계산 성격 | scatter |
|------|-----------|-----------|---------|
| `godunov_fine` (비분할 하이드로) | `godunov_fine.kjhan.f90:1151` | stencil gather → 리만/flux → scatter-reduce | L, L−1 셀로 (per-thread buffer + serial merge) |
| `force_gradient_hybrid` (φ 기울기 → 중력·제5힘) | `force_fine.kjhan.f90:249` | φ stencil gather → 중심차분 | 없음 (셀 로컬 쓰기) |
| `courant_fine_hybrid` (CFL 시간간격) | `courant_fine.kjhan.f90:129` | 셀별 리덕션(min dt, 보존량) | 리덕션만 |
| `synchro_hydro_hybrid` (중력 소스항 가스속도) | `synchro_hydro_fine.kjhan.f90:90` | 셀 로컬 v 업데이트 | 없음 |
| `upload_fine_hybrid` (프롤롱게이션/보간) | `interpol_hydro.kjhan.f90:110` | 부모→자식 보간 | 없음 |

공통점: 모두 `active(ilevel)%ngrid`를 `nvector` 배치로 도는 구조 +
셀/그리드 로컬 연산 또는 관리 가능한 scatter.

---

## 3. 판정 기준

전략 A를 새 루틴에 적용할 수 있는지는 아래로 판단한다.

1. **배치 루프 구조** — `do igrid=1,ncache,nvector` 형태로 grid(또는 입자)
   묶음을 도는가. (superbatch로 묶어 H2D 1회에 여러 배치를 보내는 게 이득의 핵심)
2. **산술강도 대비 전송량** — 배치당 계산량이 H2D/D2H 바이트를 상쇄할 만큼
   큰가. 리덕션만 하는 저-강도 루틴은 전송이 지배해 이득이 작다.
3. **scatter 충돌** — 출력이 셀-로컬이면 최선. 입자→메쉬처럼 겹쳐 쓰면
   `godunov`의 per-thread buffer + serial merge 패턴을 재사용해야 한다.
4. **메쉬 상주 재활용** — 이미 GPU에 올라간 `uold/f/son`에서 gather만 하면
   추가 업로드 비용이 거의 없다.
5. **분기/발산·테이블 룩업** — warp 내 스레드가 서로 다른 반복횟수/분기를
   타면 divergence로 GPU 효율이 급락한다.
6. **비트-동일 재현성** — 본 프로젝트는 결과 비트-동일을 중시한다. 리덕션
   순서·부동소수 누적을 결정론적으로 유지할 수 있어야 한다.

---

## 4. 후보 루틴 분석

`amr_step.jaehyun.f90`의 레벨별 단계 순서를 따라, 아직 전략 A가 없는
루틴을 판정한다.

### 4.1 강력 후보 (High)

#### `move_fine` — 입자 이동 (drift + force gather)
- 파일: `move_fine.kjhan.f90:146` (`move1`). CIC로 메쉬 `f`를 입자 위치에서
  보간해 `ff`를 만들고(`x0`, `vol`, `indp` 이웃 셀), `new_vp`, `new_xp`를
  입자별로 갱신한다 (`move_fine.kjhan.f90:176`).
- **gather-only, 입자-로컬 출력, scatter 없음.** 구조적으로 가장 깨끗하다.
- 산술강도: 입자당 CIC 8셀 gather + 몇 개 flop. 중간 정도지만 DM 우세
  시뮬레이션에서 입자 수가 압도적이라 총량이 크다.
- 필요 작업: `f`가 이 시점에 GPU 상주여야 한다. 현재 `godunov` 후
  `cuda_mesh_free_c()`로 반납되므로(`amr_step.jaehyun.f90:762`), `f`를
  `move_fine`/`synchro_fine`까지 상주 연장하거나 `move` 직전에 재업로드해야
  한다. `force_gpu_gather_batch`의 stencil 코드를 거의 그대로 재활용 가능.
- **판정: High.** 입자 우세 워크로드에서 가장 큰 이득이 기대된다.

#### `synchro_fine` — 입자 속도 kick
- 파일: `synchro_fine.kjhan.f90:143` (`sync`). `move1`의 부분집합으로,
  메쉬 `f` gather 후 `vp`만 갱신한다. scatter 없음.
- `move_fine`과 GPU 상주 `f`를 공유하므로 함께 구현하는 게 자연스럽다.
- **판정: High** (move_fine와 묶어서).

### 4.2 중간 후보 (Medium)

#### `rho_fine` — 입자→메쉬 밀도 할당 (CIC/TSC)
- 파일: `rho_fine.kjhan.f90:509` (`cic_amr`), `:1600` (`tsc_amr`). 입자를
  이웃 셀에 뿌리는 **scatter**다.
- scatter 충돌이 있어 `godunov`의 per-thread scatter-buffer + serial merge
  패턴이 필수. 패턴 자체는 이미 검증돼 있어 재사용 가능.
- 하위 단계 중 `cic_from_multipole`/`tsc_from_multipole`
  (`rho_fine.kjhan.f90:1258`, `:1954`)와 `multipole_fine`(`:905`)은
  **메쉬→메쉬** 규칙 연산이라 scatter 충돌이 덜하다 → 먼저 오프로드하기 좋다.
- **판정: Medium(–High).** 입자 deposit은 이득이 크지만 scatter 처리
  난이도가 있고, multipole 계열부터 시작하는 게 안전하다.

#### `phi_fine_cg` — Conjugate-Gradient Poisson
- 파일: `amr_step.jaehyun.f90:466`. MG의 대안 solver. matvec(라플라시안
  적용)을 오프로드할 수 있으나, MG가 이미 전략 B로 GPU화돼 있어 **중복
  투자**가 된다.
- **판정: Medium** (MG를 안 쓰는 설정에서만 의미).

### 4.3 약한/조건부 후보 (Low)

#### `cooling_fine` — 냉각/가열 소스항
- 파일: `cooling_fine.kjhan.f90`. 셀별 독립·고산술강도라 겉보기엔 최적
  후보지만, **개발자가 이미 GPU 미구현을 의도적으로 결정**했다
  (`cooling_fine.kjhan.f90:6`):
  > "GPU dispatch NOT implemented because solve_cooling uses iterative
  > Newton-Raphson with variable iteration count per cell (while-loop with
  > active cell tracking), causing severe GPU thread divergence."
- 셀마다 반복 횟수가 달라 warp divergence가 심하다. persistent-thread +
  active-cell 재큐잉 커널로는 가능하나 난이도가 높은 별도 R&D 사안이다.
- **판정: Low (현 구조) / 재설계 시 Medium.**

#### `newdt_fine` — 입자 시간간격
- 파일: `amr_step.jaehyun.f90:700`. 입자별 dt 리덕션. 계산이 가볍고
  리덕션 위주라 전송이 지배한다.
- **판정: Low.** (`courant_fine`가 하이드로용으로 이미 있지만, 그건 셀당
  보존량 계산이 더 무거워서 성립.)

#### `flag_fine` / `refine_fine` — 세분화 표시
- 정수·분기 위주의 불규칙 연산, 트리 포인터 추적. GPU 부적합.
- **판정: Low.**

#### `feedback` (SN/PBH), `sidm_scatter`
- feedback: 별/PBH 입자만 메쉬로 sparse scatter → GPU 점유율 낮음.
- SIDM: 셀 내 쌍(pair) 상호작용 + RNG + 분기, 게다가 OMP 관련 미해결
  이슈 존재. 오프로드 이득 대비 복잡도 과다.
- **판정: Low.**

### 4.4 sink particle 루틴 — 이미 GPU화됨(전략 B), 전략 A는 부적합

결론부터: **sink의 계산-무거운 단계(AGN feedback)는 이미 GPU 가속돼 있다.**
단 본 문서의 전략 A(공유 풀 동적분배)가 아니라 **전략 B(전용 스트림 +
전체 오프로드 + 자동튜닝 토글)** 방식이다.

이미 구현된 것:
- `sink_cuda_kernels.cu`의 `kernel_average_AGN`(`:47`),
  `kernel_AGN_blast`(`:162`) 두 커널. Fortran 진입점은
  `sink_cuda_interface.f90`의 `cuda_sink_average_agn_c`,
  `cuda_sink_agn_blast_c`.
- **전용 스트림** `g_sink_stream`(`sink_cuda_kernels.cu:41`) 사용 —
  공유 풀 `cuda_acquire_stream`을 쓰지 않는다.
- **`gpu_sink` 자동튜닝 토글**(`amr_step.jaehyun.f90:297`) — Phase 0 CPU /
  Phase 1 GPU 벤치 후 Phase 2에서 빠른 쪽 선택. MG의 `use_mg_gpu`와 같은
  전략 B 구조.

핵심은 **병렬화 축**이다. `kernel_AGN_blast`는 sink가 아니라 **셀 단위로
병렬화**한다(`if(tid >= ncells) return`, `sink_cuda_kernels.cu:194`). AGN
데이터는 SoA로 올리고 공간 bin(`bin_head`/`agn_next`)으로 각 셀에 영향을 주는
AGN을 찾는다. 즉 "sink 수가 적다"는 근본 제약을 **sink가 아니라 영향받는
셀들로 병렬화**해 우회했다.

이 때문에 sink에 전략 A는 잘 맞지 않는다:

- **occupancy 문제**: `nsink`은 보통 수십~수천으로 작다. `active(ilevel)%ngrid`를
  `nvector`로 도는 전략 A의 부하분산은 배치가 수천~수백만 개일 때 성립한다.
  sink 작업을 sink 축으로 나누면 GPU가 놀고, 이미 셀 축으로 나눈 AGN 커널이
  올바른 선택이다.
- **국소성**: sink 작업은 전 레벨을 균일하게 쓸지 않고 각 sink 주변 cloud
  셀에 국한된다. 전략 A의 "active grid 균일 sweep" 모델과 맞지 않는다.

나머지 sink 단계 판정:

| 단계 | 파일:라인 | 성격 | 판정 |
|------|-----------|------|------|
| AGN feedback | `sink_cuda_kernels.cu` | 셀-병렬 scatter | ✅ 이미 GPU(전략 B) |
| `create_sink`/`kjhan_make_sink` | `sink_particle.kjhan.f90:105` | 밀도피크 스캔(`do igrid` 있음) + MPI + FoF 그룹핑 | Low (분기·통신 과다, 드물게 발동) |
| `grow_bondi`/`bondi_hoyle` | `amr_step.jaehyun.f90:659,858` | sink별 cloud 셀 강착률, 작은 N | Low (희소 → 필요시 전략 B가 적합) |
| `merge_sink` | `sink_particle.kjhan.f90:1026` | FoF 트리(`Do_FoF_Tree`) | 비대상 (그래프·직렬) |

**권장**: sink에 전략 A를 새로 붙이지 말 것. 무거운 AGN feedback은 이미 셀-병렬
전략 B로 처리돼 있고, 그게 sink 특성에 맞다. sink GPU 커버리지를 더 넓히려면
`grow_bondi`/`bondi_hoyle`를 (전략 A가 아니라) AGN처럼 **셀-병렬 전략 B**로
확장하는 방향이 옳다. 다만 DE/nGR 우주론 실행에서 sink/AGN 비중은 작아
우선순위는 §5의 `move_fine`/`rho_fine`보다 낮다.

### 4.5 비대상

`make_virtual_*`/`make_reverse_*`(MPI 통신), load balance, 트리 빌드, I/O는
계산이 아니라 통신·불규칙 포인터 연산이라 전략 A 대상이 아니다.

---

## 5. 권장 우선순위 (roadmap)

1. **`move_fine` + `synchro_fine`** — 입자 gather. 구조적으로 가장 깨끗
   (scatter 없음), `force_gpu_gather_batch` 재활용, DM-우세 시뮬에서 최대 이득.
   전제작업은 GPU-상주 `f`의 수명을 `move`/`synchro`까지 연장하는 것.
2. **`rho_fine`의 multipole → cic/tsc_from_multipole** — 메쉬→메쉬 규칙
   연산부터 시작해 scatter 파이프라인을 검증한 뒤 입자 deposit으로 확장.
3. **`cooling_fine`** — persistent-thread 재큐잉 커널 R&D. 이득은 크나 별도 과제.

## 6. 엔지니어링 체크리스트 (재사용 자산)

- **스트림 풀**: `cuda_acquire_stream_c` / `cuda_release_stream_c`를 그대로 사용.
  스트림 수만큼 스레드가 GPU 워커가 되고 나머지는 CPU. 튜닝 노브는 `n_streams`.
- **GPU-상주 메쉬**: gather 계열은 `cuda_mesh_upload`가 올린 `uold/f/son`을
  재활용. 입자 단계용으로 `f`의 반납 시점(`cuda_mesh_free_c`)을 뒤로 미룰지 검토.
- **scatter**: 겹쳐 쓰는 출력은 per-thread buffer에 모았다가 병렬영역 밖에서
  serial merge (`godunov_fine.kjhan.f90:1213`의 L−1 merge 패턴).
- **폴백**: 모든 신규 GPU 경로는 `stream_slot < 0`에서 기존 CPU 루틴을 그대로
  호출해야 한다 (GPU 부재/OOM 시 우아한 강등 보장).
- **비트-동일**: 리덕션 순서·부동소수 누적을 결정론적으로 유지.
- **superbatch**: 배치를 `HYBRID_SUPER_SIZE`처럼 묶어 H2D 1회에 여러 배치를
  전송해야 전송 오버헤드가 상쇄된다.

---

## 7. 실측 CPU 점유율 (2026-08-06 벤치마크)

우선순위를 실측으로 검증했다. 세 워크로드 모두 **CPU 빌드**(GPU 미사용)의
내장 타이머(`amr/timer.f90`) 결과다.

- ① **hydro 벤치**: zoom CDM+hydro, 64³→lv9, 8랭크×2OMP, 10 coarse step,
  timer 합계 375 s (`tests/dynamic_particle_cdm.nml`, nstepmax=10).
- ② **DMO 벤치**: ①에서 `hydro=.false.`만 바꾼 동일 조건, 13.5 s.
  DE/nGR 프로덕션(DMO) 워크로드에 해당.
- ③ **FDM 프로덕션 로그**: 128랭크, 9294 s
  (`~/BACKUP/FDM/omp_optimization/runs/r128_t1_401983`).

loadbalance는 테스트 설정(`nremap=1`, 매 스텝 재분배)의 아티팩트라 제외한
**compute-only 점유율**:

| 타이머 (루틴) | ① hydro | ② DMO | ③ FDM | GPU 현황 |
|---|---|---|---|---|
| poisson - mg AMR | **67%** | **41%** | **77%** | ✅ 전략 B (`use_mg_gpu`) |
| poisson - mg base | 1% | **40%** | — | ✅ 전략 B / FFT 대체 |
| hydro - godunov | **14%** | — | 6% | ✅ 전략 A |
| poisson (rho_fine+force_fine) | 7% | 6% | 3% | force ✅ A / rho 후보 |
| refine | 7% | 6% | 0.2% | 비대상 |
| particles (move+synchro) | 0.9% | **4%** | 0.9% | 후보 |
| fdm-drift-fd | — | — | 9% | FDM 전용 |
| cooling | off | — | 0.3% | 후보(Low) |
| sinks / feedback / courant / flag 등 | 각 ≤1% | 각 ≤2% | 각 ≤1% | — |

②의 mg base가 큰 것은 zoom 구성이라 base FFT 경로를 못 타서다.
full-box 주기 경계의 프로덕션 DMO에서는 FFT base가 훨씬 싸져 AMR MG가 남는다.

### 실측이 주는 결론 (§5 우선순위의 수정)

1. **모든 워크로드에서 Poisson MG가 지배** (compute의 41~81%). 이미 전략 B로
   GPU화된 부분이 정확히 급소다. **최우선 행동은 새 전략 A 이식이 아니라,
   프로덕션 실행이 cuRamses GPU 빌드(`use_mg_gpu`)를 실제로 쓰게 하는 것.**
2. 전략 A 후보들의 상한(Amdahl): DMO에서 particles+rho_fine을 전부
   오프로드해도 compute의 ~10%, 전체 속도향상 ≤1.1×. hydro 실행에서
   particles는 ~1%로 더 작다. §5의 `move_fine` 1순위는 구조 기준으로는
   유효하나 **이득 크기는 제한적**이다.
3. godunov(hydro의 14%)는 이미 전략 A 적용됨 — 두 번째 급소도 처리돼 있다.
4. 남는 실질 과제는 (a) MG GPU 경로의 프로덕션 검증, (b) refine/loadbalance
   같은 비대상 오버헤드의 알고리즘적 절감이지, 새 오프로드 대상 발굴이 아니다.

### GPU MG의 물리 모델 커버리지 (FDM / DE / nGR)

GPU MG 커널(`cuda_mg_gauss_seidel`)은 **표준 라플라시안만** 푼다. 모델 의존성은
전부 CPU에서 만든 RHS로 들어오므로 커버리지는 다음과 같이 갈린다.

- **FDM: 커버됨.** FDM 중력은 표준 Poisson 솔버(`multigrid_fine`)를 그대로
  쓴다. GPU MG가 모델 무관하게 적용된다. FDM 고유의 drift/HJM 운동항
  (실측 9%)은 CPU+OMP 전용으로 GPU 미구현.
- **배경 DE(CPL, quintessence, k-essence, Chaplygin, RVM 등): 자동 커버.**
  DE는 Poisson 방정식의 4πG 계수(`cosmo_poisson_fourpi`)와 RHS
  `f(:,2)=fourpi*(rho-rho_tot)` (`multigrid_fine_commons.f90:1738`)로만
  들어온다. cs2≈0 클러스터링 boost도 k-독립 스칼라 곱
  (`multigrid_fine_commons.f90:1698-1707`)이고, k-의존 보정은 base FFT의
  Green 함수에서 처리(`:2777`)되므로 GPU MG 커널은 수정 없이 유효하다.
- **nGR 스칼라장 솔버(f(R), nDGP, symmetron, dilaton, galileon): 미커버.**
  `force_fine.kjhan.f90`의 별도 비선형 Newton–Gauss-Seidel 솔버들로,
  OMP 병렬은 있으나(symmetron 7 / dilaton 6 / galileon 11 pragma)
  **CUDA 경로가 전혀 없다.** nGR 프로덕션에서는 Poisson과 대등한 규모의
  solve가 추가되므로(타이머 fR-solve/nDGP-solve/...), 실측상 유일하게
  남은 대형 GPU 기회다. red–black GS 구조가 표준 MG와 같아 전략 B
  인프라(`cuda_mg_*`)를 비선형 항으로 확장하는 방식이 적합하다.

실무 전제: GPU MG는 cuRamses 빌드(HYDRO_CUDA) + `gpu_poisson` 네임리스트 +
CUDA 장비에서만 작동한다. 지금까지의 FDM/DE 프로덕션 실행은 lagRamses
CPU 빌드였다.

## 8. nGR 스칼라 솔버 GPU 구현 (구현 완료 2026-08-06)

grammar 로그인 노드(A10 3장, nvcc 13.0.2)에서 **구현·빌드·검증 완료**.
결과는 §8.5 참조. 아래 조사 노트의 설계가 그대로 구현되었다.

### 솔버 인벤토리 (`patch/cuRamses/force_fine.kjhan.f90`)

| 솔버 | solve_level | gauss_seidel | 스텐실 |
|------|-------------|--------------|--------|
| f(R) Hu–Sawicki | `:2417` | `:2665` | 6점 face |
| nDGP | `:2887` | `:3094` | **18점** (face 6 + 대각 12) |
| symmetron | `:3303` | `:3496` | 6점 face |
| dilaton | `:3750` | `:3945` | 6점 face |
| galileon | `:4188` | `:4388` | **18점** + `galileon_tracker` 변형 |

### 공통 골격 (solve_level)

warm-start 재스케일 → seed → (`level_fft_ok`면 FFT 전처리) → Newton-GS 루프
{ `*_gauss_seidel`(red-black, OMP) → `make_virtual_fine_dp(scalar_gr)` →
ALLREDUCE 수렴판정 } → save_old → `compute_fifth_force`(기울기, 계수).

- red-black 컬러: `popcnt(ind-1)` 홀짝 — 표준 MG GPU 커널과 동일 구조.
- 수렴량: res_max/src_max의 **max-리덕션** (Poisson MG의 norm2 합-리덕션과 다름).
- 6점 솔버의 이웃: `vain_face_grid`(active grid×6, `vain_prepare_uniform_cache`가
  사전계산) — 그대로 GPU 업로드 가능.
- 18점 솔버의 이웃: `morton_nbor_grid` + `scalar_sample_offset`(Morton 해시
  룩업) — GPU에서 해시 대신 **CPU 사전계산 이웃-셀 인덱스 테이블** 업로드로 대체.

### 핵심 설계 결정 (조사에서 확정)

1. **coarse-fine 경계는 레벨 solve 동안 상수다.** `scalar_sample_offset`의
   Dirichlet 값은 부모(ilevel−1) CIC 보간(가중치 0.75/0.25)이고 부모 레벨은
   이 solve 동안 고정 → **CPU에서 사전계산해 한 번 업로드**하면 GPU 커널은
   해시/트리 순회가 전혀 필요 없다. (부모 셀 미존재 시 fallback=u_c
   zero-gradient 예외만 플래그로 전달.)
2. **cuda_mg 인프라 재사용**: 스트림/업로드/halo 교환
   (`make_virtual_fine_dp_gpu`, `cuda_mg_halo_*`)과 `gpu_poisson` 토글 패턴을
   `d_scalar` 배열용으로 확장. halo 셀 리스트는 Poisson MG가 쓰는 것과 동일한
   virtual-boundary 구조.
3. **rho(소스)는 read-only** — solve당 1회 업로드.
4. **비트-동일 주의**: res_max/src_max max-리덕션은 순서 무관이라 안전.
   Newton 업데이트 자체는 스윕 내 갱신 순서에 의존(GS) — red-black이 같으면
   색 내부 순서는 수학적으로 무관(같은 색끼리 독립)이므로 GPU 병렬화 유효.
5. **구현 순서**: 6점 3종(fR/symmetron/dilaton, 공통 커널 골격 + 모델별
   소스/야코비안 함수) → 18점 2종(nDGP/galileon) 확장.

### 8.5 구현 및 검증 결과 (2026-08-06)

**신규/수정 파일** (커밋 전, 작업 트리):

| 파일 | 내용 |
|------|------|
| `patch/cuRamses/scalar_cuda_kernels.cu` (신규, ~620줄) | 5모델 공용 Newton-GS 커널(`scal_gs_kernel`, 스레드=grid, 색당 4셀), 비트연산 이웃 해석, halo gather/scatter 커널, atomicMax res/src 리덕션 |
| `patch/cuRamses/scalar_cuda_interface.f90` (신규) | bind(C) 인터페이스 + `scalar_gpu_commons`(호스트 테이블·halo 리스트) |
| `patch/cuRamses/force_fine.kjhan.f90` | 헬퍼 5개(`scalar_gpu_begin/sweep_halo/end`, `build_scalar_halo_indices`, `make_virtual_scalar_gpu`, `scalar_lookup_icell`) + 5개 solve_level 배선(`gscal_ok` 분기, CPU 폴백 유지) |
| `patch/{lagRamses,cuRamses}/{amr_parameters,read_params}.jaehyun.f90` | `gpu_scalar` 네임리스트 파라미터(&RUN_PARAMS, 기본 .false.) |
| `patch/cuRamses/adaptive_loop.jaehyun.f90` | CUDA 풀 조기 초기화 게이트에 `gpu_scalar` 추가 |
| `bin/Makefile`, `patch/cuRamses/cuda_stream_pool.cu` | 오브젝트 추가, finalize 훅 |

빌드: CPU(`make USE_FFTW=1`)·GPU(`make USE_CUDA=1 USE_FFTW=1`) 모두 클린 통과
(직렬 빌드 필수 — Makefile은 `-j` 모듈 레이스가 있음). 바이너리:
`bin/ramses_gpuscal3d`(GPU), `ramses_final3d`(=CPU, 동일 소스).

**검증** (64³→lv9 zoom DMO, 8랭크×2OMP, 4스텝, step4 입자 스냅샷 비교):

| 항목 | fR (6점) | nDGP (18점 Vainshtein) |
|------|----------|------------------------|
| 수렴 이력 (iters, res) | 사실상 동일 (다수 잔차 완전 일치) | FFT base 완전 일치, GS 레벨 미수렴 잔차까지 동일(1.364E+02) |
| 입자 위치 max\|dx\| | **2.2e-16** (머신 정밀도) | 1.4e-6 |
| 입자 속도 상대오차 | **2.0e-15** | 1.1e-3 |
| 벽시계 (CPU→GPU) | 4.7→6.0 s (문제 과소, GPU 오버헤드) | **20.7→6.6 s (3.1×)** |

해석: 6점 솔버(fR/symmetron/dilaton)는 red-black이 정확한 컬러링이라 GPU가
CPU와 머신 정밀도로 일치한다. 18점 Vainshtein(nDGP/galileon)은 같은 색
대각 이웃의 갱신 순서가 구현마다 달라(CPU 순차 vs GPU 병렬) 반복 상태가
1e-3 수준에서 갈리며, 이는 수렴 허용오차가 지배하는 영역 안이다(두 구현의
미수렴 잔차가 동일 값으로 포화). 이 순서 의존성은 CPU OMP 버전에도
이미 존재한다.

**5모델 전수 A/B 완료** (symmetron/dilaton은 A10, galileon은 A100 syn101):

| 모델 | max\|dx\| | 상대 dv | 비고 |
|------|-----------|---------|------|
| symmetron | 4.8e-15 | 3.3e-13 | 머신 정밀도 |
| dilaton | 1.1e-16 | 1.7e-14 | 머신 정밀도 |
| galileon | 4.2e-13 | 2.6e-10 | 18점, 정상 |

galileon 첫 시도에서 관측된 큰 차이(rel dv~9)는 커널 결함이 아니라 legacy
모드 z≈48의 병리적 영역(|coeff|~6e6, tracker가 solve 자체를 거부하는 구간)에서
CPU/GPU 모두 미수렴했기 때문이다. c3=1e-9로 온건화하면 위 표대로 일치한다.

**A100 강제반복(eps=1e-30, 100 sweeps) 스케일링**:

| 솔버 타이머 | CPU | GPU | 속도향상 |
|------------|-----|-----|---------|
| fR-solve (6점) | 2.95 s | 1.62 s | 1.8× |
| nDGP-solve (18점) | 31.09 s | 2.88 s | **10.8×** |
| galileon 전체 실행 | 289 s | 21 s | **13.7×** |

스텐실이 무거운 18점 솔버에서 이득이 압도적이다. nGR 프로덕션의 비용은
대부분 이쪽이므로 실용적 의미가 크다. 6점 솔버의 1.8×는 64³ 장난감 크기 +
경량 스텐실이라 전송 오버헤드 비중이 큰 탓이다.

미해결 관찰: a40 파티션에서 강제반복 GPU 벤치가 SIGKILL로 죽었다
(host MaxRSS 2.4 GB로 OOM 아님). syn101(A100) 수동 실행에서는 재현되지
않아 코드가 아니라 a40 환경(공유 GPU 경합 또는 slurm 스텝 제한) 쪽으로
보인다. a40에서 프로덕션을 돌릴 때 재확인할 것.

사용법: `&RUN_PARAMS`에 `gpu_scalar=.true.` (USE_CUDA 빌드에서만 유효,
아니면 경고 후 무시). GPU 부재/OOM/랭크 불일치 시 자동 CPU 폴백.

## 9. move_fine/synchro_fine 전략 A 구현 (2026-08-06, 검증 진행 중)

§5의 1순위 항목을 구현했다. **진짜 전략 A** — 스레드가 스트림을 잡으면 GPU
워커, 못 잡으면 CPU 워커(`move1`/`sync` 그대로), `schedule(dynamic)`으로
실시간 부하분산.

| 파일 | 내용 |
|------|------|
| `patch/cuRamses/particle_cuda_kernels.cu` (신규) | 통합 CIC 커널(스레드=입자, move/sync 모드): ilevel CIC → 8그리드 존재 검사 → 미존재 시 ilevel−1 CIC 폴백 → f gather → kick/drift, `OUTPUT_PARTICLE_POTENTIAL` 지원; 메쉬(f/son/phi) call당 1회 업로드; 슬롯별 pinned superbatch(16384 입자) |
| `patch/cuRamses/particle_cuda_interface.f90` (신규) | bind(C) + `pm_gpu_commons`(슬롯 버퍼) + `pm_gpu_dispatch`(append/flush) |
| `move_fine`/`synchro_fine` | 그리드 벡터화 + 하이브리드 dispatch (기존 CPU-전용 경로는 `#ifndef HYDRO_CUDA`로 보존) |

설계 결정:
- **트리 작업은 CPU에 유지**: `get3cubefather`/x0는 배치 패킹 단계(CPU)에서
  수행해 커널은 순수 산술 + 인덱스 읽기만 한다.
- **연산 순서 재현**: 입자별 힘 누적(ind 외측, dim 내측)과 CIC 산술을 CPU와
  동일 순서로 구현 — LCDM(coupled-DE off)에서는 비트-동일이 기대치.
- **안전장치**: CUDA 실패 시 append된 배치 기록(ind_grid/ind_part)을 그대로
  `move1`/`sync`로 재실행(replay)하고 GPU 경로를 비활성화. sync 모드는
  levelp 복원 후 재실행.
- **게이트**: `gpu_particle` 네임리스트; sink/tracer 실행은 자동 CPU 경로.
  move는 랭크-로컬이라 MPI 합의 불필요(랭크별 독립 선택 가능).

### 9.1 rho_fine deposit (입자→메쉬 CIC)

같은 인프라를 재사용해 `cic_amr` deposit도 GPU화했다. scatter 충돌은
**device 누산기 + atomicAdd**로 처리한다.

- `pm_deposit_kernel`: 입자당 8개 cloud 셀에 `atomicAdd(rho, m*vol/vol_loc)`
  와 `atomicAdd(phiw, vol)`. cloud 그리드가 없는 셀은 CPU `cic_amr`과
  동일하게 건너뛴다(coarse 폴백 없음).
- GPU 워커는 device 누산기에, CPU 워커는 기존 `cic_amr`(critical 안에서)
  호스트 배열에 각각 누적하고, dispatch 루프 종료 후 `pm_rho_merge`가
  device 누산기를 호스트 rho/phi에 더한다.
- **multipole 이중계산 방지**: 레벨 최저에서 pack 단계가 multipole을
  누적하므로, CPU replay 경로는 해당 배치 몫을 먼저 빼고 `cic_amr`을 부른다.
- 게이트: star/sink 실행과 `cic_levelmax` 특수 레벨, TSC 빌드는 CPU 경로.
  (그 조건들에서 CPU 코드가 rho_top/phi에 추가 분기를 타기 때문.)
- 기존 CPU 전용 구현은 `#if defined(_OPENMP) && !defined(HYDRO_CUDA)`로 보존.

### 9.2 검증 결과 (A100 syn101, LCDM DMO 10스텝, `gpu_particle=.true.`)

move + synchro + rho를 **동시에** 켠 A/B:

| 스냅샷 | max\|dx\| | max\|dv\| | 판정 |
|--------|-----------|-----------|------|
| output_00001 (초기) | **0** | **0** | 비트-동일 |
| output_00002 (10스텝 후) | 4.4e-16 | 7.1e-18 (상대 5.1e-15) | 머신 정밀도 |

`gather` 계열(move/synchro)은 연산 순서를 CPU와 동일하게 맞춰 비트-동일이
나온다. 10스텝 후의 1 ulp 차이는 rho deposit의 atomicAdd 누적 순서에서
비롯되며, 예고한 트레이드오프 그대로다(비트 재현이 필요한 검증 실행은
`gpu_particle=.false.`로 돌릴 것).

### 9.3 최초 측정에서 드러난 고정비용 문제

첫 측정(최적화 전): 벽시계 CPU 10.4 s vs GPU 21.8 s.

타이머로 분해하면 원인이 갈린다.

| 구간 | CPU | GPU | 성격 |
|------|-----|-----|------|
| 벽시계 | 10.4 s | 21.8 s | — |
| 타이머 총합(메인 루프) | 10.7 s | 14.2 s | — |
| 차이(=미계측) | −0.3 s | **+7.6 s** | CUDA 컨텍스트 생성·대용량 할당·종료, **1회성** |
| particles 단계 | 0.47 s | 1.78 s | 호출당 메쉬 업로드 |

7.6 s는 실행 길이와 무관한 1회성이라 프로덕션에서는 무시된다. 실제 문제는
**호출당 고정비용**이다.

```
ngridmax = 500000  →  ncell ≈ 4,000,000 셀
  f 96 MB + phi 32 MB + son 16 MB  =  호출당 144 MB
그런데 실제 사용 grid는 랭크당 5,000~8,000개 (로그 own=5194)
  → 업로드의 약 98%가 미사용 패딩
```

`schedule(dynamic)`은 **배치당** 비용만 흡수한다. GPU가 느리면 배치를 덜
가져갈 뿐이다. 그러나 메쉬 업로드는 병렬영역 진입 **전에**, 배치를 하나도
가져가지 않아도 지불되므로 동적 분배가 손댈 수 없다. 전략 A가 "최악의 경우
CPU와 동등"을 보장하려면 이런 고정비용이 없어야 한다.

### 9.4 최적화 (2026-08-06)

1. **고수위(high-water) 슬라이스 업로드** — 셀 인덱스는
   `ncoarse + (ind-1)*ngridmax + igrid` 구조다. `ngridmax`는 할당 상한일 뿐이므로
   coarse 블록 + 각 oct 슬롯의 `[1..hw]`만 복사한다. `hw`는 모든 레벨·cpu의
   그리드 리스트(+경계)를 한 번 훑어 구한 **정확한** 최대 인덱스이며, 비용은
   O(사용 그리드)로 무시할 만하다. 업로드가 약 1/50로 줄어든다.
2. **작업량 문턱값** — `pm_gpu_min_part`(&RUN_PARAMS, 기본 20000). 해당 레벨의
   입자가 이보다 적으면 업로드조차 하지 않고 CPU 경로로 간다. 업로드가
   상환되지 않는 구간에서 GPU를 강제하지 않는다.
3. **필요한 것만** — synchro는 phi를 올리지 않고, rho deposit은 son만 올린다.

### 9.5 최적화 후 재측정 (A100 syn101, 동일 조건 2회 반복)

particles 단계(move+synchro)의 CPU 대비 초과분:

| | CPU | GPU | 초과 |
|---|---|---|---|
| 최적화 전 | 0.468 s | 1.783 s | **+1.32 s** |
| 최적화 후 rep1 | 0.382 s | 0.574 s | +0.19 s |
| 최적화 후 rep2 | 0.587 s | 0.629 s | +0.04 s |

**호출당 오버헤드가 1.32 s → 0.04~0.19 s로 약 7~30배 줄었다.** 고정비용이
동적 분배의 사각지대였다는 진단이 맞았고, 슬라이스 업로드가 그것을 제거했다.

노드 경합에 대한 주의: `poisson - mg AMR`은 `gpu_poisson=F`라 이 경로와 무관한데도
실행마다 4.17/5.53/6.29/5.95 s로 크게 흔들린다(같은 CPU 실행끼리도 4.17 vs 6.29).
syn101을 다른 사용자가 동시에 쓰고 있어서다. 따라서 **총합이나 벽시계로
CPU/GPU를 비교하면 안 되고**, 단계별 타이머를 봐야 한다. rep1에서 mg AMR을
제외한 메인 루프는 CPU 4.40 s vs GPU 4.40 s로 사실상 동일했다.

벽시계는 여전히 GPU가 5~6 s 길다(13.1/15.9 vs 18.6/19.9). 이는 §9.3의
1회성 CUDA 컨텍스트·할당·해제 비용이며 실행 길이와 무관하므로 프로덕션에서는
상환된다.

정확성 재검증(슬라이스 업로드가 셀을 누락하지 않았는지 확인하는 의미도 있다):
초기 스냅샷 **비트-동일**, 10스텝 후 4.4e-16 — 최적화 전과 동일한 수준이다.

### 9.6 기본값 정책

고정비용을 제거한 뒤 `gpu_particle`의 기본값을 **`.true.`로 전환**했다.
프로덕션이 대형 시뮬레이션을 향하고 있고, 크기가 작을 때는 문턱값이 알아서
CPU 경로로 보내기 때문이다. 판단은 두 층으로 이뤄진다.

- `pm_gpu_min_part`(기본 100000): 레벨의 입자가 이보다 적으면 업로드조차 하지
  않는다. A100에서 레벨당 랭크당 2만~4만에서 겨우 손익분기였으므로 기본값을
  그보다 넉넉히 위에 뒀다. 프로덕션 구성에서 교차점을 측정해 조정할 것.
- USE_CUDA 없이 빌드하면 경고 없이 조용히 꺼진다(기본값이 켜져 있으므로
  경고를 내면 모든 CPU 실행이 시끄러워진다).

즉 작은 실행은 문턱값에서 걸러지고, 큰 실행에서만 GPU가 관여한다.

### 9.6b 교차점 계측 결과와 `pm_gpu_min_part` 근거 (A100 syn101)

`cuda_pm_report`(실행 종료 요약)가 찍은 랭크 1의 비용 분해:

```
uploads : 60 calls, 0.264 s, 4.394 ms/call, 2.9 MB/call
flushes : 80 calls, 0.041 s, 673064 particles
per particle: 61.7 ns (GPU 경로)
```

**슬라이싱 효과 확인**: 호출당 업로드가 144 MB → **2.9 MB (약 1/50)**로 줄었다.
예측대로다.

**그런데 GPU 경로는 여전히 진다** (particles 타이머 CPU 0.416 s vs GPU 0.726 s).
분해하면 이유가 분명하다.

| 항목 | 시간 |
|------|------|
| 업로드 | 0.264 s |
| 커널+전송(flush) | 0.041 s |
| 나머지(=CPU측 트리 작업 + 패킹) | 0.421 s |

CPU 경로의 전체가 0.416 s인데 GPU 경로의 **CPU측 잔여 작업만으로 이미 0.421 s**다.
즉 **`get3cubefather` 트리 순회가 이 단계의 지배적 비용이고, GPU로 옮긴 CIC
산술은 작은 몫**이다. 산술을 아무리 빠르게 해도 상한이 낮다.

구조적 이유: 이 구성은 그리드당 입자가 약 6.7개(랭크당 입자 5.7만 / 그리드
8533)로 매우 희박하다. 트리 작업은 **그리드 수**에, 산술은 **입자 수**에
비례하므로, 그리드당 입자가 적으면 트리가 이긴다. 따라서 문턱값은 원리적으로
"레벨당 입자 수"보다 **"그리드당 입자 수"**가 더 옳은 지표다(현재는 전자를
프록시로 쓴다).

업로드 2.9 MB에 4.394 ms는 **대역폭이 아니라 지연** 때문이었다(실효 660 MB/s).
성분마다 coarse 블록 + 8개 oct 슬라이스를 따로 복사해 호출당 memcpy가 45회였다.
8개 슬롯이 `ngridmax`로 균등 스트라이드되므로 `cudaMemcpy2D` 한 번으로 묶어
호출당 10회로 줄였다.

| 업로드 | 바이트/call | 시간/call |
|--------|-------------|-----------|
| 최적화 전(full ncell) | 144 MB | — |
| 슬라이싱 후(memcpy 45회) | 2.9 MB | 4.394 ms |
| **+ `cudaMemcpy2D`(10회)** | 2.9 MB | **1.358 ms** |

전송 비용은 이렇게 사실상 제거됐다(합계 0.264 s → 0.082 s). 그런데도 particles
단계는 여전히 CPU가 빠르다(0.409 s vs 0.727 s). 위 분해가 말해주듯 **병목이
전송이 아니라 CPU측 트리 순회**이기 때문이며, 전송을 더 줄여도 이 결론은
바뀌지 않는다. 이 경로가 이기려면 그리드당 입자가 훨씬 조밀해지거나
`get3cubefather` 자체를 손봐야 한다.

주의: 이 노드는 다른 사용자와 공유 중이라 단계 벽시계가 실행마다 크게 흔들린다
(동일 구성 CPU `particles`가 0.382~0.800 s). 위 결론은 흔들리지 않는 계측값
(업로드/커널 격리 측정)에 근거한다.

기본값 `pm_gpu_min_part=100000`은 이 구성(랭크당 5.7만)을 걸러낸다. 위 계측
도구가 실행마다 출력되므로 프로덕션에서 같은 방식으로 교차점을 재서 조정할 것.

### 9.7 GPU 스위치 기본값 정리

| 스위치 | 기본값 | 근거 / 안전장치 |
|--------|--------|-----------------|
| `gpu_hydro` (godunov, 전략 A) | **`.true.`** | `gpu_auto_tune`이 CPU/GPU를 벤치마크해 빠른 쪽 채택 |
| `gpu_sink` (AGN, 전략 B) | **`.true.`** | 동일한 auto-tune |
| `gpu_scalar` (nGR 솔버, 전략 B) | **`.true.`** | 1.8×(f(R)) ~ 13.7×(Galileon) 측정, GPU 없으면 CPU 폴백 |
| `gpu_particle` (CIC, 전략 A) | **`.true.`** | `pm_gpu_min_part` 문턱값 |
| `gpu_poisson` (fine AMR MG, 전략 B) | `.false.` | 비중이 가장 크고(compute의 41~81%) 구조도 유리하다(업로드 1회를 V-cycle의 수십~수백 스윕이 상환). 남은 과제는 §9.8의 업로드 낭비와 `mg_merged_rb` 재현성 검증 |
| `gpu_fft` (base cuFFT) | `.false.` | **분산 FFT는 MPI all-to-all 전치가 지배**한다. 연산을 GPU로 옮겨도 통신이 병목이라 이득이 없다(실측 확인). base/coarse는 CPU FFTW3가 담당한다 |

### 9.8 Poisson 솔버 역할 분담 (오해 방지)

- **base/coarse 레벨**: FFT 직접 해법이 담당한다(우선순위 FFTW3(CPU) > cuFFT >
  MG V-cycle). 따라서 `recursive_multigrid_coarse`가 CPU에 남는 것은 실질적
  쟁점이 아니다. 그리고 위 이유로 이 FFT는 CPU에 두는 것이 맞다.
- **fine AMR 레벨**: MG V-cycle이 담당하고, `gpu_poisson`이 이 부분을 통째로
  오프로드한다. `phi/f(:,1..3)/flag2`를 레벨 solve 시작에 한 번 올려 V-cycle
  내내 device에 상주시키고, 스윕 사이에는 halo 셀만 호스트를 경유한다
  (전체 배열 왕복은 `use_ri_gpu`가 꺼진 경우에만 발생하며 정상 경로에서는
  `use_ri_gpu = use_mg_gpu`라 타지 않는다).
- **남은 최적화**: `cuda_mg_upload`가 여전히 full `ncell`을 복사한다
  (phi+f1+f2+f3+flag2 = 셀당 36 B → 테스트 구성에서 레벨 solve당 144 MB).
  실사용 grid가 5,000~8,000개뿐이므로 §9.4의 고수위 슬라이싱을 그대로 적용할
  수 있다. MG는 상환 구조라 particle 경로만큼 치명적이지는 않다.

플래그가 아예 없는 전략 A 루틴들(`courant_fine`, `synchro_hydro_fine`,
`interpol_hydro`, `force_gradient_hybrid`)은 CUDA 빌드에서 항상 활성이며
스트림 획득 실패 시 CPU로 폴백한다.

USE_CUDA 없이 빌드하면 위 스위치는 모두 조용히 꺼진다.

## 요약 표

| 루틴 | 성격 | scatter | 판정 |
|------|------|---------|------|
| godunov_fine | 하이드로 stencil | 있음(관리됨) | ✅ 적용됨 |
| force_fine (gradient_phi) | 중력/제5힘 gather | 없음 | ✅ 적용됨 |
| courant_fine | CFL 리덕션 | 리덕션 | ✅ 적용됨 |
| synchro_hydro_fine | 가스속도 소스항 | 없음 | ✅ 적용됨 |
| upload_fine | 보간 | 없음 | ✅ 적용됨 |
| **move_fine** | 입자 force gather+drift | 없음 | **High** |
| **synchro_fine** | 입자 속도 kick | 없음 | **High** |
| rho_fine (multipole/from_multipole) | 메쉬→메쉬 | 적음 | Medium |
| rho_fine (cic/tsc_amr deposit) | 입자→메쉬 | 있음 | Medium–High |
| phi_fine_cg | CG Poisson | — | Medium(중복) |
| cooling_fine | 셀별 Newton-Raphson | 없음 | Low(발산) |
| newdt_fine | 입자 dt 리덕션 | 리덕션 | Low(저강도) |
| flag/refine | 세분화 | — | Low(불규칙) |
| feedback / sidm | sparse/pairwise | 있음 | Low |
| **sink: AGN feedback** | 셀-병렬 scatter | 있음 | ✅ 이미 GPU(전략 B) |
| sink: create/bondi/merge | 작은 N·분기·그래프 | 있음 | Low (전략 A 부적합) |
