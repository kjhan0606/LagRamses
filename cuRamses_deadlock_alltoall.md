# Multigrid `MPI_ALLTOALL` Deadlock in `build_parent_comms_mg`

> VoidSim base-12 두-유체 void zoom 프로덕션 런이 진화 초반(z≈45)에서
> **영구 정지(hang)** 한 사건의 원인 분석. OOM이 아니라 **dense collective 데드락**.

---

## 1. 데드락 위치 (file / line)

**컴파일된 소스:** `patch/lagRamses/multigrid_fine_commons.f90`
(VPATH 최우선이라 `patch/cuRamses/multigrid_fine_commons.f90`가 아닌 lagRamses 사본이 링크됨)

**루틴:** `subroutine build_parent_comms_mg(active_f_comm, ifinelevel)` — 1058~1542행
(multigrid Poisson 솔버가 코어스 레벨 부모 통신자를 구성하는 단계)

**정지한 호출:** 카운트 발견(count-discovery) 단계의 `MPI_ALLTOALL` 2곳 중 하나

```fortran
! --- STAGE 2: 코어스 그리드 MG 활성화 요청 (1175~) ---
1177   if(ordering=='ksection') then
1178      ...                                  ! sparse k-section 트리 교환
1193      call ksection_exchange_dp(sbuf_ksec, ntotal_ksec, dcpu_ksec, 2, rbuf_ksec, nrecv_ksec)
1199   else
1200      call MPI_ALLTOALL(nreq, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, &   ! ← 데드락
1201         & MPI_COMM_WORLD, info)
1202   end if
```

동일 패턴이 **STAGE 5**(새 수신 그리드 요청 공유)에도 반복:

```fortran
1356   else
1357      call MPI_ALLTOALL(nreq2, 1, MPI_INTEGER, recvbuf2, 1, MPI_INTEGER, &   ! ← 동일 위험
1358         & MPI_COMM_WORLD, info)
1359   end if
```

두 Alltoall 모두 `nreq`/`nreq2`(rank별 요청 건수) 정수 1개씩을 **전 rank가 전 rank와**
교환하는 dense global collective (`MPI_COMM_WORLD` 전체).

---

## 2. 호출 경로 (stack trace, `eu-stack -p <pid>`)

> grammar 컴퓨트 노드에는 gdb/gstack/perf가 없음. `eu-stack`은 있음.

```
#0  0x...2e6
#1  MPIDI_Progress_test            ← Intel MPI busy-wait (유저공간 spin)
#2  MPID_Progress_wait
#3  MPIR_Wait_impl
#4  MPIC_Wait
#5  MPIC_Sendrecv
#6  MPIR_Alltoall_intra_brucks     ← Bruck 알고리즘(소형 메시지 Alltoall)
#10 PMPI_Alltoall
#11 pmpi_alltoall_
#12 build_parent_comms_mg_         ← 데드락 루틴
#13 multigrid_fine_                ← Poisson V-cycle (multigrid_fine_commons.f90:134/137에서 호출)
#14 amr_step_                      ← 레벨별 재귀 (8단 이상 깊이)
```

`multigrid_fine`은 `amr_step`이 레벨마다 호출하는 Poisson 솔버이고, 그 안에서
finer 레벨 1회(134행) + 코어스 레벨 루프(136~138행)로 `build_parent_comms_mg`를
여러 번 부른다. 각 호출이 위 Alltoall을 친다.

---

## 3. 증상 (어떻게 보였나)

| 항목 | 관측값 | 정상 기댓값 |
|---|---|---|
| SLURM 상태 | `RUNNING` 1일 18시간+ | 정상 진행 또는 종료 |
| 로그 mtime | **~40시간 동결** (Fine step 190에서 멈춤) | `ncontrol=1`이면 매 fine step 출력 |
| `.err` | **빈 파일** (0 byte) | 크래시면 에러 메시지 |
| rank 프로세스 | 6/노드 전부 `Rl`, `wchan=-`, **99.3% CPU** | — |
| `/proc/<pid>/syscall` | `running` (커널 대기 아님) | — |
| 노드 CPULoad | 40~52 / 64 (OMP 스레드도 churn) | — |
| 마지막 로그 | `Fine step=190  Main step=9  level=16  a=2.195E-02` | — |

**판정 근거:** 24-rank가 정수 1개씩 주고받는 `MPI_ALLTOALL`은 정상이면 마이크로초.
1.5일째 100% CPU로 진전 0 = **collective에 전 rank가 모이지 못한 데드락**.
Intel MPI는 collective 대기 시 유저공간 busy-spin이라 "바빠 보이지만 막힌" 상태
(데드락인데 CPU 99%, `wchan=-`)로 위장됨 — OOM/segfault와 혼동 주의.

OOM과의 구별: OOM이면 `refine_utils.f90`가 `No more free memory  32  <numbf>`를
찍고 죽음. 이번엔 그 메시지 **없음** → 메모리 문제 아님.

---

## 4. 왜 데드락이 발생하나 (메커니즘)

`MPI_ALLTOALL`은 **모든 rank가 동시에 진입해야** 완료되는 동기 collective다.
`MPI_COMM_WORLD` 전체가 참여하므로, **단 하나의 rank라도 늦게 도착하면
나머지 전원이 그 rank를 무한정 기다린다.**

void zoom의 부하 불균형이 그 straggler를 만든다:

1. `ivar_refine=6` 마스크가 AMR을 void footprint(레벨-9 윈도의 약 **1.25%**)에만 가둠.
2. 기본 `ordering='hilbert'` 도메인 분할에서, 이 작은 footprint가 **소수 rank**에만 매핑됨.
3. 그 rank들은 레벨 16까지 refine된 거대한 작업을 떠안고, 나머지 다수 rank는
   해당 레벨에서 셀이 거의/전혀 없음(빈 도메인).
4. multigrid V-cycle의 `build_parent_comms_mg`에 도달하는 시점이 rank마다 크게 어긋남
   → 무거운 rank가 도착하기 전, 가벼운 rank 23개가 1200행 Alltoall에서 영구 대기.
5. `nremap=10`(10 coarse step마다 부하 재분배)인데, **첫 remap(step 10) 직전인
   Main step 9에서 멈춤** → 재분배가 한 번도 작동하기 전에 불균형이 치명화.

요약: **dense global collective + 극단적 부하 불균형(소수 rank에 footprint 집중) =
straggler 데드락.** 통신 자체의 데이터량은 1 정수로 사소함 — 문제는 *동기화 지점*.

---

## 5. autotune이 막지 못한 이유

`xchg_autotune_update` (`patch/cuRamses/virtual_boundaries.kjhan.f90:1981`)는
**AMR virtual boundary 교환**(icomp 1~7: `fine_dp/fine_int/reverse_dp/reverse_int/
pair_int/bulk_dp/bulk_rev_dp`)에 대해 `exchange_method=='auto'`일 때 P2P↔ksection을
호출별 `MPI_WTIME` 측정으로 적응 선택한다(Phase 0 P2P 측정 → Phase 1 ksection 측정 후
결정 → Phase 2 EMA 추적 → Phase 3 주기 probe, `XCHG_SWITCH_MARGIN` 넘으면 스위치).
두 방법 모두 sparse라 dense collective가 없다.

그러나 데드락 난 `build_parent_comms_mg`는 **multigrid 솔버 내부의 독립 통신 경로**로
**autotune 범위 밖**이다. 이 카운트 발견은 여전히
`if(ordering=='ksection') … else MPI_ALLTOALL` 분기에 남아 있어,
hilbert ordering에선 보호 없는 dense Alltoall로 떨어진다.

| 통신 사이트 | 보호 | 비고 |
|---|---|---|
| AMR virtual boundary (icomp 1~7) | ✅ autotune P2P/ksection | sparse, dense 없음 |
| multigrid `build_parent_comms_mg` 카운트 발견 | ❌ **무방비** | 1200/1357행 dense Alltoall |

---

## 6. 시뮬레이션 상세 특징 (재현 조건)

- **종류:** genetIC 역-클러스터 nested-zoom void, 두-유체(DM+baryon) 하이드로, lagRamses/cuRamses
- **그리드:** `levelmin=8`, `levelmax=19`; grafic 5단계 입력
  `void_N.grafic_{256,512,1024,2048,4096}` = 레벨 8~12
  (256³ 전 박스 + 64³/128³/256³/512³ 윈도, 윈도 = 박스의 1/8/변)
- **refine:** `ivar_refine=6` + `var_cut_refine=0.01` → `ic_pvar_00001` footprint 마스크
  (윈도의 ~1.25%만 AMR). `m_refine=20*8.`, `err_grad_d=0.3`
- **분할:** 기본 `ordering='hilbert'` (← 데드락 트리거). `nremap=10`
- **자원:** 4 노드 × 6 rank = **24 rank**, OMP=10, grammar `normal` 파티션 (64core/515GB)
- **하이드로:** `slope_type=2`(MC), `scheme='muscl'`, `riemann='hllc'`, `pressure_fix`
- **정지 시점:** `Fine step=190`, `Main step=9`, AMR `level=16`,
  `a=2.195E-02` → **z ≈ 44.6** (목표 z=0). mem 6.9%(여유 충분 — OOM 아님 재확인)
- **사건 job:** 316173 (kill 후 320372로 `ordering='ksection'` 재투입)

---

## 7. 조치

- **임시(검증용):** `&RUN_PARAMS`에 `ordering='ksection'` 추가.
  → `build_parent_comms_mg`가 else 분기를 피해 1193/1351행
  `ksection_exchange_dp`(O(log_k ncpu) point-to-point 트리)로 라우팅 → dense Alltoall 제거.
  (단, 도메인 분할 자체가 k-section으로 바뀜 — `ordering='ksection'`(DD)는
  네트워크 ksection 교환과 별개 개념임에 유의.)
- **정공법(권장):** 1200/1357행의 `MPI_ALLTOALL`을 autotune된 P2P/ksection 선택으로
  교체하여, DD ordering과 무관하게 multigrid 카운트 발견도 sparse 경로로.

---

## 8. 진단 재현 절차 (체크리스트)

```bash
# 1) RUNNING인데 로그가 멈췄나? (.err 비었고 mtime 동결)
squeue -u <user>; stat -c '%y' <log>; date

# 2) 프로세스가 spin인가? (R, wchan=-, 99% CPU)
ssh <node> 'ps -L -u <user> -o pid,psr,pcpu,stat,wchan:20,comm | grep ramses'

# 3) 어디서 막혔나? (gdb/perf 없으면 eu-stack)
ssh <node> 'eu-stack -p <pid>'        # → MPI_Alltoall ← build_parent_comms_mg 확인

# 4) OOM 아님 확인
grep "No more free memory" <log>      # → 매치 없어야 함
```
