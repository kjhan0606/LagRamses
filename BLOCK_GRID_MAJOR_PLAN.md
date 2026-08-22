# Block grid-major 동적 메모리 개편 계획

- 상태: **CLOSED** (2026-08-22, 사용자 승인 기본 범위)
- 최종 결론: [`BLOCK_GRID_MAJOR_CONCLUSION.md`](BLOCK_GRID_MAJOR_CONCLUSION.md)
- 작성일: 2026-08-11
- 대상 브랜치: `main` (2026-08-18 단일 브랜치 통합에 따라 갱신; 원문은 fdm-dev)
- 구현 위치: CPU/AMR source winner는 `patch/lagRamses/`, CUDA kernel/interface는
  audited VPATH에 따라 `patch/cuRamses/` 사용

## 1. 목적

현재 lagRamses는 `ngridmax`와 `npartmax`를 실행 전에 정하고, 그 크기에
맞추어 AMR grid, cell, particle 배열을 한 번에 할당한다. 실제 사용량이 예상보다
작으면 메모리를 낭비하고, 예상보다 크면 계산 중 `Increase ngridmax` 또는
`Increase npartmax`와 함께 중단된다.

이 개편의 목적은 다음과 같다.

1. `ngridmax`와 `npartmax`를 필수 입력값에서 제외한다.
2. grid와 particle capacity를 실행 중 안전하게 확장한다.
3. 기존 checkpoint를 그대로 읽는다.
4. output grid format은 namelist에서 선택하며, 기본값은 기존(legacy) format으로 한다.
5. 기존 solver의 수치 결과와 성능을 보존한다.

이 문서는 개편 당시의 계획과 gate 설계를 보존하는 역사 기록이다. 최종 승인 범위,
실행 결과와 알려진 제한은 `BLOCK_GRID_MAJOR_CONCLUSION.md`를 정본으로 삼는다.
아래에 남은 후속 항목은 폐문 뒤의 필수 작업을 뜻하지 않는다.
이번 CLOSED 판정은 아래 §10/§13의 원래 production qualification matrix를 모두
충족했다는 뜻이 아니다. 사용자가 축소·확정한 기본 범위에 대한 폐문이며, 원래
조건 중 미충족 항목은 승인 비주장 또는 알려진 제한으로 남는다.

## 2. 현재 구조의 제약

기존 fine-cell index는 대략 다음과 같다.

```text
icell_legacy = ncoarse + (ichild - 1) * ngridmax + igrid
```

즉, 모든 grid의 같은 child가 하나의 큰 평면을 이룬다. 실행 도중
`ngridmax`를 바꾸면 기존 fine-cell index가 모두 바뀐다. `father`, `son`,
`nbor`, particle-to-cell link와 MPI buffer에 저장된 index도 함께 무효화되므로,
단순한 Fortran `reallocate` 또는 C `realloc`만으로는 capacity를 늘릴 수 없다.

2026-08-11의 예비 정적 검색에서는 다음 범위가 확인되었다.

- `ngridmax`를 직접 참조하는 파일: 99개
- `npartmax`를 직접 참조하는 파일: 22개
- `ncoarse`와 `ngridmax`를 조합한 cell-index 산술을 포함하는 파일: 88개

따라서 이 작업은 allocator만 교체하는 국소 수정이 아니다. index 생성·분해를
공통 API로 모은 뒤 call site를 단계적으로 교체해야 한다.

## 3. 제안하는 block grid-major index

`B`를 한 block에 들어가는 grid 수, `C=2^ndim`을 grid당 child-cell 수라고
하자. 1-based `igrid`에 대해 다음 값을 정의한다.

```text
iblock = (igrid - 1) / B
islot  = mod(igrid - 1, B) + 1

icell_block = ncoarse
            + iblock * (C * B)
            + (ichild - 1) * B
            + islot
```

한 block 안에서는 같은 child 위치가 연속적이므로 현재 vector loop의 메모리
접근 특성을 상당 부분 유지한다. 새 block을 뒤에 추가해도 기존 block의 index는
변하지 않는다. 이 불변성이 실행 중 capacity 확장을 가능하게 한다.

초기 성능 시험에서는 `B=64`와 `B=128`을 비교한다. 최종 기본값은 hydro,
gravity, refinement, particle deposition, MPI packing의 실측 결과로 정한다.
block 크기는 checkpoint header에 기록하여 restart 시 해석이 모호하지 않게 한다.

## 4. index API

직접적인 `ncoarse+(ichild-1)*ngridmax+igrid` 계산을 금지하고 다음 역할을 하는
작은 공통 함수를 둔다.

```text
cell_index(igrid, ichild) -> icell
cell_grid(icell)          -> igrid
cell_child(icell)         -> ichild
cell_block(icell)         -> iblock
cell_is_fine(icell)       -> logical
```

성능이 중요한 loop에서는 이 함수들이 compiler에 의해 inline되어야 한다.
가능하면 scalar와 vector 형태를 모두 제공한다. debug build에는 범위 검사와
encode/decode round-trip 검사를 넣고, production build에서는 산술식만 남긴다.

공통 API로 우선 교체할 대상은 다음과 같다.

- AMR tree: `father`, `son`, `nbor`, free/active grid list
- refinement와 derefinement
- load balancing과 MPI packing/unpacking
- hydro, gravity, RT, FDM cell loops
- particle deposition/interpolation 및 particle-to-cell link
- binary/HDF5 checkpoint reader와 writer
- movie, light-cone, clump, sink 및 analysis output

## 5. 동적 capacity 관리

### 5.1 Grid capacity

grid capacity는 `B`의 정수배로만 증가시킨다. cell field의 capacity는 자동으로
`ncoarse + C * grid_capacity`가 된다. 성장량은 다음 조건의 큰 값을 사용한다.

- 최소 한 block
- 현재 capacity의 일정 비율(초기 후보 25--50%)
- 바로 필요한 grid 수에 안전 여유분을 더한 값

growth는 refinement나 load balance가 새 grid를 요구하기 전에 예측하고,
OpenMP parallel region 밖의 명시적인 safe point에서 수행한다. 미완료
nonblocking MPI request가 있는 동안에는 재할당하지 않는다.

### 5.2 Particle capacity

**`nparttot`은 `ngridtot`과 동등한 1급 목표이며, block layout 작업에 종속되지
않는다** (2026-08-19 사용자 지시). 근거는 구조적이다. cell index 식
`ncoarse+(ichild-1)*ngridmax+igrid`에 `npartmax`가 들어가지 않고, 트리 전체
정적 검색에서도 `npartmax`와 `ncoarse`가 함께 나타나는 식은 존재하지 않는다.
즉 particle capacity 확장은 index layout 변경과 **완전히 직교**하며, Phase 1의
index API나 Phase 2의 block layout을 기다릴 이유가 없다.

따라서 구현 순서에서 particle growth(Phase 3)는 grid growth(Phase 4)의 선행
조건이 아니라 **병렬 트랙**으로 취급한다. Phase 1이 끝나는 즉시 착수할 수 있고,
Phase 2와 동시에 진행해도 충돌하지 않는다.

particle 배열은 같은 capacity를 공유하는 하나의 bundle로 관리하며, star
formation, sink/feedback, tracer 생성 또는 load balance 전에 필요한 수를
확인한다. bundle의 일부 배열만 확장된 상태는 허용하지 않는다. 확장 시
`move_alloc` 복사가 순간적으로 구·신 배열을 동시에 점유하므로, admission
control은 정상 상태가 아니라 **transient peak**를 기준으로 판정하거나 bundle을
배열 단위로 staggered copy하여 peak를 최대 단일 배열 크기로 묶는다.

현재 두 capacity는 `read_params.jaehyun.f90:1290-1301`에서 대칭적으로 유도된다.
`ngridmax==0`이면 `ngridtot/ncpu`, `npartmax==0`이면 `nparttot/ncpu`이다. 다만
`ngridtot==0`은 오류로 처리되는 반면 `nparttot==0`은 조용히 `npartmax=0`을
남긴다. 자동 모드에서는 두 경로 모두 "미지정 = 자동"으로 동일하게 해석하고,
어느 쪽도 침묵하는 0을 남기지 않는다.

### 5.3 메모리 backend

우선 두 backend를 검토한다.

1. 표준 Fortran의 새 배열 할당 + 복사 + `move_alloc`
2. `ISO_C_BINDING` 기반 C-managed storage와 `realloc`, Linux에서 선택적으로
   `mremap`

첫 번째 방식은 이식성과 소유권이 명확하지만 확장 때 전체 복사가 필요하다.
두 번째 방식은 in-place 확장 가능성이 있으나 주소가 바뀔 수 있으므로 모든
Fortran pointer view를 다시 연결해야 한다. `mremap`은 Linux 전용 최적화이며
정확성을 위한 필수 조건으로 사용하지 않는다. 어떤 backend에서도 raw address를
장기간 보관하거나 MPI/OpenMP 작업을 가로질러 유지하지 않는다.

## 6. 자동 메모리 예산

"남은 메모리를 모두 사용"하는 방식은 MPI rank 간 경쟁과 OS OOM 위험이 있어
사용하지 않는다. node 단위 예산을 먼저 정하고 각 rank가 이를 공유한다.

우선순위는 다음과 같다.

1. Slurm/cgroup에 설정된 node 또는 job memory limit
2. scheduler가 제공하는 명시적 memory 정보
3. `/proc/meminfo`를 이용한 보수적 fallback

예산에서는 MPI buffer, OpenMP temporary, HDF5, runtime, page cache를 위한
reserve를 제외한다. 같은 node의 rank는 `MPI_Comm_split_type`으로 node-local
communicator를 만들고, 합산 capacity가 node budget을 넘지 않도록 합의한다.
어느 한 rank라도 필요한 growth를 확보하지 못하면 전체 rank가 일관된 진단을
출력하고 checkpoint 후 종료해야 한다. 개별 rank만 allocator failure로 죽는
상태는 허용하지 않는다.

기존 `ngridmax`와 `npartmax` namelist key는 당분간 호환성을 위해 받아들이되,
상한이 아닌 초기 capacity hint로 해석하는 방안을 사용한다. key가 없거나 0이면
자동 모드가 기본값이다. 의미 변경은 release note와 시작 로그에 명시한다.

## 7. Checkpoint format

### 7.1 Output 선택

다음 namelist key를 추가하는 설계를 사용한다.

```fortran
amr_output_layout = 'legacy'            ! default
! or
amr_output_layout = 'block_grid_major'
```

- `legacy`: 저장 직전에 index와 cell field를 기존 layout으로 변환한다.
- `block_grid_major`: 새 index layout을 기록한다.

기본 output은 `legacy`이다 (2026-08-18 사용자 결정). 기존 분석 도구
(utils/f90 reader, HOP 파이프라인 등)와 구버전 lagRamses가 모든 checkpoint를
그대로 읽어야 하고, format 이행은 도구 정비가 끝난 뒤 별도 결정으로 한다.
새 layout의 이점(변환 비용 없는 저장)이 필요한 대규모 실행에서만
`block_grid_major`를 명시적으로 선택한다.

### 7.2 Input 자동 판별

input에는 별도의 layout namelist를 두지 않고 파일에서 자동 판별한다.

Binary AMR 파일의 설계는 다음과 같다.

- 기존 파일: 첫 integer record가 양의 `ncpu`이면 legacy로 판정
- 새 파일: 예약된 음의 format tag로 시작
- tag 다음 record: format version, `B`, dimension, 필요한 index-width 정보
- 알 수 없는 tag/version: allocation 전에 명시적으로 중단

실제 magic tag 값은 기존 checkpoint corpus를 조사한 뒤 충돌하지 않는 값으로
확정한다. 값의 범위만 보고 cell index layout을 추측하지 않는다.

HDF5에서는 AMR group attribute에 다음 정보를 둔다.

```text
amr_layout = "block_grid_major"
amr_layout_version = 1
amr_block_size = B
```

attribute가 없으면 legacy HDF5로 판정한다. 한 checkpoint를 구성하는 모든 rank
파일은 layout, version, `B`가 같아야 한다. AMR header의 판정을 hydro,
Poisson, particle, RT, FDM companion file에도 일관되게 적용한다.

### 7.3 변환 경계

solver 내부에서는 block grid-major만 사용한다. 변환은 I/O 경계에 제한한다.

- legacy input -> block grid-major internal
- block input -> block grid-major internal
- block internal -> legacy output (기본)
- block internal -> block output (선택 시)

변환 대상에는 field 배열뿐 아니라 `father`, `son`, `nbor`, particle cell link 등
cell index를 값으로 저장하는 모든 record가 포함된다. legacy output 변환은
`O(Ncell + Nlink)` 시간과 packing buffer를 요구한다. 실제 비용은 checkpoint
시간 대비 비율로 benchmark하여 문서화한다.

## 8. 구현 위치와 변경 규칙

모든 수정본은 반드시 다음 디렉터리에 둔다.

```text
/home/kjhan/BACKUP/lagRamses/patch/lagRamses/
```

`amr/`, `pm/`, `hydro/`, `poisson/`, `rt/` 등의 원본 파일을 직접 수정하지
않는다. 필요한 파일은 먼저 `patch/lagRamses/`로 복사한 뒤 그 사본만 수정한다.
build system에서 patch 파일이 원본보다 우선하도록 한다. 이 문서 작성 시점에는
어떠한 구현 파일도 변경하지 않는다.

## 9. 단계별 구현 순서

### Phase 0: 기준선 고정

- 현재 branch와 compiler/MPI 조합 기록
- small hydro, gravity, AMR, particle, RT/FDM 기준 결과 보존
- 기존 binary/HDF5 checkpoint corpus 수집
- memory high-water mark와 주요 kernel wall time 기록

### Phase 1: index abstraction

- encode/decode API 추가
- 기존 layout을 그대로 사용한 상태에서 직접 index 산술을 API로 교체
- bitwise 또는 허용 오차 내 기준선 일치 확인

이 단계에서는 메모리 layout을 바꾸지 않는다. API 전환 자체와 layout 변경의
오류를 분리하기 위한 단계다.

### Phase 5 재설계: 마커 추가가 아니라 변환 경계 (2026-08-19 운전자 판단)

파손 #1(재시작)과 #3(defrag/dbl2sng)의 해법을 "헤더에 amr_block_size 레코드를
추가"로 잡았으나, 헤더가 순차 unformatted 레코드 열이라 레코드 추가는 그 자체로
포맷 변경이고 모든 외부 reader(amr2map, amr2cell, io_ramses, pario)를 깬다.

사용자가 이미 **출력 기본값을 legacy로** 확정했으므로(§7.1) 더 나은 해법이
나온다. 계획 §7.3의 변환 경계를 실제로 구현하면 된다.

- 쓰기: solver 내부는 block, 저장 직전 father/nbor/son 등 **셀 인덱스를 값으로
  담는 모든 레코드**를 legacy로 변환해 기록한다. 파일은 기존 포맷 그대로다.
- 읽기: legacy 파일을 읽어 block으로 변환해 메모리에 올린다.

이러면 #1은 사라진다. 파일이 항상 legacy 의미이므로 remap의 legacy 디코드
가정이 맞고, `ngridmax.ne.ngridmax2` 가드도 원래 의도(용량 변경 감지)대로
동작한다. #3도 사라진다. defrag/dbl2sng와 모든 외부 reader가 손대지 않고
계속 동작한다. 마커도, 포맷 버전도, 하위호환 분기도 필요 없다.

비용은 저장·복원 시 O(Ncell+Nlink)의 인덱스 변환 한 번이며, 이는 §7.3이 이미
예산에 넣은 항목이다. block 포맷 출력(`amr_output_layout='block_grid_major'`)은
변환 비용이 실측으로 문제가 될 때 켜는 선택지로 남긴다.

따라서 Phase 5의 첫 작업은 마커가 아니라 **변환 경계 구현**이며, 이것이
#1과 #3을 동시에 닫는 최소 경로다. #2(CUDA)는 파일과 무관하므로 별개로
nGR-GPU 트랙과 함께 처리한다.

### Phase 2 이후 남은 파손 목록 (2026-08-19 Fable 감사 + 실측, 사전등록)

B=64 활성화(main 55a38bd)는 CPU·fresh-start·바이너리 출력 구성에서 물리량
비트 동일로 검증됐다. 그 밖의 경로는 아직 legacy contiguous를 가정한다.
위험순 목록과 각각의 처리 방침:

1. **바이너리 체크포인트 재시작.** 헤더(output_amr.kjhan.f90:373)에 layout
   표시가 없고, init_amr의 remap은 파일이 legacy stride라고 가정해 디코드한다.
   양방향 파손이며, 특히 ngridmax가 이미 B의 배수면(예: 1000000=64x15625)
   `ngridmax.ne.ngridmax2` 가드가 remap을 건너뛰어 legacy 파일을 block으로
   조용히 오독한다. **수정 방향**: 헤더에 amr_block_size를 기록하고, remap을
   파일의 (B2, ngridmax2) 쌍으로 디코드하도록 바꾼다. 이는 외부 reader와
   defrag/dbl2sng에도 포맷 판별 수단을 준다. Phase 5의 핵심 작업.
2. **CUDA 커널 3종 — 코드 교정 완료(2026-08-22, `23e2942`)**: poisson 8곳,
   particle 3곳(그중 1곳은 live-cell upload prefix), scalar 2곳의 legacy stride를
   공용 device helper와 명시적 B/C ABI로 교체했다. Poisson의 `flag2/ngridmax`
   packing 2곳은 cell stride가 아니므로 whitelist로 유지한다. active VPATH source
   winner와 cuRamses shadow ABI를 함께 맞췄고, host/device unit, static census,
   전처리 reachability, full `USE_CUDA=1 HDF5=1 USE_FFTW=1` build가 PASS했다.
   Fable의 코드·빌드 감리는 승인됐다. Phase 6 job 322420에서는 CUDA 기능 경로
   실행 증거를 확보했지만 batch/strict numerical-parity gate는 FAIL했다.
3. **체크포인트 재작성 유틸리티**: utils/f90/defrag.f90, dbl2sng.f90이 각각
   46곳에서 raw stride로 father/nbor를 디코드·재인코딩한다. B=64 파일을 조용히
   손상시킨다. 1번의 헤더 마커가 들어가면 이들은 "미지원 layout이면 중단"으로
   먼저 막고, 이후 매크로와 동일한 식으로 옮긴다.
4. 미컴파일 디렉토리(rt/, pario/, mhd/, rhd/, aton/ 및 VPATH 그림자 사본)는
   전부 legacy이며 센서스 스캔 대상 밖이다. 해당 SOLVER/타깃을 되살릴 때
   반드시 함께 변환한다.

**게이트 원칙 확정**: layout을 바꾸는 청크는 amr_* 파일이 설계상 달라지므로
`tests/phase2_layout_gate.sh`(hydro/part/grav 비트 동일 + 에너지·타임스텝
일치)로 판정한다. 용량이 다르면 부하분산 이력이 달라져 물리량이 반올림
수준으로 흔들리므로, 비교는 **반드시 같은 용량**에서 한다.

### Phase 2 chunk 2 설계 (2026-08-19 운전자 판단, 사전등록)

인벤토리 결과 37곳이 `ncoarse+twotondim*ngridmax`를 전체 셀 배열 크기로
쓴다. block 크기 B가 ngridmax를 나누면 nblocks=ngridmax/B 이고
nblocks*C*B = C*ngridmax 이므로 배열 크기가 불변이다. 최대 셀 인덱스도
ncoarse+nblocks*C*B = ncoarse+C*ngridmax 로 기존 배열 경계에 정확히 맞는다.
따라서 chunk 2는 (a) setup에서 ngridmax를 amr_block_size의 배수로 올림,
(b) amr_block_size를 B(기본 64)로 설정, 이 둘뿐이며 37곳 할당·루프 변경이
불필요하다. 셀 값은 동일하고 메모리 슬롯 순서만 바뀌므로, 전체 셀 배열을
메모리 순서로 순회하는 감산/영기화 루프의 FP 합산 순서가 달라져 결과가
반올림 수준으로 흔들린다. 게이트는 bitwise가 아니라 **보존량 tolerance**:
완주 + mcons/econs/총질량/총운동량의 legacy(B=ngridmax) 대비 상대차 < 1e-8,
그리고 halo 구조가 tolerance 내 일치. 결과 diff는 Fable 검수를 거친다.

### Phase 2 진입 판정 (2026-08-19 운전자 판단, 사전등록)

Phase 2의 첫 청크는 layout을 실제로 바꾸지 않고, index 식만 block
일반식으로 재작성한 뒤 **block 크기 B를 ngridmax로 설정**한다. 이때 block
식은 legacy 식과 항등이다: iblock=(igrid-1)/B, igrid<=ngridmax<=B 이므로
iblock=0, islot=igrid-1 이 되어

  icell = ncoarse + iblock*(C*B) + (ichild-1)*B + islot + 1
        = ncoarse + (ichild-1)*ngridmax + igrid

로 legacy와 완전히 일치한다. 역변환(IGRID_OF/ICHILD_OF)도 같은 대입으로
항등임을 확인했다. 따라서 이 청크의 게이트는 **Phase 1과 동일한 -O3 bitwise**
이며, block 인프라(B 파라미터, 일반식, 배열 정렬 규칙)를 legacy 결과를 한
비트도 바꾸지 않고 트리에 들일 수 있다. 통과하면 layout 코드는 검증된 채로
자리를 잡고, 다음 청크에서 비로소 B를 64/128로 줄여 실제 block 배치를
켠다. B를 줄이는 청크부터는 셀의 메모리 순회 순서가 바뀌어 FP 합산 순서가
달라질 수 있으므로 게이트를 **보존량 tolerance + -O0 재현**으로 전환하고,
그 설계는 Fable 검수를 거친다.

### Phase 2: 고정 capacity block layout

- block grid-major 배열 배치 도입
- `B=64`, `B=128` 정확성 시험
- hydro/gravity/refinement/MPI kernel 성능 비교
- 아직 실행 중 growth는 사용하지 않음

### Phase 3 병렬 착수 (2026-08-19 사용자 지시로 Phase 2와 동시 시작)

Phase 2와 Phase 3는 npartmax가 cell index 식에 없다는 점에서 직교하므로,
각각 detached worktree(`lagRamses-p2`, `lagRamses-p3`)에서 병렬로 진행하고
완료 단계마다 main에 순차 머지한다. 주로 건드리는 파일이 겹치지 않는다
(Phase 2: amr_index.h, amr_commons, 배열 할당; Phase 3: read_params, pm
배열, star/sink/feedback 경로). Phase 3의 첫 청크는 read_params의 nparttot
대칭화이며, 그 자체로 동작을 바꾸지 않는 정리이므로 게이트는 -O3 bitwise.

### Phase 3: particle dynamic growth (Phase 2와 병렬 가능)

`npartmax`는 cell index 식에 없으므로 이 단계는 block layout과 독립이다.
Phase 1 완료 직후 착수하며 Phase 2를 기다리지 않는다.

- particle bundle allocator 도입
- creation, deletion, load balance, restart 시험
- `npartmax`와 `nparttot`을 optional initial hint로 전환
- 확장 시 transient peak 메모리를 admission control에 반영

### Phase 4: grid dynamic growth

- grid와 cell bundle allocator 도입
- refinement/load balance safe point에서 collective growth
- `ngridmax`를 optional initial hint로 전환
- allocator failure 시 collective checkpoint/abort 경로 추가

### Phase 5: I/O 호환성

- legacy binary/HDF5 auto detection
- 새 block format/version header
- 두 output layout 구현
- legacy <-> block restart matrix 검증

2026-08-22 legacy AMR canonical comparator를 완료했다.
`tests/compare_amr_canonical.py`는 각 파일의 `ngridmax`로 legacy cell index를
해독하고 전역 dyadic geometry로 father/nbor/son/flag1을 비교한다. capacity,
free/list 상태, local grid 번호와 rank decomposition은 제외한다. Job 321554와
321942의 output 3개씩(총 6 pair)이 모두 PASS했으며 JSON 보고서는 각 QA job
디렉터리에 보존한다. 이 도구는 topology equivalence용이며 HDF5/restart matrix를
대체하지 않는다.

### Phase 6 범위 확정 (2026-08-22 사용자 지시): CUDA/nGR 소형 결정적 게이트

VoidSim과 대형 production 실행은 이 프로젝트의 선행조건·완료조건·테스트
자산에서 제외한다. Phase 6의 코드 교정은 commit `23e2942`에서 완료됐다.
공용 helper는 Fortran `amr_index::icell_of`와 같은 block grid-major 식을 쓰며,
device측 13곳(poisson 8, particle 3, scalar 2), B/C ABI, active rho particle
deposit dispatch를 함께 교정했다. `tests/cuda_block_index_check.sh`와 full CUDA
build가 PASS했고 Fable은 코드 commit/main 반영을 승인했다. 작은 고정 zoom IC의
runtime gate는 job 322420에서 실행했으며, child run 완주와 기능 marker는 확인했지만
전체 batch와 strict numerical-parity gate는 FAIL했다.

최소 실행은 1 node, 2 MPI ranks, OMP=1, A10 1개로 한다. CPU와 GPU는 같은
USE_CUDA binary·IC·물리 입력을 쓰고 accelerator flag만 달라야 한다. IC는 QA
아래 project-local snapshot과 SHA manifest로 고정한다. GPU run은 flag만 켜는
것으로 부족하며 Poisson, scalar, particle CUDA 경로가 실제 실행됐다는 양의
counter/marker를 각각 남겨야 한다.

필수 양성 증거는 `[CUDA_MG]`의 gs/residual/restrict/interp가 모두 0보다 큼,
`[CUDA_NGR]`의 uploads/scalar_sweeps가 모두 0보다 큼, 그리고
`[CUDA_PM_GATHER]`·`[CUDA_PM_DEPOSIT]`의 B=64, C=8, count>0이다. CPU control은
`gpu_particle=.false.`를 포함해 모든 accelerator flag를 끄고 기존 fallback을
검증한다. warning/fallback/CUDA error/OOM/NaN/fatal은 모두 0이어야 한다.

nGR GPU는 색 내부 계산 순서가 달라질 수 있으므로 CPU↔GPU bitwise를 요구하지
않는다. particle ID/count와 topology는 exact, post-step 위치는
`max|dx| <= 2e-6`, 상대 속도 오차는 `<= 2e-3`을 1차 사전등록 한계로 둔다.
Poisson residual 정상 수렴, CUDA error/OOM/fallback/NaN/fatal 0도 필수다.
scalar/force L2 한계는 첫 characterization 뒤 baseline의 2배 이내로 고정한다.

반복 growth, MPI rank 수가 다른 restart, MPI+OpenMP, 작은 2-node smoke는 각각
독립된 짧은 gate로 수행한다. 장시간·production-size 실행으로 이 검증을 대신하지
않는다. HDF5/restart gate는 사용자 승인 뒤 job 322426/322427에서 실행했으며,
restore 경로는 관측했지만 one-step exact continuity는 FAIL했다.

## 10. 필수 검증 행렬

| Input | Internal | Output | Reader | 기대 결과 |
|---|---|---|---|---|
| legacy binary | block | block | new | 정상 restart |
| legacy binary | block | legacy | old/new | 모두 정상 restart |
| block binary | block | block | new | 정상 restart |
| block binary | block | legacy | old/new | 모두 정상 restart |
| legacy HDF5 | block | block HDF5 | new | 정상 restart |
| block HDF5 | block | legacy HDF5 | old/new | 지원 범위 확인 |

각 조합에서 다음을 검사한다.

- grid/cell/particle 개수와 level별 분포
- AMR parent-child-neighbor topology
- 총 질량, 운동량, 에너지 및 particle ID
- hydro/gravity/RT/FDM field checksum 또는 허용 오차
- 같은 MPI rank와 다른 MPI rank 수에서의 restart
- 한 번 및 여러 번 capacity가 증가한 뒤의 checkpoint
- 비정상/부분 checkpoint에 대한 명시적 오류 처리

## 11. 성능 평가

기준은 기존 fixed-capacity legacy layout이다. 다음을 별도로 측정한다.

- 주요 kernel wall time과 vectorization report
- MPI packing/unpacking 및 collective 시간
- capacity growth 횟수, 복사량, pause 시간
- binary/HDF5 checkpoint 시간
- block -> legacy 변환 시간과 peak scratch memory
- rank별 및 node별 resident/high-water memory

고정 capacity block layout 자체가 주요 kernel에서 유의미한 회귀를 보이면 동적
growth 구현으로 넘어가기 전에 block 크기나 loop ordering을 다시 설계한다.
평균 성능뿐 아니라 load-imbalanced rank의 tail latency를 함께 본다.

## 12. 주요 위험과 대응

1. **직접 index 산술 누락**: 정적 검색, debug bounds check, topology validator로
   검출한다.
2. **재할당 후 stale pointer**: 영구 pointer alias를 금지하고 growth 뒤 view를
   중앙에서 재연결한다.
3. **MPI 진행 중 재할당**: safe point와 request-completion assertion을 둔다.
4. **OpenMP data race**: master/serial 구간에서만 growth하고 barrier 뒤 공개한다.
5. **legacy output의 과도한 scratch memory**: level/rank 단위 streaming packing을
   우선 설계한다.
6. **외부 분석 도구 호환성**: legacy output을 유지하고 새 format specification을
   문서화한다.
7. **OOM killer**: cgroup-aware reserve와 node-local collective admission을 둔다.
8. **format 오판**: 명시적인 magic/version만 사용하고 heuristic index 판정을
   금지한다.

## 13. 원래 production qualification 조건 (이번 폐문 범위에서는 superseded)

다음 조건은 원래 production 기본값 채택을 위해 작성한 행렬이다. 2026-08-22의
사용자 승인 기본 범위 폐문이 이 조건을 모두 충족했다는 뜻은 아니다.

- `ngridtot`/`ngridmax`, `nparttot`/`npartmax` 넷 다 없이 시작하고, grid와
  particle 각각 최소 두 번의 growth를 거쳐 완주
- 기존 checkpoint 자동 인식 및 restart 성공
- 기본값인 legacy output 파일을 기존 reader가 읽음
- `block_grid_major` output 선택 시 same-rank/different-rank restart 성공
- 기준 계산과 보존량 및 과학 결과가 합의된 허용 오차 안에서 일치
- MPI, OpenMP, MPI+OpenMP 시험 통과
- memory limit 접근 시 OOM kill 대신 예측 가능한 진단과 종료
- 성능 및 I/O 변환 비용이 측정되어 문서화됨
- CPU allocator/layout 변경은 `patch/lagRamses/`의 audited source winner에,
  CUDA 전용 kernel/interface 변경은 `patch/cuRamses/`의 audited source winner에
  존재하며, `lagRamses-DE`와 다른 프로젝트는 수정하지 않음
- 작은 고정 IC의 CUDA/nGR CPU↔GPU gate에서 Poisson·scalar·particle GPU 실제
  실행 증거, 정상 수렴, exact topology/ID 및 사전등록 수치 허용오차를 모두 통과
- VoidSim 또는 다른 프로젝트의 production 자산에 의존하지 않음

## 14. 현재 결정 사항과 미결 사항

확정된 사항:

- 내부 layout은 block grid-major 방향으로 개발한다.
- output layout은 namelist로 선택한다.
- output 기본값은 `legacy`(기존 format)이다. (2026-08-18 변경)
- input layout은 자동 판별한다.
- 기존 checkpoint 호환성을 유지한다.
- CPU/AMR 구현은 `patch/lagRamses/`, CUDA 전용 구현은 실제 build VPATH의
  `patch/cuRamses/`에 둔다. source winner와 shadow ABI는 static census로 확인한다.
- HDF5/restart는 사용자 승인 뒤 job 322426/322427에서 실행했다. restore 경로는
  확인했지만 one-step exact continuity는 FAIL로 남겼다.

구현 전에 결정하거나 측정할 사항:

- 기본 block 크기: `B=64` 또는 `B=128`
- particle bundle 확장 backend: `move_alloc` 전체 복사(peak 2x)를 감수할지,
  배열 단위 staggered copy로 peak를 낮출지
- 성장 계수와 최소 headroom
- C-managed backend를 기본으로 할지, Fortran backend를 fallback으로 둘지
- binary magic tag의 실제 값과 version header 상세 규격
- HDF5 legacy writer가 지원해야 하는 정확한 외부 reader 범위
- 허용 가능한 kernel 성능 회귀와 legacy 변환 비용의 정량 기준
