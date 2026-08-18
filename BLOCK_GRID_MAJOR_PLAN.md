# Block grid-major 동적 메모리 개편 계획

- 상태: 설계 단계, 구현 미승인
- 작성일: 2026-08-11
- 대상 브랜치: `main` (2026-08-18 단일 브랜치 통합에 따라 갱신; 원문은 fdm-dev)
- 구현 위치: `patch/lagRamses/`만 사용

## 1. 목적

현재 lagRamses는 `ngridmax`와 `npartmax`를 실행 전에 정하고, 그 크기에
맞추어 AMR grid, cell, particle 배열을 한 번에 할당한다. 실제 사용량이 예상보다
작으면 메모리를 낭비하고, 예상보다 크면 계산 중 `Increase ngridmax` 또는
`Increase npartmax`와 함께 중단된다.

이 개편의 목적은 다음과 같다.

1. `ngridmax`와 `npartmax`를 필수 입력값에서 제외한다.
2. grid와 particle capacity를 실행 중 안전하게 확장한다.
3. 기존 checkpoint를 그대로 읽는다.
4. output grid format은 namelist에서 선택하며, 기본값은 block grid-major로 한다.
5. 기존 solver의 수치 결과와 성능을 보존한다.

이 문서는 구현 계획이다. 승인 전에는 lagRamses 소스코드를 변경하지 않는다.

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

particle index는 `npartmax`가 cell index 식에 포함되는 구조가 아니므로 grid보다
독립적으로 전환할 수 있다. particle 배열은 같은 capacity를 공유하는 하나의
bundle로 관리하며, star formation, sink/feedback, tracer 생성 또는 load
balance 전에 필요한 수를 확인한다. bundle의 일부 배열만 확장된 상태는 허용하지
않는다.

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
amr_output_layout = 'block_grid_major'  ! default
! or
amr_output_layout = 'legacy'
```

- `block_grid_major`: 새 index layout을 기록한다.
- `legacy`: 저장 직전에 index와 cell field를 기존 layout으로 변환한다.

기본 output은 `block_grid_major`이다. 기존 외부 도구나 구버전 lagRamses로
읽어야 할 때만 `legacy`를 선택한다.

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
- block internal -> legacy output (선택 시)
- block internal -> block output (기본)

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

### Phase 2: 고정 capacity block layout

- block grid-major 배열 배치 도입
- `B=64`, `B=128` 정확성 시험
- hydro/gravity/refinement/MPI kernel 성능 비교
- 아직 실행 중 growth는 사용하지 않음

### Phase 3: particle dynamic growth

- particle bundle allocator 도입
- creation, deletion, load balance, restart 시험
- `npartmax`를 optional initial hint로 전환

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

### Phase 6: 대규모 성능 및 안정성 검증

- MPI rank 수를 바꾼 restart
- MPI+OpenMP 혼합 실행
- 반복 growth와 장시간 실행
- node memory limit 근접 시험
- production-size load imbalance 시험

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

## 13. 완료 조건

다음 조건을 모두 만족해야 production 기본값으로 채택한다.

- `ngridmax`, `npartmax` 없이 시작하고 최소 두 번의 growth를 거쳐 완주
- 기존 checkpoint 자동 인식 및 restart 성공
- `amr_output_layout='legacy'` 파일을 기존 reader가 읽음
- 기본 block output의 same-rank/different-rank restart 성공
- 기준 계산과 보존량 및 과학 결과가 합의된 허용 오차 안에서 일치
- MPI, OpenMP, MPI+OpenMP 시험 통과
- memory limit 접근 시 OOM kill 대신 예측 가능한 진단과 종료
- 성능 및 I/O 변환 비용이 측정되어 문서화됨
- 변경 파일이 모두 `patch/lagRamses/`에만 존재

## 14. 현재 결정 사항과 미결 사항

확정된 사항:

- 내부 layout은 block grid-major 방향으로 개발한다.
- output layout은 namelist로 선택한다.
- output 기본값은 `block_grid_major`이다.
- input layout은 자동 판별한다.
- 기존 checkpoint 호환성을 유지한다.
- 구현 파일은 모두 `patch/lagRamses/`에 둔다.

구현 전에 결정하거나 측정할 사항:

- 기본 block 크기: `B=64` 또는 `B=128`
- 성장 계수와 최소 headroom
- C-managed backend를 기본으로 할지, Fortran backend를 fallback으로 둘지
- binary magic tag의 실제 값과 version header 상세 규격
- HDF5 legacy writer가 지원해야 하는 정확한 외부 reader 범위
- 허용 가능한 kernel 성능 회귀와 legacy 변환 비용의 정량 기준

