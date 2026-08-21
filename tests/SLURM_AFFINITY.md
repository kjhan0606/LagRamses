# Slurm MPI rank affinity 운용 메모

이 문서는 `syn` 계산 노드에서 한 노드/32 MPI rank CPU 작업을 다시 제출할 때
재사용하는 기록이다. Phase 4 전용 스크립트는 `phase4_gridgrow.sbatch`이며,
공용 도구는 다음 세 파일이다.

- `expand_cpuset.sh`: Linux `Cpus_allowed_list`를 오름차순 CPU CSV로 확장한다.
- `slurm_rank_bind.sh`: `SLURM_LOCALID`번째 CPU에 현재 rank를 `taskset`으로
  고정한 뒤 실제 프로그램을 `exec`한다.
- `validate_affinity_map.sh`: MPI_Init 이후의 rank/localid/CPU map을 예정된
  mapping과 대조한다.

## 이 방식이 필요한 이유

2026-08-20 job 321015에서는 Intel MPI pinning과 Slurm binding을 함께 사용해
32 rank가 모두 CPU 0으로 축소되었다. 정상 실행보다 약 37배 느려졌고 2시간 뒤
timeout되었다. 이 결과는 구현 성능이나 block grid-major의 단점으로 평가하지
않는다.

2026-08-21 job 321121에서는 `I_MPI_PIN=0`으로 바꿨지만 모든 rank가 동일한
allocation-wide mask `1-4,16-19,32-55`를 물려받아 strict probe에서 중단됐다.
이 사이트의 Slurm은 `TaskPlugin=task/cgroup`이고 `task/affinity`가 없다.
따라서 cgroup은 job의 CPU 경계만 격리하며 `srun --cpu-bind=cores`가 rank별
mask를 만든다고 기대할 수 없다. 공유 cpuset 자체는 기능상 가능하지만, 재현성
있는 검증 런에는 명시적인 one-rank/one-CPU mapping을 사용한다.

## 재사용 규칙

1. Intel MPI의 pinning은 끈다.

       export I_MPI_PIN=0
       unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_LIST I_MPI_PIN_RESPECT_CPUSET
       export OMP_PROC_BIND=false
       unset KMP_AFFINITY KMP_HW_SUBSET GOMP_CPU_AFFINITY

2. batch process의 `/proc/self/status`에서 `Cpus_allowed_list`를 읽고 다음처럼
   정확한 task 수만큼 확장한다. malformed/reverse/duplicate/overlap/개수 오류는
   모두 실패해야 한다.

       allocation_cpu_mask=...
       PHASE4_CPU_LIST=$(tests/expand_cpuset.sh "$allocation_cpu_mask" 32)
       export PHASE4_CPU_LIST

3. `srun`은 launcher 역할만 하며 per-rank binding은 helper 한 곳에서 맡는다.
   probe와 본 계산에 반드시 같은 명령 형태를 사용한다.

       srun --kill-on-bad-exit=1 --mpi=pmi2 -n 32 -c 1 --cpu-bind=none \
         tests/slurm_rank_bind.sh PROGRAM ARG...

4. MPI probe는 `MPI_Init` 이후 `MPI rank`, `SLURM_LOCALID`, 실제
   `Cpus_allowed_list`를 기록해야 한다. 본 계산은 다음을 모두 확인한 뒤에만
   시작한다.

   - 파일과 행이 각각 32개
   - MPI rank 0..31과 localid 0..31이 각각 한 번씩 존재
   - 실제 mask가 단일 CPU이고 32개 모두 고유
   - `actual_cpu == PHASE4_CPU_LIST[SLURM_LOCALID]`

5. provenance에 원래 allocation mask, 확장된 CPU 목록,
   `binding=taskset_by_SLURM_LOCALID`를 남긴다.

노드나 Slurm 설정이 바뀌어 `task/affinity`가 활성화되면 이 우회가 계속 필요한지
다시 감사한다. 그 전에는 Intel MPI pinning, Slurm `--cpu-bind`, `taskset` 중
둘 이상을 동시에 rank binding authority로 사용하지 않는다.

## 제출 전 확인

다음을 모두 통과시키고 변경 파일을 같은 commit에 포함한다.

    bash -n tests/phase4_gridgrow.sbatch tests/expand_cpuset.sh \
      tests/slurm_rank_bind.sh tests/validate_affinity_map.sh
    git diff --check
    mpiifx -syntax-only tests/mpi_affinity_probe.f90

파서는 최소한 아래 사례를 검사한다.

- `1-4,16-19,32-55`, expected 32: PASS
- malformed, 역순 range, 중복/겹침, 31개, 33개: FAIL

제출 후에는 affinity probe PASS와 초기 coarse step 시간이 과거 수준인 수십 초로
회복됐는지 먼저 본다. 약 930초라면 affinity 밖의 원인을 다시 조사하며 결과를
Phase 4 구현 평가에 사용하지 않는다.

## Slurm 작업 안전

제출에서 반환된 exact job ID만 추적한다. 취소가 필요하면 먼저 다음 명령으로
job name과 작업 디렉터리가 모두 이 프로젝트인지 확증한 뒤 그 ID 하나만 취소한다.

    sacct -j JOB_ID -o JobName,WorkDir
    scancel JOB_ID

`squeue -u kjhan`은 여러 프로젝트가 공유하는 계정 전체 목록이다. 모르는 작업은
남의 작업으로 취급하며 `scancel -u kjhan` 같은 계정 단위 취소는 금지한다.
