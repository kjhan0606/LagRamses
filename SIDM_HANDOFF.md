# SIDM 프로젝트 업무인수인계서

작성 시각: 2026-08-24 KST

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

- Monte Carlo event ratio는 `a=0.3`에서 `1.07 +/- 0.06`이고 `a=0.6`에서
  `1.006 +/- 0.044`이다.
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
