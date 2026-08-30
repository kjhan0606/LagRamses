구현과 빌드를 완료했습니다.

- 정방향 `make_virtual_fine_dp2` 추가
  - P2P: re/imag를 `%u(:,1:2)`에 패킹해 한 번에 교환
  - ksection: `2*twotondim+2` properties로 한 번에 교환
  - component 8을 사용해 P2P/ksection/auto 선택
  - 새 드리프트 경로에서만 Hilbert ksection 후보 허용
- 역방향 `make_virtual_reverse_dp2` 추가
  - 기존 reverse P2P/ksection 구현과 component 3 재사용
  - 누적 순서 차이에 따른 roundoff 가능성 문서화
- 기본 `.false.`인 `fdm_ghost2`, `fdm_ghost2_rev` 독립 토글 추가
- 지정된 CN 호출부만 변경했으며 다른 ghost 호출부는 그대로입니다.
- 기존 바이너리를 `bin/ramses_final3d.pre_ksecghost`로 백업하고 동일성을 확인한 뒤 빌드했습니다.
- `make HDF5=1 USE_FFTW=1`: `EXIT 0`
- 컴파일러 경고/오류 없음
- 지정된 `setvars.sh` 경로는 이 호스트에 없어 source 경고가 발생했지만, 기존 환경의 `mpiifx`로 컴파일과 링크는 성공했습니다.
- 요청대로 solver 실행이나 A/B 비교는 하지 않았습니다.

변경 사항은 [virtual_boundaries.kjhan.f90](/home/kjhan/BACKUP/lagRamses/patch/cuRamses/virtual_boundaries.kjhan.f90:564), [fdm_step.f90](/home/kjhan/BACKUP/lagRamses/patch/lagRamses/fdm_step.f90:1024), [fdm_commons.f90](/home/kjhan/BACKUP/lagRamses/patch/cuRamses/fdm_commons.f90:13)에 있습니다.

상세 보고서: [CODEX_A_KSECGHOST.md](/home/kjhan/BACKUP/lagRamses/docs/CODEX_A_KSECGHOST.md)