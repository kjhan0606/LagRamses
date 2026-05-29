# cuRAMSES Source File Classification

All source files reside in `../patch/cuRamses/`.
Build: `cd ../bin && make -f Makefile.cuRamses HDF5=1 USE_CUDA=1`

## Summary

| Category | Count |
|----------|-------|
| New (cuRAMSES 신규 작성) | 12 |
| Modified (기존 파일 수정) | 31 |
| Original (원본 그대로) | 57 |
| **Total** | **100** |

---

## 1. NEW — cuRAMSES 신규 작성 파일 (12)

기존 RAMSES/Horizon5에 대응하는 원본이 없는 완전 신규 파일.

### K-Section Domain Decomposition
| File | Description |
|------|-------------|
| `ksection.f90` | K-section tree 구축 + hierarchical MPI exchange |

### Morton Key Octree
| File | Description |
|------|-------------|
| `morton_keys.f90` | 64-bit Morton key 계산 (21 bits/dim) |
| `morton_hash.f90` | Per-level open-addressing hash table + neighbor lookup |
| `morton_init.f90` | Hash table 초기화 + rebuild |

### CUDA GPU Acceleration
| File | Description |
|------|-------------|
| `cuda_commons.f90` | CUDA Fortran commons module |
| `cuda_stream_pool.cu` | CUDA stream pool (C/CUDA) |
| `cuda_stream_pool.h` | Stream pool header (N_STREAMS 정의) |
| `hydro_cuda_interface.f90` | Fortran-CUDA interface (ISO_C_BINDING) |
| `hydro_cuda_kernels.cu` | 5-kernel hydro pipeline + scatter-reduce kernel |

### HDF5 Parallel I/O
| File | Description |
|------|-------------|
| `ramses_hdf5_io.f90` | HDF5 I/O module (read/write helpers) |
| `backup_hdf5.f90` | HDF5 checkpoint 출력 |
| `restore_hdf5.f90` | HDF5 checkpoint 복원 (variable-ncpu 지원) |

---

## 2. MODIFIED — 기존 파일 수정 (31)

원본 RAMSES/Horizon5 파일을 override하여 수정한 파일.

### patch/cuda → base RAMSES override (17)

| File | Original | Modification |
|------|----------|-------------|
| `amr_parameters.jaehyun.f90` | amr/ | memory_balance, mem_weight 파라미터 추가 |
| `amr_commons.kjhan.f90` | amr/ | ksec_cpumin/cpumax, grid_level 배열 추가 |
| `bisection.f90` | amr/ | 64-bit histogram, memory-weighted cost, nc_in 파라미터 |
| `read_params.jaehyun.f90` | amr/ | ksection/memory_balance namelist 읽기 |
| `init_amr.f90` | amr/ | Morton hash build, hilbert_key 최소 할당, HDF5 init |
| `adaptive_loop.jaehyun.f90` | amr/ | writemem_minmax, bulk virtual exchange 호출 |
| `update_time.f90` | amr/ | timer_m module, cmpmem 조건부 카운팅 |
| `output_amr.kjhan.f90` | amr/ | Morton-computed nbor 출력, HDF5 output |
| `virtual_boundaries.kjhan.f90` | amr/ | ksection ghost zone exchange (dp/int/reverse/bulk) |
| `load_balance.kjhan.f90` | amr/ | ksection load balance, numbp sync, profiling |
| `cooling_fine.kjhan.f90` | hydro/ | GPU dispatch (HYDRO_CUDA) |
| `force_fine.kjhan.f90` | poisson/ | GPU dispatch + hybrid CPU/GPU |
| `multigrid_fine_commons.f90` | poisson/ | ksection MG communication (4 ksec functions) |
| `courant_fine.kjhan.f90` | hydro/ | GPU dispatch |
| `godunov_fine.kjhan.f90` | hydro/ | Hybrid CPU/GPU dispatch, scatter-reduce, super-batch |
| `interpol_hydro.kjhan.f90` | hydro/ | GPU-compatible data layout |
| `synchro_hydro_fine.kjhan.f90` | hydro/ | GPU dispatch |

### patch/oct_tree → base RAMSES override (2)

| File | Original | Modification |
|------|----------|-------------|
| `refine_utils.f90` | amr/ | Morton hash 유지 (make_grid_coarse/fine/kill_grid) |
| `nbors_utils.kjhan.f90` | amr/ | Morton lookup 기반 neighbor 탐색 |

### patch/Horizon5-master-2 수정 (10)

| File | Modification |
|------|-------------|
| `amr_step.jaehyun.f90` | Bulk virtual exchange 호출 (uold/unew 5곳) |
| `feedback.kjhan3.f90` | SNII 27-bin spatial hash (35x speedup) |
| `flag_utils.kjhan.f90` | sink_refine bug fix (sink particle 직접 체크) |
| `init_flow_fine.f90` | Stream access IC reading (GRAFIC2) |
| `init_part.f90` | Stream access IC + MPI_ALLTOALL→ksection |
| `particle_tree.kjhan.f90` | MPI_ALLTOALL→ksection exchange |
| `multigrid_coarse.kjhan.f90` | MG coarse: merged red-black exchange, neighbor cache |
| `synchro_fine.kjhan.f90` | GPU dispatch (HYDRO_CUDA) |
| `newdt_fine.kjhan.f90` | GPU dispatch (HYDRO_CUDA) |
| `rho_fine.kjhan.f90` | GPU dispatch (HYDRO_CUDA) |

### poisson/ 수정 (2)

| File | Modification |
|------|-------------|
| `multigrid_fine_fine.kjhan.f90` | MG fine: nbor_grid_fine 사전계산, fused residual+norm |
| `multigrid_fine_coarse.kjhan.f90` | MG fine-coarse: merged red-black exchange |

---

## 3. ORIGINAL — 원본 그대로 (57)

수정 없이 사용하는 원본 RAMSES/Horizon5 파일.

### From amr/ (7)
`random.f90`, `hilbert.f90`, `title.f90`, `sort.f90`, `units.f90`, `write_gitinfo.f90`, `ramses.f90`, `init_refine.f90`

### From hydro/ (6)
`hydro_commons.f90`, `cooling_module.f90`, `init_hydro.f90`, `write_screen.f90`, `output_hydro.f90`, `condinit.f90`, `hydro_boundary.f90`, `boundana.f90`

### From pm/ (5)
`sparse_mat.f90`, `clfind_commons.f90`, `gadgetreadfile.f90`, `remove_list.f90`, `clump_finder.f90`, `clump_merger.f90`, `flag_formation_sites.f90`

### From poisson/ (4)
`poisson_parameters.f90`, `phi_fine_cg.kjhan.f90`, `gravana.f90`, `boundary_potential.kjhan.f90`, `output_poisson.f90`

### From patch/Horizon5-master-2 (27)
`pm_parameters.f90`, `pm_commons.f90`, `poisson_commons.f90`, `hydro_parameters.f90`,
`init_time.f90`, `physical_boundaries.kjhan.f90`, `light_cone.part.f90`,
`light_cone.hydro2.f90`, `light_cone.sink2.f90`, `movie_mod.yonghwi_org.f90`,
`kjhan.f90`, `output_sphere_hydro.f90`, `output_sphere_part.f90`,
`output_sphere_sink.f90`, `output_part.f90`, `move_fine.f90`,
`add_list.f90`, `star_formation.kjhan.f90`, `sink_particle.kjhan.f90`,
`init_sink.f90`, `output_sink.f90`, `init_poisson.f90`,
`interpol_phi.kjhan.f90`, `godunov_utils.kjhan.f90`,
`hydro_flag.kjhan.f90`, `uplmde.kjhan.f90`, `umuscl.kjhan.f90`,
`read_hydro_params.f90`

---

## Documentation

`../patch/cuRamses/docs/` contains the LaTeX manual:
- `main.tex` — main document (Legrand Orange Book style)
- `structure.tex` — chapter/section structure
- `main.pdf` — compiled PDF (81 pages)
- `*.png`, `*.pdf` — figures

---

## Build Instructions

```bash
cd /home/kjhan/BACKUP/cuRamses/bin

# CPU only
make -f Makefile.cuRamses HDF5=1

# CPU + GPU (CUDA)
make -f Makefile.cuRamses HDF5=1 USE_CUDA=1

# Clean
make -f Makefile.cuRamses clean
```

Compiler: Intel ifx (mpiifx) + NVIDIA nvcc
Dependencies: HDF5, CUDA 13.0+, Intel MPI

---

Created: 2026-02-16
Source: `git@github.com:kjhan0606/cuRAMSES-kjhan.git` (master branch)
