# M5 regrid optimization validation — 2026-10-30

## Source and build identity

- Repository: `kjhan0606/LagRamses`
- Source commit: `4e4c323450ef01baeb8848011d286b58ae2a7628`
- Change: validate M5 realizability once before transactional regrid
  publication, then use the prevalidated restore path.
- CPU build job `1243944`: passed, `CHIMES_CUDA=0`; binary:
  `.m5-ramses-momentum-build-1243944/bin/ramses_m5_ramses_momentum_closure_v13d`.
- CUDA build on Lageunha Ada: passed, `CHIMES_CUDA=1`; binary and full build
  record are retained under
  `.m5-ramses-momentum-build-lageunha-20261030-m5/`.

## Integrated CPU run

Job `1244047` ran on grammar `normal`, node `grammar016`, with four CPUs,
`MPI=1`, `OMP=4`, and no GPU allocation. It used the existing
`effective_population` reference-control profile, CPU CHIMES, M5 order 5,
live dust, and `M5_AMR_INTERFACE_CHECK=1` for four coarse steps. This exercises
the coupled coarse/fine M5 path on the binary built from the changed source.

- Effective input: `.m5-effective-amr-1244047/effective.nml`
- Runtime log: `.m5-effective-amr-1244047/run.log`
- Acceptance record: `.m5-effective-amr-1244047/result.txt`
- Slurm record: `.m5-source-dust-chimes-1244047.log`
- Result: `COMPLETED 0:0`; elapsed `00:02:45`, simulation step `00:02:38`,
  `MaxRSS=781236K`.
- Runner result:
  `M5_COUPLED_AMR_INTERFACE_PASS coarse_leaf_cells=56 fine_leaf_cells=64 MPI=1 OMP=4 raw_outputs=0 reference_controls=NONPRODUCTION`.
- Four coarse steps committed M5 order 5 at levels 2 and 3. The level-3 mesh
  had eight grids; dispatch counters show CPU work and zero GPU work. Dust
  material batches likewise report CPU=1/GPU=0. The runner's finite-state,
  dust, chemistry, source, momentum, convergence, and no-dump checks passed.
- No `output_NNNNN` raw simulation directory was produced. Inputs, logs,
  effective namelist, and compact result remain for reproducibility.

This is a bounded CPU integration check, not a performance comparison or a
physical-production qualification. The selected population and material
contracts remain reference controls. It does not establish CUDA runtime
correctness: the separate CUDA build passed, but the Lageunha runtime smoke
was refused by its preflight because only 6701 MiB of Ada memory was free,
below the required 8192 MiB. No CUDA smoke simulation was launched.

The earlier CPU native regrid contract job `1243990` also completed with exit
0 and emitted `SNRT_NATIVE_REGRID_PASS`, `SNRT_NATIVE_DEFRAG_PASS`, and
`SNRT_M5_NATIVE_REGRID_PASS`; its compact result is under
`.m5-ramses-momentum-build-1243944/regrid-native-1243990/`.
