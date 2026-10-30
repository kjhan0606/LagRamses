# M5 AGN-only CPU integration evidence — 2026-10-30

## Result

Slurm job `1244059` completed successfully on `grammar029` (`normal`): 1 MPI
rank, 8 OpenMP threads, CPU CHIMES, 16 GiB requested, 1,368,528 KiB peak RSS,
55:04 elapsed. The runner reported
`M5_SOURCE_ACTIVE_DUST_CHIMES_REFERENCE_PASS` after three coarse steps, with
M5 order 5, live two-size dust, CHIMES, an actual accreting BH sink, and stellar
SED disabled to isolate the AGN radiation source. The profile uses reference
controls and is explicitly nonproduction.

The effective namelist is retained at
`.m5-source-dust-chimes-1244059/effective.nml`; input and executable digests
are in `.m5-source-dust-chimes-1244059/inputs.sha256`. The CPU executable was
built by job `1243944` from source commit
`4e4c323450ef01baeb8848011d286b58ae2a7628`; its build manifest is in
`.m5-ramses-momentum-build-1243944/build-result.txt`.

## Coupling evidence

- The AGN sidecar records positive BH accretion and positive thermal luminosity
  at coarse step 2, satisfying the source-attribution condition. No stellar
  radiation population was loaded.
- The runtime reached `Run completed`; the wrapper result records 16,896
  momentum-assertion cells and `raw_outputs=0`.
- At the source-active M5 commit, photon-number residual was
  `1.79e-17` and radiation-energy residual was `9.98e-18` in the printed
  ledger convention. At the following commit, the energy residual was
  `-5.62e-18`. No M5 accounting rejection occurred.
- The run also records `Main step=3`; the output schedule remained in the
  future and no `output_NNNNN` directories were created.

## Runtime and limits

Total reported runtime was 3,291.5 s. RAMSES attributed 3,286.2 s (99.8%) to
its `cooling` bucket. This bucket includes the coupled M5/CHIMES path and is
not a standalone CHIMES timing. It identifies that coupled path as the next
runtime comparison target.

The three level-3, order-5 M5 commits took 449.1 s (16 substeps), 2,361.0 s
(16 substeps, source-active), and 437.2 s (1 substep), totaling 3,247.3 s or
98.7% of the run. For the source-active 16-substep call, existing tile timers
sum to 848.1 s of material preparation and 499.0 s of staged commit work;
about 1,013.8 s remains outside those tile intervals and is not split further
by current timers. The summed per-cell photo durations (5,942.2 s) overlap
across OpenMP workers and are not wall time. These measurements identify the
M5+material path as the target but do not justify attributing all its cost to
CHIMES.

The printed RAMSES `econs=9.98e-1` is computed from kinetic energy, potential
energy, and the integrated potential-work term in `amr/update_time.f90`; it
does not include the gas internal, dust, radiation, or explicit feedback-source
terms. It is therefore not a closed-system energy audit for this run. The M5
photon/radiation ledger closes to the residuals above, but a total
gas+radiation+dust+AGN energy budget was not measured here.

This establishes CPU-only coupled execution and AGN source attribution under
the declared reference profile. It does not establish production physical
calibration or GPU/hybrid runtime behavior. H100 CUDA build job `441960`
passed, and hybrid integrated runtime job `442061` is pending Slurm GRES
capacity at the time of this record.
