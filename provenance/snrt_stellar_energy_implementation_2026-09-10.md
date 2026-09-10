# Native stellar Q/E coupling — bounded completion

Outcome: **PASS for the independent BPASS H/He band comparison**. This is
actual Fortran source/runtime implementation, not only a Python validator.
No new default, global namelist flag, checkpoint tensor/version, audit gate,
or common-population/production-completeness claim.

## Scope and review

Workdir `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
[Plan](snrt_stellar_energy_plan_2026-09-10.md) reviewed once by
[Fable](snrt_stellar_energy_fable_2026-09-10.txt): Q-GOAL yes, Q-LEAN mostly
lean, approve. Applied the lean recommendation to use log1p without a series
branch and to reuse the checked identity dataset without a format bump.
The driver performed this end evaluation, per the operator's latest cadence.
The review's qualitative assertion about all ionizing means decreasing with
age is not a completion criterion or a measured result here: the actual
tabulated Q/E, not a presumed monotonic trend, determines the source.

## Physical source and implementation

Pinned source:
`external/bpass_v2.2.1/SSP_Spectra_BPASSv2.2.1_bin-imf135_300.hdf5`, SHA256
`b53d7bf4e8c50ae0a02458eae9d5b6dff5d5782f9f2b23dd40514ad85716c9b3`.
The [pinned upstream converter](https://github.com/galacticusorg/galacticus/blob/e8d9c46113eb2639515ecddb8fb29dea98a3989b/scripts/ssps/convertBPASSv2.2.1SSPsToGalacticus.py)
divides the original burst by 1e6 and converts to Lsun/Hz. Thus the staged
spectra are already per initial Msun. No second burst normalization or IMF
reweighting is applied. HDF5 units/axes and pinned checksum were checked.

- `simulation/snrt/tools/build_bpass_native_sed.py`: opt-in
  `--energy-moments` emits v3. Q integrates the piecewise-linear positive
  photon integrand f(lambda)=Lnu*Lsun/(h*lambda). E analytically integrates
  the **same** f times hc/lambda with positive endpoint weights and log1p.
  Both use exactly clipped source coverage and identical group boundaries;
  no boundary double-counting or extrapolated tail. Escape fraction applies
  once to both moments. 0--1 Myr duplicates the first spectrum for both.
- `patch/lagRamses/snrt_stellar_source.f90`: v3 loader, positivity/finite/
  band support checks and first-spectrum hold, explicit energy semantics,
  separate linear age/Z moment interpolation and interval integration.
  Missing energy output is rejected for v3 instead of silently discarded.
  v1/v2 cannot carry hidden unbound energy fields. MPI/HDF5 identity includes
  the complete E array and semantics, with old identities unchanged.
- `snrt_agn_source.f90`: optional actual source mean is distinct from Eref.
  State remains `Eref*N + shift`. Actual injected energy charges the shift;
  the source quantization error bound uses actual injected E. Default paired
  arithmetic and legacy number-only arithmetic remain unchanged. Invalid
  late-group inputs leave the entire source/state transaction untouched.
- `snrt_ramses_driver.f90`: requests interval Q/E when v3 is selected,
  computes guarded E/Q per source and passes it into paired deposition.
  This feeds existing H/He band transport/absorption/chemistry. AGN injection
  is unchanged. No RAMSES namelist or generator/mkrun change is necessary;
  existing `SNRT_STELLAR_SED` selects this external source asset.

New native asset `simulation/snrt/config/snrt_stellar_sed_bpass_independent_v3.nml`:
SHA256 `45cd3cb496832cc8061ebe817bd2576a5badc979eb1a9735ef32e235cd1ca99d`.
Sidecar `simulation/snrt/data/snrt_stellar_sed_bpass_independent_v3.json`
records units, physical provenance and limitations. 9 groups, 52 ages,
13 metallicities. Rebuilt default v2 is byte-identical to the retained v2:
SHA256 `ebd5e097e0e9eda50ca029747cb12d89f0a1d89ac02dda665f1381fcccbab840`.

## Verification and retained evidence

All task artifacts are under `.rt-stellar-energy.kzyggZ/`.

Checked CPU build:

```sh
make -f ../bin/Makefile -j1 NVECTOR=32 NENER=1 NVAR=19 SNRT=1 \
  DUST_LIVE=0 HDF5=1 USE_FFTW=0 FDMDEBUG=1 EXEC=ramses_stellar_energy
```

Binary `ramses_stellar_energy3d`, SHA256
`a0996eeaf2b3c6225c2cd9a614f953eb4fddf6fdf2369fc7eefb4631c592bdd7`.
`build.log` records ifx bounds checking, MPI, OpenMP and link inputs. Existing
MPI symbol-length warnings occur; no build errors. A subsequent source-header
comment clarification does not change runtime semantics. Makefile/VPATH was
not changed by this bundle.

- `source.log`: actual below/above-reference energy injection, unchanged
  FP32 photon/default paired results, rollback and actual-energy rounding
  bound pass. Includes a case that would wrongly pass if the rounding bound
  used Eref instead of actual E, plus negligible/entirely-lost source tails.
- `stellar.log`: analytic energy integrals including narrow wavelengths,
  exact v2/v3 photon parity, native split-interval Q/E, interior-metallicity
  Q/E interpolation, mean-of-means counterexample, legacy behavior, 9 old
  and 8 new mode/moment rejection cases pass. Existing native smoke/tests
  extended; no new test framework.
- `completed/physical.nml`, `run.env.sh`: 4^3 noncosmo hydro/stellar feedback/
  advective CR plus BPASS v3 and H/He band RT. No AGN/sinks, dust, CHIMES or
  ordinary cooling. MPI2/OMP2, `SNRT_BACKEND=auto` selects CPU/OpenMP in this
  CPU build. Four steps, noutput1/aout2/tout1e30/foutput2/fbackup1e6.
  Continuous wall time 13.9482 s. Source log explicitly confirms BPASS Q/E.
  Nonzero-source transaction iterations 20/19/12; no failed RT transactions.
- `completed-restart/physical.nml`: restart from step2 to step4, wall6.19544s.
  `evaluation.txt`: **all 82 non-header physical datasets exactly match**:
  AMR12, coarse3, domain1, gravity8, hydro38, particles17, SNRT3.
  Header times/step counters and SNRT attributes match. HDF5 format remains15,
  band model `hhe_maxent64_v1`, cell width1444. Arrays are finite, N/E are
  nonnegative, zero N has no hidden E, every mean stays within its band,
  density and gas internal energy remain positive. Step4 minimum density
  0.000999823395567; minimum gas internal energy6.21625289524e-7 (code units).
- `reject-v2/run.log`: v3 dump with v2 SED rejected at checked dataset extent
  (code14). `reject-energy/run.log`: same Q but one physically admissible
  E coefficient changed by0.1%, rejected by identity values (code1).
  Both MPI exit10 expected; neither evolves nor writes output00002.
  The modified E input independently passes the normal native source loader,
  so this specifically tests restart identity, not malformed-source rejection.

The two reader checks use the same binary/physics, MPI2/OMP2,
nrestart1/nstepmax0, noutput1/aout2/tout1e30/foutput=fbackup1e6. Launch settings
and storage policy were reported before every RAMSES launch. Free /gpfs
space was about59TiB; bounded raw output was under6MiB.

## Limits and handoff

This fills the actual stellar spectral **moment injection** gap in approved
group8. Group1 is still incomplete: BPASS binary imf135_300 and the feedback
LC18/Monash population do not become one population. The v3 source remains
`reference_control`, not production approved. The opt-in H/He closure does
not reproduce a full spectrum and still excludes DUST_LIVE and CHIMES;
node-resolved secondary electrons and stellar/dust spectral coupling remain
in their existing scope. No GPU performance or resolution convergence claim.
Do not silently enable v3 for existing fixed/dust runs or rewrite defaults.

Driver end evaluation complete. Retain inputs, logs, native test executables,
build identity, evaluator and compact numbers. Remove only the three evaluated
raw HDF5 dumps; see `snrt_stellar_energy_raw_cleanup_2026-09-10.md`.
No commit or push was requested in this continuation, and none was performed.
