# Project memo

## Paper split and TeX source locations

- Paper I: SIDM
  - Overleaf/source directory: `/home/kjhan/paper_sidm_overleaf`
  - Main TeX source: `/home/kjhan/paper_sidm_overleaf/main.tex`
  - Overleaf git remote: `https://git@git.overleaf.com/6a64afa4af99ca2536a1f4c6`

- Paper II: DE/nGR
  - Overleaf/source directory: `/home/kjhan/paper_nonstd_overleaf`
  - Main TeX source: `/home/kjhan/paper_nonstd_overleaf/main.tex`

## DE/nGR DMO validation status

- Main uniform-resolution validation campaign:
  - `/gpfs/kjhan/Hydro/DE_nonstd/DMO_lagCAMB_resolution_ladder_20260726/phase_anchor8_all`
  - Box/seed: `64 h^-1 Mpc`, seed `20260726`
  - Phase anchor: level 8, shared by L6/L7/L8
  - Models: `lcdm`, `f5`, `f6`, `n1`, `n5`, `sym_a`
  - Resolutions: L6=`64^3`, L7=`128^3`, L8=`256^3`
  - All models have z=0 outputs and measured CIC spectra; `lcdm` and `n1`
    in `phase_anchor8_all` are symlinked from `phase_anchor8_n1`.
- Final audit result as of 2026-07-29:
  - Every model uses the common primordial amplitude `A_s=2.1e-9`.
    Model-specific lagCAMB spectra are passed to LagMUSIC with
    `force_pnorm`; `sigma_8` is retained as a derived model prediction.
  - All 18 required L6/L7/L8 model-resolution combinations have an exact
    z=0 output with `match_aout=.true.`.
  - Every non-standard model passes the 1% L7->L8 convergence target over
    the full measured range `k <= 0.5 h Mpc^-1`.
  - z=0 L7->L8 low-k target is 0.1% for `k <= 0.2 h Mpc^-1`.
  - Passing low-k 0.1%: `f6`, `n1`, `n5`, `sym_a`.
  - F5 initially missed the low-k diagnostic at L7->L8 with a maximum
    residual of `0.1315%`. The direct L8->L9 extension now passes with a
    maximum low-k residual of `0.09336%`.
  - Full `k <= 0.5 h Mpc^-1` 0.1% is stricter and still fails for several
    models; use it as diagnostic, not the immediate low-k gate.
- F5 targeted extension:
  - `/gpfs/kjhan/Hydro/DE_nonstd/DMO_lagCAMB_resolution_ladder_20260726/phase_anchor9_f5_0p1`
  - Models: `lcdm`, `f5`
  - Resolutions: L8=`256^3`, L9=`512^3`
  - Phase anchor: level 9
  - F5 tolerance: `fR_eps=1e-5`; this was A/B-tested against `1e-6` at L8
    and changes z=0 raw P(k) by at most `0.03056%`, within the 0.1% budget.
  - Grammar Slurm job `389321` completed normally on 2026-07-29 with
    exact z=0 `output_00003` and exit code 0.
  - Analysis job `389324` produced all requested spectra and reports.
    Slurm labelled the job `FAILED` only because its wrapper also enforced
    the optional 0.1% threshold over the full `k <= 0.5 h Mpc^-1` range.
    The authoritative audit exited 0 with `complete=true` and no failures.

## Completed resolution-validation objective

The uniform global-spectrum gate is complete. The machine-readable result is

`/gpfs/kjhan/Hydro/DE_nonstd/DMO_lagCAMB_resolution_ladder_20260726/UNIFORM_LADDER_AUDIT.json`.

AMR should not be added to improve the global power spectrum. Add AMR only
when the next selected science analysis requires halo-scale screening or
internal halo structure.

## DE/nGR production plan

- Main volume: `256 h^-1 Mpc` with `512^3` particles and a uniform
  `512^3` force mesh.
- Every model uses the same stored Gaussian white-noise realization.
  Individual Fourier amplitudes remain Gaussian and are not fixed. No
  phase-reversed partner is used.
- Every model uses common `A_s` and `n_s` with a model-specific lagCAMB
  transfer function.
- One simulation is run for each of the following 16 models:
  1. LCDM
  2. smooth CPL
  3. clustering CPL
  4. coupled quintessence
  5. parametrized Horndeski gravity
  6. F5
  7. F6
  8. N1
  9. N5
  10. Symmetron A
  11. massive-neutrino LCDM
  12. F6 with the same massive-neutrino sector
  13. early dark energy
  14. running vacuum
  15. generalized Chaplygin gas
  16. Ratra--Peebles \(\phi\)CDM tracker (`phicdm_a01`, \(\alpha=0.1\))
- Planned outputs are `z=2`, `1`, `0.5` and `0`.
- Production starts only after the background, linear growth, transfer
  function and low-resolution diagnostic checks pass for every selected
  benchmark.

### phiCDM production gate

- `phicdm_a01` is the physical inverse-power-law scalar-field benchmark:
  \(V=A\phi^{-0.1}\), with the matter-era tracker initial condition and
  \(A\) shot to the present dark-energy density.
- The old `q1` case starts from a frozen, user-selected `phi_ini`; it remains
  a regression test but is not the production phiCDM model.
- Required checks are: present-day closure, the matter-era tracker limit,
  lagCAMB/lagRamses background agreement over the N-body interval, finite
  model-specific transfer functions at common \(A_s\), low-resolution
  completion to \(z=0\), and measured \(P(k)\) against lagCAMB linear theory.
- Gate completed on 2026-07-29:
  `/gpfs/kjhan/Hydro/DE_nonstd/PHICDM_validation_20260729/resolution_anchor6/L6_64/phicdm_a01_validation.json`.
  All eleven checks pass. The largest lagRamses/lagCAMB background-density
  residual is 0.1103% at z=49 (CAMB radiation versus the standard matter+DE
  N-body background); the equation-of-state residual is \(3.69\times10^{-4}\).
  The common-phase initial power-ratio residual is 0.0203%, and the
  normal-versus-doubled-accuracy lagCAMB linear-ratio residual is 0.0288%.
- An independent production-equivalent FFTW smoke run also completed through
  \(z=0\), with the base-level FFT solver explicitly active, and passed the
  same eleven checks:
  `/gpfs/kjhan/Hydro/DE_nonstd/PHICDM_validation_20260729/fftw_prodcheck_L5/phicdm_a01_validation.json`.
- The phase-anchored L5->L6 z=0 test passes the strict 0.1% convergence gate
  over \(k\leq0.2\,h\,{\rm Mpc}^{-1}\) (maximum 0.0412%) and the 1% gate over
  \(k\leq0.5\,h\,{\rm Mpc}^{-1}\) (maximum 0.6003%). This is a low-resolution
  code gate, not a replacement for the production-volume resolution audit.

### Coupled-DE velocity gate

- LagMUSIC DMO 2LPT now offers `dmo_velocity_source=transfer|density_2lpt`.
  The default `transfer` path uses lagCAMB's mass-weighted `vtotal(k)` for
  the linear velocity and adds MUSIC's 2LPT velocity term; the legacy
  density-derived path remains user-selectable.
- `cde10` has a 1.7467% scale dependence in its velocity-growth correction,
  so it must use the default transfer path rather than one scalar
  `vfact_scale`.
- The phase-anchored L6 gate completed through \(z=0\) with all eight checks:
  `/gpfs/kjhan/Hydro/DE_nonstd/CDE_velocity_validation_20260729/resolution_anchor6/L6_64/cde10_velocity_validation.json`.
  The maximum lagMUSIC/lagCAMB velocity-ratio residual is 0.0433%, and the
  initial density-transfer residual is 0.0324%.
- The L5->L6 \(z=0\) ratio comparison is not yet a resolution certificate:
  its maximum residual is 0.230% over \(k\leq0.2\,h\,{\rm Mpc}^{-1}\) and
  7.29% over \(k\leq0.5\,h\,{\rm Mpc}^{-1}\). This does not invalidate the
  velocity kernel, which independently passes the 0.1% linear gate, but a
  production-scale resolution audit remains required.
