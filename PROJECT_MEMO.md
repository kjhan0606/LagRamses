# Project memo

## NEXT TASK: regenerate and relaunch the corrected 1024^3 CPL production

The next DE/nGR production task is to replace the mixed-background CPL IC
without changing its Gaussian realization:

1. Work on Grammar and use a 515-GB normal node. The 258-GB debug node is
   sufficient for the validated `512^3` run but has no safe memory margin for
   the production `1024^3` IC or simulation.
2. Reuse the existing level-10 white-noise file
   `/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729/wnoise_0010.bin`.
   Do not regenerate from the integer seed alone, because that would not
   guarantee the exact production realization.
3. Generate the CPL IC into a new directory such as
   `ics_cpl_m09_p02_corrected`. Do not overwrite or delete the old
   `ics_cpl_m09_p02` until the replacement has passed all checks.
4. The LagMUSIC configuration must contain `w0=-0.9`, `wa=0.2`,
   `dmo_velocity_source=transfer`, the existing common-\(A_s\)
   `force_pnorm`, and the stored white-noise path. Record the LagMUSIC binary
   checksum and require the final log to report the same CPL values.
5. Validate all seven non-empty GRAFIC level-10 components, the particle
   count and box metadata, the shared-white-noise provenance, finite
   displacement and velocity fields, and the expected CPL velocity factor.
6. Point only the CPL production namelist at the corrected IC. Keep
   `dump_pk=.false.` and retain exact outputs at `z=2,1,0.5,0`.
7. Submit CPL followed by F6, N1, Symmetron A and Chaplygin through `afterok`
   dependencies. Exclude at least
   `grammar007,grammar012,grammar023,grammar081,grammar112`, which have
   failed the production startup or are already known bad nodes.
8. After the CPL run reaches its first output, verify the distributed FFT,
   per-rank memory, timestamp progress and empty error log before treating
   the remaining Lane-B dependencies as healthy.
9. Measure production spectra from snapshots with the external CIC
   post-processor. Do not re-enable the replicated in-situ `1024^3`
   estimator.

The implementation and validation record for the hand-off fix is
`/gpfs/kjhan/Hydro/DE_nonstd/CPL_debug_L512_N512_20260731/VALIDATION.md`.
The corrected `512^3` CPL run reached exact `z=0` with exit code 0.

## ALWAYS FIRST: identify the current server

- The user runs Codex concurrently on several shared-filesystem servers.
- At the beginning of every session and before any server/process/Slurm
  status report, run `hostname -f` locally and state the current server at
  the very beginning of the report.
- Do not infer the active server from the working-directory path, shared
  `/home`, shared `/gpfs`, or previous-session context.
- Keep the cluster roles distinct:
  - `LagEunha`: manual runs, shared `/home` and `/gpfs`, 1 TiB memory.
  - `grammar`: CPU Slurm cluster.
  - GPU cluster/login nodes are not Grammar compute nodes.

## NEXT SESSION: automatic simulation status report

At the beginning of the next session, report this campaign's status before
starting new work:

- Campaign:
  `/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729`
- Grammar IC job `394194` completed successfully.
- The original simulation chain `394195`--`394205` was retired after
  `grammar112` failed and was drained (`Kill task failed`).
- Replacement Lane A: `395189` LCDM -> `395190` F5 -> `395191` N5 ->
  `395192` coupled DE -> `395193` RVM -> `395194` phiCDM
- Replacement Lane B: `395195` smooth CPL -> `395196` F6 -> `395197` N1 ->
  `395198` Symmetron A -> `395199` Chaplygin
- Every replacement job excludes `grammar112` and uses `afterok`
  dependencies. At most two simulations, one per lane, run concurrently.
- The two replacement lane heads, `395189` (LCDM) and `395195` (smooth
  CPL), failed on `grammar081` during the initial in-situ 1024^3 power
  spectrum.  The estimator replicated two 1024^3 real arrays on every MPI
  rank and exhausted node memory; this was not a cosmological-model failure.
- Production namelists now set `dump_pk=.false.`.  Spectra will be measured
  from snapshots with distributed post-processing.  The source estimator
  also refuses replicated FFT grids larger than 512^3.
- Current replacement chains (submitted 2026-07-30; excluding
  `grammar012,grammar081,grammar112`):
  - Lane A: `395979` LCDM -> `395980` F5 -> `395981` N5 ->
    `395982` coupled DE -> `395983` RVM -> `395984` phiCDM
  - Lane B: `395985` smooth CPL -> `395986` F6 -> `395987` N1 ->
    `395988` Symmetron A -> `395989` Chaplygin
- At submission both lane heads were pending for priority; no repeated
  monitoring loop was started.

The startup report must check grammar directly and include:

1. job state, elapsed time, node and pending/failure reason from `squeue` or
   `sacct`;
2. for IC generation, the model currently being generated, number and disk
   size of completed `ics_*` directories, and the tail of the current log;
3. for each running simulation, its latest dump/output number and the latest
   scale factor/redshift found in its log or dump metadata;
4. non-empty error logs, abnormal exits, memory failures, stalled timestamps,
   or unsatisfied dependencies;
5. which jobs are running and which job is next in each lane.

This is a one-time startup report. Do not start a repeated polling loop or
send periodic duplicate reports unless the user explicitly requests
monitoring.

## Paper split and TeX source locations

- Paper I: SIDM
  - Overleaf/source directory: `/home/kjhan/paper_sidm_overleaf`
  - Main TeX source: `/home/kjhan/paper_sidm_overleaf/main.tex`
  - Overleaf git remote: `https://git@git.overleaf.com/6a64afa4af99ca2536a1f4c6`

- Paper II: DE/nGR
  - Overleaf/source directory: `/home/kjhan/paper_nonstd_overleaf`
  - Main TeX source: `/home/kjhan/paper_nonstd_overleaf/main.tex`

## Paper I SIDM figure plan (agreed 2026-07-31)

The intended main-text figure set is:

1. SIDM interaction models and Monte Carlo verification (retain the current
   implementation-validation figure).
2. The collisionless \(\Lambda\)CDM \(z=0\) full-volume slab and zoom-region
   density map (retain the current figure, and explicitly call the reference
   model "collisionless \(\Lambda\)CDM").
3. At a common \(z=0\), matched central projected-density maps for
   \(\Lambda\)CDM, SIDM1 and SIDM3.  Use the same centre, line-of-sight
   thickness and colour range, a field of view of roughly
   \(100\)--\(300\,h^{-1}{\rm kpc}\), and include a lower row showing
   \(\Sigma_{\rm SIDM}/\Sigma_{\rm CDM}\) or the corresponding difference.
4. Common-\(z=0\) radial profiles of dark-matter density, SIDM/CDM density
   ratio, one-dimensional velocity dispersion and, if practical, circular
   velocity.  This should replace the provisional \(z\simeq4\) production
   profile figure once SIDM1 and SIDM3 reach \(z=0\).
5. Redshift evolution of the core and cumulative scattering: use snapshots
   such as \(z=4,3,2,1,0\), and show either the evolving density profiles or
   \(r_{\rm core}(z)\) and \(\rho_{\rm core}(z)\), together with the mean or
   cumulative number of scatterings per particle versus radius.
6. Numerical convergence: compare density/core measurements across force or
   mass resolution, mark the adopted convergence radius, and demonstrate
   that the measured core is not set by the cell size or particle sampling.

Appendix/supporting figure candidates, in priority order, are:

- SIDM timestep/probability convergence under changes to `sidm_courant` or
  the \(P_{\max}\) limit.
- Matched-halo and mass-conservation diagnostics, including convergence of
  outer profile ratios to unity.
- Load-balancing/performance results: coarse-step wall time and rank-load
  distributions for memory/work/minimax modes.
- ADM cooling/heating curves, since the paper describes the ADM module.
- iSIDM elastic, up-scattering and down-scattering rates versus relative
  velocity, including the inelastic threshold.

Additional science figures are optional rather than required for the first
submission: radial phase space or \(Q=\rho/\sigma^3\), velocity anisotropy
\(\beta(r)\), mass-assembly history versus core growth, and subhalo mass or
radial distributions.  Large-Mpc SIDM density maps should not be added merely
as morphology panels because they will be nearly indistinguishable from CDM;
SIDM comparison maps should focus on the resolved halo core.

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

- Submitted main volume: `512 h^-1 Mpc` with `1024^3` particles and a
  uniform `1024^3` force mesh.
- Campaign path:
  `/gpfs/kjhan/Hydro/DE_nonstd/DMO_production_L512_N1024_20260729`.
- Submitted on grammar on 2026-07-29 as one IC job followed by two
  dependency lanes. Eleven parameter-identical lagCAMB/lagRamses models are
  currently in the submitted set.
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
- The five cases not yet submitted are clustering CPL, parametrized
  Horndeski, massive-neutrino LCDM, F6+neutrino and EDE. Their exact
  transfer/perturbation input gates must be completed before submission;
  they must not silently use an LCDM fallback transfer.

### CPL IC hand-off fix and debug validation (2026-07-31)

- The CPL lagCAMB density and velocity transfers were correct, but
  `patch/lagRamses/aux/dmo_benchmark_setup.py` hard-coded `w0=-1`, `wa=0`
  in the LagMUSIC configuration. LagMUSIC itself already reads and uses
  `w0` and `wa` in `H(a)`, growth, and `vfact`.
- The generator now passes `w0=-0.9`, `wa=0.2` for `cpl_m09_p02` and
  `w0=-0.9`, `wa=0` for `w09`. The old production `1024^3` CPL IC must be
  regenerated. Its velocity-factor error is small, 0.00525% at z=49, but
  the mixed-background IC is not production-consistent.
- Grammar debug jobs `397240` and `397241` passed the full `32^3` evolution
  and the `256^3` distributed-FFTW startup test.
- Corrected production-box test:
  `/gpfs/kjhan/Hydro/DE_nonstd/CPL_debug_L512_N512_20260731`.
  IC job `397247` and `512^3` run `397248` both completed with exit code 0.
  The simulation reached exact `z=2,1,0.5,0`, used at most 1.6 GB/rank,
  and retained zero mass drift. CIC analysis job `397260` also completed.
- Full details and the growth-spectrum caveat are in
  `/gpfs/kjhan/Hydro/DE_nonstd/CPL_debug_L512_N512_20260731/VALIDATION.md`.
- The successful `512^3` run does not reproduce the production crash.
  Evidence points to the node/memory path during the first `1024^3`
  distributed FFT, not to the CPL background equations.

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

## Lageunha SIDM load-balance restart handoff (2026-07-30 06:31 KST)

- The two active Paper-I SIDM zoom runs are manual Lageunha jobs:
  - `Zoom0/zoom_run_sidm1`: 32 MPI ranks x 2 OpenMP threads on NUMA1,
    \(\sigma/m=1\).
  - `Zoom0/zoom_run_sidm3`: 24 MPI ranks x 2 OpenMP threads on NUMA0,
    \(\sigma/m=3\).
- Both were asked to checkpoint and stop gracefully at the next coarse-step
  boundary by writing `0 -1` to each run's `jobcontrol.txt`.  At the request
  time SIDM1 was inside a load-balance operation after fine step 37987
  (\(a=0.3013\)); SIDM3 was in main step 655 near fine step 39428
  (\(a=0.2471\)).  Do not send a kill signal while the graceful stop/output
  is still in progress.
- Before restart, verify all of the following on Lageunha: both launcher PIDs
  have exited, their `.rc` files contain `0`, the logs contain the job-control
  stop message, and the new output directory is complete on every MPI rank.
- The load-cost fix is currently an uncommitted source change in:
  - `patch/lagRamses/amr_parameters.jaehyun.f90`
  - `patch/lagRamses/load_balance.kjhan.f90`
  - `patch/lagRamses/amr_step.jaehyun.f90`
  - `patch/cuRamses/bisection.f90`
  - unit test `tests/load_balance/test_domain_leaf_cost.f90`
- The fix gives Hilbert, bisection/ksection, and the imbalance trigger the
  same per-leaf cost.  Memory mode no longer multiplies allocated-memory cost
  by AMR subcycles; work mode uses the common grid/particle proxy and applies
  subcycle weighting.
- Four 64-cubed, two-rank integration runs passed for
  Hilbert/ksection x memory/work, and the unit test passed.  The newly linked
  executable is
  `/home/kjhan/BACKUP/lagRamses-de-nonstd/bin/ramses_final3d`
  (built with `NVAR=11`, `LONGINT`, `QUADHILBERT`, and linked to FFTW).
- On Lageunha, preserve the old `ramses_zoom3d` before replacing it.  SIDM3's
  executable resolves to SIDM1's file, so replace the shared target only once.
  Set each namelist's `nrestart` to the newly completed output number, clear
  `jobcontrol.txt`, and relaunch with the existing manual scripts so the
  original MPI/OpenMP affinity is retained.
- For the speed comparison, exclude checkpoint I/O, restart initialization,
  and the first corrective remap.  Compare at least three subsequent complete
  coarse steps against the pre-fix reference of roughly 432 s (SIDM1) and
  396 s (SIDM3), and separately report load-balance wall time and whether
  repeated remapping ceased.

## DE/nGR validation Level-up (2026-08-25)

- Every Paper-I DE/nGR entry formerly at Level 1 or 2 now has Level 3 or
  stronger evidence. The gate root is
  `/gpfs/kjhan/Hydro/DE_nonstd/DE_level3_solver_gate_20260825`.
- Massive-neutrino response and F6+neutrino are Level 3 diagnostics. Future
  tables must be generated with CAMB
  `accurate_massive_neutrino_transfers=True`; the independent CLASS gate has
  maximum per-redshift median `R_nu` error 0.765% and maximum applied-source
  error 0.00349% over the paper domain.
- QUMOND is Level 3 on a uniform isolated spherical two-Poisson test. Dilaton
  and tracker Galileon are Level 4 uniform diagnostics after isolated
  screening tests and residual-qualified 32-cubed end-to-end runs.
- Machine records are `patch/cuRamses/aux/DE_NGR_LEVEL3_SOLVER_AUDIT.json`,
  `COSMOLOGICAL_SOLVER_GATE_AUDIT.json`, and
  `NEUTRINO_CLASS_LEVEL3_AUDIT.json`. These levels do not certify AMR
  coarse--fine behaviour or agreement with an independent nonlinear code.
