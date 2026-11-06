# BIMOND quasi-static implementation design for lagRamses

**Status:** uniform CPU two-sector field and short particle-evolution paths validated  
**Date:** 2026-11-06  
**Revision:** 2026-11-07 after user-authorized Astra plan review  
**Project:** RAR comparison of standard Lambda CDM, QUMOND, and later FDM

## Q-GOAL

Add a BIMOND model to lagRamses while preserving the user's first scientific comparison: identical prescribed stellar and gas distributions, a CDM halo for standard Lambda CDM, and the full QUMOND equation for the baryon-only MOND case. Add the BIMOND branch without relabeling the existing QUMOND solver as a complete relativistic BIMOND cosmology. Preserve the later FDM extension and the requested hydro contract of stellar and AGN feedback, metal enrichment, and cooling/heating, with dust and ray tracing excluded.

## Selected equation-level branch

Use the BIMOND nonrelativistic branch with `alpha + beta = 0` and fix `beta = 1`. Milgrom identifies the branch as having QUMOND as its nonrelativistic limit. For two nonrelativistic matter sectors, the branch has a differential potential and two sector potentials. Write `rho` for ordinary-sector density and `rho_hat` for twin-matter density. Define `Phi_* = Phi - Phi_hat` and `M'(z) = nu(sqrt(z)) - 1`, where `z = |grad Phi_*|^2/a0^2`.

The physical-coordinate equations are

```text
Laplacian(Phi_*) = 4 pi G (rho - rho_hat)
Laplacian(Phi)   = 4 pi G rho      + div[M'(z) grad(Phi_*)]
Laplacian(Phi_hat) = 4 pi G rho_hat + div[M'(z) grad(Phi_*)]
```

Ordinary matter follows acceleration `-grad(Phi)`. Twin matter follows acceleration `-grad(Phi_hat)`. The shared phantom source is

```text
rho_ph = (4 pi G)^(-1) div[M'(z) grad(Phi_*)].
```

For `rho_hat = 0`, `Phi_*` is the ordinary Newtonian potential and `Phi` obeys the full QUMOND equation. The current `mond_type=1` source is therefore the ordinary-sector, no-twin-matter limit of the selected BIMOND branch. For nonzero twin density, the current single-potential QUMOND solver is insufficient because the ordinary and twin sectors feel different potentials. [Milgrom 2013, Sec. IV.4](https://arxiv.org/html/1310.3373v1), [Milgrom 2010, Sec. V.2](https://arxiv.org/html/1001.4444)

For expanding coordinates `x`, physical peculiar potentials, and an externally specified scale factor `a(t)`, the derived quasi-static equations are

```text
Laplacian_x(Phi_*) = 4 pi G a^2 (delta_rho - delta_rho_hat)
Laplacian_x(Phi) = 4 pi G a^2 delta_rho
                     + div_x[M'(z) grad_x(Phi_*)]
Laplacian_x(Phi_hat) = 4 pi G a^2 delta_rho_hat
                         + div_x[M'(z) grad_x(Phi_*)]
z = |grad_x(Phi_*)|^2 / (a a0)^2.
```

The comoving form follows from applying `grad_r = a^(-1) grad_x` to the physical nonrelativistic equations. The implementation must verify the code-unit factors against the existing `cosmo_poisson_fourpi` convention before any scientific run. The quasi-static equations do not determine the BIMOND Friedmann background or the relativistic transfer functions.

## Scope and model labels

The lagRamses implementation targets the nonrelativistic, quasi-static two-sector field equations on an explicitly supplied background. The first RAR comparison uses no twin matter, so the ordinary force must match the existing full QUMOND result. Two-sector runs require a declared twin-matter density and velocity realization. The user has not selected twin-sector abundances or primordial phases.

The current cosmological BIMOND perturbation calculation located in the literature treats a different restricted branch with `alpha = beta` and equal background densities. The published calculation reports nonlinear growth of sector-difference perturbations at arbitrarily small amplitude and requires twin-sector initial conditions. The `alpha = beta` transfer/growth equations cannot supply initial conditions for the selected `alpha + beta = 0` branch. [Milgrom 2010, cosmological fluctuation growth](https://arxiv.org/abs/1006.3809), [Milgrom 2013, Sec. V.1](https://arxiv.org/html/1310.3373v1)

The first implementation must be labeled **BIMOND nonrelativistic quasi-static branch on a prescribed background**. The implementation must not be labeled a relativistic BIMOND cosmology, a validated BIMOND Boltzmann solution, or self-consistent cosmological initial conditions. A branch-specific relativistic action, cosmological background, twin-sector contents, stable perturbation theory, and theory-matched transfer functions remain separate requirements. Earlier bimetric MOND formulations have a published generic ghost-instability warning. The warning does not exclude every BIMOND action, but requires stability analysis of the selected action before a relativistic cosmological claim. [D'Ambrosio et al. 2020](https://arxiv.org/abs/2004.00888)

The first `PTYPE_TWIN` implementation is gravity-only and collisionless. Twin particles receive the twin-sector acceleration and contribute only to `rho_hat`. Ordinary-DM microphysics, including SIDM scattering, ADM interactions, and dark cooling, must not act on `PTYPE_TWIN` by default. At startup, every enabled non-gravity module using a generic DM-like predicate must either explicitly exclude `PTYPE_TWIN` or trigger a fatal compatibility error. A later twin-sector microphysics model requires its own physical definition and code path.

## lagRamses implementation outline

1. Add an explicit `use_bimond` selection and validate mutual exclusion with other modified-gravity modes. Fix `beta=1` for the first implementation. Use the existing simple or standard MOND interpolation choice through `nu`; record the exact mapping `M'(z) = nu(sqrt(z)) - 1` in namelist and run metadata.
2. Add `PTYPE_TWIN` as a distinct value in the existing one-byte `ptypep` species code. Current code values are `PTYPE_DM=0`, `PTYPE_STAR=1`, `PTYPE_SINK=2`, `PTYPE_ISIDM_EX1=10`, `PTYPE_ISIDM_EX2=11`, and `PTYPE_ADM=20`. Use `PTYPE_TWIN=30` after confirming no collision on the implementation branch. Preserve `ptypep` through IC loading, particle migration, HDF5 restart, and output. Reject a nonzero-twin run when the IC contains no explicitly tagged twin particles. Keep ordinary gas and stars in the ordinary sector. Do not reinterpret ordinary CDM particles as twin matter.
3. Add separate ordinary and twin density deposits from the same particle and gas epoch. The CPU deposit receives `ptypep`; the CUDA CIC kernel must receive and apply a matching species mask. Check ordinary, twin, and total deposited masses independently.
4. Solve the signed differential Poisson source from the two density contrasts, compute the shared phantom density from the differential force, then solve and retain separate ordinary and twin potentials and force fields. Preserve source arrays across the existing in-place `rho` and `phi` solver calls. Use a common mesh epoch for both sector densities.
5. Apply the ordinary force to ordinary particles and gas. Apply the twin force only to `PTYPE_TWIN` particles. Keep sector potentials in transient mesh arrays and recompute potentials after restart. Preserve `ptypep`, positions, velocities, and masses in restart state. Make potential-field output optional for diagnostics.
6. Add CPU and GPU paths with the same equations and controls. The current `code/lagRamses` checkout is three commits behind `origin/main`. Base implementation on an isolated branch from the reviewed upstream revision, without altering active runs or existing outputs.
7. Keep the cosmological background external and explicitly named. Do not emit a model-specific CAMB transfer file until a valid relativistic perturbation theory and initial mode prescription are selected.

## Initial no-twin implementation and validation (2026-11-07)

The no-twin status in the present historical subsection is superseded by the native two-sector implementation subsection.

The active source branch is `codex/bimond-nr`, based on reviewed `origin/main` revision `a733cd1`. The source currently exposes `use_bimond=.true.` in `&MOND_PARAMS`. Selecting it activates `mond_type=1` and prints the selected `alpha + beta = 0`, `beta = 1` branch. This is the exact `rho_hat = 0` limit of the branch and uses the existing full QUMOND phantom-density solver.

The active `init_part` path rejects massive nonstellar/non-sink particles when this no-twin mode is selected, so a CDM particle load cannot silently be described as baryon-only BIMOND. Stars, sink particles, and gas remain in the ordinary sector. No `PTYPE_TWIN`, twin-density deposition, separate sector potentials, or species-specific twin force path has been implemented. Nonzero twin-matter runs therefore remain unsupported and must not be launched.

Example selector:

```fortran
&MOND_PARAMS
  use_bimond = .true.
  a0_mond = 1.2d-8
  mond_mu_type = 1
/
```

The full two-sector design remains the approved follow-on target; the no-twin selector does not constitute a new force prediction beyond QUMOND. It also does not add a relativistic action, cosmological transfer functions, or theory-matched ICs.

The continuation audit found that the base revision lacked the previously validated blocked-cell QUMOND neighbor-index correction. The correction from `248eeef` is now applied: both neighbor cell addresses use `ICELL_OF(actual_neighbor_grid, child)`. The no-twin particle check is also present in the alternate cuRamses initializer. The selector requires three dimensions and a solved potential (`gravity_type<=0`), and rejects FDM, ADM, SIDM and neutrino-response sources for this baryon-only mode.

The isolated CPU build source is `build/bimond_no_twin_20261107` under the RAR root, based on `a733cd16c31b0026aa2a8c3542efa791c40882f2` plus patch SHA-256 `cdf0ea1db7ba29c2a1df7a07de4dd384ddc94a53e8f459c6d29b5f575d33b6f3`. Grammar build job `1273823` completed in 2 min 10 s with exit `0:0`. Dependent verification job `1273825` completed with exit `0:0`. Both selectors ran the same static uniform `64^3` source; the BIMOND input set `use_mond=.false.` explicitly to exercise automatic selection of the full solver. The force identity has zero RMS and maximum difference. Against the independent discrete periodic reference, normalized RMS and maximum vector errors are `7.41184e-14` and `7.01191e-14`, respectively, with equal-volume weights over all active cells. All retained run-manifest checksums pass.

The verification reserved 1 GiB for two approximately 56 MB snapshots. The effective inputs, logs, snapshots, patch and result remain under `/gpfs/kjhan/RAR_QUMOND/runs/bimond_limit_1273825/`; compact records are copied to `reports/bimond_limit_1273825/` under the RAR root. Particle admission, AMR evolution, GPU parity and nonzero twin matter are not covered by this force comparison. The successful build and no-twin test do not complete the two-sector implementation.

## Native two-sector implementation (2026-11-07)

The earlier review label was incorrect. The initial reviewer was a Codex subagent and did not invoke Fable. On 2026-11-07, the user explicitly approved Astra as the replacement reviewer. Astra approved the two-solve nonrelativistic design after reviewing Q-GOAL and Q-LEAN. The user also specified driver self-review whenever Astra serves as driver. The review does not approve a relativistic cosmological model.

The new source adds `bimond_twin_matter` and particle family `PTYPE_TWIN=30`. The ordinary and twin CIC deposits use the same positions and weights. Each sector density is accumulated independently to preserve exact cancellation for identical deposits. The signed source produces `Phi_*`. A second Poisson solve uses the preserved ordinary density plus the shared phantom source. The linear force operator gives `f_twin=f_ordinary-f_*`, and the twin potential is `Phi_twin=Phi_ordinary-Phi_*`. No third solve is required.

Both particle half-kicks select the sector force. The timestep receives an acceleration bound from both fields. The ASCII particle input requires eight columns, `x y z vx vy vz mass ptype`, with positive masses and tags 1 or 30. Coordinates follow the existing centered ASCII convention. The reader rejects nonfinite states, unsupported tags, and a requested twin sector without tagged massive twins. The default no-twin selector remains available.

The first admitted two-sector configuration uses a uniform periodic grid, CPU FFTW, CIC particles, and a nonexpanding box. AMR, cosmological evolution, GPU execution, grid growth, remapping, restart, and hydro/stellar/sink/PBH evolution are rejected for the new two-sector path. The restrictions bound the initial verification stage. The final hydro objective still includes stellar/AGN feedback, metals, and cooling/heating, without dust or ray tracing.

Optional `bimond_output_fields` adds the sector densities, phantom density, differential and twin potentials, and twin forces to the normal HDF5 dump schedule. The generic RAMSES Newtonian potential-energy diagnostic is not a conserved BIMOND energy. Validation uses the signed force and momentum relations.

The first compile attempt, Grammar job `1273847`, failed because a new local variable declaration was inside an unrelated SNRT preprocessor block. The declaration was moved outside the block. The corrected frozen source is `build/bimond_two_sector_20261107_v2` under the RAR root, based on `a733cd16c31b0026aa2a8c3542efa791c40882f2` and patch SHA-256 `82470dbddf31074db6cfaf7721ed90fe96aa5f2448a70cad15d35e3734c2e069`. Build `1273852` and dependent force verification `1273854` were submitted to Grammar. Results must be recorded after completion.

### Completed CPU validation

Build `1273852` completed in 2 min 8 s. Trial `1273854` rejected an empty twin inventory because the test generator supplied an incorrectly indexed `initfile` entry. RAMSES shifts the input list from slot 1 to `levelmin`. The corrected trial `1273856` completed in 15 s and verified all retained checksums.

| Calculation | Measured result |
| --- | --- |
| Identical deposited sector densities, `32^3` | Differential density and phantom density cancel exactly. Both sector forces recover the same Newtonian field. |
| Displaced distributions, `32^3` and `64^3` | Both force fields match the independent discrete Poisson/phantom reference. The largest common-force-normalized RMS error is `5.02389e-14`. |
| Off-center CIC deposition | Sector density RMS errors are at most `1.20441e-16`. Both sector masses agree with the input inventory within `1e-12`. |
| Static particle force interpolation | Both species receive the independently interpolated field with zero reported RMS difference. |
| Two MPI ranks, `32^3` | Force errors against the same independent reference are below `2.65e-15`. Particle tags and masses are preserved in the output. |
| Four particle steps, `16^3` | Both half-kicks and drift agree with the reconstructed update to approximately `1.12e-16`. Coincident ordinary/twin particle states remain identical. |
| Shifted distributions, signed force balance | Relative residual decreases from `9.39729e-5` at `32^3` to `7.01050e-6` at `64^3`. |

The MPI, evolution, and shifted `64^3` results come from Grammar job `1273902`. The unshifted symmetric source has a signed-force residual near roundoff at both resolutions, so a strict decrease in the unshifted residual is not a useful convergence test. The translated source supplies the nontrivial resolution comparison. The four-step displaced-particle run has a maximum signed-momentum residual of `1.60311e-3`, normalized by the sum of particle momentum magnitudes. The result does not establish long-time or timestep convergence.

The evolution audit identified an input-interface limitation. `HYDRO_PARAMS` is not read when `hydro=.false.`, so the initial runs used the default `courant_factor=0.5` despite a test entry of `0.05`. The reported kick/drift checks use the actual saved timesteps. The revised source accepts `courant_factor` in `MOND_PARAMS`, checks `0<courant_factor<=1`, and reports the effective value. The revised source also excludes SIDM baryon drag, tracer evolution, and SNRT builds from the initial twin path.

The latest frozen source is `build/bimond_two_sector_20261107_v3`, with patch SHA-256 `d1949d90ec4c47944a73910ec21065f7660751ad66322b4f60751fa0b48032df` on the same base commit. Grammar build `1273945` completed in 2 min 9 s. Evolution job `1273953` completed with the logged coefficient `0.05`; the four-step displaced-particle signed-momentum residual remains below `1.5470e-5`. The smaller coefficient was tested over a shorter physical interval, so the two runs do not form a fixed-time convergence comparison. Regression `1273954` reproduces exact no-twin/QUMOND force identity and an independent-reference RMS error of `7.41184e-14`. All jobs exited `0:0` and all run manifests pass checksum verification.

Compact records are stored under `reports/bimond_two_sector_1273856`, `reports/bimond_evolve_1273902`, `reports/bimond_courant_1273953`, and `reports/bimond_limit_1273954` in the RAR root. Effective inputs, logs, source patches, and snapshots remain in corresponding directories beneath `/gpfs/kjhan/RAR_QUMOND/runs/`.

Example gravity selection for an otherwise admitted uniform CPU setup follows. A complete input also requires the explicit particle ICs, fixed mesh capacity, noncosmological periodic boundaries, and an audited output schedule.

```fortran
&MOND_PARAMS
  use_bimond = .true.
  bimond_twin_matter = .true.
  mond_type = 1
  a0_mond = 1.2d-8
  mond_mu_type = 1
  courant_factor = 0.05
  bimond_output_fields = .true.
/
```

Restart, GPU parity, evolving MPI migration, AMR boundaries, high-acceleration limits for both interpolation functions, and coupled hydro remain unvalidated. The startup restrictions remain in place. The verified CPU module is a nonrelativistic two-sector gravity implementation and does not supply a relativistic BIMOND background or initial transfer functions.

## Required equation-level acceptance conditions

- With `rho_hat = 0`, the ordinary potential and force reproduce `mond_type=1` QUMOND on a matched mesh. Require normalized RMS force difference below `1e-12` and maximum normalized force difference below `1e-10`.
- With equal ordinary and twin density fields, `Phi_*` and the shared phantom source vanish to a normalized RMS residual below `1e-12`; each sector recovers its own Newtonian source for `beta=1`.
- In a high-acceleration field with `|grad(Phi_*)|/a0 >= 1e6`, each sector recovers its own Newtonian force to within `1e-5` fractional error for both supported interpolation functions.
- A displaced ordinary/twin source pair reproduces the asymmetric branch equations. Ordinary matter can attract twin matter while twin matter repels ordinary matter. The symmetric-branch mutual-repulsion statement does not apply. Check the conserved momentum `P = integral(rho v - rho_hat v_hat)` and the closed-system force residual `||F - F_hat||/(||F|| + ||F_hat||) <= 1e-2` at the target resolution. Require the residual to decrease under a factor-of-two mesh refinement.
- The signed-source periodic zero mode is treated consistently. Ordinary and twin mean densities belong to the background and must not be silently mixed into the peculiar-field source.
- Ordinary and twin particle counts and deposited masses match input inventories to relative error below `1e-12`.
- On a matched uniform periodic mesh, CPU and GPU phantom sources and both sector force fields agree to normalized RMS error below `1e-6`. Single-rank and multi-rank runs meet the same force-field tolerance across MPI boundaries.
- A checkpoint/restart preserves the exact `PTYPE_TWIN` particle count and produces force fields within normalized RMS error `1e-6` of an uninterrupted run.
- Add an AMR coarse-fine boundary check for both sector potentials and forces. The current QUMOND source stencil skips cells without same-level face neighbors. Preserve the documented limitation; do not claim production AMR accuracy until the BIMOND shared-source operator has a validated composite-boundary treatment.

Use reproducible normalizations in the acceptance report. For cell fields on the same active-leaf mask, define `E_RMS(A,B) = sqrt(sum(V_cell*|A-B|^2)/sum(V_cell*|B|^2))` and `E_max(A,B) = max(|A-B|)/max(|B|)`. For particle accelerations, replace `V_cell` with particle mass. Evaluate zero-reference residuals by dividing by the RMS norm of the matched nonzero Newtonian reference field. Define the closed-system force residual exactly as `||F - F_hat||/(||F|| + ||F_hat||)`. Record the mask, volume/mass weights, and reference field with every reported error.

## Explicitly unresolved inputs

- Twin-matter abundance, species composition, and density/velocity initial conditions.
- A cosmological background derived from the selected relativistic BIMOND action.
- A stable perturbation theory for the selected `alpha + beta = 0` action.
- A matching relativistic transfer solution in lagCAMB and a validated transition to the nonlinear quasi-static RAMSES equations.
- Twin-sector IC generation and any relativistic transfer prescription for the selected branch.
- A production AMR solution for the shared phantom source across coarse-fine boundaries.

Until the unresolved relativistic inputs are supplied by a validated action, the code target is the nonrelativistic BIMOND field solver. Cosmological runs using an imposed background remain phenomenological. The fixed-source RAR comparison can proceed independently, with the present no-twin BIMOND limit numerically equal to QUMOND.
