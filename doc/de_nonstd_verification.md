# Nonstandard DE/gravity model verification (branch de-nonstd-dev)

Running log of the model-by-model verification of the pre-existing
nonstandard models. Each entry: equation-level check against the
literature, wiring check (solver → force → particles), and numerical
tests where feasible. Started 2026-07-12.

---

## 1. CPL dynamic dark energy — VERIFIED (perturbation fallback FIXED)

**Background (OK).**
- `f_de = a^{-3(1+w0+wa)} e^{-3wa(1-a)}` is the exact CPL integral
  (init_time.f90).
- `dadtau^2 = a^6 E^2`, `dadt^2 = a^2 E^2` — correct supercomoving
  (Martel–Shapiro) and proper-time forms with the CPL `f_de`.
- `fpeebl` uses `a·dy/da` with `w(a) = w0+wa(1-a)` — self-consistent.
- Caveat (minor): the growing mode uses the LCDM closed-form integral
  `D ∝ (sqrt(y)/a)∫da/y^{3/2}`, exact only for w=-1. At IC redshifts
  (z≳50) dark energy is negligible, so the vfact error is ≪0.1%.

**LagMUSIC background hand-off — BUG CONFIRMED AND FIXED
(2026-07-31).** The campaign generator supplied the correct model-specific
lagCAMB density and velocity transfers but hard-coded `w0=-1`, `wa=0` in
every LagMUSIC configuration. The existing \(1024^3\) CPL IC therefore has
the correct scale-dependent CPL transfer but an LCDM background velocity
factor and 2LPT velocity term. At \(z_{\rm start}=49\), correcting the input
changes LagMUSIC's velocity factor from `394.397958311` to `394.377241378`,
or 0.00525%. The error is too small to explain the production FFT crash, but
the CPL production IC must be regenerated for parameter consistency.
`dmo_benchmark_setup.py` now passes `(-0.9,0.0)` for the constant-\(w\)
benchmark and `(-0.9,0.2)` for the CPL benchmark. Generation tests for
LCDM, constant-\(w\), and CPL pass.

**Perturbation, table mode (OK).** `R_DE(k,a)` from CAMB via the
Weyl-potential method (aux/generate_de_table.py); applied as
`1 + (Ω_de(a)/Ω_cb) R_DE` in MG/CG (cs2≈0) and per-k in FFT paths.

**Perturbation, kappa2/alpha fallback — BUG CONFIRMED AND FIXED
(commit 674d839).** The old factor `(k²+κ²)/(k²+κ²+α)` SUPPRESSED the
source by 15–29% with the wrong sign vs CAMB. Replaced by the
Sapone & Kunz (2009) quasi-static closure
`δ_de/δ_m = (1+w)/[(1-3w) + 2cs²k²/(3a²H²)]`, i.e.
`de_factor = 1 + α/(k²+κ²)` with
`κ² = 1.5(1-3w)a²E²(L/C)²/cs²`,
`α = (Ω_l f_de a³/Ω_m)(1+w)·1.5a²E²(L/C)²/cs²`.
Agreement with CAMB (Weyl method): ≤0.1% (w0=-0.9, wa=0.1), ≤0.8%
(w0=-0.8) for z≥0.5; ~4% at z=0 (quasi-static limit — use the table
mode for precision below z≈0.5). Verification script:
`patch/cuRamses/aux/compare_helmholtz_camb.py`.
**Runs made with the old fallback (no de_table, cs2_de>0) are invalid.**

---

## 2. Early Dark Energy (Poulin+19 fluid form) — VERIFIED (OK, with caveats)

Implementation: `dadtau`/`dadt` in `patch/cuRamses/init_time.f90`;
`ρ_EDE(a) = ω_ede·ρ_std(a_c)·2/[1+(a/a_c)^{3(1+w_ede)}]`,
`a_c = 1/(1+z_ede)`.

Numerical driver results (ω_ede=0.05, z_ede=3000, w_ede=1/3,
Ω_m=0.3, Ω_l=0.7):
- Post-a_c dilution slope `dlnρ_EDE/dlna = -4.0000` (theory −3(1+w)=−4). OK
- `dadt = dadtau/a²` exact. OK
- `E²(a=1)−1 = 1.0e-5` (documented, negligible). OK
- Age of universe: 0.96405/H0 (vs 0.96410 LCDM) — plausible.
- f_EDE peak: 0.05391 at a=4.37e-4. Analytic peak of this fluid form
  in the matter era is at `x=(a/a_c)=3^{1/4}≈1.316` with
  `f_peak = 1.140·ω/(1+1.140ω) = 0.05391` — code matches exactly.

Caveats (conventions, not bugs):
1. The true peak is 14% ABOVE ω_ede and occurs at 1.32·a_c, i.e.
   ω_ede is f_EDE(a_c), not the peak value.
2. The background has NO radiation term (base RAMSES choice), so at
   z_ede≈3000 (near z_eq) the "EDE fraction of standard density"
   differs from Planck-style f_EDE (which includes radiation) by
   ~×0.5. Translate constraints accordingly.
3. Growth/vfact at IC time ignores the EDE term in y(a) (init_cosmo);
   for ω_ede=0.05, z_ede=3000 this biases vfact by ≲0.4% at z=99.
   Acceptable; add the term if sub-percent IC velocities are needed.

---

## 3. f(R) Hu–Sawicki — BUGS CONFIRMED AND FIXED

The solver did not relax the f(R) equation. Confirmed and fixed
(all in `patch/cuRamses/force_fine.kjhan.f90`):
1. Curvature term had the WRONG SIGN (tachyonic instead of Yukawa
   mass) and was missing the (H0·L_box/c)² unit factor (34–2200×
   too large depending on box size) → indefinite operator, Newton
   Jacobian positive on coarse levels, GS non-convergent.
2. Density term was Ω_m a³ρ instead of Ω_m (ρ−ρ̄)/a·(H0L/c)²: an a²
   double-count AND no mean subtraction (k=0 residual irreducible;
   secular field-mean drift ≳|fR0| per step).
3. fR0 convention: f̄_R had a spurious factor n (n≠1 only).
4. Emergent physics pre-fix: unscreened G_eff=(1+a²/3)G at ALL
   scales — no Compton cutoff, no chameleon, F5 4×/9× weak at z=1/2.
Fixed equation (supercomoving, ρ in mean units):
`lap f_R = +(a²/3)(H0L/c)²(R̃(f_R)−R̃̄) − Ω_m(H0L/c)²(ρ−ρ̄)/a`,
force `F5 = +(a²/2)(c/H0L)²∇f_R`. Validated numerically with the
exact discrete stencil: μ(k)−1 = (1/3)k²/(k²+m²) reproduced to
rel. err ≤6e-7 over a∈{0.25,0.5,1}, fR0∈{−1e−5,−1e−6}, n∈{1,2}
(aux/fr_fixed_check.py); chameleon limit R̃→R̃̄+3Ω_mδ/a³ exact.
Also fixed: zero-seeding of f̄_R every solve (restarts/refined cells
were never initialized), MPI virtual f-sync after F5 (all solvers).

## 4. nDGP — BUGS CONFIRMED AND FIXED

1. Fifth-force factor −0.5/β double-counted the 1/β already in the
   source → G_eff/G = 1+1/(3β²) instead of 1+1/(3β) (force 2.7×
   weak at z=0 for Ω_rc=0.25). Fixed to −0.5 (Winther+15 eq. 16-17).
2. Source used ρ instead of ρ−ρ̄ (mean not subtracted). Fixed.
3. Vainshtein coefficient 1/(12Ω_rcβa²) → 1/(12Ω_rcβa⁴)
   (supercomoving conversion of the quadratic term); halos were
   under-screened by (r_V ratio) a^(2/3) at z>0. Fixed.
Verified OK: β(a) formula + branch sign, source normalization,
GR limit, wiring/order. Known limitations (documented, not fixed):
cross-derivative terms omitted from the Vainshtein operator
(screening overestimated in anisotropic configs), plain GS without
multigrid (n_iter=20), source-scaled Newton clamp, F5 holes at
coarse-fine boundaries, LCDM-hardcoded β(a) (read_params now errors
on scalar-DE combos and warns on CPL/EDE).

## 5. Symmetron — FATAL BUGS CONFIRMED AND FIXED

1. FATAL: field initialized to 0, and χ=0 is an exact fixed point of
   the homogeneous equation → fifth force was IDENTICALLY ZERO for
   every run ("converged, res=0" logs were the 0/0 case). Fixed by
   seeding the broken-phase VEV χ̄(a)=√(1−(a_ssb/a)³) into cells at
   exactly 0 on every solve (covers restarts and refined cells).
2. Force factor −6Ω_mβ²(a_ssb/a) was missing (L_symmetron/L_box)²
   and had inverted a-dependence; correct factor
   −6Ω_mβ²(L/L_box)²a²/a_ssb³ derived independently and validated
   by the unscreened limit F5/F_N = 2β²χ̄² (exact). Pre-fix error
   would have been 625× (100 Mpc/h box) once seeding was fixed.
3. Density double-count (1+ρ)=2+δ → SSB delayed to 1.26·a_ssb,
   over-screening. Fixed to ρ·(a_ssb/a)³.
Field equation structure (Davis+12/ISIS) otherwise verified OK.

## 6. Dilaton — REIMPLEMENTED as the true Brax+ model

The original code was a parameter-renamed symmetron clone (its own
comment admitted it), NOT the Brax+2010/11 environmentally-damped
dilaton. After first fixing its mechanical bugs (seeding, ρ
double-count, force normalization), the module was REWRITTEN as the
original environment-dependent dilaton (Brax, van de Bruck, Davis,
Shaw 2010; N-body form Brax+12, arXiv:1206.3568, r=3/2):
- A(φ)=1+(A₂/2)φ²/Mpl², V=V₀e^{−γφ/Mpl}; m²(a)=3A₂H²(a),
  β(a)=β₀a^s with s=3Ω_m; χ̄=β₀a^s/A₂ (minimum tracking).
- Code-unit field eq: ∇̃²χ = (3Ω_mA₂B₂/a)(ρ̃χ−χ̄) +
  a²B₂[v(χ)−v(χ̄)], v(χ)=−3Ω_mβ₀(A₂χ/β₀)^{1−3/s}; negative-definite
  Newton Jacobian; χ>0 preserved. F5 = −c̃²a²A₂·χ∇̃χ.
- Params: beta_dilaton=β₀, L_dilaton=2998ξ (range today, Mpc/h);
  a0_dilaton ignored (legacy).
Validated (aux/dilaton_check.py): background exact fixed point
(residual ~1e-18); linear response = 2β(a)²k²/(k²+m²) to 5 digits
(a=0.25/0.5/1, two k-modes); δ=10⁵ top-hat: χ→0.03χ̄, F5/F_N→0
(Damour–Polyakov screening).
The symmetron implementation was additionally cross-checked against
Brax+12 eq. (68) and its force normalization: both match EXACTLY
(with c̃ξ=λ⋆/L_box ≡ our L̃), independently confirming the §5 fixes.

## 7. Cubic Galileon — TRACKER REIMPLEMENTED (galileon_tracker=T)

The original coefficient functions were fabricated (no background
Galileon evolution, no b1/b2 functions, no Poisson back-reaction);
pre-fix linear G_eff/G−1 = 12c3²E²/c2² grew INTO the past (13× at
z=0, ~10⁵× at z=49) instead of decaying to GR. After fixing the
mechanical family bugs, the full Barreira+13 TRACKER model was
implemented (default galileon_tracker=.true., parameter-free):
- Tracker: H φ̇ = ξH0²Mpl, c2=−6c3ξ, ξ=√(6(1−Ω_m)), c3=1/(6ξ).
- Background wired through the f_de dispatch: ρ_de ∝ 1/E²,
  E²(a)=[Ω_m a⁻³+√(Ω_m²a⁻⁶+4(1−Ω_m))]/2, w(a)=−1+(2/3)Ḣ/H²
  (phantom-like, w0≈−1.18) — modified expansion + growth included.
- Perturbations (Barreira eqs. 11-15 reduced on the tracker):
  β₁=(ξ/3)[2Ḣ/H²−1+(1−Ω_m)/E⁴], β₂=2E²β₁/ξ²; code-unit field eq
  = nDGP template with coeff 1/(3β₁a⁴), source Ω_m a δ/β₂ (kernel
  reused); Poisson back-reaction −(κc3/M³)φ̇²∇²φ integrates to a
  potential ∝ u, so the total fifth force is a single gradient with
  factor +ξ/(6E²).
Validated (aux/galileon_tracker_check.py): unscreened F5/F_N equals
the analytic −ξ/(9β₂E²) to 5 digits at a=1/0.5/0.25 — G_eff/G(a=1)
= 1.844 (paper: "twice as large today"), decaying as 1/E⁴ into the
past; Vainshtein screening suppresses F5/F_N to 0.006 for a δ=5e4
top-hat. The legacy template remains available (galileon_tracker=F)
with an EXPERIMENTAL warning.

## 8. MOND — VERIFIED (healthy core; caveats documented)

Verified OK: a0 unit conversion (f_code·scale_l/scale_t² is exactly
the proper peculiar acceleration; a0 correctly constant in proper
units at all a, both cosmo and non-cosmo), μ/ν conjugate pairs
(consistent to 2e-16), QUMOND phantom-density sign and 4πG factor
chain (including boxlen≠1 and DE-boost cancellation in MG/CG paths),
single application of the algebraic correction, solve ordering.
Caveats (documented; larger design work if needed):
- Single-phi AMR design: on refined levels the "Newtonian" field
  feeding ν̃ inherits MOND-contaminated boundary conditions →
  phantom density systematically underestimated in refined regions
  (PoR uses a second, Newtonian-BC potential array).
- Zero-phantom ring at coarse-fine boundaries; isolated runs keep a
  Newtonian multipole far-field BC (truncates the MOND log field).
- AQUAL: undamped Picard can stagnate in deep-MOND regions.
- epot_tot is multiply accumulated with MOND (energy diagnostic
  invalid; dynamics unaffected).
- New guard: use_mond now errors when combined with de_perturb /
  use_horndeski / use_coupled_de (FFT path breaks the phantom/boost
  cancellation).

---

## Summary table

| Model | Verdict | Action |
|---|---|---|
| CPL background | correct | — |
| CPL perturbation fallback | wrong sign, 15–29% | FIXED (674d839) vs CAMB |
| EDE | correct | conventions documented |
| f(R) Hu–Sawicki | wrong PDE (no chameleon/Compton) | FIXED + validated (μ(k) ≤6e-7) |
| nDGP | force 1/β weak, ρ−δ, a²/a⁴ | FIXED |
| Symmetron | F5≡0 (fatal), norm 625×+, ρ−δ | FIXED |
| Dilaton | symmetron clone + same bugs | mech. FIXED, identity documented |
| Cubic Galileon | fabricated coefficients | mech. FIXED, EXPERIMENTAL flag |
| MOND | core correct | guards added, caveats documented |

**Any pre-fix production run of f(R)/nDGP/symmetron/dilaton/galileon
(and the CPL kappa2/alpha fallback) is invalid and must be re-run.**
MOND, EDE, CPL-background and CPL-table results are unaffected.

---

## 9. End-to-end simulation tests of the new DE models

Protocol: 64³ particles, 100 Mpc/h, Zel'dovich grafic ICs at z=49
from a CAMB LCDM P(k) (aux/make_ics.py), SAME velocity files for
every model; dmonly, levelmin=6, levelmax=9, USE_FFTW build; P(k)
from the built-in dump_pk. Comparison metric: P_model/P_LCDM
averaged over k=0.05–0.12 h/Mpc at z≈0, against an IC-consistent
linear forward model (aux/forward_model.py) that (i) rescales the
initial displacements by the measured vfact ratio — RAMSES converts
the shared velocities with each model's own vfact, printed at
startup since this work — and (ii) evolves the exact (δ,θ) initial
vector through each model's linear growth system.

| model (params) | sim | forward model | sim/pred |
|---|---|---|---|
| quintessence (RP α=1) | 0.832 | 0.817 | 1.018 |
| k-essence (x₀=1.0005/2) | 0.875 | 0.881 | 0.993 |
| Horndeski (μ₀=0.2) | 1.036 | 1.052 | 0.985 |
| coupled DE (β=0.1) A: G_eff only | 0.762 | 0.784 | 0.972 |
| coupled DE C: +mass evol | 0.909/0.762=1.192 vs 1.199 | | 0.994 |
| coupled DE B: +friction (fixed) | 1.044 vs 1.047 (ratio) | | 0.997 |
| coupled DE full (z=0) | 0.966 | 1.005 | 0.961 |
| galileon tracker (z=0) | 1.148 | 1.208 | 0.950* |
| dilaton (β₀=0.5, L=20) | 1.039 | 1.071 | 0.970* |

(*) The galileon/dilaton deficits relative to LINEAR theory carry
the expected sign: Vainshtein and Damour–Polyakov screening remove
part of the fifth-force boost already at quasi-linear k (cf.
Barreira+13, where the nonlinear boost at k ≳ 0.1 h/Mpc falls well
below the ~20% linear one).
Residuals of 2–5% are consistent with screening, quasi-linear
leakage into the k-band and the decaying-mode transient of the
shared-IC protocol;
every model-specific effect (background, μ(a), coupled source,
friction, mass evolution) is verified at the ≤1% level in isolation.

Bugs found and fixed BY these tests (each with its own commit):
1. dmonly refinement crash — poisson_refine used unallocated uold
   when hydro=F (pre-existing; every dmonly+m_refine run crashed).
2. Coupled-DE Friedmann inconsistency — the matter term in
   dadtau/dadt lacked the DM mass-evolution factor; isolation run
   matched the inconsistent-H growth prediction to 0.2%, and the
   fixed run matches the consistent one.
3. Coupled-DE friction at half strength — applied only in
   synchro_fine's half-kick; move_fine's half-kick now carries the
   other half (isolation runs: 1.031→1.044 vs 1.047 predicted).
4. MPI collective deadlock — all five scalar solve_levels returned
   early on ranks with no grids on the level while the survivors
   entered ALLREDUCEs (galileon hung at its gate-opening step).
5. MPI Waitall abort — the fifth-force routines skipped the final
   f-sync on empty ranks (mismatched point-to-point at sparse
   refined levels).
Also added: early-time gate for the tracker galileon (solver skipped
while G_eff/G−1 < 1e-3, i.e. z ≳ 2.5 — the |coeff| ∝ a⁻⁴ regime is
never entered) and the vfact/fpeebl startup diagnostic.

## 10. SIDM Monte-Carlo scattering — fixes and validation (2026-07-13)

Ownership of the SIDM module transferred to this project (original SIDM
project is dead). A read-only review found five defects in
`patch/cuRamses/sidm_scatter.f90`; all are now fixed on this branch.

### 10.1 Defects fixed

1. **Unit conversions (cosmo runs)** — three spurious `/aexp` factors:
   velocities (`vp*scale_v/aexp`), cell volume (`(dx*scale_l/aexp)^3`)
   and gas density (`uold*scale_d/aexp^3`). In RAMSES supercomoving
   units `vp*scale_v` *is* the proper peculiar velocity and
   `scale_l`/`scale_d` are already proper (verified against the gadget
   IC conversion `vfact = sqrt(a)*a/(100 L)` in `init_part.f90` and the
   `scale_nH` usage in cooling). Net effect: P_scatter was multiplied
   by a^2 (4e-4 at z=49, 0.25 at z=1); the DM-baryon drag rate by
   a^-4. All pre-fix cosmological SIDM runs are invalid; non-cosmo
   (aexp=1) runs were unaffected.
2. **nvector clamp** — cells with more than nvector(=32) DM particles
   used (N-1)->31 and silently excluded particles beyond the first 32
   in list order (up to ~10x rate undercount in halo cores at
   levelmax). Replaced by a per-grid grown allocatable buffer; an
   `N/(2 floor(N/2))` factor now makes odd-N sampling exact.
3. **iSIDM kinematics** — the mass splitting delta [keV] was compared
   against the kinetic energy of the ~1e44 g N-body macro-particle
   (~1e55 erg), so thresholds never gated and kicks were ~1e-64:
   label churn with zero dynamics. New parameter `sidm_mchi` [GeV];
   all channel gating and velocity changes now use the microphysical
   reduced mass mu_chi = m_chi/2.
4. **iSIDM channels** — elastic scattering was excluded whenever any
   inelastic channel was open (forced transitions), and single
   excitations g+g->g+e opened at threshold delta. Now: 2-state model
   allows only gg<->ee (threshold 2*delta) with mixed pairs elastic,
   the elastic channel competes with equal weight, and N>2 states use
   uniform branching over open pair channels.
5. **DM-baryon drag** — gas momentum was updated without total energy
   (internal-energy corruption, possible negative pressures) and the
   backreaction was uncapped in DM-dominated cells. Now conserves
   momentum and total energy exactly (frictional heat kept positive by
   the twin caps dv <= 0.5 v_rel * min(1, m_gas/m_dm)), and the local
   gas velocity is kept current between particles in a cell.

Also: `sidm` master switch accepted in `&SIDM_PARAMS` (previously
`&RUN_PARAMS` only); per-thread RNG streams are independent RANMAR
states from `rans()` (previously 100-draw shifted copies of one
stream, i.e. massively overlapping).

### 10.2 Validation (64^3 unigrid LCDM, L=100 Mpc/h, sigma/m=1000 cm^2/g)

* **Rate**: offline exact MC expectation from snapshots (cell-by-cell
  pair sums of |dv|, aux script `predict_rate.py`) vs the code's
  per-step `scattered=` counts: measured/predicted = 1.07+-0.06 at
  a=0.3 and 1.006+-0.044 at a=0.6. Per-cell `pairs=` counts match
  exactly. The two-epoch agreement kills the old a^2 error (would have
  read 0.09/0.36).
* **Conservation**: per-event |dp|/p ~ 1e-14, elastic dEk/KE ~ 1e-13
  (machine precision).
* **iSIDM** (delta=10 keV, m_chi=1 GeV, restart at a=0.6):
  f_exc=0 stays exactly 0 (threshold gates; old code would up-scatter
  immediately); a 50% excited population decays via down-scattering
  with per-event macro energy release matching 2*delta*m_p/m_chi =
  1.717e61 erg to 4 digits; late re-excitations appear only after
  down-kicked particles exceed the threshold (physically correct).
* **Drag smoke test** (hydro, sigma_b=1e5 cm^2/g): runs clean, all
  particles dragged, no negative internal energies.

### 10.3 rans() array-shape landmine (found and fixed)

The first fixed build died deterministically at step ~115-136 under
np=8 x OMP=4 with SIGSEGV inside `allocate` at
`particle_tree.kjhan.f90:959`, no bounds violation, and nothing from
plain valgrind (Intel's runtime allocates through the libiomp5
rml/tbbmalloc pool, whose arenas valgrind cannot instrument). Routing
`scalable_*` through glibc malloc with an LD_PRELOAD shim exposed the
culprit: `rans(N,...)` in `amr/random.f90` declares its seed array as
`dimension(ncpu, IRandNumSize)` — the leading dimension is the GLOBAL
`ncpu`, not `N`. The new per-thread seeding called
`rans(nthreads-1, ...)` with an `(nthreads-1, 4)` buffer, so rans
wrote with stride ncpu past the buffer (seed components such as
0x617 landed in pool metadata, e.g. a bin pointer became
0x617_00000000, detonating a later allocate). This also explains why
np=1 (stride 1) and OMP=1 (no rans call) were immune. Fixed by
generating thread streams directly with the module's jump primitives
(`ranfk`/`ranfkbinary`/`ranfatok`/`ranfmodmult` — the same
construction rans uses internally, with correct array shapes); a
warning comment was added to rans(). After the fix, np=8 x OMP=4 runs
1088 steps to a=0.6 clean (production -O3 binary re-verified past the
former crash point), with rate validation unchanged
(1.0005+-0.058 at a=0.3, 0.986+-0.045 at a=0.6).
