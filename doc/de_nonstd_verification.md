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

## 6. Dilaton — NOT THE ADVERTISED MODEL; mechanical bugs FIXED

The "dilaton" is a parameter-renamed symmetron clone (its own
comment admits it), NOT the Brax+2010/11 environmentally-damped
dilaton (no exponential potential, no A₂ coupling, spurious phase
transition at a0_dilaton). The same three symmetron bugs were fixed
(seeding, ρ double-count, force factor −2β² → dimensionally correct
−6Ω_mβ²(L_dilaton/L_box)²a²/a0³), and the model-identity caveat is
now stated in the code. A true Brax+ dilaton would be a separate
implementation task.

## 7. Cubic Galileon — NOT Barreira+13; marked EXPERIMENTAL

CONFIRMED: the coefficient functions are fabricated (no background
Galileon evolution, no b1/b2 functions, no Poisson back-reaction);
pre-fix linear G_eff/G−1 = 12c3²E²/c2² grew INTO the past (13× at
z=0, ~10⁵× at z=49) instead of decaying to GR. The mechanical
family bugs were fixed (ρ−ρ̄ source, single 1/β_G coupling → −0.5
force factor, a⁴ Vainshtein power, c2≠0 guard), which makes the
PDE machinery self-consistent, but the model remains a simplified
nDGP-template — read_params now prints an EXPERIMENTAL warning.
Implementing the true tracker cubic Galileon (Barreira+13 b1/b2 +
back-reaction) is future work.

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
