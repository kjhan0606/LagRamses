# SNRT Radiation-Transfer Subsystem — Independent Scientific & Architecture Audit

**Repository root:** `/gpfs/kjhan/LRD_JWST` · **Audit date:** 2026-09-01 · **Mode:** read-only; working-tree state as found (dirty tree assessed, not git HEAD)

**Model identifier:** `claude-opus-5`, as declared to me by the runtime environment block. Runtime identity was not independently attested; I report the exposed string, not a verified attestation.

**Evidence standard.** Every finding is tagged **[OBSERVED]** (read directly at the cited line, or from a recorded artifact in the repo) or **[INFERRED]** (my own analysis/arithmetic from observed inputs). I executed no code, so every claim about *runtime behaviour* not read from a recorded artifact is marked as inference and paired with a confirming test. `provenance/rt_architecture_audit.md` and `provenance/production_publication_readiness_plan.md` were read only as claims; none of their conclusions is assumed.

---

## 1. Executive summary and gate verdict

### Verdict: **BLOCK** for publication of any radiation-driven ionization, temperature, or dust observable. **APPROVE WITH CONDITIONS** for continued development of the transport prototype.

The transport operator is one of the better-constructed parts of this repository: the upwind discrete-ordinates scheme is genuinely in conservative flux form, the local source/absorption sub-step is solved *analytically* rather than explicitly, the positivity condition is correctly identified and correctly enforced by construction in every driver, and the boundary-escape term in the photon ledger is algebraically identical to the discrete flux the operator actually applies. These are real, verifiable strengths (§8).

The blocking problems are not in the transport kernel. They are in **the coupling between radiation and gas thermodynamics**, and in **which solver the science pilots actually run**:

1. **The thermal sector is not coupled to the RT ionization state at all.** The cooling term is a 4-D table indexed by `(a, n_H, Z, T)` with no ionization argument (`thermochemistry.py:91`), while heating comes from a non-equilibrium H/He network. Heating and cooling are computed under mutually inconsistent ionization assumptions. This is, in my assessment, the direct cause of the recorded `mean T = 1.03835e7 K` (`P5_THERMOCHEMISTRY_VALIDATION.md:42`), and of the otherwise inexplicable fact that this mean is identical to six significant figures across two runs whose mean ionization differs by 60%.

2. **The validated solver and the production solver are different code paths.** The only genuine analytic benchmark (`radius_ratio = 0.9559` against the time-dependent Strömgren solution) validates `conservative_hydrogen.py`/`conservative_primordial.py`. The P4/P5 science pilots run `multiphysics.py`/`thermochemistry.py` — a different solver, with dust and X-ray secondaries, which **has no test anywhere that exercises it with more than a single direction or against any reference solution**.

3. **The primary science observable is set by numerical knobs, not physics.** From the repository's own ledger: halving the timestep moves mean `x_HII` from `0.106431` to `0.065732`; changing the source-cell limiter moves it to `0.0223`; switching deposition mode moves it to `0.00884`. The repo correctly declines to promote this, but that gate must stay closed.

Two further defects are material for the LRD/AGN case specifically: the He III case-B recombination coefficient is **37–42% too small** (M1), and the dominant hard AGN group's photoelectron energy (55.27 eV) falls *below* the hard 100 eV activation threshold of the secondary-ionization closure, so the largest photoelectron-energy channel deposits **100% heat and zero secondary ionizations** (M3).

Finally, a caution on how conservation is reported: **the photon-ledger residual is an algebraic identity of the discretization, not an independent test.** The `1.74e-16`/`3.89e-16` "photon ledger error" values presented as PASS gates (`P4_TRANSPORT_CONSERVATION_VALIDATION.md:49-50`) measure float64 round-off. They are correct and worth keeping as regression guards, but they certify arithmetic, not physics (M5).

---

## 2. Architecture as implemented (not as documented)

There are **two independent RT+chemistry solvers**, with sharply different validation status. Confirmed by exhaustive call-site enumeration.

| | **Solver A — "multiphysics"** | **Solver B — "conservative"** |
|---|---|---|
| Core | `multiphysics.py:248`, `thermochemistry.py:170` | `conservative_primordial.py:58`, `conservative_hydrogen.py:29` |
| Opacity closure | old-state opacity + optional time-average | C2-Ray-style time-averaged H opacity, fixed point |
| Chemistry | photon transfer + cap, then split recombination | analytic H relaxation + backward-Euler He |
| Dust | **yes** | **no** |
| X-ray secondaries | **yes** (on by default) | optional, **off** in recorded runs |
| Driven by | `p4_run_transport_pilot.py:180`, `p5_run_thermochemical_pilot.py:216` | `p6_run_conservative_thermochemical_pilot.py:171`, three `validate_conservative_*` tools |
| Analytic benchmark | **none** | `radius_ratio = 0.9559` vs analytic Strömgren |
| Timestep convergence | **fails** (38% shift on halving dt) | passes (`x_hii` L1 coarse→fine `2.79e-4`) |
| Transport-coupled test | **none** | yes (32³, S4, 116 steps) |

**[OBSERVED]** Every invocation of Solver A outside the pilots uses a degenerate angular grid: `tests/p2_p3_validation.py:26` sets `directions = jnp.zeros((1, 3))`, making `_transport_rhs` (`transport.py:56-63`) identically zero — these are 0-D local-operator tests with transport switched off. The single exception is `tests/dust_opacity.py:59`, a 2-cell, one-direction, CFL = 1 case where upwind advection is exact; it validates dust momentum/energy wiring but is not a physics benchmark.

---

## 3. Mathematical algorithm — assessment

**3.1 Transport discretization — sound** [OBSERVED + INFERRED]. `transport.py:49-63` implements first-order upwind discrete ordinates, selecting the backward difference when `v_a ≥ 0` and forward otherwise, with `rhs = -ĉ v_a·gradient`. Because `v_a` is constant per ordinate, this is algebraically identical to the conservative flux form `F_{i+1/2} = ĉ v_a I_i`, so the scheme **conserves photon number exactly up to boundary flux** — a property of the discretization, not an assumption. Boundaries (`transport.py:32-46`) pad with zeros; I verified both signs — at the high-`x` face with `v_x < 0`, `gradient = (0-I)/dx` gives `rhs = ĉ v_x I/dx < 0`, i.e. the cell drains outward with no inflow. **Correct vacuum outflow on all six faces.**

**3.2 Local source/absorption sub-step — a genuine strength** [OBSERVED + INFERRED]. `transport.py:102-151` solves the local problem *analytically*:

```
I^{n+1} = I_transported·e^{-τ} + dt·S·φ(τ),   φ(τ) = (1 - e^{-τ})/τ
absorbed = I_transported·(1 - e^{-τ}) + dt·S·(1 - φ(τ))
```

Their sum is exactly `I_transported + dt·S`. **The local photon budget closes identically**, which is why positivity requires only the directional transport CFL, exactly as `transport.py:116-117` claims. `expm1` (line 130) removes cancellation at small τ; the quadratic series fallback (lines 139-144, switched at τ < 1e-4) is a correct belt-and-braces guard, and the stated float32 overflow rationale for keeping it quadratic is sound.

**3.3 Positivity / CFL — correct criterion, enforced by construction only** [OBSERVED]. `cfl_number` (`transport.py:66-71`) computes `ĉ·dt·max_d Σ_a |n_{d,a}|/dx_a` — exactly the diagonal-dominance condition. It is the right criterion. **It is never called** (exhaustive grep: only its definition and the `__init__.py` export). Drivers instead construct `dt` to satisfy it (`p5_run_thermochemical_pilot.py:131-132` and four others), safe for the cubic cells they all use (`:130`); `build_thermochemical_step` further subdivides by `total_subcycles` (`thermochemistry.py:203-208`), which only tightens it. See **Mo7**.

**3.4 Angular quadrature — correct normalization** [OBSERVED + INFERRED]. `quadrature.py:31-35` normalizes `Σw = 1`. The source is broadcast isotropically across ordinates (`transport.py:126`), so the zeroth moment gains exactly `dt·emissivity = dt·L/V`, consistent with `sources.py:48` producing `L/V`. I checked this chain end-to-end; it is correct. Carlson S4 satisfies `Σ|μ_a| = 1.5` exactly, matching the `directional_extent` used by all drivers.

**3.5 Reduced speed of light — correct steady-state scaling** [INFERRED]. The solved system `∂I/∂t = -ĉ n·∇I - ĉκI + S` gives `I ∝ 1/ĉ` in steady state, so the absorption rate `ĉκn_γ` and flux `ĉΣwnI` are **ĉ-independent** — the correct RSLA property. `photoionization_rates` (`primordial.py:279`) correctly uses `ĉ`, and `absorbed_dust_momentum_rate` (`dust.py:175`) correctly uses the **physical** `c` for photon momentum while transport uses `ĉ` — a subtle choice that is right, and documented at `dust.py:152-153`. The *validity condition* (`v_front ≪ ĉ`) is separate and is violated in one benchmark — **M4**.

**3.6 Group treatment / SED closure — well designed** [OBSERVED]. `sed_weighted_group_closure` (`primordial.py:125-203`) computes `σ̄[g,s] = ∫N_E σ_s dE / ∫N_E dE` and the **cross-section-weighted** photoelectron excess energy — the correct pair of weightings, genuinely better than centre-energy sampling. `group_spectral_closure_from_metadata` (`:206-245`) validates shape, finiteness, positivity, and monotonicity before use. One edge defect: **Mo2**.

---

## 4. BLOCKER findings

### B1 — Cooling is blind to the RT ionization state; heating and cooling are mutually inconsistent

**BLOCKER** — [OBSERVED] code + artifact; [INFERRED] causal attribution

`thermochemistry.py:90-92`:
```python
def residual(temperature_k):
    background = net_rate(atlas, scale_factor, temperature_k, chemistry.n_hydrogen, metallicity_solar)
    return heat_capacity * temperature_k - thermal.internal_energy_density - dt * (photoheating_rate + background)
```

`net_rate` (`jax_thermal_atlas.py:87-96`) interpolates a table whose only axes are `(scale_factor, log n_H, log Z, log T)` (`:18-22`). **There is no ionization-state argument.** Meanwhile `photoheating_rate` comes from a fully non-equilibrium network tracking `x_HII`, `x_HeII`, `x_HeIII` per cell.

`P4_THERMAL_ATLAS.md:11` states the intent: *"`net_rate` is a background UVB/metal thermal term."* That intent does not survive contact with what a Grackle equilibrium table contains. Either:
- **the table includes primordial H/He cooling** — evaluated at the table's own equilibrium ionization, not the network's, *and* its UVB photoheating is added on top of the RT photoheating (double counting); or
- **the table excludes primordial cooling** — the coupled system then has **no H/He cooling channel at all**: no collisional excitation, no collisional ionization cooling, no recombination cooling, no bremsstrahlung.

**Physical impact.** The recorded result is `mean T = 1.03835e7 K` (`P5_THERMOCHEMISTRY_VALIDATION.md:42-43`) at `mean n_H = 0.1456 cm⁻³` (`P4_TRANSPORT_CONSERVATION_VALIDATION.md:18`). The decisive symptom is in the recorded table itself: runs `C0.4_N4` and `C0.1_N1` report **identical mean T to six significant figures (`1.03835e7`)**, while `C0.05_N1` reports `1.06261e7` with `x_HII` falling 38% to 0.0657. A temperature this insensitive to a 38% ionization change is the signature of a cooling function that does not know the ionization state. **[INFERRED]**

Strong supporting evidence: the *same* `_implicit_thermal_update` driven by Solver B with three soft groups (18/30/80 eV) and no secondaries yields `temperature_mean = 13760 K`, `temperature_max = 52313 K` (`data/validate_conservative_primordial_thermal.json:11-13`). Physical. The pathology appears only when hard groups with large photoelectron excess are added and the missing cooling channels can no longer be masked.

**Provenance aggravation.** `write_thermal_atlas` (`thermal_atlas.py:232-246`) writes only `format` and `format_version`. **No Grackle version, no UVB model, no chemistry flags, no data-file identity.** The generator is an external binary invoked by path (`build_grackle_atlas.py:48-52`) that is **not present in this repository** (filesystem search finds no C sources outside the vendored venv). The physical content of the dominant cooling term is therefore **unauditable from this repository**; a referee cannot determine whether the UV background is on. This contrasts unflatteringly with the dust sidecar loader (`dust.py:67-79`), which *requires* non-empty `reference_mixture`, `opacity_source`, `spectral_weighting` and cross-checks group edges. The repository knows how to do this.

**Secondary defect, same routine.** `atlas.mean_molecular_weight` is loaded (`jax_thermal_atlas.py:35`), exposed via `mean_mu` (`:99`), and **never called at runtime**; `_implicit_thermal_update` derives heat capacity from the network's `particle_number_density` (`thermochemistry.py:86`). At 10⁷ K the atlas assumes μ ≈ 0.6 while the network's `x_HII ≈ 0.02` implies μ ≈ 1.22 — so `net_rate(T)` corresponds to a different gas than the one being heated. **[INFERRED]**

**Acceptance test:**
1. Regenerate two atlases with `UVbackground = 0/1`, recording all Grackle flags into required HDF5 attributes; re-run `C0.1_N1`. A material T difference proves double counting.
2. One-zone: fix `n_H = 0.146`, `Z = 1e-6`, `T = 10⁷ K`, sweep `x_HII ∈ {0.02, 0.5, 1.0}`. Assert total net rate changes by > 10×. Currently it changes by exactly 0.
3. Add explicit non-equilibrium primordial cooling (collisional excitation + collisional ionization + recombination + free-free at the network's `x_HII`, `x_HeII`, `x_HeIII`, `n_e`); restrict the atlas to metals. Gate: `C0.1_N1` mean T below 10⁵ K, or a quantitative energy budget explaining otherwise.
4. Make `read_thermal_atlas` **reject** atlases lacking generator provenance, as `read_dust_opacity_metadata` already does.

### B2 — The production solver has no transport-coupled validation

**BLOCKER** — [OBSERVED], exhaustive call-site enumeration

The P4/P5 pilots call `build_multiphysics_radiation_step` (`p4_run_transport_pilot.py:180`) and `build_thermochemical_step` (`p5_run_thermochemical_pilot.py:216`). The complete set of non-pilot call sites:

- `tests/p2_p3_validation.py:31, 44, 56, 83` — all four use `directions = jnp.zeros((1,3))` (line 26): transport RHS identically zero, 0-D local-operator tests.
- `tests/dust_opacity.py:58` — 2 cells, 1 direction, CFL exactly 1.0.

**No test in the repository runs the production solver with a real angular quadrature against any reference — analytic, converged, or differential.** Conversely, the single genuine analytic validation (`tools/validate_conservative_hydrogen_stromgren.py:115`, artifact `data/validate_conservative_hydrogen_stromgren.json`) exercises `build_conservative_hydrogen_step`, a solver with **no dust and no secondaries**, which the pilots do not use.

Every distinctive feature of the production path — photon-conservative cap, dust opacity, X-ray secondaries, time-averaged opacity, thermal-atlas coupling — is unvalidated in the presence of spatial transport. Ray effects, front propagation, shadowing, and the cap/upwind-stencil interaction are all untested for this path.

**Acceptance test:** (1) Port the Strömgren validator to Solver A (32³/S4/`ĉ = 3e-3 c`/0.5 `t_rec`, zero dust, secondaries off, one 18 eV group); require `radius_ratio` within **5%** of `R_S(1-e^{-t/t_rec})^{1/3}`. (2) **A-vs-B differential:** in that configuration require L1 `|x_HII^A - x_HII^B| < 1e-3`. (3) Re-enable dust and secondaries one at a time, recording each shift as a documented physics delta. (4) Add a shadow test for Solver A (S8 vs 192-direction), as `p1_validation.py:34-45` already does for transport-only.

### B3 — The primary science observable is controlled by numerical parameters the repository itself declares non-physical

**BLOCKER (for publication of ionization/temperature fields)** — [OBSERVED], from the repository's own recorded ledger

From `P5_THERMOCHEMISTRY_VALIDATION.md`, mean `x_HII` over the identical 6.37 Myr problem:

| Configuration | mean `x_HII` | line |
|---|---:|---|
| `C0.1_N1` (outer CFL 0.1) | 0.106431 | :43 |
| `C0.05_N1` (outer CFL 0.05) | 0.065732 | :44 |
| `source_limit1_C0.1_N1` | 0.027458 | :81 |
| `source_limit0.25_C0.1_N1` | 0.0223196 | :100 |
| `compact3_C0.1_N1` | 0.00882267 | :125 |

Halving the timestep alone changes the answer by **38%** (`:54`). Across the numerical controls the answer spans a factor of **12**. `compact3` is explicitly *"a numerical control, not a physical source-size model"* (`sources.py:35-36`), yet it changes volume-mean ionization by a factor of 2.5.

**Credit where due:** the repository states this plainly and refuses promotion — *"The P5 state is therefore not promoted to a science result"* (`:59-60`) — and marks the `compact3` and refined-mesh rows `PASS*` with an explicit footnote that only the internal conservation gate passed (`:128-136`, `:153-162`). This is exemplary honesty and is the main reason my verdict is not harsher. But an open gate is still a closed door: **no `x_HII`, `T`, or derived quantity from this configuration may be published.**

**Root cause [INFERRED]:** (i) the source cell is optically thick by τ ≈ 3×10³ (**Mo3**), so the sub-cell I-front is unresolved and the cell average is set by the limiter/deposition choice; (ii) the photon-conservative cap (**M6**) returns un-absorbable photons to the radiation field rather than advancing the front, making the result explicitly dt-dependent.

**Acceptance test:** declare in advance a science observable and threshold (e.g. volume-mean `x_HII` stable to < 2% under simultaneous halving of dt **and** doubling of linear resolution, limiter fixed), then demonstrate it. `tests/p5_source_cell_convergence.py` and `tests/p5_refined_mesh_convergence.py` are the right harnesses; they currently check artifact pairs rather than enforcing a declared threshold.

---

## 5. MAJOR findings

### M1 — He III case-B recombination coefficient is 37–42% too small, and misattributed to Cen (1992)

**MAJOR** — [OBSERVED] code; [INFERRED] arithmetic

`primordial.py:323-329` returns the He III → He II coefficient as `4.0 * hui_gnedin_case_b_hydrogen(temperature)`. The correct hydrogenic scaling for charge Z is `α_Z(T) = Z·α_1(T/Z²)`, i.e. `α_HeIII,B(T) = 2·α_H,B(T/4)`. Evaluating the Hui & Gnedin fit at `primordial.py:316-320`:

| T (K) | code `4·α_H,B(T)` | correct `2·α_H,B(T/4)` | ratio |
|---:|---:|---:|---:|
| 2×10⁴ | 5.711×10⁻¹³ | 9.083×10⁻¹³ | **0.629** |
| 4×10⁴ | 3.024×10⁻¹³ | 5.184×10⁻¹³ | **0.583** |

The docstring at `primordial.py:5-6` claims the helium closure *"uses the standard Cen (1992) radiative plus dielectronic approximation."* That is true of the He II term (lines 326-328 match Cen 1992), but the He III coefficient is **not Cen 1992**; Cen gives `3.36e-10·T^-0.5·(T/10³)^-0.2/(1+(T/10⁶)^0.7)`. This is an undocumented ad-hoc substitution.

**Physical impact.** Too-small α means He III recombines too slowly → **He III systematically over-abundant** in AGN-irradiated gas. Not peripheral for this project: He II λ1640/λ4686 recombination lines are direct observational diagnostics for LRD/AGN spectra, and He III contributes 2 electrons per ion to `n_e` (`primordial.py:250-252`), feeding back on hydrogen recombination throughout the network.

**Acceptance test:** one-zone at fixed `n_e`, `Γ_HeII = 0`, `x_HeIII(0) = 1`, several recombination times at T ∈ {10⁴, 2×10⁴, 4×10⁴, 10⁵} K; assert the decay matches `2·α_H,B(T/4)` to < 2%. Then re-run the P6 pilot and record the `x_HeIII` shift as a physics delta.

### M2 — Case A helium mixed with case B hydrogen

**MAJOR** — [OBSERVED]

`primordial.py:326` uses Cen's **case A** He II radiative coefficient (`1.5e-10·T^-0.6353`) while hydrogen uses **case B** (`:316`). The two cases embody contradictory assumptions about recombination-photon fate. In an optically thick nebula — every configuration here, with τ_cell ≈ 10³ (**Mo3**) — case B is appropriate for both. Case A > case B, so He II→He I is over-recombined → **He I over-abundant → He I opacity overestimated**. Combined with **M1** (He III over-abundant), the helium ionization structure carries errors in *both* directions, which will not cancel.

**Acceptance test:** substitute a case-B He II coefficient (Hummer & Storey 1998 tabulated, or Verner & Ferland 1996); record `x_HeII`/`x_HeI` profiles before and after. Add a unit test asserting `α_HeII,B(T) < α_HeII,A(T)` across the tabulated range.

### M3 — The dominant hard AGN group's photoelectron energy falls below the 100 eV secondary threshold

**MAJOR** — [OBSERVED] code + production ledger; [INFERRED] energy budget

`secondary.py:36`: `active = energy >= 100.0`, with `heating = where(active, heating, 1.0)` and all ionization fractions zeroed below. The tested quantity is the *photoelectron* excess energy — physically the correct variable — but the cutoff is a hard step.

The production AGN ledger (`data/p4_pilot_agn_photon_ledger.json:40-49, 91-112`) gives, for group 3 (`[54.42, 500] eV`):

| | H I | He I | He II |
|---|---:|---:|---:|
| photoelectron excess (eV) | **55.27** | **46.59** | **15.52** |

**All below 100 eV.** Group 4 (`[500, 2000] eV`) has excesses 612.2/603.1/576.0 eV — all above. The closure is active only for group 4.

**Energy budget [INFERRED]** from the ledger's own `total_photon_rate_s`:

| group | photons s⁻¹ | H I excess (eV) | photoelectron energy (eV s⁻¹) | secondaries? |
|---|---:|---:|---:|---|
| 3 `[54.42,500]` | 7.151×10⁵³ | 55.27 | **3.95×10⁵⁵** | **no** |
| 4 `[500,2000]` | 1.609×10⁵² | 612.2 | 9.85×10⁵⁴ | yes |

**Group 3 carries 4× more photoelectron energy than group 4, and 100% of it is deposited as heat with zero secondary ionization**, solely because its cross-section-weighted excess sits 45 eV below a hard switch. Physically, a 55 eV photoelectron in weakly ionized gas does produce secondaries. `use_secondary_ionization=True` is hard-coded at `p4_run_transport_pilot.py:188` and is the P5 default (`thermochemistry.py:185`), so this is the production configuration.

**Knife-edge sensitivity.** At `x_e = 10⁻⁴`, just above threshold the SvS heating fraction is ≈ 0.111; just below it is 1.0. A marginally harder SED, or a shifted 500 eV edge, flips the dominant heating channel by **~9×**. A production configuration should not sit 45 eV from a factor-9 discontinuity. This is a plausible contributing cause of B1's 10⁷ K.

`P2_MULTIPHYSICS.md:19-21` documents the 100 eV limitation honestly, but **does not report that the production SED places its dominant channel below the threshold** — the gap between a disclosed limitation and a quantified impact.

**Acceptance test:** (1) Replace the step with the tabulated Furlanetto & Stoever interpolation already named as the plan at `P2_MULTIPHYSICS.md:21`, or at minimum a continuous low-energy extension. (2) Regression assert `|f_heat(99.9 eV) - f_heat(100.1 eV)| < 0.05` at `x_e = 10⁻⁴`; the current code fails by ~0.89. (3) Re-run P5 `C0.1_N1` with group-3 excess forced above threshold; if ΔT or Δ`x_HII` exceeds 20%, the current result is threshold-artifact dominated.

### M4 — RSLA validity is violated in the P1 benchmark, and its docstring claim is contradicted

**MAJOR** — [OBSERVED] config + table; [INFERRED] arithmetic

`ionization_front.py:53-55` claims: *"The reduced speed of light is chosen to remain above the characteristic Strömgren-front speed."* Default `reduced_light_speed_fraction = 1.0e-2` (line 45).

Reconstructing the P1 configuration from `p1_validation.py:20-22` and `ionization_front.py:38-67` (Q = 10⁴⁹ s⁻¹, n_H = 1, T = 10⁴ K, 128 pc box, S4 with `Σ|μ| = 1.5`, courant 0.2, `n_steps = 4·size`):

- integration time = **7.027×10¹¹ s = 22 268 yr**, independent of grid size — matching `P1_CONVERGENCE.md:23` ("22,265.6 yr") to 4 digits, confirming the reconstruction.
- `α_B(10⁴ K) = 2.5919×10⁻¹³`, `t_rec = 3.858×10¹² s`, so `t/t_rec = 0.1822`.
- Analytic: `R/R_S = (1 - e^{-0.1822})^{1/3} = **0.5502**`.
- Reported: **0.4369 / 0.4399 / 0.4393** (`P1_CONVERGENCE.md:27-29`).

**The code is ~20% low in radius — ~50% deficit in ionized volume — and converges to that wrong value across all three grids.** [INFERRED]

The cause is almost certainly the RSLA, not a discretization bug: light travel to the final radius (30 pc) at `ĉ = 0.01c` takes `3.1×10¹¹ s`, i.e. **44% of the total run time**; even at the end `v_front/ĉ = 0.166`. The docstring's claim is false in the early R-type phase and gives only a factor-6 margin at the end.

**Compounding:** `P1_CONVERGENCE.md:22-31` presents this as "B01 spatial convergence" with a 0.70% spread — and **never compares to the analytic value**. The assertion (`p1_validation.py:28-29`) is `assert spatial_spread < 0.02` — grid-to-grid self-consistency only. **A converged wrong answer passes.**

**Contrast — done right elsewhere.** `validate_conservative_hydrogen_stromgren.py:115` *does* compare to `R_S(1-e^{-t/t_rec})^{1/3}`, uses `ĉ = 3e-3 c` (20× front-speed margin), and records `radius_ratio = 0.9559`. That artifact also contains the diagnostic of the residual error: `photons_in_domain/emitted = 3.122e62/1.929e63 = **16.2%** still in flight`, against a volume deficit `1 - 0.9559³ = 12.7%`. **The residual 4.4% radius error is quantitatively attributable to RSLA photon storage.** [INFERRED]

**Acceptance test:** (1) Correct the docstring, or raise `ĉ` until `v_front,max/ĉ < 0.05` throughout. (2) Add the missing assertion `|R_ion/R_analytic - 1| < 0.10` to `p1_validation.py`; it will currently **fail** at 0.44/0.55 — that is the point. (3) Run the conservative validator at `ĉ/c ∈ {1e-3, 3e-3, 1e-2, 3e-2}` and show `radius_ratio → 1`, converting an unquantified approximation into a documented error bar that justifies the production `ĉ = 0.01c`.

### M5 — The photon-ledger residual is an algebraic identity, not a physics test

**MAJOR (reporting/claims)** — [OBSERVED] code; [INFERRED] proof

`ledger.py:130`: `residual = final - initial - emitted + absorbed + escaped`. I verified this closes **identically** in exact arithmetic:
- Local identity (§3.2): `I^{n+1} + absorbed = I_transported + dt·S`.
- Interior upwind fluxes telescope: `Σ_cells I_transported·V = initial - dt·(boundary outflow)`.
- `_face_escape` (`ledger.py:27-40`) evaluates exactly `ĉ·A·Σ_d w_d I_d·max(±n_d,0)`, which I confirmed is precisely the discrete upwind boundary flux the operator applies.
- Substituting: `residual ≡ 0`.

The ledger is *correctly constructed* — a genuine strength as a round-off and wiring-mismatch guard, and the deliberate provision of `photon_ledger_from_absorbed` (`ledger.py:101`) for the capped/iterated case shows real care. **But it cannot detect a wrong absorption coefficient, SED, cross section, front speed, or ĉ.**

The reported values are round-off, not validation: `P4_TRANSPORT_CONSERVATION_VALIDATION.md:49-50` reports `1.74e-16`/`3.89e-16` as PASS, with the doc noting float64 (`:36`). The `unallocated primary = 3.18e-17` column is likewise near-zero *by construction* — the photon-conservative cap (`multiphysics.py:329-365`) exists to drive it to zero, and `_advance_species` clamps again (`:165-172`). Two self-fulfilling metrics presented side by side as gates. The prior audit (`rt_architecture_audit.md:28-29`) lists these as *"useful evidence"* without the distinction.

**Remediation:** relabel as "arithmetic closure / regression guards" in all docs and gate tables; keep them. Add a genuinely independent check: compare total absorbed photons against `∫Q dt·(1 - escaped)` from an **independent** ray-traced optical depth along a few sightlines, not from the same operator.

### M6 — The photon-conservative cap trades I-front accuracy for robustness

**MAJOR** — [OBSERVED]

`multiphysics.py:329-365` computes a per-cell `scale` limiting gas absorption to available neutral H/He, multiplies `absorbed_intensity` by `dust_fraction + scale·gas_fraction`, and **returns un-absorbed photons to the radiation field** (line 365). This is positivity-preserving and drives the "unallocated photon" diagnostic to zero. But physically, a cell that would be over-ionized within one step should **advance the ionization front**, as a photon-conserving C2-Ray scheme does via the time-averaged optical depth. Returning photons to the local field makes the front position explicitly dt-dependent — exactly the pathology in **B3**.

State inconsistency in the same block: `channels` is built from `opacity_state` (line 313) but the availability cap reads `neutral_number_densities(state)` (line 330) — the un-averaged initial state. These differ when `time_averaged_absorption_iterations > 0`.

**Dead branch, same routine.** `multiphysics.py:374-395`: the `if time_averaged_absorption_iterations:` and `else:` branches contain **byte-identical** calls to `_advance_species` with the same arguments and the same `mean_hii` computation. An intended distinction was evidently lost. In the stiffest part of the solver, that is a signal worth investigating, not a cosmetic duplication.

**Acceptance test:** in the Solver-A Strömgren port (**B2**), record the fraction of cells where `scale < 1` per step; require < 10⁻³ at the production timestep, or replace the cap with the C2-Ray time-averaged-optical-depth update Solver B already implements correctly (`conservative_primordial.py:190-199`). Resolve or delete the dead branch with a comment recording intent.

### M7 — Dust energy leaves the system unaccounted; no dust temperature, re-emission, or scattering

**MAJOR (for LRD claims)** — [OBSERVED]

`thermochemistry.py:260-268` passes only `radiation.gas_heating_rate` to the thermal update; `dust_heating_rate` is accumulated as a diagnostic (`:281`) and never returned to the gas or re-emitted. No dust temperature is computed. `dust.py:152-153`, `P4_DUST_OPACITY.md:42-53`, and `P2_MULTIPHYSICS.md:15` all declare this out of scope with a correct technical reason (scattering cannot be treated as absorption). **The scope gate is honestly declared.**

The gap: **no magnitude bound is reported.** For "little red dots" — a class *defined* by a red, dust-affected continuum — the fraction of the photon energy budget silently removed by dust is the central quantity. The repo records `cumulative_dust_absorbed_photons` and `cumulative_dust_heating_energy` (`P4_DUST_OPACITY.md:60-64`) but no run reports what fraction of `L_bol` that represents.

**Acceptance test:** in the next dusty P5 run, report `Σ dust_heating_energy / Σ emitted photon energy` as a headline number. Above ~10%, no SED-adjacent claim can be made without an IR re-emission closure, and that must appear in the paper's systematics section.

### M8 — Helium opacity is evaluated at the end state while hydrogen is time-averaged

**MAJOR** — [OBSERVED]

In `conservative_primordial.py`, hydrogen uses analytic relaxation and its **time-averaged** neutral fraction for opacity (`:119`, `:193-199`, `:277-278`) — the correct C2-Ray construction, and what earns the 0.9559 Strömgren result. Helium opacity is taken from the backward-Euler **end state** (`:120-121`, `:279-282`); the docstring is candid (`:75-77`). The species are integrated at different temporal orders, so the He I-front will systematically lead the H front (end-state opacity under-estimates the step-averaged He optical depth). For AGN-driven He II→He III fronts — again directly relevant to LRD He II diagnostics — this is a first-order error at the production timestep.

**Acceptance test:** apply time averaging to He; the three-state backward-Euler solve (`implicit.py:70-106`) is algebraically correct — I verified the elimination term by term — and can be wrapped in the same mean-value factor. Gate: He I-front radius stable to < 2% on halving dt, matching the H criterion.

### M9 — The AGN ledger silently drops the 2–10 keV group defined in the P0 configuration

**MAJOR** — [OBSERVED]

`config/p0_photon_group_edges_ev.txt` defines **10 edges / 9 groups**, extending to **10 000 eV**. The production ledger (`data/p4_pilot_agn_photon_ledger.json:9-59`) defines **5 groups**, 11.2 eV to 2000 eV. Its stated rationale (`:119`) — *"Photons above 2 keV are excluded because the current P0 group layout ends at 2 keV"* — **is false with respect to the config file in the repository**, which has an edge at 10 000 eV.

**Physical impact.** For the Sazonov et al. (2004) continuum with `νL_ν ∝ ν⁰` above 1 keV (photon number `∝ ν⁻²`, per `:3`), the 2–10 keV band carries ~27% as many photons as 0.5–2 keV, each 3–4× more energetic with cross sections 4–64× smaller. These are the **most penetrating** photons — the ones that set diffuse pre-ionization and pre-heating far from the source, the physically interesting AGN signature. Omitting them truncates the hard tail while the retained 500–2000 eV group still triggers the full SvS closure (**M3**), further skewing the heating/ionization split. This also means the P0 9-group contract is unmet, consistent with the prior audit's no-go (`rt_architecture_audit.md:88`) — but the *reason recorded in the ledger is wrong*, a provenance defect in itself.

**Acceptance test:** regenerate over all 9 configured groups; add a hard check that ledger group edges equal the configured edges — the dust loader already does exactly this (`dust.py:105-108`); reuse the pattern. Report the change in `x_HII` and T at large radius.

### M10 — Failed thermal brackets snap the temperature to the table floor or ceiling

**MAJOR** — [OBSERVED]

`thermochemistry.py:145-146, 159`: when the 32-sample log-spaced search finds no sign change, `has_crossing` is False and the cell is set to `bound_temperature` — the atlas floor or ceiling. `cumulative_thermal_bound_hits` records these (`:163-166`), and `P5_THERMOCHEMISTRY_VALIDATION.md:21-22` correctly declares bound hits a hard failure, with accepted runs reporting zero. **Good instrumentation, correctly enforced gate.**

Residual risk is structural: the fallback is a silent, unphysical jump of up to 8 decades (bounds `[10, 10⁹] K`, `data/validate_conservative_primordial_thermal.json:38-41`), and the bracketing exists because the cooling curve is non-monotone (`:104-107`). With B1 unresolved — heating against an ionization-blind curve — nothing in the current tests bounds the probability of entering a regime with no local root. Zero bound hits in three runs is encouraging, not sufficient.

**Acceptance test:** make a non-zero bound hit raise, not just record. Add a one-zone stress sweep over `(n_H, Z, T₀, Γ)` spanning the production range; assert zero bound hits across the grid.

---

## 6. MODERATE findings

**Mo1 — Secondary-ionization ionized fraction uses the wrong denominator.** [OBSERVED] `multiphysics.py:106` and `conservative_primordial.py:147-150` define `electron_fraction = n_e/(n_H + n_He)`. The Shull & van Steenberg parameter is the hydrogen ionized fraction (`n_e/n_H` in standard usage). With `n_He/n_H = 0.079` this under-estimates `x` by ~7.3% uniformly; the fits are steep near `x → 0`. *Test:* switch to `n_H`; assert `f_heat` changes < 5% at `x_e = 0.1` and record the shift at `x_e = 10⁻⁴`.

**Mo2 — Closed group intervals leak sub-threshold cross sections across group edges.** [OBSERVED] `primordial.py:175` selects `(grid >= lower) & (grid <= upper)` — both edges inclusive — so a threshold on a group's upper edge contributes a sliver to the trapezoid. Three independent confirmations in the production ledger: `σ_HI = 6.54e-22` for the **sub-Lyman-limit** group `[11.2, 13.6]`; `σ_HeI = 5.26e-22` for `[13.6, 24.59]`; `σ_HeII = 9.15e-23` for `[24.59, 54.42]` (`:70, 78, 86`). Each must be **exactly zero**. Magnitude is negligible (I estimate ~0.01% of total H ionizations), but non-ionizing photons formally ionize hydrogen. *Test:* make the upper edge exclusive; assert `σ_s[g] == 0` for every group entirely below threshold `s`.

**Mo3 — The ionization front is unresolved by ~3 orders of magnitude in optical depth.** [OBSERVED + INFERRED] `P4_TRANSPORT_CONSERVATION_VALIDATION.md:16-18`: cell width `3.707×10²¹ cm` (1.2015 kpc), mean `n_H = 0.1456 cm⁻³`. At the Lyman limit with neutral gas, `τ_cell ≈ 6.3e-18 × 0.1456 × 3.707e21 ≈ 3.4×10³`. A single cell is opaque by three orders of magnitude; sub-cell I-front structure is unmodelled. The reported `mean x_HII = 0.539` (`:49`) is a cell-average where the sub-cell profile dominates. This is the physical root of **B3**: volume-mean ionization is not a converged observable at this resolution regardless of timestep. *Test:* report `τ_cell` in every run's HDF5 attributes and in the docs; gate any ionization claim on `τ_cell < 1` in the cells that carry it.

**Mo4 — Two inconsistent mean-molecular-weight definitions coexist.** See B1; `atlas.mean_mu` is never called at runtime while `thermal.particle_number_density` is.

**Mo5 — The crossing-beam benchmark is vacuous for S_N.** [OBSERVED] `crossing_beams.py:75-80` checks `|F|/(cN) → 0` at the symmetry point between two equal sources. Discrete ordinates satisfies this **by construction** — it is the diagnostic that distinguishes S_N from M1/flux-limited diffusion, and S_N cannot fail it. Useful as a contrast with moment methods; not evidence of S_N accuracy, and should not be listed among validation evidence without that caveat (cf. `rt_architecture_audit.md:28`).

**Mo6 — Explicit recombination branch lags the electron density.** [OBSERVED] `multiphysics.py:208` and `photon_coupling.py:69` use `electron_number_density(state)` — the **pre-photoionization** state — in the recombination exponent. First-order in dt, mitigated when `implicit_recombination_iterations > 0` (the pilots set 24, `p4_run_transport_pilot.py:189`), but `multiphysics.py:256` defaults it to **0**, so any new caller silently gets the lagged branch. *Test:* change the default, or require it explicitly.

**Mo7 — `cfl_number` is exported but never called.** [OBSERVED] Exhaustive grep finds it only at `transport.py:66` and `__init__.py:13, 30`. CFL is enforced by construction in each driver — correct for the cubic cells they all use — but there is no runtime assertion and no protection for a future anisotropic-cell or non-cubic AMR caller. *Test:* assert `cfl_number(config, directions) < 1` inside `build_explicit_step` and in every driver.

---

## 7. MINOR findings

- **Mi1** `TransportConfig.cell_width` is annotated `tuple[float, float, float]` (`transport.py:16`) but receives a `jnp.ndarray` from `ionization_front.py:62` and `crossing_beams.py:33`. Works (closed over as a constant) but defeats the frozen-dataclass hashability the annotation implies. [OBSERVED]
- **Mi2** `residual` is rebound from a function to an array at `thermochemistry.py:162`, shadowing the closure at `:90`. Harmless (after last use). [OBSERVED]
- **Mi3** `benchmarks.py:43-89` (`make_stromgren_problem`, the `coupling.py`-based B01) is superseded by `ionization_front.py` and referenced by no test — a dead path that could mislead a reader into thinking it is the validated benchmark. [OBSERVED]

---

## 8. Evidenced strengths

Load-bearing; preserve through remediation.

1. **Exact local source/absorption integration** (`transport.py:102-151`) with a closing photon identity, `expm1` cancellation control, and a documented float32 overflow rationale. Better than the explicit treatment most codes at this maturity use.
2. **Conservative-form upwind transport with correct vacuum boundaries** on all six faces (verified for both signs of every direction cosine).
3. **A correct RSLA implementation**, including the non-obvious and correct use of the physical `c` for dust photon momentum while transport uses `ĉ` (`dust.py:175`, documented `:152-153`).
4. **A real analytic benchmark with a recorded, good result**: `radius_ratio = 0.9559`, with a well-chosen configuration (`ĉ = 3e-3 c` → 20× front-speed margin, cell = `R_S/8`, 0.5 `t_rec`).
5. **Genuine timestep convergence for Solver B**: `x_hii` L1 coarse→fine `2.79e-4`, `x_heii` `6.53e-4`, mean `|Δlog₁₀T|` `0.0057`.
6. **An analytic implicit-solver check**: `p2_p3_validation.py:79-81` verifies the backward-Euler bisection against the exact quadratic root `(√(1+4a)-1)/(2a)`.
7. **Correct, verified algebra in the helium backward-Euler elimination** (`implicit.py:90-103`) — I derived the 3-state system independently; it matches term for term, including the `α_HeIII` back-coupling.
8. **Structural positivity of ion fractions**: I traced every path and confirmed `x_HeII + x_HeIII ≤ 1` is preserved by construction in the transfer caps (`multiphysics.py:165-193`), the explicit branch (`:220-222`), and the implicit solve (`implicit.py:152-155`) — so `n_HeI` can never go negative and produce negative opacity.
9. **Correct bisection bracketing** in `implicit_case_b_recombination_with_recombinations`: `f(0) < 0 ≤ f(n_e,initial)` always brackets the root, and reported recombination counts are *exactly* consistent with the fraction update (not an independent estimate) — which is why those ledger residuals are meaningful fixed-point diagnostics.
10. **The SED closure** (`primordial.py:125-203`): correct photon-number weighting for cross sections, cross-section weighting for photoelectron excess, full input validation, offline evaluation so arrays are static for XLA.
11. **The dust sidecar contract** (`dust.py:48-109`): mandatory provenance strings, unit/dimension checks, in-band weighted-energy check, exact group-edge agreement with the photon ledger. The model the thermal atlas should follow.
12. **Deliberate design of `photon_ledger_from_absorbed`** (`ledger.py:101-131`) to avoid reconstructing absorption when the kernel has capped or iterated it — real care about ledger/operator consistency.
13. **Documentation honesty.** `P5_THERMOCHEMISTRY_VALIDATION.md` reports its own non-convergence, marks rows `PASS*` with explicit footnotes, and declines promotion; `P4_TRANSPORT_CONSERVATION_VALIDATION.md:22-27` retains and labels invalid controls. Better practice than most production codes.

---

## 9. Publication claim classification

**VERIFIED**
- Explicit S_N transport conserves photon number to round-off up to boundary escape, with a closing finite-volume ledger. *(Algebraic identity — see M5 for correct framing.)*
- Positivity holds under the stated directional CFL, enforced by construction in every driver.
- The conservative H-only closure reproduces the time-dependent analytic Strömgren radius to **4.4%** at 32³/S4/`ĉ = 3e-3 c`, with 16.2% of photons still in flight accounting for the deficit.
- The conservative H/He and H/He+thermal closures are timestep-converged at the recorded settings.
- Backward-Euler recombination matches the exact quadratic root at `α n_e dt ≫ 1`.
- Angular convergence: S8 within 0.44% of a 192-direction product quadrature on the transport-only clump-transmission test.
- Sharded and unsharded transport agree to 10⁻⁶ across 2 devices.
- Dust absorption removes the correct photon number and energy and deposits `E/c` momentum along the correct ordinate (2-cell, CFL = 1 exact-advection test).

**CONDITIONALLY SUPPORTED**
- *"The AGN photon budget is auditably wired from sink diagnostic → rate ledger → photon ledger → deposition."* — the chain is real and the SED closure is serialized and validated, **but** the ledger omits the configured 2–10 keV group and misstates why (**M9**), and carries three spurious sub-threshold cross sections (**Mo2**).
- *"S8 is the minimum production quadrature."* — supported for the transport-only shadow test; never established for the coupled production solver (**B2**).
- *"Photon-conservative absorption resolves the unallocated-photon problem."* — it does, by construction, at the cost of dt-dependent front speed (**M5**, **M6**).
- *"Dust is wired end-to-end."* — loader, opacity, heating, momentum diagnostic are correct and provenance-gated; but no dust temperature, re-emission, or scattering, and no reported bound on energy removed (**M7**).

**UNSUPPORTED**
- Any quantitative `x_HII`, `x_HeII`, `x_HeIII`, or `T` field from the P4/P5 pilots (**B1**, **B3**).
- Any claim of correct photoheating or thermal structure (**B1**, **M3**).
- *"B01 spatial convergence"* as evidence of transport correctness — the assertion tests grid self-consistency only, and the value is ~20% from analytic (**M4**).
- Any helium ionization-structure or He II recombination-line result (**M1**, **M2**, **M8**).
- Any X-ray secondary-ionization result in the production configuration (**M3**).
- Full P0 9-group RT (**M9**).
- Any result depending on the thermal atlas's physical content, which is not auditable from this repository (**B1**).

**OUT OF SCOPE** (correctly declared, not defects)
- Dust scattering and IR re-emission; radiation-pressure/hydro coupling; H₂/Lyman-Werner chemistry (group `[11.2, 13.6]` interacts only with dust); live RAMSES/lagRamses RT→accretion→feedback loop (AGN sources are post-processed snapshots); non-equilibrium metal chemistry; AMR (the solver is static-grid by construction); implicit transport/DSA (`implicit.py` is local chemistry only, correctly noted at `rt_architecture_audit.md:61-63`).

---

## 10. Wiring diagram (as implemented)

```
 AGN sink diagnostic ──▶ p4_build_agn_rate_ledger.py ──▶ p4_build_agn_photon_ledger.py
                                                                    │
                     [SOS04 continuum, 5 groups 11.2 eV–2 keV]      │  ✗ M9: config defines
                     [σ̄, ⟨E⟩, photoelectron excess per group]       │     9 groups to 10 keV
                                                                    ▼
 config/p0_photon_group_edges_ev.txt ····(9 groups, NOT enforced)···┤
                                                                    │
 static RT input (32³, dx = 1.2 kpc, τ_cell ≈ 3e3) ─┐               │
                                                    ▼               ▼
                                        sources.deposit_point_sources
                                        (point | compact3)  ✗ B3: mode changes answer 2.5×
                                                    │
                                                    ▼  emissivity[g,x,y,z] = L/V
     ┌──────────────────────────────────────────────────────────────────────┐
     │  transport.advance_with_absorption                                   │
     │    upwind S_N (conservative flux form) + vacuum boundaries  ✓        │
     │    exact local  I·e^{-τ} + dt·S·φ(τ)  ✓                              │
     │    ĉ = 0.01 c   ✗ M4: RSLA margin unquantified in production         │
     └───────────────┬───────────────────────────────┬──────────────────────┘
                     │ absorbed_intensity            │ next_intensity
                     ▼                               ▼
   ┌─────────── SOLVER A (P4/P5 pilots) ───────────┐  ┌──── SOLVER B (P6 + validators) ────┐
   │ multiphysics.py / thermochemistry.py          │  │ conservative_primordial.py         │
   │ ✗ B2: NO transport-coupled test               │  │ ✓ analytic Strömgren 0.9559        │
   │                                               │  │ ✓ dt-converged                     │
   │ partition by κ_HI:κ_HeI:κ_HeII:κ_dust         │  │ time-averaged H opacity  ✓         │
   │   ├─▶ dust ──▶ heating + E/c momentum   ✓     │  │ end-state He opacity  ✗ M8         │
   │   │        ✗ M7: energy exits system          │  │ no dust                            │
   │   └─▶ gas ──▶ SvS secondaries                 │  │                                    │
   │            ✗ M3: 100 eV cut kills group 3     │  │                                    │
   │            ✗ Mo1: x_e denominator             │  │                                    │
   │        photon-conservative cap  ✗ M6          │  │                                    │
   │        ──▶ H/He network                       │  │ H analytic relax + He backward-    │
   │            ✗ M1: α_HeIII 37–42% low           │  │ Euler (algebra verified ✓)         │
   │            ✗ M2: case A He + case B H         │  │ same α defects M1/M2               │
   └───────────────────────┬───────────────────────┘  └──────────┬─────────────────────────┘
                           │ gas_heating_rate                    │ gas_photoheating_rate
                           ▼                                     ▼
              ┌────────────────────────────────────────────────────────────┐
              │ thermochemistry._implicit_thermal_update                   │
              │   C·T = u + dt·(photoheating + net_rate(a, n_H, Z, T))     │
              │   ✗✗ B1: net_rate has NO ionization argument               │
              │   ✗  B1: atlas provenance not recorded; generator external │
              │   ✗  M10: no-bracket → snap to table floor/ceiling         │
              └────────────────────────┬───────────────────────────────────┘
                                       ▼
                       ledger.photon_ledger_from_absorbed
                       residual ≡ 0 identically  ✗ M5: not a physics test
                                       ▼
                        HDF5 artifact + validation gate tables
```

---

## 11. Ranked remediation sequence

Ordered by *what unblocks the most downstream work per unit effort*.

1. **Close the thermal-coupling defect (B1).** Add explicit non-equilibrium primordial cooling evaluated at the network's ion fractions and `n_e`; restrict the atlas to metals; make `read_thermal_atlas` reject atlases without generator provenance (mirroring `dust.py:67-79`). *Nothing temperature-dependent can be trusted until this is done.* Confirm with the one-zone `x_HII` sweep and the UVB on/off differential.
2. **Fix the three cheap physics errors in parallel** — each is localized with an unambiguous correct value: `α_HeIII = 2·α_H,B(T/4)` (**M1**, `primordial.py:329`) plus docstring correction; case-B He II (**M2**, `:326`); continuous secondary closure across 100 eV (**M3**, `secondary.py:36`). Each gets a one-zone regression test with a closed-form or tabulated reference.
3. **Validate the production solver (B2).** Port the Strömgren validator to Solver A; add the A-vs-B differential test. The single highest-value test addition in the repository, and likely to surface further findings.
4. **Quantify and bound the RSLA (M4).** `radius_ratio` vs `ĉ/c ∈ {1e-3 … 3e-2}`; add the missing analytic assertion to `p1_validation.py:29`; correct the `ionization_front.py:53-55` docstring.
5. **Regenerate the photon ledger over all 9 configured groups (M9)**; add the ledger-vs-config edge assertion reusing the `dust.py:105-108` pattern; fix the closed-interval edge leak (**Mo2**) in the same change.
6. **Replace the photon-conservative cap with the C2-Ray time-averaged update in Solver A (M6)** — Solver B already contains the correct construction (`conservative_primordial.py:190-199`). Resolve the dead branch at `multiphysics.py:374-395`. Most likely route to closing **B3**.
7. **Establish the resolution gate (B3, Mo3).** Declare the observable and threshold *before* the next run; report `τ_cell` in every artifact; do not promote ionization results from cells with `τ_cell ≫ 1`.
8. **Time-average the helium opacity (M8)**; extend the He convergence gate to match H.
9. **Report the dust energy fraction (M7)**; decide the IR scope before any SED-adjacent claim.
10. **Housekeeping:** assert `cfl_number < 1` at runtime (**Mo7**); fix the SvS denominator (**Mo1**); use `atlas.mean_mu` or delete it (**Mo4**); relabel ledger residuals as arithmetic guards in all gate tables (**M5**); remove the dead `benchmarks.make_stromgren_problem` path (**Mi3**).

---

## 12. Gate verdict

**BLOCK** — for publication or production interpretation of any radiation-driven ionization, temperature, dust, or helium-line observable.

**APPROVE WITH CONDITIONS** — for continued development, and for publication of the *transport operator and conservative H/He closure* as a methods contribution, conditional on items 1–4 of §11 and on this claim discipline:

- Report ledger residuals as arithmetic closure, never as physics validation.
- Report `radius_ratio` against the analytic solution, with its `ĉ` dependence, wherever Strömgren validation is cited.
- State that the P4/P5 configuration is not timestep- or resolution-converged, as the repository's own P5 ledger already does.

The blocking issues are concentrated in the radiation↔thermodynamics coupling and in test coverage, not in the transport mathematics. The transport kernel, the conservative H/He closure, the SED closure, and the dust provenance contract are of publishable quality. The distance from BLOCK to APPROVE is measured in the four remediation items above, not in a rewrite.

**Model identifier:** `claude-opus-5`, as declared by the runtime environment. Runtime identity was not independently attested to me.
