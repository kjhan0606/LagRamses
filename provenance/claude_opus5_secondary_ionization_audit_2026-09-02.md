The plan file tool is disabled in this session, so I'm delivering the audit inline.

# Stage-2 independent audit — FS2010 secondary ionization

**Verdict: CONDITIONAL PASS**

I re-ran all five named gates. Every one passes and reproduces its documented number exactly (`continuity_delta=0.000335457`, `host_error=1.11022e-16`, `electron_closure=3.47804e-12`, `delta_mean_xhii=2.24653e-08`, `A_B_L1=6.48652e-09`). The reported numbers are honest.

## Independently verified as correct

Measured from the vendored tables and source, not taken from the docs:

- **Column mapping.** Header is `Photon Energy | f_ion | f_heat | f_exc | n_Lya | n_{ion,HI} | n_{ion,HeI} | n_{ion,HeII} | Shull heating`; `secondary.py:98-103` takes columns 1,2,3,5,6,7. Correct.
- **Grid/file pairing.** `secondary.py:16-50` reproduces the 21cmFAST `x_int_XHII` array and file order, including the non-obvious `2.318e-4`/`4.677e-4` entries whose values don't equal their `-3.6`/`-3.3` filename labels.
- **Renormalization is defensible, not a fudge.** Over all 14×258 rows, `max |f_ion+f_heat+f_exc − 1| = 3.9e-5` (worst: `log_xi_-4.0`, 339.493 eV). The all-channel renormalization at `secondary.py:227-234` is a ≤4e-5 correction.
- **The ionization split is the right decomposition.** `Σ n_i E_i^th / (f_ion·E) ∈ [1.000140, 1.000404]` over every row with `f_ion·E > 1e-6`. Counts and `f_ion` are mutually consistent to ≤0.04%, so weighting by `n_i·E_i^th` (`secondary.py:215-225`) reproduces per-species fractions to that accuracy, and the round trip back to counts (`multiphysics.py:132`) recovers `n_i`.
- **No pathological rows.** Zero rows with `f_ion>0` and all counts zero; zero with `f_ion=0` and counts>0; no negatives. The `safe_ionization_weight` fallback is never exercised.
- **The 100 eV discontinuity is genuinely gone, with no new one.** Measured 99.9→100.1 eV max channel delta 3.35e-4. The 10 eV floor is *exactly* continuous — all 14 files have `f_ion=0, f_heat=1, f_exc=0` at E=10, matching the sub-floor branch. The high-E clamp is sound; the table is flat there (f_ion 0.387618→0.387748 over 9034→9937 eV).
- **H II is the correct coordinate.** Headers give `fHI` per file, so the tables are parameterized by `1−x_HII`. The old code passed an electron fraction into an `x_HII`-parameterized fit — this is a real physics correction, and the composition caveat is disclosed in four places.
- **Single-charging holds** in both families (`multiphysics.py:217-228`+`265-294`; `conservative_primordial.py:205-214`+`296-320`), with ledger residuals algebraically exact at the fixed point.
- **The B2 composition change is correct.** `n_He/n_H = 0.079` is the primordial value for Y≈0.24 (`Y/(4(1−Y))=0.0789`); the tables are only defined for primordial H/He, so the old H-only fixture was outside the table domain. The band moves are physically consistent (adding He roughly doubles absorbed photons at 200 eV) and centred on the measurement.

## Findings, severity ordered

### B1 — BLOCKER: FS2010 runs on a zero-helium cell, violating the He ledger

`tests/p2_p3_validation.py:25` sets `n_helium=zeros`; `:44-52` builds a 1 keV multiphysics step where `use_secondary_ionization` takes its **default `True`** (`multiphysics.py:305`). I reproduced it:

```
x_HII      = 0.09066587686538696   <- the number P2_MULTIPHYSICS.md now quotes
sec_HeI    = 0.0058762431          <- helium secondaries in a gas with no helium
sec_HeII   = 4.204407e-05
HeI ledger = 0.0058762431          <- 100% of the spurious channel, unconserved
He I photoionization rate = 8.42e+35 s^-1
backward-Euler denominator product = 8.43e+67   (float32 max 3.40e+38)
```

`minimum_n_he = max(1e-12·n_helium, tiny)` (`multiphysics.py:215`) collapses to `tiny`, so the rate at `:221-224` becomes ~1e36 and the denominator in `implicit.py:97-102` overflows to `inf`. The solve survives *only* because `inf` lands in the denominator, not the numerator.

This is the same defect class the B2 fixture change repaired — `B2_PRODUCTION_SOLVER_VALIDATION.md` explicitly cites preserved "zero-He underflow failures" — fixed in one fixture, left in another, with the unrepaired one's output quoted in a canonical doc. The existing assertion at `:53` passes either way. Acceptance item 6 forbids leaving this unrepaired and undeclared.

**Remedy.** Route He secondary fractions into heat when helium is negligible (preserving energy closure), or reject `n_He == 0` at build time; set the P2 fixture to `n_He/n_H = 0.079`; add a He-ledger assertion; regenerate the quoted number.

### C1 — No in-solver photoelectron-energy ledger on the multiphysics path
`conservative_primordial.py:45,196-202` exports `photoelectron_energy_ledger_residual`. `ThermochemicalStepResult` (`multiphysics.py:31-58`) has no equivalent. The unit test reconstructs closure externally (`tests/secondary_furlanetto_stoever.py:242-251`) only because that fixture is single-group with uniform 200 eV excess. **P5 runs on multiphysics** with a multi-group AGN ledger and has no runtime photoelectron-energy closure diagnostic — `thermal_energy_closure_relative_error` covers the thermal solver, not the photoelectron split. Items 2 and 3 ask for this in both families.

### C2 — B2's Solver A/B cross-check lost chemistry independence
`conservative_hydrogen.py` now calls the same `coupled_photo_collisional_hhe_update` as multiphysics, replacing the distinct `hydrogen_photoionization_relaxation`. Evidence: A/B mean |ΔxHII| fell `1.15422e-6 → 6.48652e-9`, and `maximum_fixed_point_residual` is now **bit-identical** (`2.396106719970703e-05`) in both blocks of the B2 JSON. The doc still calls Solver B "the independent H-only conservative solver" and presents the error drop without explanation.

### C3 — "Independent interpolation" is independent of implementation, not of rule
`host_interpolate` (`tests/secondary_furlanetto_stoever.py:79-134`) is a line-for-line NumPy transcription of `_bilinear_interpolate` plus the identical clip/split/renormalize sequence. The 1.11e-16 agreement tests JAX-vs-NumPy execution; it cannot detect a wrong *rule*. Nothing checks against 21cmFAST's actual `interp_*` routines or a published FS2010 value, yet `secondary.py:164-165` and two docs claim the interpolation "match[es] the authors' public 21cmFAST implementation". I could not verify that claim offline.

### C4 — The nearest-root repair is sound in principle but its policy is untested
The mechanism checks out: bracket endpoints are guaranteed (`residual(0) = −n_e0 ≤ 0` searching down; `residual(n_max) ≥ 0` searching up, since all fractions are clipped), the sign conventions at `implicit.py:225-226` match the bisection update at `:236-237`, and the `stationary` shortcut at `:176` is mathematically exact — `n_e = 0` really is the solution for a hot, exactly neutral, unilluminated cell.

What is untested is **root selection itself**. `tests/coupled_photo_collisional_hhe.py` checks a neutral cell stays neutral (`:35-48`) and electron self-consistency to 3.48e-12 (`:74`) — but self-consistency holds at *every* root, including the spurious one, so it cannot discriminate. Also untested: the 65-point quartic sampling (`:189-198`) resolves only ~6% of the interval at its far end and can step over a *pair* of roots; and if `crossed` is all-False, `argmax` returns 0 and the routine silently returns the incoming density (`:213-224`) with no diagnostic. Separately, "nearest root" biases toward the incoming state, so for `dt ≫ t_ion` it under-ionizes — a defensible conservative choice that belongs in the docs as a stage-5 timestep caveat.

### C5 — The P5 control is near-vacuous by construction and the reason is unstated
The pair ran at `<x_HII> = 0.978`, `<T> = 5.58e6 K`. Measured from the tables at ~199 eV: `f_ion = 0.3756` at `x_HII=1e-4` but **`0.0022` at `0.9` and `0.0008` at `0.99`**. FS2010 is suppressed ~2 orders of magnitude in this regime *by construction* — which is why `Δ<x_HII> = 2.25e-8`. The docs say the effect "is small for this particular hot, UV-dominated control" and "not evidence that secondary ionization is generally negligible" — honest, but they never give the mechanism a reader needs to judge the claim.

The bands are also ~3.5 orders of magnitude wider than the measurement (`0 < Δ<x_HII> < 1e-4` vs 2.25e-8; `−100 K < ΔT ≤ 0` vs −0.0075 K). They do **not** hide a failure — both signs are physically correct and the strict `> 0` bounds would catch a dead path — but they have almost no regression power.

### C6 — Excitation energy is a silent sink
`fractions.excitation` is exported as a rate and never returned to the gas thermal budget in either family. Defensible (Lyα escape) and unchanged from the old code, but it is a claim boundary stated nowhere, and it directly contributes to the sign of the P5 ΔT.

### Advisory (later gates)
- **A1** — measured cost: B2 runtimes 3.83→68.23 s, 3.90→59.52 s, 5.25→79.07 s (~15–18×). Sharper risk is memory: `electron_residual` runs on a `(65, *grid)` array (`implicit.py:199-204`), so every temporary is 65× the field — fine at 32³, ~4.4 GB per temporary at 256³ float32.
- **A2** — P5 provenance (`validate_p5_secondary_ionization.py:216-230`) omits `primordial.py`, `primordial_cooling.py`, `transport.py`, `dust.py`; covered only because `b2_multiphysics_artifact.py` binds all of `snrt_core/*.py`.
- **A3** — `13.60/24.59/54.42` duplicated as literals in `multiphysics.py:132,137,141` and `conservative_primordial.py:174,178,184,197-199` while `secondary.py:12-14` defines them unused.
- **A4** — `use_secondary_ionization` defaults `True` in multiphysics/thermochemistry, `False` in conservative_primordial; the `True` default is what makes B1 reachable.
- **A5** — no assertion that the 10 eV floor is continuous; it's a property of the data, and the 99.9/100.1 eV test cannot see a floor discontinuity.
- **A6** — linear (not log) interpolation in `x_HII` across a grid with a 0.1→0.5 gap; upstream-consistent, worth listing as a limitation.

## Acceptance checklist

| # | Item | Status |
|---|------|--------|
| 1 | Replace 100 eV SvS branch with FS2010 tables + continuous interpolation | **closed** |
| 2 | Continuity, table limits, independent interpolation, finite/nonneg, energy closure | **conditional** — C1, C3, A5 |
| 3 | Carry all five channels through both families, ledgers, names, runners, sharding | **closed** |
| 4 | Pin provenance, license, commit, hashes; B2/P5 fail closed | **conditional** — A2, C3 |
| 5 | Matched OFF/ON P5 effect control passing all gates | **conditional** — C5 |
| 6 | Production-reachable numerical issues repaired with a justified test or declared | **open** — B1, C4 |

**Overall: CONDITIONAL PASS.** The core stage-2 physics is sound and independently confirmed — tables, column mapping, split, renormalization, removal of the 100 eV discontinuity, H II coordinate, single-charging, provenance chain. `BLOCK` would be disproportionate: nothing I found invalidates the shipped B2 or P5 numbers, which use real helium. `PASS` is unavailable because item 6 has a live, reproduced, production-reachable violation (B1), and item 2's energy-closure gate is absent on the very solver P5 runs (C1).

Stage-3+ items (RSLA/refinement, AGN nine-group, P5 spatial/timestep convergence, coupled H+He front, dust closure) are not implicated — this implementation does not make their current state unsafe.
