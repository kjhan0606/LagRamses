# Stage-1 helium recombination audit

**Model:** `claude-opus-5`
**Date:** 2026-09-01 · read-only · no shell available (see §5 on hashes)

## 1. He III case-B relation — PASS

`primordial.py:344-347` implements `alpha_HeIII,B(T) = 2 * alpha_HII,B(T/4)` literally, on top of the correct HG97 H fit (`primordial.py:316-320`: `2.753e-14 λ^1.5/(1+(λ/2.74)^0.407)^2.242`, `λ = 315614/T`). Hand-verified: `alpha_HeIII,B(1e4) = 1.5446e-12` vs recorded `1.5447607e-12`.

Minor: the α path uses `4 × 157807 = 631228 K` implicitly while `primordial_cooling.py:120` uses HG97's tabulated `631515 K` — a 0.045% inconsistency. Both are defensible; not blocking.

## 2. He II radiative case-B + dielectronic — PASS

- `primordial.py:330-334`: `1.26e-14 λ_He^0.75`, `λ_He = 2×285335/T` — HG97 case-B ✓
- `primordial.py:337-341`: `1.9e-3 T^-1.5 exp(-4.7e5/T)[1+0.3exp(-9.4e4/T)]` cm³ s⁻¹ — the standard Aldrovandi & Péquignot term HG97 adopts, correct units ✓
- **No omission / no double counting:** `helium_ii_dielectronic_recombination` has exactly one call site (`primordial.py:359`), inside `case_b_helium_recombination`, and no production caller invokes `hui_gnedin_case_b_helium_ii_radiative` directly. Every caller goes through `case_b_helium_recombination`.
- Semantics are right: dielectronic recombination has no case-A/case-B split, so adding it unchanged to the case-B radiative rate is correct; and He III correctly gets **no** dielectronic channel (needs ≥2 bound electrons) — the docstring at `primordial.py:351` says exactly this.
- Independent cross-check on the matched HG97 constant pair: `ζ/ξ = 1.24e-13/1.9e-3 = 6.53e-11 erg = 40.7 eV`, the expected autoionizing-level energy.

All eight recorded coefficients hand-verified to ~4 s.f. (e.g. at 1e5 K: radiative `4.6526e-14` + dielectronic `6.105e-13` = `6.570e-13` ✓ artifact line 15).

## 3. Temperature-resolved tests — PASS, weakly

Covers `1e4, 2e4, 4e4, 1e5 K` ✓ (`tests/helium_recombination.py:31`). Two substantive issues:

**(a) Two assertions are exactly the tautology you asked about.** `tests/helium_recombination.py:109-114` compares `alpha_heiii` against `2 * hui_gnedin_case_b_hydrogen(T/4)` — which *is* the body of `hui_gnedin_case_b_helium_iii` (`primordial.py:347`). `:115-120` compares the function against itself. Neither can fail.

**(b) The one-zone gate is temperature-degenerate by construction.** `tests/helium_recombination.py:71` sets `dt = 3/(alpha·n_e·512)`, so `alpha·dt = 3/512` at *every* temperature. The artifact confirms it: `helium_ii_final_fraction` is `0.05022486181495307` four times (`data/helium_case_b_recombination_validation.json:41-52`). The reported error `0.008793316446238375` is exactly the analytic backward-Euler truncation error `(1+3/512)^-512/e^-3 - 1 = 0.0087933` — so the 2% threshold has ~2.3× headroom over a number fixed by the scheme, carrying no physics information, and no refinement study makes it a convergence test.

**It is not vacuous, though.** Because `dt` comes from `primordial.py` and the α inside `helium_photoionization_backward_euler` (`implicit.py:88`) comes from the solver, a case-A/case-B swap in `implicit.py` would show up as ~84% error at 1e4 K. So it is a real cross-module consistency gate at four temperatures — just not four independent accuracy measurements. `:124` (`case_b < case_a`) is a genuine ordering check. The `REFERENCE_ALPHA_*` arrays (`:32-49`) are self-generated regression locks rather than independent references, but I verified all eight by hand, so they are correct.

## 4. Caller wiring — **FAIL (one thermal caller)**

Chemistry/transport callers are all wired to the shared coefficients, with no stale case-A anywhere in executable production source:

| Caller | Line | Status |
|---|---:|---|
| `implicit.py` (3 sites) | `:88`, `:144`, `:227` | ✓ |
| `multiphysics.py` | `:198-199`, `:256`, `:259` | ✓ |
| `conservative_primordial.py` | `:117-118`, `:260`, `:263` | ✓ |
| `photon_coupling.py` | `:70-71`, `:73-74` | ✓ |
| `primordial.py` explicit network | `:376-377` | ✓ |
| `hui_gnedin_case_a_helium_ii_radiative` | `primordial.py:323` | ✓ test control only (`tests/helium_recombination.py:122`) |

**BLOCKING — `simulation/snrt/snrt_core/primordial_cooling.py:134-140`: He III case-B recombination cooling is 4× too large.**

```python
helium_iii = 8.0 * 3.435e-30 * temperature * lambda_heiii**1.970 / (...)**3.720
#                              ^^^^^^^^^^^ must be temperature / 4.0
```

HG97 gives `beta_HeIII,B(T) = 8·beta_HII,B(T/4)`. Unlike `alpha_HII,B`, `beta_HII,B` carries an **explicit `T` prefactor**, so substituting `λ_HeIII` alone is not enough — the prefactor must also become `T/4`. Correct: `2.0 * 3.435e-30 * temperature * λ^1.970/(...)`.

Three independent confirmations:

- current `β/α` at 1e4 K = `6.47e-24/1.5448e-12 = 3.03 k_B T` per recombination. A Maxwellian electron has mean kinetic energy `1.5 k_B T` and recombination favours low energies, so `β/α < 1.5 k_B T` is a hard bound. **Unphysical.**
- current value is 3× the Cen (1992) **case-A** value (`2.11e-24` at 1e4 K); case B must be *below* case A.
- corrected `1.617e-24` gives `β/α = 0.758 k_B T`, matching the H II pattern (`0.664 k_B T` at 1e4 K).

Reachable in production: `thermochemistry.py:107` and `:182` → `primordial_net_rate`. `tests/b1_thermal_coupling.py:56-61` only asserts sign and a >10% separation, so B1 cannot catch it.

The rest of `primordial_cooling.py` is correct and now *consistent* with the chemistry: He II radiative cooling is case-B (`:127`), dielectronic cooling appears exactly once (`:128-133`, summed at `:184`), collisional/bremsstrahlung/Compton terms all check out against Cen (1992).

## 5. Artifacts, hashes, B2 evidence

**Hashes could not be recomputed** — no shell tool exists in this session, only file reads. What I *can* establish:

- The helium artifact's `primordial_sha256` (`1ca0b4f7…`) and `implicit_sha256` (`a540cb6d…`) are **byte-identical** to the entries for the same files in `b2_multiphysics_transport_validation.json:138` and `:131`. Two independently generated artifacts agree, so both runs observed the same sources.
- B2's `git_head` `3d5d1d61acaf…` matches current HEAD `3d5d1d6` ✓.
- Doc ↔ artifact agree on all eight coefficients and on `0.00879332` (`HELIUM_RECOMBINATION_VALIDATION.md:24-27,39`) ✓.
- The doc's SHA for the JSON (`:43`) is unverified.

**B2 did not hide a physics-payload regression.** `tools/validate_multiphysics_b2.py:115-118` sets `n_helium = 0`, `x_helium_ii = 0`, `x_helium_iii = 0`, and `:338-339` zeroes the He cross sections. Every He term is multiplied by zero, so the payload is *provably* invariant under the coefficient change — a stronger argument than the diff. It also means B2 supplies zero helium evidence, which `HELIUM_RECOMBINATION_VALIDATION.md:50-53` states accurately.

Stage-6 deferral (coupled H+He opacity/front, timestep convergence) does **not** invalidate the local implementation; I did not weigh it against the verdict.

## Minimum fixes for re-audit

1. **`primordial_cooling.py:134-140`** — change the prefactor to `temperature / 4.0` (or `2.0 * 3.435e-30 * temperature`). Rerun B1.
2. **`tests/b1_thermal_coupling.py`** — add `beta/alpha < 1.5 * k_B * T` for H II, He II, He III at 1e4/2e4/4e4/1e5 K. First-principles bound, so not self-referential.
3. **`tests/helium_recombination.py:109-120`** — replace the two tautologies with a closed-form NumPy re-evaluation written independently of `primordial.py`, or literature values at 1e-3 tolerance; note the `REFERENCE_ALPHA_*` provenance in the docstring.
4. **`tests/helium_recombination.py:71` + `HELIUM_RECOMBINATION_VALIDATION.md`** — state that the gate is a four-temperature *consistency* check (`alpha*dt = 3/512` by construction, hence four identical outputs), and add a fixed-`dt` variant so the temperatures differ in stiffness.

Non-blocking: `chemistry.py:19-23` holds a second H case-B convention (`2.59e-13 (T/1e4)^-0.7`) reachable only via `coupling.py:47` → `benchmarks.py`, off the production path; `primordial_cooling.py` duplicates the HG97 constants rather than importing them, which is what let finding 4 survive.

## Verdict: **CONDITIONAL PASS**

The helium recombination-coefficient work itself is correct and correctly wired through every chemistry and transport caller. Closure is conditional on the four fixes above — principally the 4× He III recombination-cooling error, which is production-reachable and is the thermal counterpart of the exact relation stage 1 asserts. This is a stage-1 verdict only, not a statement about the full RT stack.
