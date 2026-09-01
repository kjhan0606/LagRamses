# Stage-1 helium recombination re-audit

**Model:** `claude-opus-5`
**Date:** 2026-09-02 · read-only · no shell (tool set = Read/Grep/Glob only; see §6)

## 1. Condition 1 — He III recombination cooling — CLOSED

`primordial_cooling.py:150-152` now reads

```python
helium_iii = 8.0 * _hui_gnedin_case_b_hydrogen_recombination_cooling(temperature / 4.0)
```

against the single shared fit at `primordial_cooling.py:117-129`, which is also the H II path (`:138`). The `T/4` argument propagates to **both** the explicit prefactor and `lambda` — exactly what the 4× defect missed.

Hand-verified at 1e4 K: `lambda(T/4) = 126.2456`, `beta_HeIII,B = 8 × 2.0215e-25 = 1.6172e-24`, i.e. the `1.617e-24` predicted by the first audit. `beta/alpha = 1.6172e-24 / 1.5448e-12 = 0.758 k_B T`, matching the H II pattern (`0.664 k_B T`). The unphysical `3.03 k_B T` is gone.

Duplicate removal confirmed by search: `3.435e-30` occurs **exactly once** in the tree (`primordial_cooling.py:125`), and the hard-coded `631515 K` flagged as a minor 0.045% inconsistency no longer appears anywhere — alpha and beta now share `4 × 157807 = 631228 K`, so that minor finding closed incidentally.

## 2. Condition 2 — first-principles thermal bound in B1 — CLOSED

`b1_thermal_coupling.py:59` sets the four temperatures; `:76-83` asserts `0 < beta/alpha < 1.5 k_B T` over the three **radiative** channels; `:84-89` checks the dielectronic ratio separately against the matched pair `1.24e-13 / 1.9e-3`. The comment at `:73-75` states the reason correctly.

The exclusion is physically right: `beta_diel/alpha_diel = 6.53e-11 erg = 40.7 eV` is the autoionizing-level energy, about 31 times the bound at 1e4 K. The bound has teeth where it matters — reintroducing the 4× error drives He III to `3.03 k_B T` and fires the assertion.

*Observation, not a defect:* for He II radiative the check is structurally trivial, since `primordial_cooling.py:139-143` defines `beta_HeII,rad = alpha_HeII,rad · k_B · T`, giving `beta/alpha = 1.0 k_B T` identically. That identity is HG97's own prescription, so the physics is right; the assertion just cannot fail for that channel. All three channels are present as the condition required.

## 3. Condition 3 — tautologies removed — CLOSED

`helium_recombination.py:97-116` is a pure NumPy re-derivation. It does **not** import or call `hui_gnedin_case_b_hydrogen`; it recomputes `lambda_HeII = 2(285335)/T`, the `1.26e-14 lambda^0.75` term, the dielectronic term, and the `2 alpha_HII,B(T/4)` substitution from scratch, including its own `hydrogen_temperature = temperature / 4.0` at `:109`. The old self-comparison is gone; `:171-182` compares production against this separate evaluation.

`REFERENCE_SOURCE` (`:49`, HG97 Appendix A + DOI) is propagated to the artifact at `:266` and lands in the JSON at `helium_case_b_recombination_validation.json:95`.

Spot-verified by hand: `alpha_HeII,tot(1e4) = 2.6161e-13` ✓; `alpha_HeII,tot(1e5) = 4.6532e-14 + 6.1051e-13 = 6.5705e-13` ✓; `alpha_HeIII(1e4) = 1.5447e-12` ✓.

*Separate from the condition:* the arrays are 17-digit locks against the published **formulas**, not published **tabulated values** — they catch transcription and wiring drift, not a wrong choice of fit. The doc says exactly this (`HELIUM_RECOMBINATION_VALIDATION.md:36-39`).

## 4. Condition 4 — labelled gate + fixed-time run — CLOSED

- Three-`t_rec` run labelled `"cross-module consistency at a fixed dimensionless alpha*dt"` (`helium_recombination.py:238`, JSON `:37`); doc explains the degeneracy at `HELIUM_RECOMBINATION_VALIDATION.md:41-45`.
- `decay_over_fixed_time` (`:119-153`) uses the single literal `FIXED_ELAPSED_S = 2.0e12` (`:30`) with `dt = elapsed_s / 512` (`:130`) — one scalar for all four temperatures.
- Non-degeneracy asserted at `:222-223`; artifact spread is 0.559 (He II), 0.581 (He III).
- 2% limit asserted at `:224-225`, predeclared at `:260`.

I reproduced the solver semantics independently: with zero photoionization, `implicit.py:91-104` collapses to `x/(1 + alpha n_e dt)` per step for both ions, so `exp(-alpha n_e t)` is the correct reference. The truncation error `N(x - ln(1+x))` then predicts:

| quantity | predicted | recorded |
|---|---:|---:|
| He II max rel. err | 0.0016872 | `0.00168485` |
| He III max rel. err | 0.0093271 | `0.00932733` |

About five-significant-figure agreement confirms the recorded values are genuine outputs of this scheme. Per-step `alpha*dt` now spans 3.7e-4 to 6.0e-3 across the eight (ion, T) combinations — a real temperature-resolved measurement, not a relabelling.

## 5. Wiring — no regression

Every helium consumer still routes through `case_b_helium_recombination`: `implicit.py:88,144,227`; `multiphysics.py:199`; `conservative_primordial.py:118`; `photon_coupling.py:71`; `primordial.py:377`. `hui_gnedin_case_a_helium_ii_radiative` has exactly one caller, `tests/helium_recombination.py:184`, as a test control. `chemistry.py:23` / `benchmarks.py:58` retain the off-production-path `2.59e-13` convention — unchanged, non-blocking.

## 6. Provenance and B2

Cross-artifact agreement (the check that needs no shell): the helium and B2 artifacts were generated independently yet record identical hashes for `primordial_cooling.py` (`2cb0187b…`, JSON `:102` vs `b2_…json:139`), `primordial.py` (`1ca0b4f7…`), and `implicit.py` (`a540cb6d…`). Both runs observed the same **post-fix** source. The `primordial.py`/`implicit.py` hashes are byte-identical to those quoted in the 2026-09-01 audit, confirming the fix was confined to `primordial_cooling.py` and the tests.

Fail-closed coverage now spans the fix: `helium_recombination_artifact.py:54-62` re-hashes five files including `primordial_cooling.py` and `b1_thermal_coupling.py`, and `:50-51` re-asserts non-degeneracy from the stored payload. `b2_multiphysics_artifact.py:79-81` re-hashes **every** `snrt_core/*.py`, so the canonical B2 artifact could not have survived the cooling edit without a rerun.

B2 remains provably helium-free: `validate_multiphysics_b2.py:115-118` zeroes `n_helium`, `x_helium_ii`, `x_helium_iii` and `:338` zeroes the He I cross section, so every He term — including the corrected He III cooling, gated by `n_helium * x_helium_iii` at `primordial_cooling.py:167,197` — is multiplied by zero. The unchanged payload is guaranteed, not merely observed.

**Unverified:** the JSON SHA256 claimed at `HELIUM_RECOMBINATION_VALIDATION.md:62` (`9623d39b…`). No test enforces it and I have no shell. Non-blocking — the machine-enforced chain runs artifact to source, the direction that matters, and it is complete.

## 7. Future improvements (explicitly not stage-1 defects)

1. Independent tabulated cross-check for `REFERENCE_ALPHA_*` (§3).
2. He II radiative bound is structurally trivial (§2).
3. Doc-level JSON hash is unenforced (§6).
4. Import HG97 constants into `primordial_cooling.py` rather than restating `2 × 157807` — now cosmetic, since the shared-helper refactor removed the duplication that let the original defect survive.
5. Stage-6 coupled H+He front/convergence — deliberately out of scope, not weighed in the verdict.

## Verdict: **PASS**

**All four original conditions are closed.** The 4× He III recombination-cooling error is fixed at the shared-function level with the prefactor correctly evaluated at `T/4`; B1 gates it with a first-principles bound that has teeth; the two self-referential assertions are replaced by a genuinely separate NumPy evaluation with named provenance; and the degenerate one-zone gate is both honestly relabelled and supplemented by a non-degenerate fixed-time run whose recorded errors I reproduced analytically to five significant figures. No new blocker. Stage-1 verdict only, not full RT production or publication readiness.
