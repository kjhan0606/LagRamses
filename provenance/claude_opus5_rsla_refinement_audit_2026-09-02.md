The Write tool is disabled in this session, so I can't create the plan file. The audit itself is read-only and complete — here it is.

---

# Independent stage-3 audit — RSLA and refinement

Read-only on the dirty worktree at HEAD `ca90a39`. Every quantity below was re-derived from the physics rather than taken from the artifacts.

## Reproduction (all four run, all pass)

```
RSLA_REFINEMENT_ARTIFACT_OK ratios=0.845542,0.955902,0.991949,1.00199 production=0.991892
                            error_envelope=0.0125762 mesh=0.00425058 angular=0.00131851
P1_VALIDATION_OK spatial_spread=0.0096 rs_ratios=0.5144022,0.5114129,0.5094979
                 analytic_ratios=0.9351,0.9296,0.9261
B2_MULTIPHYSICS_ARTIFACT_OK  radius=0.955685 A_B_L1=1.52855e-05 shadow=0.00439426
P5_SECONDARY_IONIZATION_ARTIFACT_OK  delta_mean_xhii=2.24653e-08
```

## Independently confirmed

- **Analytic reference** `R(t)=R_S(1−e^{−t/t_rec})^{1/3}`, `t_rec=1/(α_B n_H)`, `R_S=(3Q/4πα_B n_H²)^{1/3}` — reproduced to 7e-8.
- **Fixed-duration construction** — all **9** step counts reproduced exactly (39/116/385/1153/385/231/134/769/444), `validate_rsla_refinement.py:103-114`.
- **Source normalization** — `quadrature.py:35` normalizes weights to 1, `transport.py:126,146` adds `dt·emissivity` per direction ⇒ `emitted = Q·t` exactly.
- **Radius estimator** `(3Σx_HII·V_cell/4π)^{1/3}` is the *exactly* right choice: it equals the analytic ODE's `N_HII/n_H`. Verified `N_HII/N_analytic = r³` to machine precision on all 9 runs.
- **Photon-storage attribution is real, not asserted.** Measured storage vs the free-streaming estimate `R/(ĉt)` gives ratio **0.88–0.93 across all 9 runs over a 26× range in storage**. This independently validates the report's central physical claim.
- RSLA `3e-03` reproduces the historical B2 `solver_b` radius ratio to **2.3e-7** — the "old 4.41%" is provably the same number.
- Solver A/B: mean |Δx_HII| = 4.35e-6 (0.017% of mean x_HII), Δr = 5.7e-5.
- Hash closure: four documented SHA256s match on disk; the artifact tests recompute validator + 32 `snrt_core` + cross-artifact hashes. B2 regeneration changed **only** runtimes/`git_head`/`ionization_front` hash — every physics value bit-identical.
- P1 duration is exactly 22,265.6 yr, grid-independent, and numerically identical to the old config (`0.03c×12/cell ≡ 0.01c×4/cell`); analytic ratio 0.5501305 reproduced exactly.

## Findings

**F1 (MAJOR) — P1's residual deficit is attributed to a helium correction that is exactly zero.** `P1_CONVERGENCE.md:24-26` and the surviving sentence in `ionization_front.py:55-57` claim helium supplies "the expected correction". At P1's **20 eV** group, `σ_HeI = 0.0000e+00` (Verner threshold 24.59 eV, `primordial.py:34`); `photon_coupling.py:36-38` therefore gives helium **zero opacity**, and with no collisional channel `x_HeII ≡ x_HeIII ≡ 0` gives it **zero electrons**. The correction is identically zero, not conservative. This is the same docstring M4 flagged for a false claim — one false claim was removed and a second left in, then promoted in `P1_CONVERGENCE.md` to explain the residual ~7%.

**F2 (MAJOR) — the envelope's RSLA term uses `0.03c`, which retains ~0.5% of its own deficit.** Storage scales as `K/ĉ` (`f·ĉ/c` = 4.36e-4 → 5.05e-4) and `r` is linear in `f` near 0; fits on the two finest pairs give `r(ĉ→∞)` = 1.00813 / 1.00707 (quadratic: 1.00695).

| | as declared | with `r(∞)≈1.0075` |
|---|---:|---:|
| RSLA gap of production `0.01c` | 0.010081 | **0.015591** |
| Linear envelope | 0.012576 (63% of gate) | **0.018086** (90% of gate) |

Still passes 2%, but "deliberately conservative" is only half earned — the linear sum is conservative, the reference is not. Corollary (supports the report's own hedge): `r(∞)≈1.0075` puts converged numerics ~0.75% *above* analytic, consistent with the partial-ionization skin (numerical recombination ∝ Σx² vs analytic ∝ Σx; measured 0.899× analytic at `0.03c`).

**F3 (MINOR) — the tightest gate is the one missing from the documented contract.** `mesh_refinement_does_not_worsen_analytic_error` (`:522-525`) sits at **85% of its 0.005 allowance** (measured 0.0042506), the tightest of 19 criteria (next: 63%), and is absent from the seven-bullet contract at `RSLA_REFINEMENT_VALIDATION.md:15-25`, with no derivation for 0.005. Mitigating: mtimes (validator 05:27:14, JSON 05:37:02, runtimes summing to 9m07s) are consistent with one run started right after the validator was finalized.

**F4 (MINOR) — stale HEAD in the rebound B2 report.** `B2_PRODUCTION_SOLVER_VALIDATION.md:110-112` still says HEAD `438aff4…`; the regenerated artifact records `ca90a391…`. The JSON SHA line was updated in the same edit, the HEAD line was not.

**Informational** (no action): Solver B is independent in *chemistry closure only* — `multiphysics.py:33,459` calls the same `transport.advance_with_absorption`; S4→S8 also changes Δt (`Σ|n|` 1.5→1.732), so that term bounds a combined angular+temporal change; P1 and the RSLA validator use different radius estimators; `p1_validation.py:24` omits `n_H` from `t_rec` (correct only because `n_H = 1`); the documented `--output` commands hit the by-design `FileExistsError` (repo-wide convention, pre-existing); the float32 `α_B` gives the analytic radius ~7e-8 error.

Later hydro/source-variability and coupled-species limits are correctly declared out of scope and are **not** counted against this gate.

## M4 closure table

| # | M4 action | Status | Evidence |
|---|---|---|---|
| 1 | Correct the false front-speed docstring | **Closed with defect** | Front-speed claim removed; a second false claim (helium correction) survives — **F1** |
| 2 | Add analytic radius assertion to P1 | **Closed** | `p1_validation.py:32` `abs(ratio−1) < 0.10`; passes at 0.9351/0.9296/0.9261 |
| 3 | Run `ĉ/c = {1e-3,3e-3,1e-2,3e-2}`, show convergence | **Closed** | Monotone 0.845542→1.001993; step counts and storage mechanism independently reproduced |
| 4 | Mesh/angular refinement for the 4.41% | **Closed** | Two same-`ĉ` families; correctly labelled sensitivity, not Richardson |
| 5 | Defensible production `0.01c` error range | **Closed with caveat** | All three constructions (0.81% direct / 1.26% declared / 1.81% extrapolated) pass 2%; declared value understates the RSLA term — **F2** |

# Verdict: **CONDITIONAL PASS**

The numerics are correct, reproducible, and independently confirmed — the RSLA matrix, the mesh/angular wiring, the P1 correction, and the hash closure all hold, and the production `0.01c` bound passes 2% under every defensible reference I could construct. It is not a PASS because two stage-3 findings remain open: a false physical attribution in the very docstring M4 targeted (**F1**), and an error bound labelled "conservative" whose leading term is referenced to a non-converged `0.03c` value (**F2**), plus two minor governance/provenance items (**F3**, **F4**). All four are documentation-only — no solver, validator, or JSON change is needed, so remediation requires no artifact regeneration or hash churn.
