The plan-file tool is disabled in this session, so I'm delivering the audit inline.

# Stage-2 independent re-audit — FS2010 secondary ionization

**Verdict: CONDITIONAL PASS**

All six reproduction commands re-run and reproduce their documented numbers exactly (`continuity_delta=0.000335457`, `host_error=1.11022e-16`, `floor_delta=0.00227216`, `reference_error=0`, `electron_closure=3.47804e-12`, `rk4_absolute_error=5.54127e-05`, `xray_xHII=0.090666`, `A_B_L1=1.52855e-05`, `delta_mean_xhii=2.24653e-08`, `delta_mean_temperature_k=-0.00752274`). Every JSON/report SHA256 I checked matches. Evidence below is measured from source and raw data, not read from the docs.

## Verified closed

**B1 — zero-helium violation.** Reproduced the original scenario in float32 (the dtype where it overflowed). At `n_He = 0`: `secHeI = secHeII = 0`, both He ledgers exactly `0`, all fields finite, root bracketed. The mask at `multiphysics.py:149-153` uses `state` (start-of-step), not `opacity_state`, so it cannot switch inside the fixed point; unavailable ionization energy is routed to heat at `:182-205` and the ledger closes (`pe_ledger/pe = 7.9e-8` in float32).

One bounded residual: the mask is *relative* (`1e-12 × n_helium`), so it detects "helium ionized away", not "helium absent". At `n_He = 1e-20` or `1e-35` the He I channel is fully active at the same magnitude as at 0.079, with `helium_i_ledger_residual = 0.00587624` — 100% unconserved, the original B1 signature. It is not production-reachable: `p4_coeval_static_rt_input.h5` has `n_He/n_H = 0.0789474` in all 32768 cells, zero anomalous cells; and P4 (`:424-435`) / P5 (`:479-488`) gate all three species ledgers and raise. Worth a sentence in the docs; not a blocker.

**C1.** `ThermochemicalStepResult` carries `photoelectron_energy` + residual (`multiphysics.py:43-44`), accumulated through `thermochemistry.py`, gated fail-closed at 1e-12/1e-5. P5 ON `1.3120e-17`, OFF exactly `0.0`; B2 ON `6.9768e-9`.

**C2 — now stronger than the pre-stage-2 baseline.** `conservative_hydrogen.py` is clean at HEAD and uses the H-only analytic closure with mean-state `n_e` and no collisional term; Solver A uses the bisected coupled update. I re-ran both: **1612/32768 cells differ, max |Δ| = 2.12e-3, L1 = 1.52855e-5**, reproducing the canonical value. The pre-change code shared `hydrogen_neutral_relaxation` between both solvers, which is why L1 was only 1.15e-6 then. The bit-identical `maximum_fixed_point_residual` the first audit cited is `402 × 2⁻²⁴` — a float32 quantization coincidence at the max, demonstrably not shared state.

**C3.** I re-derived the pinned 6-vector independently from raw `log_xi_-1.0.dat` (x=0.1 is a grid point, so only energy interpolation applies) and got **bitwise 0.0 difference**. It is a genuine table pin. The channel-*split* rule is still applied by the test with production's own rule, so it remains internally validated only — acceptable given it is self-consistent to 0.04% and correctly scoped in the docstring as the project's construction.

**C5.** Verified independently from the tables: `f_ion` at 198.9 eV = `0.37561` (x=1e-4) → `0.00225` (0.9) → `0.00083` (0.99), matching the doc. Control confirmed at `<xHII>=0.9781`, `<T>=5.5788e6 K`. Bands now carry ~2× margin each side. All three channels nonzero ON (H 0.1115, He I 7.06e-3, He II 4.30e-4), all zero OFF.

**C6.** Traced every use of `excitation_rate`: produced and only *reported* (P4:299); never added to `gas_heating_rate` or the thermal update in either family. Policy recorded at P4:464, P5:549 and asserted by the artifact test.

**A2 / A3 / A5.** Provenance binds the full `snrt_core/*.py` glob plus 14 tables and inputs; hashes independently verified, both P5 HDF5 files hash-match. Thresholds now imported from `secondary.py:12-14` with no literals left in either family. The `2.272165e-3` floor delta reproduces exactly and traces to `log_xi_-4.0.dat` — it is the **Lyα excitation channel turning on at the 10.2 eV threshold**, real physics smeared by linear interpolation, not a floor discontinuity. The 10 eV floor itself is exactly continuous.

**A4.** All five library defaults are now `False`.

## C4 — the one item blocking PASS

The mechanism is sound: bracket endpoints are guaranteed (`residual(0) ≤ 0`, `residual(n_max) ≥ 0` from the clips), sign conventions match the bisection update, and `root_bracket_found` (`implicit.py:279`) is hard-gated in B2/P4/P5, all measuring zero.

But the shipped RK4 fixture does not test what it is said to test:

- I enumerated the roots of its residual: **exactly one sign change**, at `n_e = 0.000323770982` — the solver's answer — with `residual(1) = +0.68`. There is no second root to discriminate against, so this fixture exercises bisection, not nearest-root selection.
- It is **not a resolved step**. `k_coll(1e5 K) = 3.824e-9` is self-amplifying: RK4 moves x from 1e-4 to 2.684e-4 (2.7×) while the solver gives 3.238e-4. The 5.54127e-5 discrepancy is **33% of the step's physical change**, and the 1e-4 tolerance is 60% of it. Describing this as "a resolved-step H-only update agrees … to 5.55e-5" (`SECONDARY_IONIZATION_VALIDATION.md:70`) is inaccurate.

I probed four genuinely stiff regimes (T=4e6/dt=1e11, T=4e6/dt=1e6, T=2e4/dt=1e13, T=1e5/dt=1e10). Each has a unique root and the solver matches RK4 truth. The only true multi-root case is `x_e = 0` exactly, which the `stationary` shortcut handles exactly and which *is* tested. So the policy is safe — I found no case where it misfires — but the first audit's C4 ("sound in principle, its policy untested") has been re-labelled rather than closed, and a stage-2 numerical-justification claim is overstated in the docs.

**Remedy:** either shrink the reference step until backward-Euler error is genuinely negligible and re-state the tolerance, or re-word the claim to say the RK4 comparison bounds one-step accuracy on a single-root fixture, with nearest-root discrimination exercised only by the exact `x_e = 0` case.

## Advisory (not stage-2 defects)

- **A1.** 17 MB at 32³ float64 — correctness unaffected at the grids stage 2 uses. The real measured cost is runtime: B2 Solver A `4.20 s → 102.85 s` (~25×), recorded in the artifact but not discussed in the report. 256³ memory is a later performance gate.
- **B2 measures no helium ledgers** (`grep -c` → 0) despite its 200 eV fixture now carrying primordial helium and driving both He channels. P4/P5 gate them, so the production path is covered.
- The `5e-5` A/B bound and the widened fixture bands (`[0.50,0.75]`, `[0.015,0.025]`) were set in the same change that produced the measurements. Normal practice for a changed fixture and disclosed, but "predeclared" overstates it.
- "Secondary ionization is opt-in" (`SECONDARY_IONIZATION_VALIDATION.md:23`) is true of the library API only: P5 `--secondary-ionization` defaults to `fs2010` and P4 hard-codes `use_secondary_ionization=True` at `:205` with no switch.
- **A6.** Linear-in-x across the 0.1→0.5 gap costs ≈±25% in `f_ion` versus a power-law continuation — bounded and upstream-consistent, worth listing.
- 100% Lyα escape is aggressive for LRD-density gas and biases ΔT negative; belongs to the radiation-pressure/Lyα gate.
- The `2.47e-4 → 1.2778e-5` He II improvement: the "after" value is in the canonical JSON; the "before" value is from a discarded intermediate and is not bound anywhere.

## Closure table

| # | Finding | Status |
|---|---|---|
| B1 | zero-helium violation | **closed** (relative-mask residual: unreachable, fails closed) |
| C1 | runtime photoelectron ledger | **closed** |
| C2 | Solver A/B chemistry independence | **closed**, stronger than baseline |
| C3 | independent interpolation reference | **closed** at table-pin level; split rule internally validated |
| C4 | nearest-root policy | **open** — single-root, unresolved-step fixture presented as a resolved-step branch-selection reference |
| C5 | P5 effect strength | **closed** |
| C6 | excitation policy | **closed** as a declared boundary |
| A2 | provenance | **closed** |
| A3 | shared thresholds | **closed** |
| A4 | unsafe defaults | **closed** at library layer; runners ship ON by design |
| A5 | table-floor continuity | **closed** |

**Overall: CONDITIONAL PASS.** Everything the first audit raised is honestly closed except C4, and nothing I found invalidates a shipped B2 or P5 number — I stress-tested the root policy myself and it holds in every regime I could construct. `BLOCK` would be disproportionate. `PASS` is unavailable because C4's remediation re-labels rather than closes the first audit's finding and attaches an inaccurate "resolved step" justification to a stage-2 acceptance claim; the fix is a re-worded sentence or a smaller reference step, not new physics. Stage-3+ items (RSLA/refinement, AGN nine-group, spatial/timestep convergence, coupled H+He front, dust closure) are not implicated.
