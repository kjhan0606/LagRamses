The plan-file tool is disabled in this session, so I'm delivering the audit inline.

# Stage-2 final re-audit — C4 closure

**Verdict: PASS**

All four commands pass and reproduce their documented numbers: `electron_closure=3.47804e-12`, `rk4_absolute_error=2.38311e-10`, `rk4_step_relative_error=0.00057435`, `continuity_delta=0.000335457`, `A_B_L1=1.52855e-05`, `fixed_point=2.39611e-05`, `delta_mean_xhii=2.24653e-08`, `delta_mean_temperature_k=-0.00752274`. Everything below is measured from source, not read from the docs.

## C4 is now honestly closed

**The single-root label is correct.** I rebuilt the residual `n_e − implied(n_e)` directly from `hydrogen_neutral_relaxation` for the new fixture. Over the whole physical domain `[0, n_max=1]` at 2e6 samples there is **exactly one sign change**, at `n_e = 1.004151e-4`, with `residual(0) = −1.003e-4` and `residual(n_max) = +0.9988`. A 4e6-sample refinement of `[0, 4e-4]` also finds one.

**The step is genuinely resolved and the reference is a true value.** RK4 at 4096 substeps agrees with RK4 at 262144 to `6.78e-20` — nine orders below the measured 2.38e-10. The reference contributes nothing to the discrepancy. Contrast the old fixture: 33% of the step change, now `5.74e-4` of it.

**The arithmetic is exact.** Independently: solver `1.00415162467615e-4`, RK4 `1.00414924156117e-4`, absolute error `2.383114985866e-10`, physical change `4.149241561167e-07`, ratio `5.743495e-4`. Matches the doc's `2.38e-10` / `5.74e-4` and the prompt's `5.74350e-4`. Margins are 2.098× on the `5e-10` limit and 1.741× on the `1e-3` limit.

**The fixture measures what it says.** The error is bit-identical at 24, 32, 48 and 64 bisection iterations, so it is time-integration error from freezing `n_e`, not root-finder truncation. `hydrogen_neutral_relaxation` (`implicit.py:51-68`) is the exact analytic constant-rate solution at fixed `n_e`, and the test's RHS matches its rate structure — the right construction for isolating integration error. A no-op update would miss the relative gate by ~1000×.

**The neutral fixture really is multi-root, and the hazard is real.** At `T=4e6`, `dt=1e11`, no photons: `residual(0) = 0.0` exactly, plus a second root at `n_e ≈ 1.157770`, residual negative in between, `residual(n_max) = +2.29e-4`. A global bracket over `[0, n_max]` would land on the ionized root and fully ionize an untouched neutral cell. The `stationary` shortcut (`implicit.py:178, 259-263`) returns 0, and the test asserts it. That is a substantive check, so splitting the two claims is a real separation, not relabelling.

**The policy itself now survives stress-testing.** Across 1260 configurations (`T` 1e4–1e7, `dt` 1e5–1e13, `x0` 0–0.5, `Γ` 0–1e-10) there are **zero** non-stationary multi-root cases. A second 980-configuration sweep compared the solver's `n_e` against a brute-force outward first-root search: zero bracket-flag failures, zero cases where the answer was not a root. Two configurations first flagged as "root skipped" were an artifact of my own test — both are `T=1e4, Γ=0` near-no-op states where the residual across the 5.6e-17-wide window is smooth and monotone with exactly one crossing, and the solver's value *is* that root to 5.6e-10 relative.

**Audit #1's sub-finding (c) is closed in shipped code.** An all-False `crossed` does not silently return the incoming density; it falls back to the guaranteed bracket `[0, n_max]` (`implicit.py:230-239`) and returns `stationary | has_crossing`, hard-gated as `electron_root_bracket_failure_count == 0` in P4 (`:433`) and P5 (`:489`). `git diff` confirms the function is a new 173-line addition.

## Wording

`SECONDARY_IONIZATION_VALIDATION.md:78-82` claims "a resolved, single-root H-only step", attributes the multi-root check separately to the exact electron-free limit, and defers large-step accuracy to stage 5 — every clause verified above. `B2_PRODUCTION_SOLVER_VALIDATION.md:34-36` carries the same deferral. The old `5.54127e-05` claim survives only in the prior audit records, where it belongs. Remaining "nearest-root" mentions (`B2:32`, `B2:55`, `implicit.py:130`, the test's module docstring) describe the mechanism or the returned flag, not a validation result.

## Prior advisories — none escalates

Three are now in the docs: the trace-helium contract (`:43-45`), "opt-in in the library API" with the runner policy stated (`:23-25`), and the linear-in-x 0.1→0.5 gap (`:47-49`). Unchanged and still non-blocking: B2 reports no helium ledger residuals (grep 0 in doc, artifact and JSON) — that coverage lives in the P4/P5 gates and the one-cell fixture at `coupled_photo_collisional_hhe.py:193-195`; the unbound "before" value `2.47e-4` at `:36`; A1 memory and the 25× runtime; 100% Lyα escape.

## Remaining boundary

Single-rootedness for non-stationary states is established **empirically** (1260 + 980 configurations), not proved analytically. The "quartic sampling steps over a pair of roots" hazard therefore has no demonstrated instance and no proof of impossibility; it is bounded by the hard-gated bracket flag, which never fired. The `5e-10` / `1e-3` limits are post-hoc round numbers with under-2× margin — regression bands rather than physics tolerances, which is adequate for a deterministic float64 fixture with a converged reference. C3's channel-split rule remains internally validated only, as accepted in re-audit #1. None of these is a stage-2 defect.
