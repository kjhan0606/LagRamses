# F-P1.5-R AGN effective-efficiency convention bundle — implementation evidence

Date: 2026-09-04 (KST)  
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)  
Branch at evidence capture: `main`  
Source revision at evidence capture: `343619e` (dirty worktree)  
Operator approval: received in the current run  

Status: **implementation and focused evidence complete; Claude Opus 5 final
bundle-end re-audit PASS; F-P1.5-R engineering bundle accepted**

No runtime activation, large RAMSES job, commit, or push was performed.

## Implemented boundary

The bundle closes the AGN coefficient and mass-convention mismatch between the
coarse-state ledger and the SNRT source driver.  It is an engineering and
accounting closure only; it does not approve an AGN SED, obscuration model,
escape fraction, jet/radiation-pressure prescription, or live RT-hydro run.

- `patch/lagRamses/snrt_agn_efficiency.f90` is the single pure,
  RAMSES-independent resolver used by both consumers.  It computes
  `inflow=min(max(Bondi,0),max(Eddington,0))`, the Eddington ratio, the strict
  MAD low-state transform below `X_floor`, and the thermal/high-state value.
  Raw/resolved base efficiency is strict `(0,1)` and effective efficiency is
  `[0,1)`.  Spin-disabled `.1`, spin-enabled zero/uninitialized, non-finite or
  non-positive raw values, `raw>=1`, invalid rates, zero Eddington, and invalid
  MAD floors receive explicit status/mode information.  Non-promotable states
  are rejected by the live driver without a hidden clamp.
- `patch/lagRamses/sink_particle.kjhan.f90` calls that helper once per coarse
  writer record.  `Lbol` uses `effective_radiative_efficiency * inflow * c^2`.
  The JSON record preserves raw sink, resolved base, effective coefficient,
  status, mode, and contract flag, and remains at the pre-feedback/pre-reset
  boundary.
- `patch/lagRamses/snrt_ramses_driver.f90` calls the same helper and passes
  only the effective coefficient to `snrt_agn_photon_budget`.  The source mass
  is the increment of the stable-`idsink` cumulative supplied ledger
  `min(dMBH_coarse,dMEd_coarse)`, not retained `dMsmbh`.
- `dMsmbh` is observed separately through a retained cursor and is used only
  for the one-sided check
  `delta_retained <= (1-epsilon_eff)*delta_inflow*(1+tol)` plus the documented
  numerical tolerance.  Gas-floor clipping, saved-energy suppression, and
  AGN-blast decreases are treated as admitted slack/rebase conditions; an
  unexplained retained excess fails closed.  The accounting marker advances
  only after the all-group transaction commits.  The marker and retained
  cursor are keyed by `idsink`, not mutable array position.
  The cumulative `dMBH_coarse`/`dMEd_coarse` supply is rate-based before the
  gas-floor retention step, so the photon budget is an upper-bound supplied
  inflow convention; it is not an equality claim for actually retained BH
  mass.  A failed or skipped source transaction does not consume the retained
  cursor, so the one-sided check is re-armed on the next AMR-level call.
- `patch/lagRamses/snrt_agn_source.f90` names the first photon-budget argument
  `delta_inflow_mass_code`, rejects non-finite inputs and `radiative_efficiency
  >= 1`, and retains the existing atomic multi-group transaction.  The source
  smoke program records the supplied-inflow positional contract.
- `simulation/snrt/snrt_core/sink_diagnostic.py` and the P4 converter preserve
  raw/resolved/effective fields plus `efficiency_status` and
  `efficiency_contract_ok`.  Coarse records now require those two fields; the
  P4 coarse path emits them and refuses any false contract before opening the
  output artifact.  Non-promotable rows remain readable in the append-only
  ledger for diagnostics, but `read_agn_coarse_state` excludes them from
  `AgnCoarseState` promotion and the converter refuses a requested step that
  contains one.  The legacy `sinkprops` path leaves the status/contract
  columns blank and labels the source
  `legacy_sinkprops_mode_unresolved`, rather than claiming helper parity.  The
  reader promotes raw zero only for the explicit promotable
  `spin_disabled_default` status; spin-enabled uninitialized zero remains
  readable for diagnostics but is rejected at state promotion and conversion.
  This prevents an offline converter from promoting a known initialization or
  configuration divergence.
- `bin/Makefile` links the helper unconditionally and gives the writer and
  driver direct helper/module dependencies, including the driver's direct
  `amr_commons` and `pm_commons` edges.  The existing SNRT/CUDA runtime gate
  and transport/transaction graph are unchanged.

## Focused evidence

All commands below were run from `/gpfs/kjhan/LRD_JWST`.

### Native and arithmetic tests

```text
python3 -m py_compile simulation/snrt/snrt_core/sink_diagnostic.py \
  simulation/snrt/tools/audit_agn_coarse_ledger.py \
  simulation/snrt/tests/agn_effective_efficiency.py \
  simulation/snrt/tests/agn_ledger_transaction.py
```

Passed.  The native runner was then executed:

```text
bash simulation/snrt/tests/run_fp15_agn_efficiency.sh
```

Result:

```text
SNRT_AGN_EFFICIENCY_OK thermal=1 mad_low=0.01 boundary=high spin_init=visible invalid=fail_closed
SNRT_AGN_SOURCE_OK luminosity=  2.246888E+53 photons=  7.011986E+63
SNRT_AGN_EFFICIENCY_NATIVE_TEST_OK helper=compiled source_api=compiled runtime=disabled
```

The native cases include thermal, MAD-low/high, strict floor boundary, zero
Eddington, spin-disabled default, spin-enabled uninitialized, invalid raw
efficiency, non-finite/negative rates, and disabled MAD floor.  The source
smoke includes supplied-inflow scaling, zero-inflow, direct unity-efficiency
rejection, and all-or-nothing group deposition.

```text
python3 simulation/snrt/tests/agn_effective_efficiency.py
```

```text
AGN_EFFECTIVE_EFFICIENCY_TEST_OK algebra=thermal,mad_low,mad_high,boundary native_parity=run_fp15_smoke supplied_mass=min_bondi_edd carryover=no_duplicate epsilon=effective
```

```text
PYTHONPATH=simulation/snrt simulation/snrt/.venv/bin/python \
  simulation/snrt/tests/agn_ledger_transaction.py
```

```text
AGN_LEDGER_TRANSACTION_TEST_OK duplicates=collapsed conflict=reject idle_effective=accepted null=reject year=365.25 transaction=atomic
```

That test also accepts the explicit spin-disabled raw-zero fallback and keeps
the spin-enabled uninitialized raw-zero record readable while rejecting its
promotion.  The fixture was
refreshed to carry raw/resolved/effective/status/mode/contract fields, with a
thermal record whose effective value equals the resolved base and an explicit
idle MAD-quenched record.  Existing AGN regressions remained green:

The same test invokes `p4_build_agn_rate_ledger.py` and verifies that the
accepted coarse CSV contains `efficiency_status=0`,
`efficiency_contract_ok=true`, and the helper provenance label.  It also
constructs a valid-algebra `MAD_FLOOR_DISABLED` row with
`efficiency_contract_ok=false`; the converter rejects it and creates no CSV.

```text
AGN_PHOTON_LEDGER_TEST_OK p0_groups=9 legacy_groups=5 hard_xray_group=positive subthreshold_opacity=zero edges=exact
AGN_NINE_GROUP_ARTIFACT_OK hard_q=3.52996e+51 hard_to_soft_q=0.219425 hard_supported_sed_fraction=0.0437846 hard_bolometric_fraction=0.00794605
```

### Machine-readable ledger audit

```text
PYTHONPATH=simulation/snrt simulation/snrt/.venv/bin/python \
  simulation/snrt/tools/audit_agn_coarse_ledger.py \
  --input simulation/snrt/data/agn_coarse_state_transaction_fixture.jsonl \
  --output simulation/snrt/data/agn_coarse_ledger_transaction_audit.json \
  --helper patch/lagRamses/snrt_agn_efficiency.f90
```

Result:

```text
AGN_COARSE_LEDGER_AUDIT_PASS records=2 output=simulation/snrt/data/agn_coarse_ledger_transaction_audit.json
```

The generated `snrt_agn_coarse_ledger_audit_v2` report has `passed: true`, no
failed static criteria, two canonical fixture records, one semantic duplicate
collapsed, and `physical_closure_claim: false`.  The helper SHA256 recorded in
the report is:

```text
34915c8cafc688763e38aa11641e613662ec2bbf420f96f64448cabeaa2bcc01
```

The audit checks both helper call sites, no independent MAD formula, effective
photon coefficient, supplied-mass minimum, one-sided retained bound, stable
sink identity, atomic commit ordering, source API, JSON fields, direct
Makefile dependencies, and the CUDA runtime gate.

The latest Opus conditional audit and its F1 finding are recorded in
`provenance/opus5_fp1_5_agn_effective_efficiency_convention_bundle_conditional_audit_2026-09-05.md`.
F1 is now closed in the current worktree; F2 (variable naming) and F3
(exact-status-bit assertion strength) remain bounded non-blocking follow-ups.
F1 and C1 are closed in the current worktree.  The final Opus `PASS` is
recorded in
`provenance/opus5_fp1_5_agn_effective_efficiency_convention_bundle_final_audit_2026-09-05.md`.
The rank-local marker limitation, driver variable naming, exact-status-bit
smoke strengthening, and Python-file provenance hashes remain non-blocking
follow-ups before live activation.

### Build evidence

```text
make -C bin -B -j4 SNRT=1 USE_CUDA=1 snrt_ramses_driver.o
```

Passed.  This is an actual SNRT/CUDA production-driver object compilation;
only pre-existing Intel global-name-length warnings from the topology module
were emitted.

```text
make -C bin -B -j4 ramses
```

Passed.  The default non-SNRT `ramses_final3d` compile/link included
`snrt_agn_efficiency.o` and `sink_particle.kjhan.o`.

```text
make -C bin -n ramses
make -C bin -n SNRT=1 USE_CUDA=1 ramses
git diff --check
```

All passed.  The two dry-runs contained 8 and 34 commands respectively after
the normal generated-Makefile step.

The project venv contains the requested CPU JAX environment:

```text
JAX 0.11.1 JAXLIB 0.11.1
```

No system-Python JAX assumption is made; SNRT Python tests use
`simulation/snrt/.venv/bin/python`.

## Limitations and explicit non-claims

- No live RT-hydro/AGN run was launched and `SNRT_RT_ENABLE` was not set.
- No AGN SED shape, obscuration, escape fraction, spectral calibration,
  radiation-pressure/jet coupling, or publication-level AGN physics is
  approved.  The existing spectral split remains outside this bundle.
- The helper's single-MAD-formula claim is scoped to the coarse writer and
  SNRT driver.  The legacy RAMSES `accrete_bondi` path retains its own
  accretion transform and `AGN_blast` retains its own feedback-energy
  prescription; those physical conventions are not silently declared
  unified by this bundle.
- No durable cross-restart crash journal is implemented.  On the first driver
  observation after restart, the retained cursor is initialized to the
  observed value; a retrospective retained-mass check is therefore not
  claimed.  Cross-coarse-step deferred-`Esave` carry-over is handled by the
  in-memory cursor and has no persisted crash journal.
- The fixture is arithmetic/transaction evidence and must be refreshed from a
  real production dump before any physical claim.
- Stellar fate/yield, dust/IR, and hydro closure remain outside this bundle.

The required Claude Opus 5 bundle-end gate is complete with `PASS`.  No commit
or push has been performed.

## First-audit repair record

The first Opus bundle-end audit returned `FAIL` because the native MAD-high
smoke expected `inflow=1.0` for `Bondi=0.2, Eddington=1.0`, and because the
evidence hash predated the `pure` helper declaration.  The smoke now expects
the defined minimum `0.2`; the evidence, fixture, and generated audit report
were refreshed.  The audit tool's commit-order regex was also updated for the
safer block-form commit, and the retained cursor is now consumed only after a
successful transaction.  The native smoke and full focused suite pass after
these repairs.  The subsequent Opus conditional-pass F1 finding exposed the
missing parser/converter propagation of status and contract metadata; the
current worktree carries those arrays, emits explicit CSV columns, rejects
false coarse contracts, and tests both accepted and rejected paths.  The final
Opus re-audit then verified F1/C1 closure and accepted this engineering bundle
with the non-blocking live-activation follow-ups recorded above.
