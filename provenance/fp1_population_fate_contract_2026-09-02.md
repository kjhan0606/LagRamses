# F-P1 stellar population and fate contract — 2026-09-02

## Current decision after audit remediation

**SOURCE-CELL/ADMISSION BUNDLE IMPLEMENTED; CLAUDE OPUS 5 BUNDLED AUDIT
PENDING; SCIENTIFIC GATE REMAINS BLOCKED.**

The production implementation now makes source normalization, IMF support,
population mode, binary fraction, and remnant ownership explicit and
fail-closed.  This does not select or approve the physical stellar-evolution
and fate package.  The provisional Kroupa default and all staged yield/fate
candidates therefore remain unsuitable for a production or publication
claim.

## Implemented contract

- `yield_source_basis` distinguishes per-star, per-event, SSP-cumulative, and
  SSP-rate inputs.  The current production path accepts only
  `per_star_cumulative` and performs exactly one IMF convolution.  An already
  SSP-normalized source is rejected before integration.
- `imf_mass_min_msun`, `imf_mass_max_msun`, and `binary_fraction` are mandatory
  in channel-resolved mode.  The implemented IMF lower support is 0.08 Msun;
  every enabled wind/AGB/SNII window must lie inside the configured support.
- Salpeter, Kroupa, Chabrier, and the current provisional Pop-III shapes use
  piecewise analytic mass normalization, satisfying
  `integral(m*phi(m) dm)=1` over the configured support.
- The Chabrier high-mass power law carries the low-branch value at 1 Msun as
  its amplitude.  The former factor-3.583 discontinuity and associated
  massive-star overweighting are removed.
- `single_star_ssp` requires zero binary fraction.  `binary_ssp` requires a
  positive fraction but remains rejected by production initialization until
  period, mass-ratio, transfer, and common-envelope assumptions are selected
  and implemented.
- The production timestep driver now executes `stellar_population_ledger`
  using cumulative states already evaluated for the source increment, so no
  second IMF integration is introduced.  It aggregates returned and
  tracked/untracked ejecta once, rejects remnant mass from a non-owner, rejects
  over-return, and checks the prospective RAMSES particle mass against
  `living + remnant`.  The identity
  `initial = living + remnant + returned` defines living mass; its useful
  guards are non-negativity, channel consistency, ownership, and agreement
  with the independently advanced particle mass.
- Production startup now runs the same strict table audit as the native
  oracle: tracked ejecta may be below returned mass, while duplicate
  coordinates, incomplete Cartesian grids, decreasing cumulative material or
  energy, non-owner remnants, and nonzero age-zero cumulative fields are
  rejected.  The residual is deposited into generic metallicity together with
  all tracked metals, regardless of individual field switches.
- A successfully read legacy namelist commits its element/channel switches as
  one legacy-mode transaction instead of committing only the mode.
- The production runtime selects piecewise-constant mass-cell assignment with
  the explicit half-open source-node convention. It does not authorize
  metallicity, rotation, engine, or age interpolation.
- The population ledger records an IMF-weighted unresolved initial-mass bucket
  from the explicit F-P1 intervals. The bucket is diagnostic-only: it is not
  part of returned/living/remnant closure and is never deposited.
- The F-P1 admission sidecar hashes the fate map, resolver, source, and
  physics contracts and cross-checks unresolved intervals and approval id.
  The current zero-node/review-only sidecar therefore remains blocked.
- SNIa and PISN remain fail-closed.  They cannot fall through to ordinary IMF
  integration; their DTD/binary and metallicity/core-mass fate gates remain
  F-P2 and F-P3.

## Reproducible evidence

- `tests/run_stellar_feedback_policy_unit.sh`: PASS, including missing basis,
  invalid IMF support, inconsistent binary fractions, and a channel window
  outside IMF support.
- `tests/run_stellar_population_contract_unit.sh`:
  `STELLAR_POPULATION_CONTRACT_UNIT_OK`.  It checks all four IMF
  normalizations, support dependence, double-convolution rejection,
  terminal-remnant ownership, untracked ejecta, population mass closure, and
  mandatory timestep-path ledger execution against production sources.
- `tests/run_stellar_yield_audit_contract_unit.sh`: PASS against production
  sources.  It independently mutates residual mass, tracked over-return,
  duplicate coordinates, a missing Cartesian corner, age-zero values,
  decreasing energy, missing age-zero support, and remnant ownership.
- `tests/run_stellar_residual_deposition_unit.sh`: PASS.  All three linked
  deposition adapters accept an untracked residual, deposit the full generic
  metal budget when mapped, and reject tracked over-return.
- `tests/run_stellar_imf_jax_differential.sh`:
  `FP1_IMF_JAX_DIFFERENTIAL_OK rows=8`, JAX 0.11.1 CPU.  A JAX shape evaluated
  by independent 256-point Gauss--Legendre quadrature per branch agrees with
  the Fortran analytic normalizations for four IMFs at two support intervals;
  evidence is
  `simulation/snrt/data/fp1_imf_jax_differential.json`.
- `simulation/snrt/tests/run_g2_population_ledger.sh`:
  `G2_POPULATION_LEDGER_RUN_OK` on the native differential mirror.
- `simulation/snrt/tests/run_fp1_population_fate_contract.sh`:
  population/fate audit, resolver test, checksum-sidecar mutation test, and
  admission audit pass with `blocked_review_only` status.
- The prior `/gpfs` production-linked build/link/smoke record predates the
  current F-P1H source changes. Source parity therefore reports
  `blocked=production_linked_build_evidence`, as required. A clean rebuild and
  regenerated evidence are mandatory before promotion; the scientific F-P1
  fate gate also remains blocked by the unresolved intervals below.
- `P0_DIAGNOSTIC=1 simulation/snrt/tests/run_g1_native_contract.sh`:
  `G1_NATIVE_CONTRACT_TEST_OK`, exact six-query JAX differential, diagnostic
  marker retained because it is not itself production evidence.
- `simulation/snrt/tests/run_p04_production_negative.sh` reports
  `P04_PRODUCTION_NEGATIVE_OK baseline=3 snia_fail_closed=pass`.

## Open scientific conditions

F-P1 is not a scientific PASS until one exact population/fate model is pinned
with source/version/license/checksums and the following are approved:

- the default IMF and its published functional convention, not merely its
  normalized code shape;
- the single/binary population choice and, if binary, its complete binary
  parameter distributions and interaction prescriptions;
- metallicity-dependent lifetimes, wind/terminal phase ownership, failed
  explosions, fallback, compact-remnant prescription, and the AGB/ECSN/CCSN
  transition;
- a source sidecar whose declared per-star basis agrees with the production
  integrator and whose full mass-metallicity-age domain covers every enabled
  channel.

The provisional 40--120 Msun interval has winds but no terminal owner (about
6.8% of initial mass for the current Kroupa shape). It must be resolved by the
selected fate map in F-P1. PISN ownership is now aligned as “no remnant owner”
in the contracts, but physical PPISN/PISN eligibility and terminal semantics
remain a separate F-P3 blocker.

## First-audit disposition

| Finding | Disposition before consolidated re-audit |
|---|---|
| D1 Chabrier discontinuity/self-copying normalization check | Fixed; continuity regression plus independent numerical quadrature. |
| D2 residual ejecta contract/deposition mismatch | Fixed; strict inequality contract and generic-metal residual on production path. |
| D3 ledger linked but not executed | Fixed; mandatory timestep-driver execution and particle-mass cross-check. |
| D4 weak production audit | Fixed; production/native audit modules are byte-identical and production-specific mutation tests pass. |
| D5 legacy partial commit | Fixed; successful legacy reads commit all runtime switches together. |

The first consolidated re-audit returned an engineering `CONDITIONAL PASS`
after finding two residual implementation issues: box-mass-scaled runtime
tolerances and obsolete equality closure in three linked deposition adapters.
Both were independently reproduced and repaired.  They are recorded in
`claude_opus5_fp1_population_fate_reaudit_2026-09-02.md` and will be reviewed
with the next substantial F-P1 bundle rather than by another immediate audit.

Existing candidate audits show that no staged package currently satisfies all
of those conditions.  No value was inferred to manufacture an approval.
