# F-P1 high-mass internal-gate implementation — 2026-09-03

Project: `/gpfs/kjhan/LRD_JWST`
Gate scope: F-P1H-A--E review controls
Status: **internal controls pass; physical source admission blocked**

## Implemented controls

- **A, build identity:** production fate admission requires exact equality of
  runtime and compiled fate-map SHA256 plus approval id. Review builds compile
  blank identities and cannot be opened by namelist self-assertion.
- **B, source-node schema:** an immutable sidecar defines 84 required fields and
  preserves all 12 resolver axes. Missing values are not physical zeroes;
  failed/direct-collapse nodes require explicit wind, zero terminal ejecta, and
  remnant records. Axis reduction requires an approved freeze or population
  marginalization record.
- **C, ownership/deposition:** channel 3 is an 8--120 M☉ candidate domain that
  must eventually be filtered by source-node fate, not a universal exploding
  interval. Wind, terminal ejecta,
  remnant, scalar radial momentum, energy kind, deposition, and packet ownership
  are explicit. Runtime deposition remains disabled.
- **D, coverage/closure:** candidate mass nodes retain F23 single/binary, LC18,
  WH07, Z9.6, W18/N20, Limongi set-R, and Pop III branch identity. Flattened
  branch unions cannot define interpolation. The high-mass review checks source-
  precision mass closure, failed-node completeness, stable-wind consistency,
  radioactive reference epochs, and duplicate isotopes.
- **E, package admission:** a nine-gate, checksum-bound package contract defines
  source rights, coordinate/population coverage, fate/remnant,
  lifetime/wind, terminal closure, decay, energy/momentum/deposition, pair
  instability, and runtime reproduction. No gate may currently pass: executable
  gate-specific validators must first be implemented and code-registered. It is
  a fail-closed gate, not a source selection.

## Current measured evidence

- Sukhbold W18/N20 high-mass outcome records: 18.
- Non-positive explosion-energy outcomes: 12.
- Terminal-yield records with review mass closure: 6.
- Failed nodes with source remnant records: 0.
- Maximum rounded-source relative mass residual: 0.0033014, below the review-
  only 0.007 tolerance; exact or production closure is not claimed.
- Common W18/N20 wind nodes: 60, 80, 100, and 120 M☉.
- The 100 M☉ stable-wind difference is approximately `1.0e-4 M☉`, entirely in
  Mg in the staged records.
- Radioactive reference-epoch warnings: 12; K-40 cross-segment duplication is
  present in all four common records.
- Qualified physical packages: 0; admitted physical nodes: 0.

## Automated evidence

- `simulation/snrt/config/fp1_source_node_contract_v1.json`
- `simulation/snrt/config/fp1_terminal_deposition_contract_v1.json`
- `simulation/snrt/config/fp1_physical_package_admission_contract_v1.json`
- `simulation/snrt/data/fp1_high_mass_seam_review.json`
- `simulation/snrt/data/g2_candidate_grid_coverage_audit.json`
- `simulation/snrt/data/fp1_physical_package_admission_audit.json`
- `simulation/snrt/config/fp1_fate_admission_sidecar_v1.json`

`run_fp1_population_fate_contract.sh` passes with seven checksum-bound fate
artifacts. Candidate-grid coverage and source branch identity tests also pass.
All successful results are review/control results; no canonical physical row or
runtime terminal deposition is enabled.

## Remaining physical blockers

The next source package must provide or explicitly and quantitatively resolve
the required metallicity/population domain, all source-node fates and remnants,
age-resolved winds or an approved lumping error model, decay epoch/network and
isotope duplication, injected-energy mapping, radial momentum and deposition,
PPISN/PISN ownership, and source rights. The LC18 failed-model Wind discrepancy
requires author clarification or a corrected release. Sukhbold W18/N20 remains
validation evidence and is not a project-wide fate law.

F-P1H-F cannot start until at least one physical package passes all nine gates
and the physical-node inventory is non-empty. At that point AGY and Fable will
receive the same frozen implementation/evidence bundle for independent audit.

## Post-audit remediation

The first independent AGY/Fable audit agreed that the physics remains blocked.
Fable also identified internal-control defects that were independently
reproduced and corrected before re-audit:

- high-mass residual, failed-node, remnant, wind, radioactive-epoch, and K-40
  conclusions are computed from source records and covered by mutation tests;
- the source-node contract now validates complete 84-field records, unique
  identities, resolver axes, finite/nonnegative values, cumulative-wind
  monotonicity, direct-collapse zero-terminal/remnant semantics, and PISN
  zero-remnant semantics;
- canonical conversion requires a source-node mapping sidecar, and the
  converter and asset auditor hash and inspect the actual contract, mapping,
  and asset files rather than accepting arbitrary digest-shaped strings;
- the initial per-gate evaluator was found to trust a checksum-bound but
  unexecuted validator artifact. That route is now disabled completely: a
  candidate cannot pass any gate until gate-specific executable validators are
  implemented and code-registered;
- the F-P1 runner is now part of G2 preflight, which also executes the 121 M☉
  coverage negative and the high-mass runtime refusal test;
- the actual stellar-enrichment namelist route now carries fate policy, map
  digest, and approval identity transactionally. A runtime claim still cannot
  authorize a review build because compiled identities are blank;
- the enrichment driver explicitly refuses a channel-3 upper mass above
  40 M☉ while the compile-time source-node fate consumer capability is false;
  no physical source/deposition behavior is inferred or fabricated;
- stale source-parity wording and the phrase that overstated a configured mass
  window as physically resolved were corrected.

Bounded regression evidence after these changes is:

- `G1_NATIVE_CONTRACT_TEST_OK` and exact CPU JAX/native agreement for six
  interpolation queries;
- `G2_POPULATION_LEDGER_RUN_OK`, including explicit refusal of the high-mass
  SNII window without a source-node consumer;
- `FP1_POPULATION_FATE_CONTRACT_OK` and all converter/asset/coverage tests;
- expected fail-closed `G2_PREFLIGHT_BLOCKED`, because no physical package is
  qualified;
- expected `STELLAR_SOURCE_PARITY_BLOCKED` on missing current
  production-linked build evidence.

These remediations strengthen the review gate only. They do not populate a
single physical high-mass node, provide the missing source package, implement
the source-node runtime consumer, or authorize F-P1H-F.

## First re-audit follow-up

AGY confirmed all first-round remediations. Fable independently reproduced
three additional future-promotion weaknesses, and local probes confirmed all
three before they were changed:

- outcome-aware node validation now rejects nonzero direct-collapse terminal
  components, null failed/direct-collapse wind histories, failed nodes without
  baryonic remnants, untyped rights identifiers, component mass non-closure,
  resolver-branch cell overlap, and cells extending to 121 M☉;
- arbitrary evidence JSON plus a hashed one-line “validator” can no longer
  qualify a package. Gate evidence activation is disabled until executable,
  gate-specific validators are implemented and registered in reviewed code;
- conversion now requires the approved repository source-node contract, a
  matching approval id, and a source-node id whose mass and metallicity match
  the canonical row. The asset auditor independently repeats the contract audit,
  path, approval, node-id, and coordinate checks;
- both interval and cumulative driver entry points have explicit high-mass
  refusal regression tests;
- roadmap wording no longer describes contract-only ownership/exactly-once
  declarations as runtime implementation.

The remaining absence of gate-specific executable validators is an explicit
physical-package blocker, not a path by which declarative evidence can pass.

## Second re-audit follow-up

Round 2 confirmed the prior fixes and current fail-closed state. Fable then
found five latent promotion defects that were independently accepted and
closed:

- approved nodes now require non-null rights/provenance fields, a SHA256 package
  fingerprint, and an approval id matching the node contract; metallicity is
  finite and typed, binary/population state is explicit, and half-open mass-cell
  membership is enforced;
- F-P1H-E evidence paths must remain repository-relative and cannot escape
  through symlinks;
- converter and asset admission now require the approved repository node
  contract **and** an admitted F-P1H-E physical package, matching approval id,
  selected source-package hash, and selected source-node mapping hash;
- canonical rows are checked against their node's channel, lifetime, cumulative
  wind history, terminal outcome, ejecta, remnant, selected energy semantics,
  and source-frame momentum. A dedicated projection test covers interpolation,
  valid direct collapse, wrong channel/Z, and the exact inconsistent terminal
  payload used in the audit reproduction;
- the mapping is one-way hash-bound by the package selection; it does not embed
  the package-contract hash, avoiding a cryptographic fixed-point dependency.

Because F-P1H-E remains blocked and gate-specific executable validators do not
yet exist, successful production conversion is intentionally untestable today.
Unit tests exercise deterministic normalization and physical projection, then
require the final conversion call to stop on package admission.

## Final latent-path hardening

After round 3 confirmed F1--F5, four lower-severity future-state findings were
also closed:

- all five F-P1H-E evidence paths are pinned to their exact repository paths;
  absolute, `..`, and symlink escape remain prohibited;
- F-P1H-E executes the source-node contract audit itself, and converter/asset
  checks require the package evidence hash to equal the exact node contract they
  bind;
- an admitted package must have valid package/mapping SHA256 identities and its
  package SHA must equal every admitted node's package fingerprint as well as
  the normalized source/sidecar source hash;
- approved research-use and redistribution statuses use the closed vocabulary
  `approved|verified|permitted`, and all conditional binary-population axes are
  null, a non-empty identifier, or a structured record—never a bool or integer.

AGY and Fable independently re-audited this final state. Both marked every
latent-path hardening item `VERIFIED FIXED`, found no new concrete bypass, and
returned **PASS** for the internal fail-closed controls. Both independently
retained **BLOCK** for physical resolution, F-P1H-F, production, and publication.
