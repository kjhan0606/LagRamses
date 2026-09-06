# F-P2 SNIa DTD contract — 2026-09-02

Status: **review-only mathematical kernel; physical SNIa activation blocked**.

## Implemented

- A Fortran interval kernel evaluates a power-law DTD, including the exact
  logarithmic integral for `alpha = -1` and the analytic non-log branch.
- The kernel integrates over `[age_old, age_new]`, clips to `[t_min, t_max]`,
  returns zero outside delay support, and rejects backward/non-finite/negative
  inputs transactionally.
- A cumulative wrapper is defined as the same interval kernel from age zero.
  Adjacent intervals therefore telescope without endpoint-rate sampling or
  hidden restart state.
- A separate event ledger maps the interval-integrated expected event count to
  per-event returned mass, tracked ejecta, diagnostic net yield, energy, and
  signed source-frame momentum. It preserves the untracked residual and emits
  only channel 4, with no terminal remnant owner.
- The identical kernel is present in the native mirror and the
  `patch/lagRamses` production source order. It is compiled by the native
  contract runner; the existing production driver still rejects SNIa until
  the physical event contract is approved.

## Evidence

- `simulation/snrt/tests/run_fp2_snia_dtd_contract.sh`
- `simulation/snrt/data/fp2_snia_dtd_contract_audit.json`
- `simulation/snrt/config/fp2_snia_dtd_contract_v1.json`
- `simulation/snrt/config/fp2_snia_dtd_candidate_matrix_v1.json`
- `provenance/fp2_snia_dtd_literature_dossier_2026-09-02.md`
- `simulation/snrt/config/fp2_snia_event_yield_candidate_matrix_v1.json`
- `provenance/fp2_snia_event_yield_dossier_2026-09-03.md`
- `manifests/fp2_snia_keegans2023_review_v1.json`
- `simulation/snrt/data/fp2_snia_event_yield_asset_audit.json`
- `simulation/snrt/data/fp2_snia_keegans_format_audit.json`
- `simulation/snrt/tools/convert_snia_event_yields.py`
- `simulation/snrt/tools/audit_fp2_snia_event_yield_asset.py`
- `simulation/snrt/tools/audit_fp2_snia_keegans_format.py`
- `simulation/snrt/tests/fp2_snia_event_yield_converter.py`
- `simulation/snrt/native/phase0/stellar_snia_dtd.f90`
- `patch/lagRamses/stellar_snia_dtd.f90`
- `simulation/snrt/native/phase0/stellar_snia_event_ledger.f90`
- `patch/lagRamses/stellar_snia_event_ledger.f90`

The native test covers pre-delay zero, the exact inverse-delay integral,
zero-age and post-support cumulative limits, timestep subdivision,
restart splitting, invalid intervals/normalization, and a non-log power-law.
The event-ledger test covers expected-count scaling, channel ownership,
untracked residual preservation, deterministic repeat evaluation, and
transactional rejection of negative counts and tracked over-return.
The contract audit also verifies the four DTD candidates and two event-yield
candidates, their review-only status, the unresolved decision fields, and both
dossier paths. It also binds the clean review-only asset audit and converter
identity, and binds the source-format audit that exposes the missing H/He/C/N
project fields. It confirms the kernel is review-only and that all physical
parameters and event-source fields remain unset.
The same runner also compiles the DTD and event-ledger mirrors in
`patch/lagRamses` without linking or activating SNIa in the production binary.
The event-yield converter test additionally checks deterministic output,
source checksum admission, license/provenance admission, no-remnant semantics,
mass closure, and overwrite refusal using an isolated unit-only fixture; it
does not create a production asset.

## Intentionally unresolved

This work does not select a DTD family as the project truth. The candidate
`alpha=-1` shape is a mathematical fixture only. A physical F-P2 pass still
requires an approved population/binary model, minimum and maximum delay,
events per initial SSP mass, SNIa yield source and checksum, event energy and
momentum, composition/decay convention, conversion hash, and named approval
id. Until then `enable_snia=true` remains fail-closed in the production driver.
