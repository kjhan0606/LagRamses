# Claude Opus 5 closure audit — F-P2.4 native nine-group remediation

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

Auditor: Claude Opus 5

Date: 2026-09-05

Mode: read-only; no edits, jobs, network, Python, or RAMSES evolution.

## Verdict

**PASS**

### C1 — spectral upper bounds: closed

`validate_species_table` applies `max_cross_section_cm2 = 1.0d-17` and the
interval-derived excess-energy bound to every supported H I, He I, and He II
group. Canonical edges are validated before the species pass and all three
thresholds are exact group boundaries, so there is no unchecked straddling
branch. The largest checked-in cross section is `4.66e-18 cm^2`; reference
values remain valid. The one-decade H I transcription-slip probe is rejected.

### C2 — direct rejection and checkpoint evidence: closed

The native runner drives real loader calls for unset environment, missing file,
malformed namelist, wrong version, malformed identity, edge digest mismatch,
unknown fraction semantics, candidate status, intrinsic-fraction runtime
blocking, and reference-control opt-in. The checkpoint smoke calls the real
`snrt_state_checkpoint_write` and `_read` with a nonzero `8 × 80 × 9` payload.
It rejects a candidate contract with `ierr=4` before state mutation and then
round-trips every direction/group/slot value and the cell-slot map under the
reference contract. Writer/reader record order and declared character lengths
match.

### C3 — native science limitations: closed as an explicit boundary

`SNRT_NATIVE_GROUP_CONTRACT.md` quantifies the 7--36% emission-mean versus
absorber-weighted heating gap and names the missing secondary-ionization and
recombination channels, including the `2564.90 eV` group-9 excess. These remain
later G3/G4 science gates and are not presented as solved.

## Additional hardening verified

- `fraction_semantics` is part of the contract and checkpoint identity.
  Resolved-domain runtime accepts only `escaped`; `intrinsic` remains
  inspectable but is blocked.
- `reference_control` requires `SNRT_ALLOW_REFERENCE_CONTROL=1`.
- The overlong opt-in environment value is checked with statement-level guards
  before any substring is formed, so no out-of-bounds substring is evaluated.
- The recorded GNU Fortran and `mpiifx` native runner outputs pass, and the
  latest `make -C bin SNRT=1 USE_CUDA=1 ramses` link post-dates the final
  remediation. No stale four-group or CUDA ABI issue was found.

## Non-blocking observations

The `1e-17 cm^2` comment's stated margin is closer to 1.35x than 2x relative
to `7.4e-18 cm^2`; this is documentation wording only and not a physical or
runtime defect. Upper-bound probes use H I directly, while He coverage follows
the same shared validator structurally. The remaining F1/F2, He chemistry,
HDF5 restart, dust, transport accuracy, and 40--120 M☉ yield-seam tasks are
outside F-P2.4 and remain in the high-level roadmap.

No physical AGN/stellar SED approval, live RT+feedback run, or publication
science validation is granted by this audit.
