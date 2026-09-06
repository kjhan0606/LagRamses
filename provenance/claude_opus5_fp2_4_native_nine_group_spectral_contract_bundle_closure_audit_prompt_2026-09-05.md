# Claude Opus 5 closure audit request — F-P2.4 final C1 remediation

You are the final closure auditor for F-P2.4 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`, branch `main`.

Read-only only: do not edit files, run jobs, launch RAMSES, use the network, or
invoke Python. Inspect the actual native Fortran/CUDA source and the recorded
evidence. The project goal is production/publication-ready high-level RAMSES
RT, stellar/AGN feedback, and dust; F-P2.4 is only the native SNRT
nine-group/source-contract wiring bundle.

The first end audit and the intermediate follow-up are recorded at:

- `provenance/claude_opus5_fp2_4_native_nine_group_spectral_contract_bundle_end_audit_2026-09-05.md`
- `/home/kjhan/.claude/plans/claude-opus-5-rippling-stearns.md` (the latest
  follow-up verdict was `CONDITIONAL PASS` solely because the initial
  `1e-16 cm^2` cross-section ceiling was too loose)

The current implementation evidence is:

`provenance/fp2_4_native_nine_group_spectral_contract_bundle_implementation_evidence_2026-09-05.md`

The final remediation changed the ceiling to `1e-17 cm^2`, documented the
physical margin (largest supported threshold-adjacent primordial cross section
about `7.4e-18 cm^2`), added a one-decade H-I transcription-slip test, removed
the Fortran non-short-circuit substring hazard in the reference opt-in path,
and made the checkpoint payload test nonuniform over all directions, groups,
and slots.

Verify narrowly:

- `snrt_spectral_contract.f90` applies `1e-17 cm^2` and the interval-derived
  excess-energy upper bound to every supported H I/He I/He II group, while the
  reference values remain valid.
- The actual GNU Fortran and `mpiifx` runner evidence covers all loader
  rejection paths, candidate/intrinsic/reference runtime gates, upper-bound
  probes, and the nonzero version-4 checkpoint write/read plus pre-mutation
  identity rejection.
- The latest production `make -C bin SNRT=1 USE_CUDA=1 ramses` evidence refers
  to the current source, and no stale four-group or ABI issue was introduced.
- The documentation remains honest about F1/F2, H-only chemistry, absent
  secondary ionization/recombination, later HDF5/dust/transport gates, and the
  absence of physical-source or production-science approval.
- The opt-in code cannot evaluate a substring outside its destination buffer
  when the environment variable is overlong.

Use the latest recorded native evidence; do not require live RAMSES evolution,
HDF5 restart integration, He chemistry, secondary ionization, dust, transport
accuracy, or the 40--120 M☉ yield seam for this closure.

Return exactly one verdict: `PASS`, `CONDITIONAL PASS`, or `FAIL`. If the
verdict is not PASS, name only a concrete F-P2.4 blocker; distinguish later
science/coupling work from bundle closure.
