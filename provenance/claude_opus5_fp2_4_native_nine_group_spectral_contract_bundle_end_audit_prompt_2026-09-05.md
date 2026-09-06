# Claude Opus 5 bundle-end audit request — F-P2.4 native nine-group SNRT contract

You are the single bundled Claude Opus 5 end auditor for F-P2.4 in
`/gpfs/kjhan/LRD_JWST`, repository `kjhan0606/LagRamses`.

This is a read-only algorithm, implementation, physics-boundary, and
provenance audit. Do not edit files, launch simulations/jobs, or use external
network tools. Inspect the repository using only permitted read/search tools.
The final purpose is production-ready and publication-ready high-level
radiation transport plus stellar/AGN feedback and dust physics in the
LagRamses/RAMSES workflow. Python is not the deliverable here; it is only an
offline reference/provenance layer. Judge the native Fortran/CUDA wiring
directly.

The operator pre-approved this coherent bundle. Give one decisive end-of-
bundle verdict, not a step-by-step approval. Treat unrelated pre-existing
dirty-tree changes as out of scope unless this bundle overwrites or depends on
them incorrectly.

Read these files first:

- `provenance/fp2_4_native_nine_group_spectral_contract_bundle_plan_2026-09-05.md`
- `provenance/fp2_4_native_nine_group_spectral_contract_bundle_implementation_evidence_2026-09-05.md`
- `simulation/snrt/SNRT_NATIVE_GROUP_CONTRACT.md`
- `simulation/snrt/P4_SOURCE_LEDGER.md`
- `simulation/snrt/P4_TRANSPORT_CONSERVATION_VALIDATION.md`
- `simulation/snrt/config/snrt_group_contract_reference_control_v1.nml`
- `patch/lagRamses/snrt_spectral_contract.f90`
- `patch/lagRamses/snrt_spectral_contract_smoke.f90`
- `patch/lagRamses/snrt_state.f90`
- `patch/lagRamses/snrt_ramses_driver.f90`
- `patch/lagRamses/snrt_nlte_coupling.f90`
- `patch/lagRamses/snrt_nlte_coupling_smoke.f90`
- `patch/lagRamses/snrt_cuda_multigroup_smoke.f90`
- `patch/lagRamses/snrt_transport_step.f90`
- `bin/Makefile`
- `simulation/snrt/data/p4_pilot_agn_photon_ledger.json`
- `simulation/snrt/config/p0_photon_group_edges_ev.txt`

Assess the complete bundle, especially:

1. Does the native state actually use nine groups and the exact ten canonical
   boundaries, including the `[2000,10000] eV` (2--10 keV) group, throughout
   state allocation, transport dimensions, source transaction, CUDA ABI, and
   the optional checkpoint format?
2. Is the namelist loader genuinely fail-closed? Check missing files, malformed
   identity, candidate-only status, edge-array mismatch, edge SHA256 mismatch,
   interval-convention mismatch, out-of-band means, negative fractions, and
   H/He threshold opacity/excess rules. Is the runtime gate consistent with the
   stated distinction between reference control, candidate, and approved
   production status?
3. Is the source-energy wiring correct? Verify that group photon energy is the
   resolved radiated energy times the declared group fraction and that removal
   of the former `0.5` factor is physically justified by the ledger contract.
   Check whether the represented/unrepresented fraction is reported without
   pretending that omitted energy was transported.
4. Is the H photoheating boundary consistent with the absorber-weighted excess
   energy in the source closure? Does the retained mean-energy fallback remain
   safe for monochromatic tests? Identify any remaining H/He chemistry
   limitation without treating the validated-but-inactive He arrays as live
   physics.
5. Is checkpoint version 4 identity binding correctly serialized and checked?
   Does it prevent a state from being read under a different source/edge
   contract? Is the stated absence of RAMSES HDF5 backup/restore call-site
   integration accurately scoped rather than a hidden production claim?
6. Are the Fortran native smoke, 9-group CUDA smoke, production `mpiifx` SNRT
   object compile, and full RAMSES link enough evidence for this engineering
   bundle? Separate direct evidence from claims that still require a live
   hydro/RAMSES run or physical-source approval.
7. Look for algorithmic or wiring errors that Python-side tests would miss:
   array ordering, unit conversions, state/checkpoint dimensions, source
   fraction semantics, stale four-group assumptions, unsafe fallback behavior,
   and any double counting at the driver boundary.

Return exactly one decisive verdict at the top: `PASS`, `CONDITIONAL PASS`, or
`FAIL`. Follow it with severity-ranked findings and file/line references,
evidence supporting each finding, and an explicit disposition for every F-P2.4
acceptance gate. Distinguish a blocker for this native bundle from a later
high-level task such as live stellar photons, He chemistry, dust/radiation
pressure, SNRT HDF5 restart integration, or the 40--120 M☉ physical-yield
seam. Do not invent physical-source approval or production science validation.
