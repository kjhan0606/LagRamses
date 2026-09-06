# Opus 5 DUST-6 bundle-end audit

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Bundle: native primary-photon dust absorption and heating contract
- Auditor/model: Claude Opus 5 / `claude-opus-5`
- Mode: read-only static audit; no edits, jobs, subagents or web
- Verdict: **CONDITIONAL PASS**

The short `opus-5` alias was not present in the installed Claude model
catalog, so the audit was retried with the full model ID `claude-opus-5`.
The full-ID audit completed successfully. It found no blocking defect in the
bounded native reference boundary: the units, optical-depth attribution,
finite H/He versus non-depleting dust distinction, closure algebra, failure
nonmutation, zero-dust arithmetic branch, heating conversion and build graph
were judged sound. The bundle was judged appropriately scoped and not
over-instrumented. No live CUDA/RAMSES or physical dust-mixture approval was
claimed.

## Findings and dispositions

The audit requested the following before a future CUDA/adapter bundle:

- anchor the inner H/He inventory tolerance to the pre-partition inventory;
- keep an H/He callee numerical remainder in the `unassigned` ledger rather
  than reclassifying it as dust heating or returned photons;
- clamp the small positive-roundoff dust direct term and remove the absolute
  one floor from the outer closure scale;
- add a sequential multi-group test sharing one H/He reservoir; and
- test positive raw removal with zero total optical depth.

All five bounded repairs were applied in the current worktree. The
zero-dust statement was narrowed to the non-saturating regression regime. In
the saturating zero-dust regime, the existing H/He routine's
`err_inventory` behavior and this adapter's explicit returned-photon behavior
are documented as intentionally different. The saturation transfer is also
documented as a bounded sub-step approximation with the audit's noted mild
high bias in dust heating; it is not presented as a dust-mixture result.

## Verification after audit

`simulation/snrt/.venv/bin/python simulation/snrt/tests/dust_ir_transport.py
--native` passed after the repairs for GNU gfortran 13.2 and Intel ifx
2025.3. The smoke includes the new shared three-group reservoir and
zero-total-optical-depth cases. The original DUST-5 differential remains at
machine precision (largest reported relative energy error `1.11e-15` for
gfortran and `6.66e-16` for ifx; temperature error at most `4.44e-16`). The
default non-native DUST-5 regression also passed.

## Remaining conditions / handoff

This conditional pass is only for the native reference contract. The next
bundle must separately prove a fourth-species CUDA ABI, FP32-versus-FP64
allocation and photon conservation, live transaction wiring, AMR/MPI/restart
semantics, and production-scale behavior. Dust opacity-mixture provenance,
depletion/evolution, photoelectric gas heating, IR scattering, persistent
spectral state and cosmological/publication qualification remain deferred.
