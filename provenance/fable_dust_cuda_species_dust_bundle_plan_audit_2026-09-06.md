# Fable pre-bundle audit: DUST-7 fourth-species CUDA dust boundary

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Base: `ad935dc`
- Auditor/model: Fable / `claude-fable-5-1`
- Mode: read-only plan audit
- Verdict: **CONDITIONAL APPROVE**

## Scope decision

Deferring live RAMSES wiring to DUST-8 is correct. The current driver
re-partitions CUDA totals in FP64 and has no dust term in its optical-depth or
heating transaction. Combining that change with a new FP32 CUDA ABI would
make failures unattributable. DUST-7 must emit every ledger quantity needed by
DUST-8 so the ABI is not cut twice.

## Required conditions adopted

1. Implement a separate CUDA kernel and wrapper. Leave the legacy kernel and
   shared implementation body untouched, and rerun the legacy smoke with its
   previous values.
2. Add a zero-dust branch and compare state, total/group absorption and H/He
   inventories bitwise to the existing species ABI.
3. Preserve the legacy `0.99995` FP32 inventory guard band as returned photon
   remnant, not as dust-eligible physical excess.
4. Use `expm1f` for weak dust optical depth.
5. Return independent raw, H/He, dust, returned and assigned-total group
   ledgers; CUDA `unassigned` is identically zero. Document the memory layout
   and do not feed assigned-total to the host partition as raw input.
6. Test dust-only low-energy groups, which have no H/He opacity in the
   production nine-group contract, and leave H/He inventories untouched.
7. Validate total optical depth against the four component sum with at least
   eight FP32 epsilons relative (or derive it on device). Reject invalid
   inputs before state/output mutation.
8. Make CUDA architecture selection overridable. A missing device must print
   a distinct marker and return nonzero; compile-only is not a bundle pass.
   Record the GPU model in the evidence.

The boundary remains a native transport contract, not dust-mixture approval,
live transaction qualification, AMR/MPI/restart evidence, or production /
publication approval.
