# Fable bundle-end audit: DUST-7 fourth-species CUDA dust boundary

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Auditor/model: Fable / `claude-fable-5-1`
- Mode: read-only source audit; no edits, jobs, or web search
- Session: `0ec337fe-90db-4240-be34-b29dfb3f4b0e`
- Audit duration: approximately 233 seconds
- Verdict: **CONDITIONAL PASS**

## Auditor conclusion

Fable found the DUST-7 boundary algorithmically sound, correctly isolated
from the legacy ABI, and backed by actual GPU execution as reported. The
conditions were test-coverage and documentation gaps rather than design
defects. Fable could not execute the smoke or independently verify hashes in
its read-only audit session, so those claims were treated as evidence claims
to be checked in the project record.

The audit confirmed from source that the legacy species kernel and shared
wrapper body contain no DUST-7 references, that the new entry point uses
separate kernels, and that the driver still calls only the legacy species
ABI. It also confirmed the ledger identity, explicit guard-band return,
`expm1f` weak-dust path, dust-only groups, device validation ordering,
architecture override, and nonzero no-device behavior. The replicated
128-lane reduction and zero-dust branch make the bitwise legacy comparison
meaningful. No new validator, dataset, gate, or framework was added; the
size-overflow checks were described as mild and harmless over-engineering.

## Conditions and dispositions

1. Fable required a host-side quantitative check of the direct-dust plus
   finite-excess formula, not closure/sign checks alone. Added to the smoke:
   dust-only `dust=raw` and `returned=0` checks plus the reconstructed
   `dust_direct + excess*(1-exp(-tau_dust))` comparison.
2. Fable required a non-saturating mixed group. Added a second native call
   with a small group-5 photon packet and ample inventory, checking the
   opacity-proportional H/He shares, dust share, and zero returned packet.
   The saturating case also checks reservoir exhaustion within the retained
   FP32 guard-band tolerance and zero H/He assignment in later groups.
3. Fable required a component-total mismatch case independent of the sign
   check. Added an all-nonnegative inconsistent `tau_total` test and verified
   caller state, inventory, and sentinel outputs remain unchanged.
4. Fable required evidence caveats. The implementation evidence now records
   caller photon units, FP32 tolerances, the non-negative/CFL-respecting
   assumption behind the bitwise path and `fmaxf` clamp, the informational
   nature of the `nvidia-smi` device line, and the remaining FP32/FP64
   reconciliation work.

The four conditions were applied in the same bundle and the full native
smoke was rerun successfully with both Intel mpiifx/ifx 2025.3 and GNU
gfortran 13.2. The resulting status remains **conditional** because the
audit does not qualify live RAMSES integration, AMR/MPI/restart behavior,
opacity-mixture science, production cost, or publication convergence.

## Handoff judgement

Fable judged DUST-8 wiring justified: the new ABI emits every ledger needed
by the FP64 RAMSES transaction, the legacy ABI remains untouched, and the
evidence is actual A10 execution rather than compile-only. DUST-8 must
resolve the intentional FP32/FP64 guard/residual differences before using
the CUDA ledger as a production reference and must prove primary
photon/H/He/dust heating closure and transaction rollback across the live
AMR/MPI path.
