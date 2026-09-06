# Claude Opus 5 follow-up audit — F-P2.5 native H/He thermochemistry

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
Audit mode: read-only native Fortran/CUDA source and recorded evidence review

## Verdict

`CONDITIONAL PASS`

The first-audit HIGH findings and the claimed MODERATE findings were confirmed
closed in the current source. No further source remediation is required. The
bundle can close after the three record-only corrections listed below.

## Direct source findings

- HIGH-1 (closed): explicit absolute/relative inventory tolerances,
  tolerance-sized clipping, an optional `unassigned_absorption_code` ledger,
  and a `0.99995` FP32 guard protect the live species-aware CUDA path.
- HIGH-2 (closed): the species-aware CUDA ABI uses the Fortran
  `(leaf,group,species)` storage as CUDA `[species][group][leaf]`; positive
  opacity masks are applied; one shared H I/He I/He II inventory is consumed
  across groups and substeps; and the RAMSES driver uses this ABI. The scalar
  compatibility entry point is not on the production driver path.
- MODERATE-3 (closed): raw `fion+fheat+fexc` closure is checked over all
  258x14 table entries before the loaded flag is set.
- MODERATE-4 (closed): the hot-temperature He II smoke separates radiative and
  dielectronic case-B terms and checks a non-negligible dielectronic share.
- MODERATE-5 (closed): a saturated-H-II smoke confirms unavailable secondary
  ionization energy is routed to heat.
- MODERATE-6 (closed): authoritative H/He state and the neutral-H mirror are
  double precision.
- MODERATE-7 (closed): checkpoint version 6 serializes and validates spectral
  and FS2010 identities before state mutation.
- LOW-10 (closed): the SNRT driver prerequisites are present in the Makefile.

## Residual findings

- N1 (recorded limitation): tolerance-clipped or post-transport
  partition/chemistry-failed absorption has no thermal or ionizing receiver;
  the current `unassigned` quantity is diagnostic only.
- N2 (recorded limitation): the ascending-group partition is feasible only
  under the nested species-eligibility invariant enforced by
  `validate_species_table`.
- N3 (recorded limitation): ascending group order is a physics policy that
  prioritizes lower-numbered groups in saturated cells and can affect the
  species-resolved ionization structure.
- N4--N10 (non-blocking follow-ups): derive guard-band margins, add an old
  checkpoint rejection smoke, initialize unused CUDA ghost slots, align
  negative-intensity handling, remove or test dead compatibility surfaces,
  and avoid duplicated group-count constants.

## Acceptance gates

| Gate | Disposition |
|---|---|
| FS2010 interpolation and 100 eV continuity | MET |
| Three-species opacity, partition, and state ownership | MET |
| H/He simplex, non-negative inventory, recombination conservation | MET |
| Photoelectron closure and heating-only RAMSES receiver | MET |
| Version-6 checkpoint identity and double-precision state | MET |
| Native-source/native-smoke audit rather than Python-only evidence | MET |

The audit directly inspected source. Build, CUDA smoke, table-manifest, and
recorded numerical transcripts remain evidence claims rather than executions
performed by the auditor.

## Scope boundary

F-P2.5 does not approve a physical AGN or stellar SED, a live RT+feedback
evolution, global implicit opacity/chemistry convergence, cooling receiver
integration, dust/radiation pressure, HDF5 restart integration, publication
scale convergence, or the 40--120 M_sun yield seam.

## Closure conditions applied

1. `SNRT_NATIVE_GROUP_CONTRACT.md` now records the `unassigned` receiver
   limitation and the ascending-group priority bias.
2. The same document records the nested-eligibility invariant required by the
   ascending greedy partition.
3. This bundle's remediation transcript now includes the complete spectral
   smoke output, including `SNRT_SPECTRAL_CONTRACT_OK`.
