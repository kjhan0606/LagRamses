# Native SNRT nine-group contract

## Purpose

The RAMSES SNRT state is now dimensioned for the same nine groups used by the
P0/P4 source ledgers. The canonical boundaries are

```text
0.01, 1.0, 5.6, 11.2, 13.6, 24.59, 54.42, 500, 2000, 10000 eV
```

They use the shared `left_closed_right_open_except_final_closed` convention.
The old native four-group table (`18, 35, 70, 200 eV`) has been removed from
the state module. In particular, the `[2000,10000] eV` hard-X-ray interval is
now an actual transport group rather than being silently omitted.

## Runtime input

The native module is
[`patch/lagRamses/snrt_spectral_contract.f90`](../../patch/lagRamses/snrt_spectral_contract.f90).
It reads a strict Fortran namelist from the path in `SNRT_GROUP_CONTRACT`:

```bash
export SNRT_GROUP_CONTRACT=/gpfs/kjhan/LRD_JWST/simulation/snrt/config/snrt_group_contract_reference_control_v1.nml
# Required only when deliberately running the checked-in reference control.
export SNRT_ALLOW_REFERENCE_CONTROL=1
```

The reference-control file is
[`config/snrt_group_contract_reference_control_v1.nml`](config/snrt_group_contract_reference_control_v1.nml).
It serializes the group means, source energy fractions, H I/He I/He II
group-averaged cross sections, and absorber-weighted photoelectron excess
energies copied from the reviewed P4 pilot ledger. Its status is explicitly
`reference_control`; it is not a physical SED approval.

The process does not parse the JSON ledger. The namelist must carry:

- version `1` and the exact canonical edge array;
- the canonical edge-file SHA256
  `d28f78f1703730c6c0b9a7d183edfe0c5e6337979e737ce002a572b66fc53ff1` and
  the exact interval convention;
- an explicit `fraction_semantics` value, either `intrinsic` or `escaped`;
- finite in-band group representative energies;
- non-negative group fractions with a sum no greater than one;
- zero opacity and zero excess energy below each species threshold;
- positive, bounded opacity and excess energy in every supported threshold
  group;
- a non-empty source identifier, a 64-hex source digest, a 40-hex source
  commit binding, and an approval/control identifier.

Any failure leaves the contract unloaded and causes `SNRT_RT_ENABLE=1` to
return without changing the radiation state. A candidate explicit SED status
may be inspected by the loader but is not runtime-admissible. A
`reference_control` contract additionally requires
`SNRT_ALLOW_REFERENCE_CONTROL=1`; an `approved_production` contract does not.
The resolved-domain source boundary requires `fraction_semantics='escaped'`,
so an intrinsic contract remains inspectable but is blocked from runtime until
an explicit upstream escape conversion is implemented and approved.

## Source and chemistry boundary

The AGN source loop uses the ledger's group fraction directly. The previous
extra `0.5` multiplier was removed: the emitted photon energy in each group is
now exactly the resolved radiated energy times that group's declared fraction.
The driver reports both the represented fraction and the unrepresented
fraction; omitted source energy is not silently put into a neighboring group.
The source fraction is explicitly tagged as escaped versus intrinsic in the
contract. The current resolved-domain injection path accepts only escaped
fractions; the checked-in reference control is tagged `escaped` because its
unresolved nuclear escape factor is one.

The native chemistry boundary accepts the absorber-weighted H I, He I, and
He II photoelectron excess energies from the contract. The live RAMSES path
now computes all three primordial opacities, passes their sum to transport,
partitions the returned absorption against the start-of-step species
inventories, and updates H II/He II/He III with FS2010 secondary channels and
case-B recombination. The production CUDA adapter receives the group/species
optical-depth tensor and consumes one shared H I/He I/He II inventory in
deterministic group order across all transport substeps; redistribution is
masked by positive group opacity. The compatibility H-only routines remain
available for older isolated benchmarks, but are not the production SNRT
driver path.

The F-P2.6 fixed-point residual is the difference between successive
under-relaxed opacity predictors, not an unqualified raw chemistry-map norm.
The raw chemistry result, photon field, group absorption, and heating from the
same final trial are committed together after that predictor residual passes;
the maximum configured cost is 32 prepared transport evaluations per level,
times the transport subcycle count.

## Restart binding

SNRT state checkpoint version 6 records the spectral status, source id,
source digest, source commit binding, approval/control id, edge digest,
interval convention, and fraction semantics. A checkpoint read is rejected
unless the same runtime-admissible contract is already loaded. This prevents a
nine-group state from being interpreted with a different group ordering,
source closure, or intrinsic/escaped convention. The authoritative H II,
He II, He III, and H I mirror fractions are stored in double precision, and
the checkpoint also records the pinned FS2010 source id, upstream commit, and
manifest identity. A missing or mismatched secondary-table identity is
rejected before state mutation. The existing RAMSES HDF5 backup and restore
call sites do not yet invoke these SNRT checkpoint routines; their integration
remains a G5 restart task.

## Native science limitations carried to the next feedback gates

These are explicit science limitations, not claims that the reference-control
wiring is production physics.

1. The emission-side group representative energy is a photon-number-weighted
   mean, while the current H heating path uses the Verner-absorber-weighted
   photoelectron excess. Therefore the per-absorbed-photon energy does not
   close exactly in the number-only state. For the checked-in reference
   control, the emission/heating energies and relative gaps are approximately:

   | group | interval (eV) | emission mean (eV) | 13.6 + H excess (eV) | gap |
   | --- | ---: | ---: | ---: | ---: |
   | 5 | 13.6--24.59 | 17.66 | 16.45 | 7% |
   | 6 | 24.59--54.42 | 34.38 | 30.42 | 12% |
   | 7 | 54.42--500 | 106.64 | 68.87 | 35% |
   | 8 | 500--2000 | 869.63 | 625.81 | 28% |
   | 9 | 2000--10000 | 4023.59 | 2578.50 | 36% |

   The resolution is a later G3/G4 energy-residual or absorber-energy-state
   gate; restoring `mean - 13.6` is not an accepted fix.

2. The native driver now includes the FS2010 secondary-ionization partition,
   case-B H/He update, and the bounded F-P2.6 local opacity/chemistry fixed
   point. The fixed point is a level-local, start-inventory/frozen-temperature
   operator; a global implicit radiation/chemistry solve is not claimed. It
   also keeps recombination radiation, excitation lines, collisional
   ionization, background/Compton cooling, and metal cooling outside the
   RAMSES thermal receiver. Those are separate gates, not silently folded
   into the gas-heating term.

3. Species-inventory partition is currently a deterministic ascending-group
   greedy policy. In a saturated cell, lower-numbered groups claim the shared
   H I/He I/He II inventory first; this can therefore affect the
   species-resolved ionization structure. The greedy is feasible because the
   current spectral contract enforces nested species eligibility through
   `validate_species_table` (species opacity is zero below threshold and
   positive above it). That invariant is a prerequisite of this algorithm and
   must be revisited before adding a non-nested opacity family such as a future
   resonance-line or dust group.

4. `unassigned_absorption_code` is now a production gate inside the F-P2.6
   transaction: an above-tolerance residual fails the collective level and
   rolls back, while the tolerance-sized residual remains visible in the
   global ledger. It is deliberately not silently reassigned. A later
   receiver/energy-closure gate is still required before a live production
   RT+feedback run; the older scalar transport adapters remain compatibility
   diagnostics and are not the production driver path.

## Scope status

This bundle closes the native group-count, source-boundary, local H/He
thermochemistry, and bounded level-local RT/chemistry transaction wiring. It
does not approve the pilot AGN SED, activate SNRT in a production hydro run,
add stellar photons, radiation pressure, dust opacity,
collisional/background/metal cooling, or a global distributed implicit
RT/chemistry solve. CUDA availability and the independent physical-source
gates remain required for any live run.
