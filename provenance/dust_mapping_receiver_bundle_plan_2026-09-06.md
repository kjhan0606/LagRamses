# DUST-9: source-bound dust mapping and persistent thermal receiver

Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)

This is the next native dust bundle after DUST-8.  DUST-8 proved that the
CUDA fourth-species ledgers reach the FP64 RAMSES trial boundary, but it
intentionally kept `ZERO_SCAFFOLD`: no cell dust abundance, persistent dust
thermal state, or dust receiver existed.  This bundle closes that native
interface without pretending that the current RAMSES hydro state already
contains a dedicated dust fluid.

## Scope

1. Validate the source-bound native opacity arrays after the upstream sidecar
   loader has checked the JSON source identity, source-table hash, and group
   edge hash.  The native boundary requires the arrays, binding status, source
   identity, and all three hashes; it does not parse JSON or recompute hashes.
2. Map a cell abundance only from explicit caller-owned
   `metallicity_solar * dust_to_metal`.  A direct abundance field remains a
   valid upstream input, but no implicit redshift law, depletion law, solar
   normalization, or use of legacy `kappa_IR` is introduced here.
3. Prepare cell/group dust optical depths from explicit `n_H`, path length,
   mapped abundance, and opacity per H.  The result is FP64 and caller-owned.
4. Add a transactional persistent thermal receiver.  It consumes the dust
   absorbed-photon ledger and group mean energies, advances caller-owned dust
   thermal energy and temperature using an explicitly supplied volumetric
   heat capacity, and commits only after all cells pass validation.  The
   constant-capacity closure is an interface contract, not a claim that the
   current RAMSES build has a validated grain heat-capacity table.

## Safety and physical boundaries

- The receiver enforces `absorbed_energy = staged_energy - old_energy` and
  rejects positive absorption for a zero-abundance cell.
- No H/He inventory is changed by the dust receiver.  No gas `uold`, momentum,
  IR photon field, scattering field, or RAMSES restart payload is changed.
- The current live driver remains `ZERO_SCAFFOLD` because `hydro_commons` has
  no dedicated dust state and the native sidecar loader/thermal table is not
  yet a live RAMSES input contract.  The new module is linked into the SNRT
  production graph so the boundary is source-to-binary tested, but it is not
  silently called from the live driver.
- Dust momentum, IR re-emission, grain-size distributions, destruction and
  growth, restart/migration, and cosmological production qualification remain
  later gates.

## Evidence and exit condition

One native smoke covers valid/invalid source binding, abundance mapping,
optical-depth construction, thermal staging, exact state nonmutation on
failure, energy closure, zero-dust behavior, and commit admission.  GNU and
Intel compilers are used when available.  The SNRT production link is rerun
with the new module in its module graph.  The bundle exits as a native
interface PASS, conditional for live/publication dust physics until a
dedicated RAMSES dust state and validated thermal/opacity loader are admitted.
