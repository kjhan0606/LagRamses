# Conservative atomic photo-chemistry and native CHIMES hot-cell connection

Continuation of medium group 8, started 2026-09-10 and evaluated 2026-09-11.
This delivers a finite-inventory photo operator and an actual native
photo/nonradiative CHIMES cell connection. It does NOT enable general live
CHIMES+spectral RT: cold molecular spectra, grain competition and live
checkpoint wiring remain unfinished. No new namelist option or defaults
were introduced. Both existing live spectral/CHIMES guards remain intact.

## Implementation

- `snrt_chimes_photo.cpp`: CVODE 5.8 BDF/SPGMR evolves actual CHIMES157
  species, one survival fraction per occupied absorbing spectral node, and
  eight budget counters together. Directional N/E is reconstructed before
  angular summation; common attenuation at a fixed node preserves different
  ray spectra. There is no dense node-by-direction chemical matrix.
- Each shell reaction consumes one primary photon, its actual energy and
  its reactant; creates the mapped product and KM93 electron multiplicity.
  Species and photons use the same reaction rate, not an after-the-fact cap.
  The node/shell primary electron receives E-binding. Existing FS2010 raw
  interpolation and finite atomic-target limiting supply heat, additional
  HI/HeI/HeII ionization and excitation. No molecular electron-degradation
  model, Auger cascade or fluorescence spectrum is invented.
- `snrt_chimes_band_photo_step` returns heat, three secondary ionization
  energies, excitation, unresolved binding/cascade reservoir, total absorbed
  energy (all eV/cm3), and primary photon count (cm^-3). Reservoir energy is
  not identified with a completed species-resolved chemical binding ledger
  and is never thermalized silently. This operator has no recombination,
  nonradiative chemistry/cooling, dust or gas-temperature update by itself.
- `chimes_cell_band_hot_atomic` in the actual Fortran CHIMES module then
  adds photo heat using the changed total particle count and calls the
  existing charged CHIMES receiver with ZERO photons. Thus recombination,
  collisions and cooling occur once and do not reapply photo rates. This is
  an explicit first-order photo-then-dark split. Its ninth budget entry is
  the signed nonradiative gas thermal change, not a predicted escape spectrum.
- That combined native comparison requires dust-free, zero solid charge,
  T>100000 K before/after the photo/dark steps and no molecular carriers.
  This is the existing CHIMES hot-network threshold, not permission to delete
  H2/H2+ when it crosses the boundary. A crossing or molecular input rejects
  the entire staged split. H-/C-/O- belong to the upstream atomic network and
  are not incorrectly excluded. The stand-alone photo operator separately
  retains H2 ionization into H2+; it does not perform dissociation.
- All outputs stage privately. Species/fraction/counter nonnegativity is
  imposed by CVODE constraints and checked before publication. Nuclear and
  charge errors must be <=1e-8 (nuclear relative, charge per H); photon and
  energy budget residuals must be <=1e-7 of incoming totals, with the
  thermal/secondary/reservoir partition closing to1e-8. Solver relative
  tolerance is1e-8. Failures do not publish partial state; no final physical
  clipping or renormalization manufactures acceptance.

The shared internal bank type fixes cross-translation-unit handle ownership.
Make links the photo object and installed SPGMR library under CHIMES=1 only.
VPATH and external CHIMES ABI5 library are unchanged. No generator/GUI edit
is needed because there is no new namelist/runtime selector.

## Single selective review and disposition

[Fable plan review](snrt_chimes_spectral_receiver_fable_2026-09-10.txt),
read-only Claude CLI model `fable`, exit0, covered the operator and following
integration. Q-GOAL: useful operator, not yet full runtime. Q-LEAN: proceed
with trims. No routine end audit or new operator approval was requested.

Adopted all three executable conditions: shared bank/build wiring, explicit
solver constraints and acceptance tolerances, and H2/hot-network gating.
The reported missing source/private type reflected a mid-edit read; the
final actual Makefile build verifies both are resolved.

Also adopted the FS cache trim: copy the existing 258x14x6 raw grid and
store energy brackets/weights, not84 values at every shell/node energy.
If every node in every band is occupied, there are73900 unique above-threshold energies:
the old raw-sample cache alone would be49.66MB per cell; grid plus brackets
and five evaluated fractions is about4.31MB, excluding reaction/ray buffers
and container overhead. No runtime HDF5 reading occurs inside the solve.

Retained occupied-node survival variables rather than adding absorber-column
unknowns and a dense solver. Photon number plus integrated captures is a
linear invariant of this formulation. Matrix-free convergence is explicitly
tested through photon exhaustion below; no universal stiffness/performance
qualification is claimed. Column-integral optimization remains optional, not
a new completion gate. New tests extend the existing smoke; no framework.

## Numerical issue found and corrected

At dt=1e12s in a photon-starved HI cell, CV_NORMAL initially returned a tiny
negative final state despite accepted-step constraints (status51). The
solver could step past the requested time and interpolate back. Setting
CVodeSetStopTime to the actual dimensionless final time resolves the case
without clipping or looser tolerances. Both CV_SUCCESS and CV_TSTOP_RETURN
are handled, and the reached time is checked. The failed and corrected logs
are retained as `native-thick.log` and `native-stop.log`.

## Evidence and reproduction

Private build: `/gpfs/kjhan/LRD_JWST/.chimes-photo.oaOwoR`.
Same pinned1,036,386-byte atomic input bank as the preceding foundation;
no new physical table or radiation source was generated.

```sh
make -f ../bin/Makefile -j4 snrt_thermochemistry_smoke \
  SNRT=1 DUST_LIVE=1 NENER=1 HDF5=1 CHIMES=1 USE_FFTW=0 \
  CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.medium-pah-charge.wBiAfM/chimes \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install
source environment.sh
./snrt_thermochemistry_smoke
```

Actual Intel mpiifx/icpx, OpenMP4, SUNDIALS5.8, double-precision CHIMES:

- `native-final.log`:142 PASS, including40 added operator/split checks;
  `SNRT_NATIVE_THERMOCHEMISTRY_OK`.
- Finite HI target and photon exhaustion, H- affinity/charge, H2+ retention,
  metal/Auger products, shell-wise FS agreement with the existing native
  reference, nuclei/charge/energy budgets, zero identity and private rollback.
- Concurrent CVODE photo cells match serial species/N/E/budgets bitwise.
- Actual hot-cell split at dt=1e8,5e7,2.5e7s: consecutive state differences
  6.13522619e-6 then3.05642201e-6; temperature differences5.96678113K then
  2.97233283K. This is first-order convergence on a monochromatic test,
  not a claim of second-order splitting or general continuum accuracy.
- Zero-radiation split matches the existing nonradiative receiver. Cold
  input, preexisting H2, and a post-photo crossing of the hot-network
  boundary reject without publishing any part of the trial state.
- `legacy-chimes-final.log`:86 PASS with the spectral bank unset.
- Separate clean CHIMES=0 build `.chimes-photo-legacy.m0sNVx`:40 PASS.

Executable SHA256:
`ea5c576d4136f0b51a8dab4a2fbe88fb5820569079ad9f3268b9974812f0c6c0`.
Photo source SHA256:
`a3026b97556673149eb696a5af8b9a6fbe83bc623fb314cdb8aa0b788f1d8256`.

No RAMSES/MPI evolution was launched. No raw simulation output was generated
or deleted. Physical inputs, native binaries, logs and compact evidence stay.
No commit or push in this continuation.

## Remaining existing connection

Cold molecular spectral absorption/dissociation and its energy routing;
competing grain absorption/heat; native driver and checkpoint identity;
then one combined short live conservation/restart evaluation. Also retain
the existing continuum/threshold quadrature qualification requirement.
The new hot native comparison does not waive any of these or complete the
whole medium-term programme. No new approval wait is inserted.
