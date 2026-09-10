# Native absorbed-spectrum secondary electrons

Driver end evaluation: **PASS for the opt-in H/He node-resolved receiver**.
This closes the node-secondary item in medium group8, not its dust/CHIMES
spectral coupling and not the other medium groups.

Select `SNRT_SPECTRAL_MODEL=hhe_maxent64_fs2010_v1`. Existing `fixed` and
`hhe_maxent64_v1` retain their algorithms; the former remains the default.
No namelist field, hydro carrier, Python runtime or extra restart payload.
Makefile changes only add the actual thermochemistry dependency/link object;
VPATH order is unchanged.

## Physics and wiring

The existing positive maxent64 reconstruction and Verner H/He cross sections
produce direction-summed absorbed node counts. After each species' finite
primary cap, integrate q*(E-I)*f(E-I,xi) using the existing pinned
[Furlanetto--Stoever2010](https://arxiv.org/abs/0910.4410) interpolator through
a read-only C callback. xi remains the step-start HII fraction. Callback
failure rejects the entire staged transport transaction. All tables are
loaded before the OpenMP region.

Eight FP64 scratch entries per cell hold three accepted primary counts and
five energy channels: heat, HI/HeI/HeII ionization and excitation. They are
summed over all transport substeps, consumed by the existing chemistry
routine, and discarded. Primary reservations precede the summed secondary
caps; unspendable ionization energy becomes heat. FP64 counts pair with the
actual energy ledger rather than introducing FP32 count-rounding energy.
FP32 public counts and reconstructed excess are checked within8 FP32 eps;
the live driver separately checks total absorbed energy within2e-12 relative.
Recombination/atomic-cooling branches and source transactions remain intact.

Excitation energy is a local chemistry result treated here as escaping line
energy, NOT a transported line-radiation carrier. The selected FS2010 table
interpolation, endpoint clamping and step-start-ionization approximation
remain; this is not a nonlocal fast-electron cascade or full spectral recovery.
CHIMES/live-dust combinations and forced CUDA remain explicitly unsupported
for this spectral option. The callback does not supply their missing physics.

Exact model indices agree across MPI ranks. The native binary format is9
for this option (8 for old band,7 for fixed); HDF5 adds16 to fixed formats,
yielding23 in this stellar-source fixture. The cell width stays1444. Latched
model strings and source/table identities are checked. Cross-model rejection
is intentional reproducibility policy, not incompatible state dimensions.

## Evidence

Retained artifacts: `.rt-node-secondary.WxSPuz/`.

- Checked Intel/ifx CPU build NVAR19/NENER1, SNRT=1, DUST_LIVE=0, HDF5=1,
  USE_FFTW=0, FDMDEBUG=1. `build.log`/`build-final.log`; final executable
  `ramses_node_secondary3d`, SHA256
  `4d305e09efe96ae57f02a3af7c42d26bc72ff93c907c52b8ec805afb1402edd2`.
- Existing C++ backend tests pass Intel and GNU (GNU undefined-behavior
  sanitizer). `backend-intel.log`, `backend-gnu.log`: independent node
  quadrature, finite primary cap, exact OpenMP1/4 results, callback failure
  rollback, existing fixed/paired transport cases. CUDA comparison skipped
  because this is a CPU build, not claimed as GPU evidence.
- Existing chemistry suite passes Intel and GNU; GNU bounds checking and
  invalid/zero/overflow traps enabled. `chemistry-intel-final.log`,
  `gnu/chemistry-final.log`: monochromatic parity, actual FS2010 broad versus
  mean heating difference6.2704332645%, finite secondary inventory/heat
  rerouting, inconsistent primary/energy and negative-energy rejection.
  This numerical example is not a universal astrophysical correction factor.
- Native checkpoint tests: `checkpoint-new-final.log`, `checkpoint-fixed.log`.
  New-mode exact N/E restart; fixed<->band and mean<->node policy rejection;
  unchanged fixed version7 tests pass. The initial new-mode test expected
  literal version8 and failed that assertion only; updated it to the actual
  selected version and added the cross-secondary-model test before rerunning.
- `live/physical.nml`: MPI2/OMP2, periodic noncosmological4^3 gas, stars,
  BPASSv3 Q/E, H/He RT and advective CR, four steps. AGN/sinks, ordinary
  cooling, dust and CHIMES intentionally off. `environment.sh` retains
  effective source references. Auto backend selects OpenMP. Wall15.721s.
- `restart/physical.nml`: same binary and inputs, output1(step2) -> step4,
  wall7.920s. `evaluate.py`/`evaluation.txt`: all82 non-header physical
  datasets exactly match: AMR12, coarse3, domain1, gravity8, hydro38,
  particles17, SNRT3. SNRT identity and header times/counters match.
  All numeric datasets finite, density/internal energy positive, photon
  means within their bands and zero photons carry no hidden energy.
  Final density minimum0.0009998214592542483, gas internal minimum
  5.953222257966276e-7 code units; HII range0.75935--0.80339.
- No failed RT transaction, HDF5 diagnostic, MG failure or nonzero field-NaN
  counter. The pre-existing SFRD print is NaN in this noncosmological fixture;
  it is not accepted as a finite observable. Physical arrays are checked
  separately. Feedback-active hydro mass/energy diagnostics are not closed-box
  conservation proofs.

Output policy audited before both launches: noutput1/aout2/tout1e30,
foutput2/fbackup1e6. Three raw files total4,161,000bytes; remove only after
this evaluation, retaining logs/inputs/compact evidence. See cleanup record.
No new framework, per-helper audit, approval wait, commit or push.
