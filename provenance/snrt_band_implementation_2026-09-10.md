# Native H/He intragroup spectral closure

Bounded group-8 continuation complete: opt-in `SNRT_SPECTRAL_MODEL=hhe_maxent64_v1`
now runs through actual RAMSES source -> N/E transport -> H/He absorption ->
chemistry -> hydro -> MPI/HDF5 restart. Default is `fixed`, with unchanged
fixed-group coefficients and default checkpoint formats. This is not closure
of all approved medium-term physics groups or universal production approval.

## Single plan review and disposition

[Fable review](snrt_band_fable_2026-09-10.txt) answered Q-GOAL yes and Q-LEAN
mostly lean, approving bounded changes to the [plan](snrt_band_plan_2026-09-10.md).
Adopted direct FP64 species-column arguments instead of division by reference
cross sections; explicit FP32 edge tolerance; zero-opacity bypass; accepted
energy into chemistry; disjoint restart identity without new per-cell state;
and a non-dust integration build. Driver performed end evaluation, as directed.

Did **not** adopt a single cap shared by all species in a band: exhausting H
would then veto He absorption despite remaining He. Each species instead has
one cap constant across the band's directions/nodes. Returned counts and
energies use that SAME species cap. Native regression explicitly proves He
still absorbs after H exhaustion. Band ordering retains the existing ascending
group inventory policy; this is not an exact time-dependent competition solve.

## Physics and actual wiring

- State remains FP32 directional photon CODE density `N`, plus FP64 correction
  `shift`; actual energy is `Eref*N+shift`, in photon CODE density times eV.
- Reconstruct a positive **discrete maximum-entropy** spectrum from N and E:
  64 logarithmic energy nodes per canonical band, trapezoidal dE prior, and
  exponential tilt `p_j proportional w_j exp(beta*x_j)`. Solve the monotone
  mean-energy equation. It is a stated two-moment ansatz, not unique recovery
  of the input SED or a multi-frequency photon tensor stored per cell.
- Evaluate the ground-state H I/He I/He II Verner1996 cross sections at nodes.
  Primary source: [author's fits](https://www.pa.uky.edu/~verner/photo.html),
  [Table 1](https://www.pa.uky.edu/~verner/dima/photo//photo.dat), first three
  rows; coefficients and rounded thresholds match existing `primordial.py`.
- Driver passes `n_species_code * scale_nH * c_eff * dt` directly, also from
  the updated midpoint neutral state on nonlinear retries. The prepared
  adapter divides it by the transport substep count. Fixed reference tau is
  still the iteration/convergence proxy, NOT the spectral absorption opacity.
- The existing conservative upwind stencil transports N and E, including
  ghost/coarse fields. Node-dependent attenuation partitions H/He absorption;
  finite atom inventories are shared between bands and substeps. The survivor
  spectrum hardens. Absorbed energy/species count minus the appropriate
  threshold supplies the chemistry photoelectron energy, no old fixed excess.
- Existing FS2010 secondary ionization uses the **mean** accepted photoelectron
  energy per species/band. It is not node-integrated nonlinear secondary
  deposition. No source-population or dust/network spectral equivalence claim.
- Runtime `auto` selects existing OpenMP paired-energy transport for this mode;
  explicit CUDA remains unsupported. No spectral CUDA kernel added or claimed.
  MPI startup requires the same spectral selector on every rank.

### Numerical boundaries found and repaired in integration

FP32 N / FP64 E admits roundoff-level edge excursions. Reject outside a band
by more than `8*epsilon(FP32)`; inside that tolerance use the endpoint limit,
preserving E and changing reconstructed N only within that rounding tolerance.
Upper endpoints use left-limit opacity except the final closed interval.
Band4's upper H edge must not ionize H; explicit regression checks this.

The BPASS reference table has extraordinarily small X-ray tails. Exact FP64
source energy paired with a rounded/flush-to-zero FP32 photon increment can
create N=0,E>0, correctly rejected by the old paired source validation. For
this new mode, `quantize_source` charges the **representable FP32 increment's**
energy. The absolute source-energy rounding summed over directions/bands
must be <= `8*epsilon(FP32)` of total source energy, or the entire source
transaction rejects. Tables are unchanged; this is a numerical source
discretization, not a physical tail cutoff or exact unrounded source ledger.
Default fixed and existing moving-dust source arithmetic remain unchanged.
Tests cover negligible-tail admission and rejection when the entire source
would be lost. Tiny additions to an existing large N still accumulate in E.

Optically thick extinction exposed a second consistency issue: independent
`exp(-tau)` left a tiny energy-only survivor after `-expm1(-tau)` had already
rounded absorption to 1. Use `transmit=1-loss` from the SAME rounded loss.
When the FP64 absorption ledger has absorbed the whole packet, no spurious
subnormal survivor remains. This is roundoff-consistent extinction, not a
temperature floor or timestep reduction. Extremely small otherwise physically
non-negligible states outside FP32 representability still fail closed.

## Restart and defaults

No new per-cell tensor. Native fixed files retain v7 (v6 number-only migration
still allowed); band files use v8 with the SAME N+shift layout. The version
identifies the complete `hhe_maxent64_v1` closure including numerical policies.
Native readers reject changing model both ways before state mutation.
HDF5 band versions add 8 to current fixed format IDs and require `band_model`;
the stellar-source fixture uses format15 and cell width1444 (80 directions).
Fixed formats1--8 remain unchanged; no automatic fixed-to-band migration.
Bad moments reject at transport and checkpoint payload admission.

`DUST_LIVE=1` or `SNRT_CHIMES` builds explicitly reject band mode. There are
no justified corresponding subband dust/network opacities in this change.
No namelist parameters were added; mkrun/generator defaults were not changed.
Makefile VPATH order is unchanged. The older dust build also recompiles.

## Evidence and cost

Work root `.rt-band.3bDWXk/`; existing native smokes extended, no new framework.
Checked Intel MPI/ifx build:

```
make -f ../bin/Makefile -j1 NVECTOR=32 NENER=1 NVAR=19 SNRT=1 \
  DUST_LIVE=0 HDF5=1 USE_FFTW=0 FDMDEBUG=1 EXEC=ramses_band
```

Executable `ramses_band3d`, SHA256
`f0d8dcd13bcc9121901a42ba9ae04795f703d158506e7d7f1621407d302facb4`.
Working-tree baseline `80265d5a432053bfcd29d5b6773e1feafbc350d3-dirty`;
this turn does not commit or push prior accumulated work.

- Native `backend.log`: reconstruction/endpoints/Verner/finite inventory,
  non-starvation, full extinction, zero-opacity transport, invalid rollback
  pass; existing paired-energy and fixed/dust tests still pass. CUDA comparison
  is SKIP in this CPU build, not GPU spectral evidence.
- Representative 64-vs-256 node absorption N/E moments: maximum relative
  difference `1.026809734521e-3` (0.103%). This is numerical quadrature evidence
  for the sampled bands/columns, not a universal physical accuracy bound.
- Native conservative N ledger max relative error `4.345956261531e-9`;
  E ledger `8.881784197001e-16`. Example mean energy hardens
  16.3475 -> 19.717429866 eV with H absorption.
- `source.log`: exact original paired-source parity/rollback plus bounded
  quantized-source tail tests pass. Native fixed/band checkpoint smokes pass.
- Effective NMLs `completed/physical.nml`, `completed-restart/physical.nml`:
  4^3 grid, 80 directions, MPI2/OMP2, noncosmo periodic hydro, channel stellar
  returns, BPASS **independent-population reference** radiation, advective CR,
  H/He photoionization/heating/recombination. No sinks/AGN, live dust, CHIMES,
  external UV or separate atomic/metal cooling. The inherited dust-specific
  CIE fixture was rejected before evolution, then intentionally removed from
  the non-dust test; no admission rule was weakened.
- Source `run.env.sh`, then `SNRT_BACKEND=auto mpirun -np 2 ../ramses_band3d physical.nml`.
  The source file sets OpenMP2, reference/secondary/stellar/yield contracts;
  auto chooses OpenMP for this paired-energy implementation.
- Four steps finish at code time ~0.187520; wall time14.6798s. Nonzero-source
  nonlinear solves use 20,19,15 iterations. Existing timing labels overlap
  (the cooling wrapper includes RT); do not sum them as independent costs.
  Reconstruction work scales with cells * directions * active bands * 64
  nodes * moment-root iterations, each transport substep and nonlinear retry.
  It is substantially more expensive than fixed coefficients and stays opt-in.
- `evaluation.log`: **82 nonheader datasets exactly identical** at step4
  after step2 restart: AMR12, coarse3, domain1, gravity8, hydro38, particles17,
  SNRT3. All numeric datasets finite; physical clocks/SNRT attributes exact.
  Whole HDF5 hashes differ due to header/build/run metadata, not physical state.
- Density min `9.998245813144558e-4`, gas internal energy min
  `6.655696740385254e-7` CODE units. HII fraction0.76319--0.80748;
  HeII0.96069--0.96595; HeIII0.0003534--0.0003666. Positive N/E and band bounds
  pass; group5 mean24.3243 eV vs injection17.6629 eV demonstrates live hardening.
- `reject-fixed` is reader-only: actual spectral HDF5 -> fixed mode fails
  before a new dump. Dust320 checked build passes; its native admission test
  rejects band mode explicitly (`.rt-angular.TwRcDH/band-admission.log`).

Driver verdict: this bounded native H/He comparison capability is complete.
Full dust/CHIMES spectral coupling, independent matched stellar injection
spectral moments, node-resolved secondary deposition, cosmological redshifting,
spectral CUDA and general multi-level AMR qualification are NOT established
by this run. Do not relabel the remaining source/PAH groups as complete.

## Evaluated raw-output cleanup

Fresh/restart schedule: `noutput=1,aout=2,tout=1e30` outside this noncosmo run;
`foutput=2,fbackup=1000000`; dumps at steps2/4 and one restart step4 dump.
Measured each HDF5 file1,338,072bytes. Reader rejection uses nstepmax0 and
foutput/fbackup1000000. Pre-launch free space ~62TiB, estimate <=9MiB.
One intermediate failed evolution also retained a step2 dump until diagnosis.

After completed test AND driver evaluation, the exact four files below were
permanently removed (total **5,352,288bytes**). All
NML/env/log/build/evaluation/native-smoke files and output headers are retained. Recreate
raw fields by rerunning the retained fresh input; deleted checkpoints cannot
be used directly for restart. No unrelated outputs or production data touched.

Paths relative to `.rt-band.3bDWXk/`, pre-delete SHA256:

| File | SHA256 |
| --- | --- |
| `completed/output_00001/data_00001.h5` | `44eb31df1d7a2d31509ae91a894fd7a39601cfc63e872459d29a980a67824028` |
| `completed/output_00002/data_00002.h5` | `d5b52a29f17f52806a37443e933e5280ffe7ecb00d5c12dd9f9de58378c60a5f` |
| `completed-restart/output_00002/data_00002.h5` | `adfbfab1769382c0871f2826513f9b8b35418d5f1920805de9ac2ebeed71159f` |
| `final/output_00001/data_00001.h5` | `7e3cf08e8e0092cf5d14bb5130d6717720e43786928090344b9127b12d886628` |
