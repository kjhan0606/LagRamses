# CHIMES atomic spectral connection: native foundation

Subsequent progress: [the conservative photo operator and bounded hot-cell
CHIMES connection](snrt_chimes_photo_implementation_2026-09-11.md) now evolve
finite N/E and species together. The record below describes the preceding
moment-provider implementation, not completion of the live receiver.

## Scope and outcome

Continues the preapproved medium group 8 after the bounded static Fe
spectral connection. This implements the native atomic spectral coefficient
provider, NOT a completed live CHIMES spectral receiver. No new approval
gate, routine audit, live selector, namelist change or simulation was added.
The same existing physical cross sections and maxent128 closure are used;
no new energy-partition prescription is admitted here.

`patch/lagRamses/snrt_chimes_spectrum.cpp/.h` and the ISO C bindings in
`snrt_chimes.f90` provide an explicitly loaded immutable bank. CHIMES=1
links it through the existing Makefile. VPATH is unchanged. The original
grey reaction tables, receiver ABI, checkpoint formats and defaults are
unchanged. The existing CHIMES+spectral rejection remains in force.

## Data and physical meaning

Extended the existing `build_chimes_group_tables.py` with `--band-nodes`;
the default grey branch still requires Leiden data and is not replaced.
The new bank has nine existing photon groups, 128 nodes per group, 311
ordered atomic/H-/H2-ionization reactions and 682 partial shells. The four
categories have 8, 116, 35 and 152 rows. It preserves shell binding energies
and probability-weighted Verner95/96 cross sections with KM93 Auger
branching, instead of collapsing shell energies before a nonlinear future
secondary-electron calculation. This is NOT the 30 empirical molecular
dissociation reactions or resolved Leiden molecular lines.

Reconstruct each direction's N/E separately, then sum its node populations.
For each band/reaction, return number-opacity, photon-energy-opacity and
primary-electron-energy-opacity moments. Multiply by reduced light speed
and target number density to get rates. Do not interpret them as accepted
finite-inventory captures or as net thermal heating. Missing cascade,
fluorescence and molecular energy channels are not fabricated or assigned
to heat. The upper band endpoint uses left-limit opacity except in the
last band, as in the existing H/He spectral kernel.

Admission checks the pinned main identity, exact ordered reaction-map hash,
edges, shape, grid, positive finite cross sections and shell bindings. A
valid-range but wrong-species mapping rejects. Full table digest is returned
for future caller/restart binding. No HDF5 operations occur during cell
moment evaluation. Handle lifetime belongs to the caller: load/free outside
threaded evaluation. Rejected loads and invalid N/E publish no partial output.

## Retained inputs and reproduction

Work directory: `/gpfs/kjhan/LRD_JWST/.chimes-spectral.TrQxAG`.
Existing `.dust-optics-venv` supplies scipy, h5py and mpi4py; default Python
was missing mpi4py. No global environment was changed.

```sh
.dust-optics-venv/bin/python simulation/snrt/tools/build_chimes_group_tables.py \
  --chimes-tools .dust-extension.AOz7mU/chimes-tools \
  --main-data .dust-extension.AOz7mU/chimes-data/chimes_main_data.hdf5 \
  --band-nodes 128 --output <new-directory>
```

The retained `bank-final/manifest.json` records all raw fit, upstream tools,
wrapper and table hashes. The generated HDF5 is 1,036,386 bytes, not a
simulation output. It and the input source tables are retained.

- Main SHA256: `8bde78faacd59249fad5e810cae43311ed03ef09131c62b0a0a5c39e8407cb0a`.
- Table SHA256: `998970ed5bb4cfeda01913a72cc3fc62af8ce2fdb55dc924f2704e3abb5391de`.
- Ordered little-endian int32 reaction-map SHA256:
  `8cf424068bd05d205155e9359d07a84c8abbfb5f880b268cf33becf94148d351`.

Actual Makefile build in that private directory:

```sh
make -f ../bin/Makefile -j4 snrt_thermochemistry_smoke \
  SNRT=1 DUST_LIVE=1 NENER=1 HDF5=1 CHIMES=1 USE_FFTW=0 \
  CHIMES_DIR=/gpfs/kjhan/LRD_JWST/.medium-pah-charge.wBiAfM/chimes \
  SUNDIALS_DIR=/gpfs/kjhan/LRD_JWST/.dust-extension.AOz7mU/sundials-install
source environment.sh
./snrt_thermochemistry_smoke
```

## Evaluation

Intel mpiifx/mpiicpx, CPU, OpenMP4; no RAMSES or MPI simulation launched.
Existing smoke extended, not a new test framework. `native-final.log`:
102 PASS checks including 16 new spectral checks, reporting
`SNRT_NATIVE_THERMOCHEMISTRY_OK`. Tests cover an independent analytic HI
Verner endpoint reference, primary binding-energy subtraction, positive
energy-bounded moments, directional sum versus separate reconstructions,
non-equivalence to merged-mean opacity, OpenMP/serial bitwise equality,
threshold double-assignment exclusion, hard-photon Auger support and
transactional invalid-N/E rejection. `legacy-final.log`: 86 PASS checks
with the new bank unset. Existing neutral/charged CHIMES tests still pass.

Final smoke executable SHA256:
`c61d2b6b7aa69d7bdc9cd5161140c41860d3033a736c26824009bd9b40f47c78`.

`admission.log`: six temporary malformed HDF5 variants (schema, edge, node
shape, reaction species, negative cross section, NaN binding) reject with
handle, counts and identity unchanged. Temporary variants are removed.
Physical inputs, compact results and logs are retained. No raw simulation
output was generated and none required deletion.

## Still required within the existing work stream

Connect per-cell spectral coefficients to evolving finite photon N/E and
species in the native receiver, without global mutable reaction tables or
double-counted dust absorption. Retain shell-wise secondary partition and
track ionization, excitation and escaping energy separately from heat.
Supply a declared molecular spectral treatment; no zero-rate placeholder
or grey-rate relabeling. Qualify continuum quadrature/threshold accuracy
before live admission, then perform the combined native/live conservation
and restart comparison. These are existing coupling tasks, not newly
invented approval gates. This foundation alone does not complete group 8
or establish production/publication readiness.
