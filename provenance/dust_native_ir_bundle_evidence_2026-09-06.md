# DUST-5 native operator evidence

Project `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`; base `1a2c3c2`.
Operator preapproved continued work. Fable CONDITIONAL APPROVE; native
implementation and focused checks complete. Opus CONDITIONAL PASS, with
both documentation-only conditions applied; no code repair was requested.

## Delivered code, not live activation

- `patch/lagRamses/snrt_dust_ir.f90`: native FP64 Planck-table admission,
  bath-relative weak-excess emission, log-temperature inversion, reciprocal
  equal-width cell-set transport, analytic absorption/source, six-face vacuum
  loss, local re-emission iteration and commit-only-on-success API.
- `patch/lagRamses/snrt_dust_ir_smoke.f90`: actual native module caller with
  physical-opacity fixture, rollback, zero and weak-source checks.
- Existing `tests/dust_ir_transport.py --native` creates one temporary plain
  fixture and compares compiled results to JAX. No new Python test file,
  validator, sidecar schema, HDF5 reader, namelist or gate framework.
- Only our object-list addition and dependency line are intended in the
  already-dirty Makefile. Unrelated existing build/AGN/stellar changes remain
  untouched. `make -n -C bin SNRT=1 USE_CUDA=1 snrt_dust_ir.o` resolves the
  actual new source through the production build rule.

The module is included in the SNRT build graph, but **the live RAMSES driver
does not call it**. No primary-photon, chemistry, gas or momentum fields are
modified. This is native operator compilation/differential evidence, not
full production link, live-runtime, AMR/MPI/restart or source-science approval.

## Executed verification

`JAX_PLATFORMS=cpu .venv/bin/python tests/dust_ir_transport.py --native`

Real GNU gfortran 13.2.0 and Intel ifx 2025.3, isolated temporary build/module
directories under the project. GNU uses O2 and runtime array/bounds checks;
Intel uses production O3/OpenMP and the current SNRT preprocessing defines
(NDIM/NPRE/NVAR/NVECTOR/NENER/SOLVER/LONGINT/QUADHILBERT/particle-potential/
HYDRO_CUDA/SNRT/FFTW/stellar-enrichment). Module has no RAMSES stub dependency.
The define list exercises the compiler flag path only: this source has no
preprocessor conditionals, so those defines add no physics/branch coverage.

Numerical fixture: DUST-4's 136-frequency, S4, 2^3 manufactured cube,
dx=1e18 cm, nH*relative_dust=1e6 cm^-3, full c, T_CMB=13.1 K and identical
dense-raw-grid primary heating equivalent to 20 K without IR reabsorption.
Peak box tau approximately .759. Fortran computes its own Planck powers from
the supplied frequency quadrature; it does not import JAX's thermal powers.
Both solves use 1e-9 local/global tolerance; differential criterion is 1e-8.
The predeclared tighter-stop fallback was not needed.

| Relative difference versus JAX | GNU | Intel |
|---|---:|---:|
| Stored energy | 1.1102e-15 | 6.6613e-16 |
| Escaped energy | 0 | 4.4409e-16 |
| Reprocessed energy | 3.3307e-16 | 2.2204e-16 |
| Maximum cell temperature | 4.4409e-16 | 0 |

Both native balances: 7.407e-10. Frequency-integrated cell energy and emitted
photon arrays satisfy the 1e-8 checks as well. Three forced failures
(one-iteration nonconvergence, CFL excess, broken neighbor reciprocity)
prove bit-identical caller energy, temperature, emitted-photon and scalar
diagnostics; the check compares field bits, not uninitialized derived-type
padding. Zero source remains exactly zero; a 1e-42 source remains positive.
No numerical thresholds were relaxed.

The first fixture parse exposed a Fortran text-format detail: exponents with
three digits omit the E under plain ES25.16. The smoke output now explicitly
uses ES26.16E3; no physics correction was required. The full native check
then passed both compilers, and passed again with the complete production
preprocessing define set recorded above.

Opus end record: `opus5_dust_native_ir_end_audit_2026-09-06.md`. Conditions
were documentation-only: clarify compiler flag-path scope and finite
iteration-residual budgeting (7.407e-10 versus 1e-9 tolerance, not machine-
precision closure). Native code and the successful test results are unchanged.

`JAX_PLATFORMS=cpu .venv/bin/python tests/dust_ir_transport.py`

Default compiler-free gray/spectral suite passed unchanged, including DUST-4
analytic, weak-source, frequency and timestep comparisons. No new parameter
matrix, broad infrastructure gate or simulation job was run. `git diff
--check` passed. Temporary fixtures/binaries were automatically cleaned.

## Open handoff

Separate native energy operator now exists; persistent live IR state, its
frequency resolution/memory layout, primary dust absorption mapping and
conservative primary/gas/force transaction remain subsequent integration.
No silent dust-to-metal/depletion/evolution or physical mixture approval.
RSLA energy inventory is still not full-c LTE energy; no live force/gas
exchange should use it without a defined coupling derivation. The module's
reciprocal equal-width graph is not a coarse/fine/MPI boundary adapter.
