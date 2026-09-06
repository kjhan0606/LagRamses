# DUST-5: native conservative dust IR operator

Project `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`.
Operator preapproved continuation absent material blockers on 2026-09-06.
DUST-4 end verdict is CONDITIONAL PASS; bounded repairs and complete test
reruns are finished. Fable returned CONDITIONAL APPROVE; implementation
proceeds under the operator's preapproval. Base commit: `1a2c3c2`.

## Purpose and bounded deliverable

Move the actual dust source/reprocessing/energy-transport algorithm into
native Fortran used by the lagRamses build. Do not spend another bundle on
Python instrumentation. The final project goal remains physically defensible
production/publication RT/feedback/dust, not indefinitely expanding pilots.

The native primary state is fixed nine-group photon number (FP32); the dust
secondary field is energy with a separate frequency grid. Reusing primary
species-atom budgets or silently reinterpreting that state as energy is wrong.
Provide a distinct FP64 native operator on caller-owned arrays, with a clear
future driver boundary and a real compiled end-to-end numerical comparison.
No automatic live activation, dust abundance assumption or new namelist.

## Implementation bundle

1. `patch/lagRamses/snrt_dust_ir.f90`: native table admission from frequency
   nodes/weights/sigma and strictly increasing thermal nodes that include the
   exact CMB. Compute Planck powers in native code, their total, weak-excess
   slope integration and temperature inversion. Frequency weights enter
   emission/energy only. No JSON/HDF5 loader or new approval schema.
2. Native explicit-array cell-set step: same-level equal-width cells with
   reciprocal six-face neighbor indices and vacuum boundary sentinel zero.
   Conservative S_N upwind transport, analytic absorption/source response,
   bounded local equilibrium reprocessing from the SAME old radiation state,
   boundary escape, absorbed energy, emitted photons, temperatures and
   source/balance diagnostics. It must preserve all caller-owned state on
   failure and return a distinct error for table/state/CFL/nonconvergence.
   The cell-set API is not an AMR/MPI qualification claim. No primary photon
   fields, chemistry, gas thermal or momentum fields are changed.
3. Add the module to the SNRT build graph without rewriting other dirty
   Makefile/native changes. Expose a callable native API, not a dead script.
   Use isolated compiler directories and no active calculation outputs.
4. One native Fortran smoke driver, called from the existing dust IR test's
   new opt-in `--native` option (default remains compiler-free). Exchange
   only a compact numerical fixture in a
   temporary directory. The existing physical quadrature supplies input;
   Fortran independently builds the thermal powers and advances the radiation.
   Test raw Draine 136-node moderate-depth case against JAX using same full-c
   S4/dt/source/geometry; compare temperature, frequency-integrated radiation,
   total absorption and boundary escape (1e-8 relative). Native closure
   tolerance 1e-9. Also exercise zero/weak source and failure rollback.
   Existing Python spectral tests remain the reference, not production code.

Fable conditions: weights normalized to ONE (not the native primary 4 pi
convention), unit directions, strict thermal/total-power and reciprocal graph
admission, distinct thermal-range error, unchanged arrays on failure.
Mirror the bath-relative slopes, log-T-in-power inversion, quadratic phi
series below 1e-4, old-field escape, half-relaxation and zero-primary inventory
scale. Compare at 1e-8 with both solves at 1e-9; if a difference is caused only
by stop-iterate choice, rerun BOTH at 1e-12 without loosening the comparison.
Use 136-node state only for this test/callable operator; live spectral
resolution/memory layout remains an explicit future design variable.

## Cost and audit discipline

One coherent native physics operator, one smoke driver, one extension of the
existing test; no new validator/manifest framework, HPC jobs or re-audits of
coarse AMR trees/checkpoints/source-to-binary identity. Compile/check this
module with available GNU and Intel compilers in private directories. Do not
claim a full production link/runtime qualification from module compilation.
One Fable plan review and one Opus end review for the whole bundle.

## Deliberately not authorized by this engineering result

No implicit dust-to-metal mapping, depletion, grain evolution or dust-mixture
science approval; those choices remain explicit G4.1 physics prerequisites.
No live force/gas energy deposition without RSLA coupling derivation; no
coarse/fine interface, MPI ghost, restart or accelerated IR claims. Persistent
live IR state and a RAMSES dust-density/source mapper require a subsequent
coupling design. This bundle removes the Python-only algorithm boundary;
it does not itself make the complete dust subsystem production-ready.
