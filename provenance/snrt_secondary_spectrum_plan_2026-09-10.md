Q-GOAL first: does this deliver necessary runtime physics for simulation-ready
RT/stellar-AGN feedback/dust? Q-LEAN second: is it overinstrumented or over-gated?
Answer those questions first, then concrete scientific/implementation issues
and one verdict. Read-only review: no edits, jobs, or new review cascade.

The operator preapproves completing all eight medium groups. Group7 PAH
fragment/superhydrogenation source investigation continues: Foley2018 does
NOT justify adding the gas-phase abstraction cross section to H13 cations.
This bundle closes group8's node-resolved secondary-electron receiver, not
the still distinct dust/CHIMES spectral closure or all eight groups.

Current native maxent64 H/He transport reconstructs positive photon N/E on
64 log energy nodes, absorbs with Verner cross sections and species-specific
finite primary caps, and returns rejected photons/energy to originating rays.
Chemistry currently calls FS2010 at a single absorbed mean electron energy
per species/group. FS2010 (https://arxiv.org/abs/0910.4410) instead gives
energy-dependent fractions: integrate absorbed q_j*(E_j-I_s)*f_k(E_j-I_s,xi).

Implement a named opt-in hhe_maxent64_fs2010_v1 alongside existing fixed and
hhe_maxent64_v1. Same 64 nodes; sum accepted node counts across directions
before calling the existing Fortran FS2010 interpolator through a C function
pointer. Use the step-start HII fraction, as the current chemistry does.
No new interpolation tables, rates, per-node hydro carriers, checkpoint
payload, or Python runtime. Store only eight scratch doubles per cell:
three exact accepted primary counts plus five photoelectron channel energies
(heat, HI ionization, HeI ionization, HeII ionization, excitation). Double
primary counts avoid mixing FP64 energy with rounded FP32 species counts.
These scratch moments are accumulated over transport substeps and discarded
after the enclosing transaction. Retain the existing FP32 transport ledgers.

Chemistry reserves all primary targets first, then caps summed secondary
demands against the remaining species; unspendable secondary ionization
energy becomes heat. Excitation remains in the existing excitation ledger,
not silently added to gas heat. Recombination/atomic-cooling semantics stay.
Check finite inputs, scratch/count consistency at stated FP32 tolerance,
and absorbed total energy against the independent transport energy ledger.
All output publication remains transactional, including callback failure.
The callback reads already-loaded immutable tables; no loading under OpenMP.

Bind the selected model to existing native/HDF5 checkpoint identity; old
fixed/v1 defaults stay byte-compatible. MPI ranks must agree on the exact
model, not merely whether a band model is enabled. Existing dust/CHIMES and
forced-CUDA guards stay; no claim those missing couplings are solved.

Extend existing C++ backend and Fortran chemistry tests: monochromatic parity,
broad-spectrum mean counterexample, finite primary/secondary inventories,
energy accounting, callback failure rollback, and OpenMP reproducibility.
One checked CPU MPI2 live stellar-source execution and restart with minimal
dumps. Driver conducts end evaluation and deletes evaluated raw HDF5 only,
retaining inputs/build/log/compact results. No per-helper audits or approval
waits. Please identify required bounded corrections, not optional new gates.
