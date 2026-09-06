# Codex self end audit: DUST-9 mapping and thermal receiver

- Project: `/gpfs/kjhan/LRD_JWST` (`kjhan0606/LagRamses`)
- Auditor: Codex, independent end-of-bundle source/evidence review
- Date: 2026-09-06
- Verdict: **PASS within declared native-interface scope; CONDITIONAL for live/publication dust physics**

## Findings

1. The source-bound contract is explicit.  The native code requires finite,
   monotonic group edges, nonnegative per-H opacity, in-band representative
   energy, an admitted binding status, a source identity, and syntactically
   valid SHA-256 tokens.  It does not silently accept an unbound opacity table.
2. Cell dust abundance is not inferred from `kappa_IR`, redshift, or a hidden
   depletion law.  The only native derived path is the explicit product of
   caller-provided `metallicity_solar` and `dust_to_metal` arrays.
3. The thermal receiver is transactional at its API boundary: it computes
   staged energy and temperature from old state, rejects positive absorption
   at zero abundance, validates all candidates, and only then copies the
   candidate into persistent caller-owned arrays.  The smoke confirms failed
   calls preserve the persistent state.
4. The production build and existing consolidated SNRT gate pass after adding
   the object to the Makefile.  GNU and Intel focused native execution pass.
5. The live driver correctly remains `ZERO_SCAFFOLD`.  The current RAMSES
   hydro state has a metal passive scalar but no dedicated dust state or
   admitted native sidecar loader/thermal capacity table.  Calling the new
   receiver from that driver now would be an unjustified physical activation,
   so it is intentionally not wired into the live step.

## Residual conditions

- The constant volumetric heat capacity is an explicit receiver contract, not
  an approved grain-material model.
- The native hash arguments are identity tokens supplied by an upstream
  loader; native JSON parsing and file hashing are still outside this bundle.
- A follow-on bundle must provide a dedicated RAMSES dust state, a physical
  temperature-dependent heat-capacity/IR closure, native sidecar admission,
  and restart/MPI semantics before changing `ZERO_SCAFFOLD`.
- Dust momentum, radiation pressure, scattering, re-emission, grain
  evolution, and cosmological production convergence remain deferred.

The DUST-9 engineering objective is met and no unrelated RT, feedback, or
AMR changes were needed.
