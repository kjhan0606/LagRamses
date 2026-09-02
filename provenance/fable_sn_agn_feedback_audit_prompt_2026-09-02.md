# Independent SN/AGN feedback physics and implementation audit

Audit `/gpfs/kjhan/LRD_JWST` read-only with the actual `fable` model. Do not
edit files, create commits, submit jobs, or change external state. The current
repository HEAD is `7e6dab63d87707dc4ee1749f242d3a809191cc00`; unrelated
untracked B3 validation files may be present and are outside this audit.

This executes now the one-time audit originally requested for 2026-09-02
09:00 Asia/Seoul. Inspect the current GPFS source, especially
`patch/lagRamses`, `simulation/snrt`, and the existing provenance. Do not use a
different `/home` checkout as the implementation under review.

Return one overall verdict—`PASS`, `CONDITIONAL PASS`, or `BLOCK`—for whether
the SN/AGN feedback physics and implementation are production-ready. A
development scaffold may be judged useful while the production verdict
remains blocked. Cite file:line evidence and distinguish independently
verified facts from inferences.

## Required scope

1. **Phase-0 yield contract and runtime flow**
   - wind, AGB, SNII, SNIa, and PISN channel ownership;
   - cumulative versus differential release semantics, age/time units,
     interpolation domain, IMF/population normalization, remnant/fallback;
   - mass, tracked/untracked element, energy, and radial-momentum conservation;
   - native Fortran reader/query/deposition path and restart/idempotence.

2. **Re-audit the four earlier Fable source findings**
   - suspected years/Gyr mismatch at the runtime boundary;
   - interval direction and `[age,age+dt]` cumulative increment semantics;
   - possible missing or duplicated `1e51 erg` conversion in legacy/Phase-0
     kinetic-SN paths;
   - three-species versus 11-species metadata, He indexing, `NVAR=17`, generic
     metallicity, and delayed-cooling ownership.

   State separately whether each is closed, still open, or was a false alarm.

3. **Physical yield assets**
   - determine whether the staged AGB/CCSN/Pop-III candidates constitute one
     approved full mass-metallicity-age production grid;
   - check citations, licenses, checksums, source adapters, decay horizon,
     rotation/explosion weighting, wind/terminal ownership, SNIa DTD, and
     PISN/PPISN eligibility;
   - do not treat synthetic fixtures or zero-row fail-closed adapters as
     physical completion.

4. **AGN bookkeeping and feedback**
   - coarse-state diagnostic timing relative to accumulator reset;
   - Bondi/Eddington/inflow/retained-BH-mass convention, radiative efficiency,
     Lbol, thermal/jet energy and momentum;
   - deduplication and restart replay of `(coarse_step,sink_id)`;
   - AGN SED, escape/obscuration, source deposition, and use by SNRT;
   - identify any double counting or unclosed mass/energy/momentum paths.

5. **Connection to RT/dust/hydrodynamics**
   - distinguish post-processed source ledgers from live RAMSES feedback;
   - identify what is implemented, partial, missing, or scientifically
     unapproved for stellar SEDs, AGN emission, dust, radiation pressure,
     thermochemistry, and live RT-hydro coupling;
   - use the current nine-group source closure where relevant, without
     granting field convergence that belongs to the separate B3 gate.

6. **Comparison baseline**
   - assess the stopped output-00011 calculation only as the recorded
     `transitional_feedback_baseline`, not pure legacy feedback;
   - state what comparisons it can and cannot support.

## Required output

- executive verdict and concise scientific rationale;
- table of implemented / partial / missing / unapproved components;
- closure table for the four prior Fable findings;
- severity-ranked physics/algorithm/wiring findings, not ordinary style bugs;
- concrete acceptance tests and artifacts needed for each blocker;
- prioritized P0/P1/P2/P3 implementation sequence;
- explicit list of claims that remain prohibited for production or
  publication.

Read and test as much as needed, but remain strictly read-only. Save no files;
return the full report in the CLI response.
