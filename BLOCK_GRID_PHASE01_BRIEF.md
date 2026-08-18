# Phase 0-1 work order: index API abstraction (layout unchanged)

Companion to BLOCK_GRID_MAJOR_PLAN.md. Scope is exactly the plan's Phase 0
and Phase 1: introduce the index API and convert every direct index
arithmetic site to it, while the memory layout stays legacy and every result
stays bitwise identical. No block layout, no dynamic growth, no I/O change.

Written 2026-08-18 against main 36fbf36. Census on the 131 sources the
tracked bin/Makefile actually compiles (VPATH resolved): **59 files carry 414
index-arithmetic sites** (367 iskip assignments, 25 direct three-term forms,
21 reverse mappings, 1 mod form). Another ~300 ngridmax mentions are capacity
semantics (declarations, allocations, loop bounds, ncoarse+twotondim*ngridmax)
and are OUT of scope for Phase 1.

## Division of labour

- Phase 0 baselines: session driver (me), before the first conversion chunk.
- API module + conversions: Codex, ONE CHUNK PER INVOCATION (whole-tree jobs
  time out; measured 2026-08-07).
- After every chunk: driver rebuilds serial and runs the bitwise gate.
- Module design audit: Fable, once, when its budget refills.

## Phase 0 artifacts (driver, prerequisite)

1. Reference binary from current main, tracked bin/, serial
   `make HDF5=1 USE_FFTW=1`, sha256 recorded.
2. Reference outputs: PBH smoke A/B (pbh_prod/smoke, ~3 min) with output
   directories preserved, plus one 64^3 cosmo+hydro+AMR run to z=5 that
   exercises refinement (the smoke stops before refinement and would miss
   layout bugs in refined paths).
3. Gate script: byte-compare all output_*/ payload files (amr_*, hydro_*,
   part_*, info_*) and the filtered log lines (PBHDIAG, PBHCACHE, SGS_DT,
   Main step, mdm/mgas) between reference and candidate. Timing lines are
   excluded. PASS = zero differing bytes in payloads.

## The API module (Codex chunk 0)

New file `patch/lagRamses/amr_index.f90`, module `amr_index`. The names
cell_index/cell_grid/cell_child from the plan COLLIDE with local variables in
sink_particle.kjhan.f90:7038 and clump_finder.f90:561; use these instead
(verified collision-free across the tree):

```fortran
module amr_index
  use amr_commons, only: ncoarse, ngridmax
  implicit none
contains
  pure elemental integer function icell_of(igrid, ichild)
    integer, intent(in) :: igrid, ichild
    icell_of = ncoarse + (ichild-1)*ngridmax + igrid
  end function
  pure elemental integer function igrid_of(icell)
    integer, intent(in) :: icell
    igrid_of = mod(icell - ncoarse - 1, ngridmax) + 1
  end function
  pure elemental integer function ichild_of(icell)
    integer, intent(in) :: icell
    ichild_of = (icell - ncoarse - 1)/ngridmax + 1
  end function
  pure elemental logical function is_coarse_cell(icell)
    integer, intent(in) :: icell
    is_coarse_cell = icell <= ncoarse
  end function
end module
```

Rules for the module:
- Confirm where ncoarse/ngridmax are declared (amr_commons vs amr_parameters
  in THIS tree) and use the right module; do not duplicate the variables.
- Phase 1 bodies are exactly the legacy arithmetic above; nothing else.
- Guarded bounds checks under `#ifdef AMR_INDEX_CHECK` (1<=ichild<=twotondim,
  1<=igrid<=ngridmax, icell>ncoarse for the reverse maps); the production
  build must compile to the bare arithmetic.
- Add `amr_index.o` to bin/Makefile's module list immediately after
  `amr_commons.kjhan.o`, and a dependency line `amr_index.o: amr_commons.kjhan.o`.
- LONGINT caveat: icell fits default integer today (ncell < 2^31 at current
  ngridmax); keep default integer, do NOT widen.

## Conversion recipe (Codex chunks A..G)

Per file, mechanical, no reformatting, no logic changes:

1. `iskip = ncoarse + (ind-1)*ngridmax` followed by `... iskip + igrid`:
   DELETE the iskip assignment where it is only used for cell indexing, and
   replace every `iskip+X` use with `icell_of(X, ind)`. If iskip is also
   used for something else in the same scope (rare), keep it and convert only
   the cell-index uses. The compiler re-hoists the invariant; bitwise identity
   is the check that nothing changed.
2. Direct `ncoarse+(ind-1)*ngridmax+igrid` -> `icell_of(igrid, ind)`.
3. Reverse child `(icell-ncoarse-1)/ngridmax+1` -> `ichild_of(icell)`.
4. Reverse grid `icell-ncoarse-(ind-1)*ngridmax` -> `igrid_of(icell)`
   ONLY when ind is the child of that same icell; otherwise leave and flag.
5. `icell <= ncoarse` / `> ncoarse` tests used to mean "is coarse cell"
   MAY be converted to is_coarse_cell(icell) but this is optional polish;
   skip when in doubt.
6. DO NOT touch capacity uses: declarations, allocate bounds,
   `ncoarse+twotondim*ngridmax`, loop bounds over 1:ngridmax, MPI buffer
   sizing, headl/taill/numbl indexing, or anything in comments.
7. Add `use amr_index` (only: list) to each converted procedure or module.
8. End of chunk: for each touched object run `cd bin && make <obj>.o` and
   report compile status per file. Do not run make clean.

Two variant idioms found during the pilot (2026-08-18), add them to every
chunk's checklist:

- **Reverse child WITHOUT the -1**: `(icell-ncoarse)/ngridmax + 1` (seen in
  bisection.f90). This is a latent off-by-one in the ORIGINAL code: it returns
  child+1 whenever igrid==ngridmax, and the paired grid recovery then goes
  nonpositive. Convert it to ichild_of anyway and RECORD the site in the
  commit message as a semantic fix; the bitwise gate is unaffected because
  the last grid slot is never populated at our capacities (and bisection runs
  only under ordering='bisection', which no gate run uses).
- **Reverse grid with swapped operands**: `icell - ncoarse - ngridmax*(ind-1)`
  (the census regex only matched `(ind-1)*ngridmax`). Same quantity; convert
  to igrid_of when ind is the child of that same icell.

The census regexes in the driver's completeness check must cover both
variants.

Chunk-B lessons (2026-08-18), binding on every later chunk:

- When a converted variable's declaration is removed, its mentions in OMP
  `private` clauses must be removed in the SAME operation, keyed to that
  exact variable and unit. Never run a general "strip undeclared names from
  private lists" sweep afterwards: names in private lists may be
  host/module-associated arrays or `type()` variables a naive declaration
  scan cannot see, and stripping one (e.g. a GPU state struct) silently turns
  it shared and creates a race. The safe repair when in doubt is to take the
  pre-conversion OMP line and drop only names with zero non-OMP references
  left in the unit.
- OMP directives start with `!` and must be exempted from any comment filter,
  or the clause cleanup silently stops firing.
- A third no-minus-1 reverse-child was found at multigrid_fine_coarse:698 and
  fixed via ichild_of, same as bisection. Three independent copies of the
  same off-by-one is the strongest argument yet for the API.

Known trap for reviewers: with the legacy layout, a WRONG conversion can
still be numerically right in small tests (the aliasing degeneracy discussed
in the plan). Phase 1's shield is bitwise identity on runs that refine plus
the driver-side census re-run: after all chunks, the census regexes must
return ZERO index-arithmetic sites outside amr_index.f90.

## Chunk map (one Codex invocation each, sites from the census)

| chunk | theme | sites | files |
|---|---|---|---|
| A (pilot) | bisection, nbors_utils.kjhan, interpol_hydro.kjhan | ~13 | 3 |
| B | poisson/multigrid | 149 | 11 |
| C | amr core | ~90 | 9 |
| D | fdm | 49 | 4 |
| E | io | 39 | 9 |
| F | hydro | ~31 | 9 |
| G | pm + misc | ~43 | 13 |

The pilot chunk covers all four idioms on small files; the driver validates
the whole pipeline (convert -> build -> bitwise gate) on it before B..G run.

Full per-file census:

| chunk | sites | file | dir | direct | iskip | rev_child | rev_grid | mod |
|---|---|---|---|---|---|---|---|---|
| amr | 23 | virtual_boundaries.kjhan.f90 | cuRamses | 0 | 23 | 0 | 0 | 0 |
| amr | 21 | flag_utils.kjhan.f90 | lagRamses | 0 | 20 | 0 | 0 | 1 |
| amr | 18 | load_balance.kjhan.f90 | lagRamses | 0 | 16 | 2 | 0 | 0 |
| amr | 13 | refine_utils.f90 | cuRamses | 1 | 12 | 0 | 0 | 0 |
| amr | 6 | init_amr.f90 | lagRamses | 3 | 3 | 0 | 0 | 0 |
| amr | 5 | bisection.f90 | cuRamses | 0 | 2 | 3 | 0 | 0 |
| amr | 4 | physical_boundaries.kjhan.f90 | cuRamses | 0 | 4 | 0 | 0 | 0 |
| amr | 3 | nbors_utils.kjhan.f90 | cuRamses | 1 | 2 | 0 | 0 | 0 |
| amr | 2 | amr_step.jaehyun.f90 | lagRamses | 0 | 2 | 0 | 0 | 0 |
| amr | 1 | adaptive_loop.jaehyun.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| fdm | 36 | fdm_step.f90 | lagRamses | 0 | 36 | 0 | 0 | 0 |
| fdm | 9 | fdm_hjm.f90 | lagRamses | 0 | 9 | 0 | 0 | 0 |
| fdm | 2 | light_cone.fdm.f90 | lagRamses | 0 | 2 | 0 | 0 | 0 |
| fdm | 2 | output_fdm.f90 | lagRamses | 0 | 2 | 0 | 0 | 0 |
| hydro | 13 | godunov_fine.kjhan.f90 | cuRamses | 1 | 12 | 0 | 0 | 0 |
| hydro | 5 | interpol_hydro.kjhan.f90 | cuRamses | 0 | 5 | 0 | 0 | 0 |
| hydro | 4 | init_flow_fine.f90 | lagRamses | 0 | 4 | 0 | 0 | 0 |
| hydro | 3 | courant_fine.kjhan.f90 | cuRamses | 0 | 3 | 0 | 0 | 0 |
| hydro | 2 | cooling_fine.kjhan.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| hydro | 2 | diag_eint.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| hydro | 2 | sgs_fine.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| hydro | 2 | synchro_hydro_fine.kjhan.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| hydro | 1 | hydro_boundary.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| hydro | 1 | hydro_flag.kjhan.f90 | lagRamses | 0 | 1 | 0 | 0 | 0 |
| hydro | 1 | init_hydro.f90 | lagRamses | 0 | 1 | 0 | 0 | 0 |
| io | 17 | restore_hdf5.f90 | lagRamses | 7 | 10 | 0 | 0 | 0 |
| io | 7 | backup_hdf5.f90 | lagRamses | 0 | 7 | 0 | 0 | 0 |
| io | 4 | write_screen.f90 | cuRamses | 4 | 0 | 0 | 0 | 0 |
| io | 3 | output_amr.kjhan.f90 | lagRamses | 0 | 3 | 0 | 0 | 0 |
| io | 3 | output_sphere_hydro.f90 | cuRamses | 0 | 3 | 0 | 0 | 0 |
| io | 2 | light_cone.hydro2.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| io | 1 | movie_mod.yonghwi_org.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| io | 1 | output_hydro.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| io | 1 | power_spectrum.f90 | lagRamses | 0 | 1 | 0 | 0 | 0 |
| misc | 3 | boundary_potential.f90 | poisson | 0 | 3 | 0 | 0 | 0 |
| misc | 3 | pbh_evap_fine.f90 | lagRamses | 3 | 0 | 0 | 0 | 0 |
| misc | 2 | morton_hash.f90 | cuRamses | 2 | 0 | 0 | 0 | 0 |
| misc | 2 | sidm_scatter.f90 | cuRamses | 0 | 2 | 0 | 0 | 0 |
| misc | 1 | init_radiation.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| misc | 1 | observe.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| misc | 1 | rad_backup.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| misc | 1 | rad_step.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| pm | 15 | sink_particle.kjhan.f90 | lagRamses | 0 | 15 | 0 | 0 | 0 |
| pm | 4 | clump_finder.f90 | cuRamses | 0 | 4 | 0 | 0 | 0 |
| pm | 4 | feedback.kjhan3.f90 | cuRamses | 0 | 4 | 0 | 0 | 0 |
| pm | 4 | star_formation.kjhan.f90 | cuRamses | 0 | 4 | 0 | 0 | 0 |
| pm | 3 | init_part.f90 | lagRamses | 0 | 3 | 0 | 0 | 0 |
| pm | 1 | particle_tree.kjhan.f90 | lagRamses | 0 | 1 | 0 | 0 | 0 |
| poisson | 45 | force_fine.kjhan.f90 | cuRamses | 3 | 42 | 0 | 0 | 0 |
| poisson | 24 | rho_fine.kjhan.f90 | lagRamses | 0 | 24 | 0 | 0 | 0 |
| poisson | 19 | multigrid_fine_fine.kjhan.f90 | cuRamses | 0 | 13 | 3 | 3 | 0 |
| poisson | 17 | multigrid_coarse.kjhan.f90 | cuRamses | 0 | 17 | 0 | 0 | 0 |
| poisson | 16 | multigrid_fine_commons.f90 | lagRamses | 0 | 12 | 2 | 2 | 0 |
| poisson | 11 | multigrid_fine_coarse.kjhan.f90 | cuRamses | 0 | 5 | 3 | 3 | 0 |
| poisson | 11 | phi_fine_cg.kjhan.f90 | cuRamses | 0 | 11 | 0 | 0 | 0 |
| poisson | 3 | boundary_potential.kjhan.f90 | cuRamses | 0 | 3 | 0 | 0 | 0 |
| poisson | 1 | init_poisson.f90 | lagRamses | 0 | 1 | 0 | 0 | 0 |
| poisson | 1 | interpol_phi.kjhan.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |
| poisson | 1 | output_poisson.f90 | cuRamses | 0 | 1 | 0 | 0 | 0 |

## Acceptance (drives the Phase 1 -> Phase 2 gate)

**Amended 2026-08-18 after chunk B.** Bitwise identity at -O3 is NOT a valid
acceptance criterion for this conversion and was withdrawn. Replacing an
explicit `iskip` invariant with a function call changes what ifx hoists and
how it reassociates, so the generated code differs even when the arithmetic
is identical. Chunk B diverged from the reference at -O3 while being provably
exact: rebuilt at -O0, where inlining, vectorisation and FP reassociation are
all off, reference and candidate were bitwise identical on the refining run
(job 448708, 1h00m). The -O3 divergence pattern is consistent with codegen
and not physics: output_00001 identical before any evolution, the SGS_DT
timestep sequence identical throughout, the first difference the last ulp of
the mcons conservation diagnostic (1.76e-16 -> 3.52e-16), and a maximum
hydro relative difference of 1.9e-8 after ~20 coarse steps of gravitational
amplification.

1. Census returns zero convertible sites outside amr_index.f90.
2. **-O0 bitwise identity** between reference and candidate on the refining
   64^3 run. This is the semantic proof and is mandatory per chunk (~1 h).
3. **-O3 sanity**: identical SGS_DT timestep sequence, and conservation
   diagnostics agreeing to rounding. Payload bitwise identity is NOT required.
4. `-qopt-report` on godunov_fine.kjhan and force_fine.kjhan shows the hot
   loops still vectorised.
5. AMR_INDEX_CHECK debug build passes the same runs with zero assertions.
6. Fable reviews the chunk diff. The -O0 gate proves behavioural equivalence
   on the exercised paths; it cannot prove a site was converted with the
   right meaning where the values happen to coincide, nor find sites the
   census regexes never matched. Those are review questions.
