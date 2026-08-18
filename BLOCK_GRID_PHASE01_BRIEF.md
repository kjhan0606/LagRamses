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

## REVISION 2026-08-18: switch the API from module functions to fpp macros

Fable's review of 20d0893+45ecb5d found zero wrong-meaning conversions but
established that **`icell_of` is never inlined**: FFLAGS carries no `-ipo`, so
cross-file inlining cannot happen, and phi_fine_cg's disassembly contains 31
`callq amr_index_mp_icell_of_@PLT`. An opaque call in the Poisson hot loop
kills vectorisation, which both explains the -O3 divergence (FP reduction
order changes, last-ulp drift amplified by gravity, with -O0 identity proving
the arithmetic exact) and is a real per-step regression. Acceptance criterion
"hot loops still vectorised" currently FAILS.

Decision: **implement the mapping as fpp macros in a shared include**, keeping
the module only for the checked build.

    patch/lagRamses/amr_index.h        (new, fpp include)
    #ifdef AMR_INDEX_CHECK
    #  define ICELL_OF(g,c)   icell_of(g,c)
    #  define IGRID_OF(k)     igrid_of(k)
    #  define ICHILD_OF(k)    ichild_of(k)
    #else
    #  define ICELL_OF(g,c)   (ncoarse+((c)-1)*ngridmax+(g))
    #  define IGRID_OF(k)     (mod((k)-ncoarse-1,ngridmax)+1)
    #  define ICHILD_OF(k)    (((k)-ncoarse-1)/ngridmax+1)
    #endif

Rationale: textual substitution reproduces the original expression exactly, so
codegen is unchanged, the performance regression disappears, and **-O3 bitwise
identity is restored as the acceptance gate** — which retires the 1 h/chunk
-O0 run. Type safety is not lost outright: the AMR_INDEX_CHECK build still
routes through the existing pure elemental functions with their assertions.
`-qipo` was rejected as it perturbs codegen across the whole build.

Conversion rules for the macro form: parenthesise every macro argument (done
above), keep `use amr_index` ONLY in units compiled for the checked build,
and add `#include "amr_index.h"` after the last `use` and before
`implicit none` (the placement rule below still applies). Every source using
the macros must be compiled with `-fpp`, which this Makefile already does
globally.

Rework scope: the 15 files already converted must be re-expressed in macro
form, and the 43 sites listed below folded in at the same time.

## The 43 sites the v1 census missed (Fable, reproduced by tests/phase1_census.py)

Two blind spots: array-element child operands `(hhh(idim,1,ind)-1)*ngridmax`
and `&`-continuation-split expressions. `tests/phase1_census.py` is the
corrected matcher — it joins continuations and accepts array-element
operands, and it independently reproduces Fable's count of 43 in 7 files.
Tree-wide remaining after chunks 0/A/B: **365 sites in 58 files**.

| file | count | note |
|---|---|---|
| force_fine.kjhan.f90 | 24 | ih_left/ih_right hoists + 6 inline scalar_gr |
| multigrid_fine_fine.kjhan.f90 | 7 | 198,221,414,446 are neighbour indices in the SAME loops whose central cells were converted — the worst Phase-2 landmine, no compile error would flag the mixed layout |
| rho_fine.kjhan.f90 | 4 | `indp(j,ind)=...+igrid(j,ind)`; `igrid` can be 0 on masked paths, so a naive ICELL_OF trips AMR_INDEX_CHECK on legitimate runs |
| nbors_utils.kjhan.f90 | 4 | reverse pairs; v1 census scored this file 0 |
| multigrid_fine_coarse.kjhan.f90 | 2 | reverse pair |
| boundary_potential.kjhan.f90 | 1 | `ind_ref(ind)` operand |
| multigrid_fine_commons.f90 | 1 | continuation-split neighbour site |

## Corrections to earlier claims in this brief

- The no-minus-1 corner is **reachable**, contrary to what the pilot commit
  message says. The free list is FIFO: init builds 1..ngridmax
  (init_amr.f90:314-317), allocation pops the head, `kill_grid` appends at
  the tail (refine_utils.f90:1138-1146). Slot ngridmax is handed out once
  cumulative grid creations exceed the list length — a churn condition,
  reachable in long production runs, invisible to short gates. The
  conversions to ICHILD_OF at bisection and multigrid_fine_coarse:698 fix a
  real latent `sink_per_grid(0)` / `lookup_mg(0)` access.
- Phase 2 must separately inventory the **flag2 packing**
  (`flag2(icell)/ngridmax`, `flag2+ngridmax*scan_flag`, `flag2(ngridmax+i)`
  scratch at multigrid_fine_fine 188,400,778-780 and
  multigrid_fine_commons 1562,1618-1970,2034). These use ngridmax as a
  packing base, not a cell stride, and must NOT be converted.
- Pre-existing, not part of this work: multigrid_coarse.kjhan.f90:85,243,675
  have `do i=1,boundary(ibound,ilevel)%ngrid` outside `do ibound=...`, so the
  bound reads an undefined ibound. Harmless in periodic runs because the
  inner loop never executes, but `-check bounds` builds abort there.

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

**Confirmed 2026-08-18**: with the macro form, the -O3 gate passes bitwise
(job 449202). The regression really was the un-inlined call, so the criteria
revert to the original, cheap form:

1. Census returns zero convertible sites (tests/phase1_census.py; the
   definition site amr_index.f90 is exempt).
2. **-O3 bitwise identity** between reference and candidate on the refining
   64^3 run, ~8 min per chunk. The -O0 route is retained only as a fallback
   diagnostic if a future chunk fails this.
4. `-qopt-report` on godunov_fine.kjhan and force_fine.kjhan shows the hot
   loops still vectorised.
5. AMR_INDEX_CHECK debug build passes the same runs with zero assertions.
6. Fable reviews the chunk diff. The -O0 gate proves behavioural equivalence
   on the exercised paths; it cannot prove a site was converted with the
   right meaning where the values happen to coincide, nor find sites the
   census regexes never matched. Those are review questions.
