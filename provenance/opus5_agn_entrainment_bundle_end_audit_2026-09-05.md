# Opus 5 — AGN entrainment bundle-end audit

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Requested/reported reviewer: `claude-opus-5`
Session: `cc6b2ccd-53bf-42ec-9b91-4a77fa47a16d`
Duration: 443.6 seconds
Read-only Read/Grep/Glob; no subagents or permission denials; exit 0.
Verdict: **PASS** for the approved bundle, no required in-scope repairs.
Driver follow-up disposition: `agn_entrainment_bundle_evidence_2026-09-05.md`.
This is not full production/publication or live-MPI qualification.

## Verbatim reviewer result

## Bundle-end audit — AGN entrainment bundle (`/gpfs/kjhan/LRD_JWST`, 2026-09-05)

Read-only inspection of `agn_feedback_deposition.f90`, `agn_feedback_deposition_smoke.f90`, `sink_particle.kjhan.f90` (`AGN_feedback`/`average_AGN`/`AGN_blast`), `bin/Makefile`, `read_hydro_params.f90` layout, `run_fp15_agn_efficiency.sh`, plus the legacy `patch/cuRamses` sibling for regression provenance. (Note: `Write` is disabled in this session, so this verdict is delivered inline rather than as a plan file.)

### Verdict: **PASS**

No in-scope defect found. Five nonblocking follow-ups below.

### Event closure (verified in source, not just claimed)

- **Mass.** `agn_withdraw_cell` scales `row(1:4)` and `row(fields)` by `(1-ratio)`, `ratio=(m_load/V)/ρ`, `m_load=min(request, 0.25ρV)` — removed mass is exactly `m_load`. `agn_jet_delta` gives `drho_i = 0.5 m_load w_i/W±`, so `Σ drho_i V_i = m_load` for the same membership set used to build `W±` in `average_AGN` (`Tpsy_norm += weights*vol_loc`, same `agn_jet_geometry` call, same `rmax`). Closure is exact by construction, not by tolerance.
- **Vector momentum.** Withdrawal removes `m_load·v_donor`; deposition returns `drho1·(v_d+u a) + drho2·(v_d−u a)`, integrating to `m_load·v_donor` because each lobe integrates to `m_load/2` independently. This is the whole point of the separate normalization and it is implemented correctly.
- **Energy.** Withdrawal removes only `½m_load|v_d|²` (`staged(5)=internal+kinetic*(1-ratio)`), i.e. the approved cold-loading convention. Return is `½(drho1|v+|²+drho2|v−|²) + (1−f_ek)EAGN/vol_gas`, with `uBlast=sqrt(2 f_ek EAGN/mAGN)` using the *actual capped* `mAGN`. Sum = `½m_load|v_d|² + EAGN`; net = `EAGN`, with the cap overflow deferred to `Esave` via `agn_deposit_cell`. Correct.
- **Lobe overlap second moment.** `agn_jet_delta` sums per-lobe kinetic terms, so the midplane cell receives more energy than the bulk KE of its net momentum; `agn_deposit_cell`'s `energy_delta >= kinetic_input - tol` guard is satisfied and the excess thermalizes — exactly the approved "midplane thermalizes opposed streams" behaviour. Smoke line 46 asserts it.
- **Species.** Withdrawal exports `fractions = row(fields)/ρ` and scales donor fields by the same `(1-ratio)`, so donor abundances are invariant (asserted, smoke:101); deposition adds `drho*fractions`, integrating to `frac*m_load`. Total metal is a separate slot, not double-counted.

### Source model vs. gas advection

`AGN_feedback` selects `nelements=nelt`, or `n_stellar_elements=11` under `use_channel_resolved_feedback()`, and never touches `runtime_field_map` or `active_element`. This is the correct distinction: `stellar_ramses_runtime` zeroes `element_index` for inactive elements (production side), while AGN transports all *stored* fields. `agn_scalar_map` rejects `<=5`, `>nvar`, duplicates and reserved `idelay/ivirial/ixion/isgs`; `NENER/=0` rejected before any mutation. I confirmed the production layout arithmetic in `read_hydro_params.f90` (`nener=0, metal, delayed_cooling` → `imetal=6, idelay=7, ichem=8`, elements 8..18, `NVAR=18`), which is precisely the case asserted in the smoke (`agn_scalar_map(18,6,8,11,[7,0,0,0])`). Reservoirs are untouched and asserted untouched (`gas(7,:)` invariant).

### Ownership, sequencing, MPI/OpenMP

- `agn_contains_donor` is half-open in cell-centre coordinates, giving unique ownership across level boundaries too — a genuine improvement over the legacy closed `<= dx_loc/2` test (which could tie two cells).
- Ownership/geometry consensus (`MPI_ALLREDUCE` on `geometry_send`, `owner_send`, then `MPI_MAX` on `donor_error`) **precedes** the sequential withdrawal loop; withdrawal reads the then-current `uold` row, so shared donors compound correctly (smoke:160-163).
- Payload propagation uses the existing `nsink`-indexed `MPI_SUM` with nonowner zeros, and `AGN_blast` consumes `vloadAGN`/`cloadAGN`/`mAGN` after the broadcast — remote receivers do get the real donor payload. All five collectives are unconditional and identically ordered on every rank; the routine is entered with `nAGN=0` too (early returns are on globally uniform `hydro/ndim/nsink`).
- OpenMP: donor hits are recorded under `critical(agn_donor_location)`; cell mutation is grid-partitioned so `uold` writes are thread-disjoint; `deposition_error` uses `reduction(max:)`; deferred energy uses threadprivate `TEsaveAGN` reduced afterwards. Fallbacks are serial as planned.
- Fallback classification: `jet_fallback = ... .and. any(psy_norm<=0)` is computed from the *globally reduced* sums, all distributed deposition is skipped for those events (`AGN_blast:7085`), and the donor loop no longer requires `vol_gas==0` (`7140`). No half event; no double deposit (resolved events `cycle` the fallback loop). The SNRT/legacy exclusivity guard and strict invalid-receiver/collective-exit policy are unchanged.

### Evidence honesty

Adequate. The evidence explicitly disclaims full-link/live evolution, calls the pack/sum/unpack test serial arithmetic, states the MPI wiring is compiled + statically checked only, states cold-loading donor heating (4/3 per withdrawal, compounding, outside `T2maxAGN`), untouched delayed/virial/ionization/SGS reservoirs, no dust, and — importantly — the lack of qualification for *general overlapping mass-weighted thermal/replay sources while other jets change receiver density*. That limit is real: `average_AGN` withdraws donor mass after `vol_gas/mass_gas` were accumulated, so a thermal or replay bubble overlapping another AGN's donor cell will not integrate to exactly `EAGN`. It is pre-existing ordering, correctly named rather than hidden.

### Follow-ups (nonblocking; none required before this bundle is closed)

1. **New hard-abort path.** `owner_all/=1` → `clean_stop`. The pre-existing `3*rmax` AGN bin search can miss the donor when the sink's leaf is several levels coarser than `nlevelmax` (`dx_loc/2 > 3 rmax`); legacy silently skipped, this aborts. Fail-loud is acceptable, but add it to the evidence limits and print sink id/level in the message.
2. **`use_sgs`.** `isgs = ichem+1` (`read_hydro_params:457`), so any map with ≥2 element fields overlaps SGS and AGN now refuses the configuration. That is the approved "layout error" rejection, but it is a new refusal for a previously-runnable combination — one evidence line.
3. **Dead serial branch.** `AGN_blast`'s `#ifndef _OPENMP` path uses `nthreads`/`PEsaveAGN`, which only exist under `_OPENMP`; harmless because the Makefile always sets `-qopenmp` and the file uses `nthreads` unguarded elsewhere, but it is misleading and should be deleted or annotated.
4. **Stack sizing.** `load_send/load_all` are `(4+nscalar)*nsink` automatic arrays (12 MB at `nsink=1e5`); fine now, revisit if sink counts grow.
5. **Wording.** "Two synthetic rank contributions" is one payload plus one zero payload. It faithfully mirrors owner+zeros, but say so.

No new physics, validator, harness or infrastructure work is requested, and nothing here reopens HDF5/CPU-box/checkpoint or parked AGB items. The bundle is closed as PASS; full linked/live and production qualification remains future work as the evidence already states.
