# Opus 5 AGN native coupling bundle-end audit

Date: 2026-09-05 (Asia/Seoul)
Project: `/gpfs/kjhan/LRD_JWST`, `kjhan0606/LagRamses`
Requested/reported audit model: `claude-opus-5`
Session: `87abc57c-c6ef-456b-a7bf-918da6675a9d`
Read-only Read/Grep/Glob; exit 0, no permission denials or subagents.
Duration approximately 702 seconds. Original response is retained below;
the driver's repair disposition is in the implementation evidence record.

## Auditor response (verbatim)

I inspected the changed Fortran directly (deposition module + smoke, `AGN_feedback`/`average_AGN`/`AGN_blast`, the efficiency latch, the `read_params` preflight, the driver enable/source-accounting sections, the Makefile, and the extended runner). Here is the bundle-end verdict.

# Verdict: CONDITIONAL PASS

Every conservation/semantics claim in `agn_native_coupling_bundle_evidence_2026-09-05.md` that I could check against source is **true as written**, and the evidence's "Claim limits" section is honest — it does not present this as live radiation+jet qualification. The three blockers below are all *fail-safety* defects in newly added lines, not errors in the approved physics corrections.

## What I verified as correct

| Claim | Verified at |
|---|---|
| Energy identity `E_old+E_inj = E_new+E_deferred` exactly | `agn_feedback_deposition.f90:52-72`; `deferred=(trial_energy-staged(5))*volume` is the algebraic complement, so closure is exact, not approximate |
| Cap limits only *new* internal energy | `:62-65` — `internal_limit=max(internal_old, …)`, so an already-hot cell keeps `internal_old` and defers all input. Not a cooling operator |
| Volume-weighted jet normalization | `average_AGN` accumulates `psy*vol_loc` (`sink_particle.kjhan.f90:6651`); `agn_jet_delta` divides by that sum (`:84`) ⇒ `Σ dρ_i·V_i = m_load` |
| Kinetic budget matches allocated mass | `uBlast=sqrt(2·f_ekAGN·EAGN/mAGN)` (`:6983`) with `p_gas=(1-f_ekAGN)·EAGN/vol_gas` (`:6984`); summing `kinetic_delta·V` over the cylinder returns `f_ek·EAGN` plus the (declared) bulk/imbalance terms |
| Midpoint semantics | `:92` — `speed²=|v_bulk|²+u_jet²` with `momentum=dρ·v_bulk`: the exact first and second moments of two equal opposed half-lobes |
| Radicand clamp in **both** passes | `:6646` and `:7071`, arithmetically identical |
| Fallback uses `vol_blast`, never the stale `vol_loc` | `:7130` |
| Donor-velocity fallback closes the remove/return pair exactly | capture at `:6713` *before* the mutation at `:6722-6727`; return at `:7126-7128` restores ρ, ρv and ½ρv² bit-for-bit. Metal closes too (`ZAGN·dm` at `:6719` vs `:7159`) |
| Metal committed only after the staged row passes | `deposit_agn_cell:7157-7166` |
| MPI-uniform exclusivity, before either source | `read_params.jaehyun.f90:1690-1715` uses MIN/MAX allreduce on the latch plus rank-uniform namelist `sink`/`sink_AGN`. `read_hydro_params` is called at `:1423`, so `sink_AGN` is populated — the gate is live, not vacuous. Defensive duplicates at `sink_particle.kjhan.f90:6180` and `snrt_ramses_driver.f90:246` are also uniform |
| Coarse-step photon re-emission removed | no `nstep_coarse` key remains in the accounting block; `snrt_ramses_driver.f90:585-596` keys purely on `idsink` |
| Makefile | `agn_feedback_deposition.o` in `MODOBJ` (`bin/Makefile:169`) with direct deps at `:292-293` |
| Runner extension is proportionate | `run_fp15_agn_efficiency.sh:28-31` compiles the *production* module, not an oracle copy |

The evidence's reproduction arithmetic also checks out: the old kernel gave `dρ = m·ψ/Σψ/V_cyl` → 0.5 of 1 for two unit cells, and the old cap branch gave 15+3 vs. the correct 13.

---

## Blockers

**B1 — The inflow cursor is volatile but the ledger it reads is durable.**
`accounted_inflow` starts at `0.0d0` for any sink not previously *owned by this rank* (`snrt_ramses_driver.f90:602`), and for every rank after restart. But `dMBH_coarse`/`dMEd_coarse` are written and restored via HDF5 (`backup_hdf5.f90:792-795`, `restore_hdf5.f90:2238-2241`) and — precisely because legacy feedback is now excluded — are **never reset** (`sink_particle.kjhan.f90:6352` is unreachable when `sink_AGN=.false.`). So `supplied_mass` at `:656` is run-cumulative while the cursor is zero, and `delta_inflow` at `:665` becomes the entire accretion history, emitted in one timestep. This fires on the first RT step after any restart, and on any load-balance remap that moves sink ownership — the cursor on every non-owning rank is permanently 0, so migration is the common case, not the exceptional one.

This is not the declared "durable restart state" follow-up. The correct idiom already exists 20 lines below: `retained_initialized` (`:671-674`) rebases the *retained* cursor on first observation. The inflow cursor needs the same first-observation rebase (seed from `supplied_mass` instead of 0). Fable's condition was "carry the cursor across coarse steps **or fail that mode closed too**"; carrying is done, but neither half holds across the ownership axis.

**B2 — New fail-closed paths in `AGN_blast` are rank-local and will hang, not abort.**
`deposition_error` is strictly per-rank per-cell, yet `:7135-7138` calls `clean_stop` on it; same for `Invalid AGN injected energy` at `:6972-6975`, whose loop runs over the rank-local `nAGN`. `clean_stop` is `MPI_FINALIZE` + `stop` (`patch/cuRamses/update_time.f90:351-365`) — collective. One rank tripping it leaves the others in a hang until the Slurm wall clock kills the job. `AGN_blast` contains no other collective, and `EsaveAGN_mpi/EsaveAGN_all` (`:6856`) are already declared and currently unused, so the reduction buffer is sitting right there. The driver gets this right (`snrt_ramses_driver.f90:355-366, 373-382, 392-401` reduce the decision before stopping); `AGN_blast` does not.

**B3 — A pre-existing bad cell aborts the run.**
`agn_feedback_deposition.f90:58` rejects `internal_old < -tol`, i.e. any cell already carrying `E_tot < E_kin`. RAMSES produces such cells transiently in strongly shocked, poorly resolved regions — exactly the galactic-nucleus gas an AGN sits in — and the legacy branch tolerated it. That rejection now escalates through `deposition_error` to a whole-run stop (and via B2, a hang). This needs to separate "invalid AGN-side input" (a code bug ⇒ globally reduce, then stop) from "invalid incoming cell state" (skip the cell, roll its intended energy into the deferred bucket, count and report). The smoke covers the NaN-input rejection but not this path.

---

## Declared, but understated in the evidence

**S1 — Passive scalars beyond `imetal` are not entrained, and the fix amplified the error by ~N_cyl.** `deposit_agn_cell` transports only `imetal`; `average_AGN:6717-6721` removes only `imetal`. With `NVAR=18`/`NDIM=3` this build carries ~13 passive scalars (`hydro_parameters.f90:98-103` plus the chemistry/dust block). Because the resolved jet now returns the *full* loaded mass instead of ~1/N of it, the specific-abundance dilution in cylinder cells and concentration in the donor grows by the same factor N. Total `Σρ X·V` is untouched, so this is species mis-mixing, not mass loss — but at N× the previous magnitude. Evidence line 105 flags it qualitatively; it should say the normalization fix scales it.

**S2 — Same amplification for the BH-frame recoil.** Removal takes `dm·v_donor`; the resolved return injects `m_load·v_BH ± lobe imbalance`. The net unbalanced momentum source per jet event was previously ~1/N of the loaded mass and is now the full mass. Evidence line 104 declares this open; the "NOT calibration-neutral" line at 26 should say explicitly that the correction scales the *intended coupling and both known non-conservations together*.

---

## Minor

1. **Duplicated membership arithmetic** (`:6644-6652` vs `:7070-7073`). Bit-identical only if the compiler contracts both sites the same way; a 1-ulp difference at `drjet <= rmax` puts a cell in `psy_norm` but not in the deposit (or vice versa), leaking mass. A shared helper or a small margin removes the class.
2. **Non-OpenMP build of `AGN_blast` cannot compile**: `nthreads` is declared only inside the `#else` of `#ifndef _OPENMP` (`:6844`) but used at `:7112`; `PEsaveAGN` is allocated only under `#ifdef _OPENMP` (`:6869`) yet zeroed at `:6990`; `TEsaveAGN` would be unassociated at `:7090`. Dead in practice — `bin/Makefile:41` is `mpiifx -qopenmp`. Guard consistency only.
3. **Refactor leftovers**: `nCOM` (`:6895`), `scale_m` (`:6894`), `onethird` (`:6879`), `EsaveAGN_mpi/EsaveAGN_all` (`:6856`) are computed/declared and unused; ~12 names in the OMP `private` clause (`:7013-7016`) no longer appear in the loop body, which makes the clause misleading to the next reader.
4. **Undeclared semantic change**: `agn_eddington_ratio` turns zero-Eddington from `NaN` (both `.lt.` and `.ge.` false ⇒ AGN skipped entirely) into `0` (⇒ radio/jet branch). Unreachable in practice — `dMEd_coarse=0` implies `dMsmbh=0` implies `Mfrac<jetfrac` — but one line in the evidence would close it.
5. `agn_eddington_ratio` living in a module named `agn_feedback_deposition` is a cohesion nit.
6. Smoke gap (proportionate, not a gate): no case exercises the `f_ekAGN`/`p_gas` thermal-plus-kinetic split. The cylinder case does correctly exercise the tight `energy_delta == kinetic_input` boundary that `tol` at `:57` exists for.

---

## Not gates

I am **not** re-raising as new conditions: simultaneous RT+mechanical physical allocation, AGN SED/obscuration, spin partition, GPU/CUDA sink parity (confirmed inert — `amr_step.jaehyun.f90:424` is the only dispatch and it is CPU), AMR/MPI jet-ownership qualification under full evolution, or AGB raw-source work. Those are correctly declared and correctly deferred.

## Bottom line

The approved scope — pure staged depositor, exact energy identity, volume-weighted normalization, donor-velocity fallback, midplane second moment, cap headroom, radicand clamps, `vol_blast` in the fallback, rank-uniform pre-source exclusivity, coarse-step re-emission removal — is implemented correctly and the evidence describes it accurately without overclaiming. B1–B3 should close before either a production legacy-AGN run or any SNRT-only run; none of them require reopening the physics decisions the operator already made.
