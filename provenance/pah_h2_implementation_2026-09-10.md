# Group7: native PAH H2 vacancy refilling

Working tree `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
New optional model: `dust_pah_model='pah_h2_rehydrogenation_v1'`.
Native implementation, four-step live evolution, exact restart comparison
and driver end evaluation **PASS for this restricted comparison**.

## Physical scope and review

[Fable plan review](pah_h2_fable_2026-09-10.txt) approved this bounded
subchannel after Q-GOAL and Q-LEAN. Applied its restrictions: no additional
carriers or CHIMES ABI, no stored molecular-energy field, no old-model control
simulation, and do not describe the result as an upper bound on H2 effects.
The driver independently read the paper that Fable could not retrieve.

[Montillaud et al. 2013, section3.4/table4](https://arxiv.org/html/1301.6507v1)
provides a sensitivity model in which a PAH cation captures two H atoms per
H2 at an experimentally bounded rate5e-13 cm3/s. It is not a measured best
estimate. We retain only refilling two vacant peripheral sites:
`C24Hn+ + H2 -> C24H(n+2)+`, n=0..10. No single-vacancy abstraction,
H2 superhydrogenation, catalytic H2 formation or carbon destruction is
inferred. Restricting channels means the resulting effect is NOT an upper
bound even though the coefficient is the literature bound rate.

The existing M13/DL01 comparison's two-H PAH bond difference is8eV.
Ground-state H2 dissociation energy is4.4781eV
([NIST CODATA compilation, equation11](https://physics.nist.gov/cuu/Constants/codata.pdf)).
Thus the declared daughter excitation increment is3.5219eV plus the captured
molecule's translational kinetic energy. This uses the existing approximate
PAH bond ladder, not new molecule-specific laboratory enthalpies. Constant
k selects3kT/2. Rotational/vibrational gas-H2 states are not introduced;
gas capacity remains the existing fixed-gamma CHIMES particle-count closure.
Neutral/cation IP and normal-H optics assumptions are unchanged.

## Runtime connection

- `dust_pah_hydrogen.f90`: reuse implicit finite donor solve with optional
  stride2 and donor binding. Old stride1 arithmetic stays unchanged. After
  atomic capture, cation H2 capture debits the finite molecular donor, before
  the existing photon/IR/H-loss step. This is explicitly split, not a claim
  of simultaneous H/H2 integration. Energy balance includes the positive
  chemical-energy cost of dissociating captured H2. Failures leave the
  caller's PAH, HI/H2, electrons, gas heat and IR outputs unchanged.
- `dust_pah_live_model.f90`: append H2 coefficients and D0 to the existing
  identity only for this model, with new tags. Existing identity arrays are
  unchanged. No population shape or mass convention changes.
- `dust_pah_mixed.f90`, `snrt_dust_live.f90`: stage H2 across iterative IR
  material solves/substeps and update heat capacity per molecule. Add only
  the substep's molecular binding **change** to the material energy balance;
  do not add D0 to thermal energy or store a negative chemical-energy bath.
- `snrt_ramses_driver.f90`: pinned CHIMES157 index138 is H2, stored as one
  number-equivalent weight per molecule, not twice-H mass. Convert to cm^-3,
  pass to the native solve, and stage the result in `chemical_trial`. The
  existing PAH mass inventory supplies the corresponding bound-H change;
  there is no separate edit to the total element ledger.
- `dust_mass_physics.f90`, `read_hydro_params.f90`, `mkrun.py`, and
  `patch/cuRamses/aux/ramses_nml_generator.py`: expose the same named choice,
  restrictive-domain description and existing build/state requirements.

Selection requirements: noncosmo, coadvected, no Fe, DUST_LIVE/CHIMES ABI5,
H0--13 x2 charges x128 excitation =3584 carriers; hydro NVAR3771/NENER1.
Primary photons <=13.6eV and gas10--10000K. Both original neutral/ion optics
tables are required. Defaults remain unchanged. Makefile/VPATH is unchanged.
This model is an explicit comparison, not general production qualification.

## Tests and artifacts

Artifacts: `.pah-h2.syBrDR/`.

- Checked build `build-h2.log`: Intel MPI/ifx, `-check bounds`, CPU/OpenMP,
  SNRT=1 DUST_LIVE=1 CHIMES=1 HDF5=1, NVAR3771 NENER1 NVECTOR32,
  USE_FFTW=0 FDMDEBUG=1. Existing pinned CHIMES library at
  `.medium-pah-charge.wBiAfM/chimes`, SUNDIALS at
  `.dust-extension.AOz7mU/sundials-install`.
- Binary `ramses_pah_h23d` SHA256
  `d6ab65b5ad260f86111f5354ff83b157aaa0d92dd6f9359a8a8584b140e7e555`.
- Existing `dust_mass_smoke` extended, no new test framework. `smoke.log`:
  analytic finite-H2 donor, energy and H nuclei, charged physical optics,
  zero-donor exact old-model parity, identity preservation, missing-donor
  rejection, and post-staging excitation-overflow rollback pass. Physical
  H2 consumption1.38228766128345e-5 cm^-3; relative energy residual7.42e-15,
  H-nuclear residual4.85e-17 cm^-3. Existing Fe/neutral/fixed-H tests pass.
- `gnu/run.log`: same tests with GNU bounds checking and traps for invalid,
  division by zero and overflow. Physical energy residual about-3.54e-16.
  Unchanged dependency objects reused from prior GNU smoke build; modified
  mass/H modules and smoke were rebuilt. No independent full GNU run claim.
- One split-refinement loop with both finite H/H2 donors: L1 population
  differences versus128 substeps at4/8/16 substeps are
  0.08541785/0.05856481/0.03562564, decreasing as expected. This establishes
  a refinement trend, not full timestep/excitation convergence or a calibrated
  physical accuracy bound.
- `gui.log`:47 tests,46 passed,1 display-dependent skip; new model's generated
  namelist/env/build description and incompatible configurations covered.
- `live/physical.nml`, `environment.sh`: periodic noncosmo4^3, MPI2/OMP2,
  hydro/advection/CHIMES/IR/PAH and advective CR; no stars, AGN, MHD or incident
  FUV. Initial neutral H12 and cation H10 populations plus nonzero H2.
  Four steps complete in134.408s. This fixture checks live H2 coupling and
  IR, while physical-optics/FUV coverage comes from native tests.
- `restart/physical.nml`: step2 -> step4 with the same binary/settings.
  Both launch policies reported before execution: noutput1/aout2/tout1e30,
  foutput2/fbackup1e6; free /gpfs about59TiB, expected raw output45MiB total.

## Final evaluation

`evaluation.txt` and its one-off read-only evaluator retain the measured
results. Restart completed in77.4944s. All7542 hydro,2 RT,8 gravity,12 AMR
datasets and the11562-coefficient PAH identity dataset match exactly. Physical
clocks/step counters match. Maximum nuclear relative residual6.51e-16,
charge/electron residual1.68e-15, skeleton-count residual5.56e-16; finite,
nonnegative chemical/PAH states. IR relative balance<=1.837e-15.

Mean H per skeleton is11.6307101 at step2 and11.6761321 at step4. Gas H2
mean (native code density) is3.69899999649 and3.69899997819. The latter is the
**combined** live CHIMES+PAH evolution, not a measurement attributing all H2
change to PAHs. The active native H2 test isolates the added channel.
The chemical-binding scalar is included in local IR accounting; no global
hydro energy-closure claim across CHIMES's external cooling is inferred.

Driver end evaluation complete. Raw cleanup is recorded in
`pah_h2_raw_cleanup_2026-09-10.md`; retained inputs/logs/build/evaluation allow
reproduction. No further audit or approval gate added.

## Remaining scope

This advances group7's H2 pathways but does not close that whole item:
PAH-catalytic H2 formation, single-vacancy reactions, molecular
superhydrogenation, higher charge states and general hard-source survival
remain absent. Carbon skeleton/fragment chemistry remains unresolved for the
previously documented carrier, enthalpy and optical-domain reasons. It is
not silently replaced by this H2 channel, rejected as a goal, or shifted to
an operator-source-finding requirement. No commit/push requested or performed.
