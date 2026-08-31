# Paper III: High-Redshift SMBH and Little Red Dot Zoom-Ins

## Status and central question

This note records the initial simulation design. The project asks whether
rare, rapidly assembling overdensities at z > 14 produce the gas-obscured,
rapid SMBH growth and observable dual-AGN phases that could underlie a
subset of JWST little red dots (LRDs) at z about 5--9.

The z > 14 objects are progenitors in the simulation, not direct
counterparts of the observed LRD sample. The project must make falsifiable
predictions for LRD-like abundance, host overdensity, companion excess,
dual-AGN fraction, and their dependence on seed and feedback models.

## Simulation architecture

### Parent and DMO screening

- Use the existing full HR5 volume as the parent volume. Do not replace it
  with an arbitrary small standalone DMO box: rare z > 14 peaks and their
  cosmic context are central to the experiment.
- Build a selection field from the parent DM distribution and merger tree.
  If particle IDs or initial-condition mappings are insufficient for
  Lagrangian tracing, construct an exact-phase DMO companion for that
  purpose. It is a screening and tracing calculation, not the science
  volume.
- Smooth the z = 14 density field on R = 0.5 and 1.0 cMpc. Select
  independent peaks in the top 0.5 percent and 1 percent, respectively.
- Define targets without using a BH or AGN outcome. Initially retain z = 5
  descendants with 10^10.5 <= Mvir/Msun <= 10^11.8 and rapid z = 14--5
  assembly, including a q >= 1:4 merger or multiple massive progenitors.
- Match controls by z = 5 halo mass. Include both dense, quiet controls and
  typical-density controls.

### Zoom regions and sample

- Trace all DM particles within 3 R200 of each z = 5 descendant back to the
  initial conditions. Use their buffered Lagrangian hull as the high
  resolution region; do not impose a fixed Eulerian cube.
- Expect an initial high-resolution patch of roughly 2--5 cMpc, with
  multi-mass particle shells outside it. Reject a zoom if low-resolution
  contamination reaches 3 R200 by z = 5.
- Run 12 high-density rapid-assembly targets, 4 dense quiet controls, and
  4 mass-matched typical-environment controls: 20 fiducial hydro zoom-ins.
- Re-run 2 targets and 2 controls at the highest-resolution convergence
  level. Re-run a representative subset of 6 fields for seed and feedback
  sensitivity.

### Resolution

| Quantity | Fiducial hydro suite | Convergence suite |
| --- | --- | --- |
| Final redshift | z = 5 | z = 5 |
| DM particle mass | 1--3 x 10^4 Msun | <= 3 x 10^3 Msun |
| Initial baryon mass | 2--6 x 10^3 Msun | <= 6 x 10^2 Msun |
| Finest cell at z = 5 | 4--5 pc physical | about 1 pc physical |
| Finest comoving cell | about 24--30 cpc | about 6 cpc |
| Output cadence | 5--10 Myr, 0.25--1 Myr for close pairs | same |

For a 200 cMpc parent convention, 24 cpc corresponds approximately to
lmax = 23 and 6 cpc to lmax = 25. The actual AMR level must be set from the
chosen parent convention, not copied blindly.

Refine quasi-Lagrangian on DM and gas mass, enforce at least 8 cells per
Jeans length, and use a resolution-aware pressure floor only where necessary.

## Physics requirements

- Include non-equilibrium H/He/H2 chemistry, metal cooling, stellar
  feedback, BH accretion, and BH feedback.
- Treat light and heavy seed prescriptions as separate hypotheses, not as
  nuisance parameters. At minimum compare a light-seed branch with a
  10^4--10^5 Msun heavy-seed branch.
- Preserve BH identities and relative orbits after close approach. A merge
  at a fixed number of cells is a numerical capture event, not a physical
  SMBH coalescence. Store such systems as unresolved pairs or evolve them
  with an explicitly documented subgrid model.
- Create JWST mock observations with dust and radiative-transfer
  post-processing. Hydro luminosities alone cannot identify an LRD.

## Audit of the current lagRamses reference configuration

The reviewed reference namelist is
/home/kjhan/BACKUP/lagRamses/runs/lageunha/snapshot-2026-08-24/Agn/Run0/run_cdm/cosmo.nml.
It enables cooling, metals, density-threshold star formation, kinetic and
delayed-cooling stellar feedback, Bondi sink accretion, sink AGN feedback,
gas drag, and a 10^4 Msun seed. The source contains an optional RT path and
an ATON path, but this namelist does not enable either one; rt and aton
therefore retain their false defaults. haardt_madau is a uniform background,
not local radiation transport.

The AGN routine in sink_particle.kjhan.f90 uses Bondi/Eddington accretion
and switches between thermal blast and low-accretion jet-like feedback. It
does not propagate photons through absorbing gas.

The reference parameters are not portable without recalibration:

- n_star = 0.1 cm^-3 is a coarse-volume star-formation threshold, not a
  suitable default for a 4--5 pc nuclear zoom.
- delayed_cooling, f_w = 3, and the SN efficiencies can double count
  unresolved losses once individual dense structures are resolved.
- rAGN = 4 kpc is far larger than a pc-scale nuclear region. AGN coupling
  must instead be tied to a local resolved mass or to a radius of order a
  few finest cells, with an explicit convergence test.
- Mseed = 10^4 Msun, Bondi accretion with boost_acc = 2, and boost_drag = 2
  are one heavy-seed and unresolved-flow prescription. They cannot be used
  as the sole physical explanation for LRDs.

The current stellar/BH subgrid model is a baseline comparison, not the Paper
III fiducial model.

## Calibration policy

Feedback calibration is mandatory at this resolution. It must be called
calibration rather than fine tuning: parameters are constrained using
observables and numerical tests that are independent of the final LRD and
dual-AGN claim, then frozen before the science comparison.

Calibrate only the parameters that change the resolved physics:

- star-formation threshold and efficiency;
- SN energy/momentum coupling and any delayed-cooling prescription;
- AGN coupling radius, thermal/kinetic partition, and mode switch;
- unresolved BH accretion and drag boosts;
- seed criterion and seed mass, as separately labelled physical models.

Use a staged calibration ladder:

1. Numerical tests: Jeans resolution, isolated SN momentum/energy retention,
   BH feedback injection-radius convergence, and close-pair identity
   preservation.
2. Galaxy calibration: non-LRD z = 5--9 stellar-mass--halo-mass relation,
   UV luminosity function, size, gas fraction, and star-formation rate.
3. BH calibration: bright non-LRD AGN abundance and broad BH--host
   consistency, without using LRD counts or the dual-AGN fraction as a
   fitting target.
4. Freeze the calibrated model. Apply it unchanged to the high-density and
   matched-control zoom sample.

Do not perform an unrestricted multi-parameter search. Use a small,
pre-declared grid around physically motivated choices and report the full
parameter response. A signal that appears only after selecting one feedback
combination is not a Paper III result.

## Radiation-transport decision

### What does not require on-the-fly RT

- The 20-field hydro suite may measure halo assembly, gas inflow, BH-pair
  demographics, and the environmental component of dual-AGN statistics
  without on-the-fly RT.
- Dust and JWST photometry require radiative-transfer post-processing, but
  do not by themselves require an RHD calculation for every field.

### What requires RHD or local radiation transport

- Any claim about the physical formation of a massive seed at z > 14.
  Local Lyman-Werner and ionizing radiation regulate H2 cooling,
  fragmentation, and the direct-collapse-like pathway.
- Any claim that a local radiation-regulated gas cocoon creates the LRD
  obscuration phase.
- Any claim that AGN or stellar photoionization changes the dual-AGN duty
  cycle rather than merely its observed color.

Use a targeted RHD branch for 4 rare targets and 2 controls after the
fiducial hydro suite identifies representative assembly histories. Its
minimum photon treatment must cover Lyman-Werner and H/He ionizing groups;
the handling of AGN hard UV/X-ray radiation must be specified before the
branch is called predictive. Verify the available RT solver, photon groups,
chemistry coupling, reduced-speed-of-light choice, and GPU/MPI scaling before
production.

## Go/no-go gates

1. Confirm that HR5 initial conditions, particle IDs, and merger trees permit
   uncontaminated Lagrangian zoom generation.
2. Run one target and one matched control to z = 5 at fiducial resolution.
3. Demonstrate that changing resolution and AGN coupling radius does not
   create or erase the dual-AGN signal.
4. Only then expand to the 20-field suite and the targeted RHD branch.

