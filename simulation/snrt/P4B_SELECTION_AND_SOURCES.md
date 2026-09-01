# P4-B: high-density target, thermal table, and source catalogue

## Density target

Select the target by maximizing the mean gas density within a fixed, periodic scout-grid cube. This is intentionally not the single highest-density AMR cell: a peak cell can be a transient unresolved clump and does not define a useful RT volume. `select_high_density_region` returns the left edge and width in code coordinates, so the same physical region can be restaged at the final AMR level.

The required user-level choice is the physical target width. The scouting and final grids must represent that same width. Report the cube's mean density, maximum cell density, source count, and total gas mass before calling it the primary target.

## Mean molecular weight table

For primordial gas with hydrogen mass fraction `X` and helium fraction `Y=1-X`, the exact composition relation is:

```text
1 / mu = X (1 + x_HII) + Y / 4 (1 + x_HeII + 2 x_HeIII)
```

`PrimordialMuTable` stores `mu`, `x_HII`, `x_HeII`, and `x_HeIII` on `(log10 T, log10 n_H)` axes. It is stored as a portable `.npz` file and uses bilinear interpolation with edge clamping. During hydro staging it solves `T=P*mu(T,n_H)*m_p/(rho*k_B)` by damped fixed-point iteration and transfers the corresponding initial ionization fractions into the RT state.

The current `output_00016` has no cooling/metallicity output. Therefore the only defensible default table is neutral primordial H/He (`mu=1.219512`), which is implemented as `neutral_primordial_mu_table`. It is not a substitute for an equilibrium chemistry table. A future ionized table must be generated from an auditable H/He chemistry calculation with its UV background and redshift recorded. The requested search did not locate such a formula/table in the readable `Eunha.A1` source tree.

## Recommended source catalogue

Use one spatial catalogue with three named source scenarios, all binned onto the exact final RT grid:

1. `STAR`: every star particle, with age and metallicity interpolated through a binary stellar-population SED. Use BPASS as the baseline and record the IMF, metallicity interpolation, and photon luminosity in each RT group. Binary populations affect both the duration and hardness of ionizing stellar emission, so a single-age luminosity conversion is not adequate. [Stanway, Eldridge, and Becker (2016)](https://arxiv.org/abs/1511.03268)
2. `AGN`: every resolved sink/BH with luminosity derived from the instantaneous accretion rate, `L_bol=epsilon_r*dot(M)*c^2`, then converted through a declared AGN SED to the same RT groups. Do not derive luminosity from stored BH mass: this keeps the RT source ledger independent of the previously identified seed/accretion-mass bookkeeping risk.
3. `STAR+AGN`: the union of the first two catalogues without any global escape-fraction multiplier. Resolved gas and dust should determine the transmission. Any unresolved birth-cloud attenuation must instead be a clearly labelled, separately varied source-side nuisance parameter.

The initial paper comparison should run all three scenarios on the same selected density peak. This directly separates stellar pre-ionization from AGN-driven hard-photon escape and prevents an AGN claim from being set by an unrecorded stellar prescription.

For the registered transitional checkpoint `output_00011`, the native stellar
metadata reader is now complete: it supplies star positions, current and
initial masses, proper-time ages, birth metallicities, and yield-table
progress values. It does not assign a stellar SED or photon luminosity, so the
STAR row becomes usable only after an explicit population-synthesis asset and
group integration contract are selected.

The group integration contract is now implemented and tested in
[`tools/p4_build_stellar_photon_ledger.py`](tools/p4_build_stellar_photon_ledger.py);
see [`P4_STELLAR_SED.md`](P4_STELLAR_SED.md). The remaining blocker is the
versioned production SED/IMF asset, not the CSV-to-`q_group_N_s` wiring.

## Deferred inputs

- physical width of the selected cube and its final cell count
- auditable chemistry source for an ionized `(T, n_H)` table
- sink/BH catalogue reader that preserves instantaneous accretion rate
- stellar population-synthesis SED asset, IMF, and its approved age/metallicity
  interpolation range
- production choice of source-side escape/birth-cloud attenuation parameters
