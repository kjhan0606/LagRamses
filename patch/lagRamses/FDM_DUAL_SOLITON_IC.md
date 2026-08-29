# Dual-soliton, dual-SMBH pure-FDM seed

The opt-in `&fdm_params` switch `fdm_dual_soliton_ic=.true.` bypasses the
Grafic density/velocity transform and initializes the wave field as the
coherent sum of two moving soliton-like profiles,

```text
psi = psi_1 + psi_2
rho_i = rho0_i / [1 + c (r_i / rc_i)^2]^8 .
```

For each component provide `fdm_dual_soliton_rho0(i)`,
`fdm_dual_soliton_rc_box(i)`, `fdm_dual_soliton_center_box(i,1:3)`,
`fdm_dual_soliton_velocity(i,1:3)`, and
`fdm_dual_soliton_phase(i)`.  Centres and radii are fractions of the periodic
box; centres must be in `[0,1)` and radii below `0.5`.  Velocities are in code
coordinate units.  `fdm_dual_soliton_profile_c` defaults to `0.091`.

This is an all-wave seed.  If `fdm_use_hjm=.true.`, then
`fdm_first_wave_level` must not exceed `levelmin`; an HJM fluid/wave seam is
rejected because a coherent two-core overlap has no unique single-stream
Madelung state.  The seed is deliberately marked as requiring relaxation and
conservation validation.  It is not a stationary two-soliton solution or a
calibrated FDM decay law.

SMBHs are not created or modified by this switch.  lagRamses' established
`ic_sink` reader must supply exactly two sinks, with their initial mass,
position, velocity, angular momentum, SMBH mass, and zero CDM fraction.  The
FDM_TOY command

```bash
python scripts/materialize_dual_soliton_ic.py seed.yaml output_directory
```

writes a matching all-wave namelist fragment and two data-only `ic_sink` rows
from one typed seed manifest.  Its `box_length_code` must match lagRamses
`&AMR_PARAMS boxlen`.  It creates no run and does not replace the required
relaxation, paired-resolution, phase-replicate, wave-seam, or boundary checks.
