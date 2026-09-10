# Carbon-fragmentation input assessment

Working tree `/gpfs/kjhan/LRD_JWST`, origin `kjhan0606/LagRamses`.
This is source/model selection for approved medium-term group7, not a
completed runtime implementation or an extra validation gate.

## Existing daughter domain

`dust_pah_radiation.f90:pah_charge_optics` uses the source relation
NC=468*(a/0.001 micron)^3 and a minimum tabulated radius 0.0003548 micron.
Hence NC_min=20.902425613056. C24 -> C22 is within the table; C22 -> C20
is not. This is an optical-table bound, not a physical proof that C20
instantaneously becomes gas, nor a justification for an inert C22 reservoir.
The existing C24 H-state operator does not evolve carbon size.

## Gas products and actual photoprocessing data

The pinned CHIMES157 species enum contains neutral C2, but none of C2H,
C2H2, C2H+, C2H2+ or C2+. Existing molecular carbon cannot be relabelled
as these distinct chemical species. Source:
`.medium-pah-charge.wBiAfM/chimes/src/chimes_proto.h`.

Original Leiden files are retained at `.pah-fragments.cqIULN/`:

- [C2H2 cross sections](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/C2H2/C2H2.txt):
  12127 knots, 6.26--300 nm, source file dated 2017-03-01.
  SHA256 `ccd74ca08afc565dedf130593b0c3a89d9c5455bffd9d74c0f8545093ce522fb`.
- [C2H cross sections](https://home.strw.leidenuniv.nl/~ewine/photo/data/photo_data/cross_sections/C2H/C2H.txt):
  686 knots, 90--178.72 nm, source file dated 2022-11-16.
  SHA256 `31ab22e4d53433b960824bd393e1075dafa4113b2bd5553f37b45f0a3ec581d3`.

All supplied cross sections are finite/nonnegative. Interpolation within the
tabulated wavelength range at the existing representative group energies:

| Species | Energy (eV) | sigma_diss (cm2) | sigma_ion (cm2) |
| --- | ---: | ---: | ---: |
| C2H2 | 10.56942670578782 | 6.512681752877522e-17 | 0 |
| C2H2 | 12.2954112566775 | 0 | 3.350476394268827e-17 |
| C2H | 10.56942670578782 | 8.069241325545856e-26 | 0 |
| C2H | 12.2954112566775 | 0 | 4.83244e-17 |

These are sample values, not integrated group coefficients or universal
branching probabilities. The [Leiden documentation](https://home.strw.leidenuniv.nl/~ewine/photo/cross_sections.html)
lists C2H+H as the acetylene threshold products, and C2+H for ethynyl.
It explicitly warns that higher-energy branching can differ. Thus a two-step
neutral-only instantaneous conversion C2H2 -> C2 + 2H is not established by
these data, even inside the present <=13.6 eV domain. Photon consumption,
ion chemistry, electron energy and fragment chemical energy cannot be omitted.

## PAH fragmentation prescriptions

[Murga et al. 2020, Table1](https://arxiv.org/html/2007.06568v1) provides
microcanonical fragmentation parameters, including a normal/dehydrogenated
C2H2-loss activation energy 4.6 eV and entropy 10 cal/K/mol. The body's
`kcal` wording conflicts with the table's `cal`; 1000x scaling is not used.
Activation energy is not automatically the net reaction enthalpy between a
particular parent isomer, daughter isomer and gas fragment. The daughter
representation and energy partition must be specified before connecting a
live conservative receiver.

[Lange et al. 2025](https://doi.org/10.1051/0004-6361/202347722) treats
competing H/H2/C2H2 pathways, acknowledges uncertain pyrene extrapolation,
and stops its cascade at carbon-backbone damage. A first-damage probability
is not a complete atomization or daughter-evolution model.

No new passive tensor, chemical ABI, namelist option or destruction rate has
been activated on the strength of these source checks alone. The Fable plan
assessment addresses the least expensive physically defensible live scope.
