# P4 transport pilot

`p4_run_transport_pilot.py` runs a bounded static-grid calculation with the
P4 photon ledger and gas input. The pilot uses the Carlson S4 quadrature,
reduced light speed `0.01 c`, and a directional CFL number of `0.4` by default.
It deposits sources in physical cgs units, evolves photon-conserving H/He
chemistry with the P3 fixed-iteration recombination closure, and enables the
P2 high-energy secondary-ionization closure.

Dust is exactly zero because the staged snapshot does not provide metallicity
or dust information. Gas temperature is held fixed: gas heating is output as
a rate and is not fed back into hydrodynamics. The initial cross-epoch input
from output 00016 gas and output 00017 sources was an interface test only.
The current coeval input is resampled from complete output 00017 adaptive AMR
leaves and uses the same output 00017 source ledger.
