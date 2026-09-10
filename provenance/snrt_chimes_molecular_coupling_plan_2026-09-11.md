# Molecular / competing-grain spectral receiver continuation

Goal: simulation-ready RT/feedback/dust; finish the existing spectral coupling,
not another general validation framework. Operator approved this bundle.

Implementation decision after the review and retained-table scan: **do not
raise Tmol_K to 1e9**. The candidate below was rejected, not silently adopted.
The new native cold API enforces the molecular cooling grid's actual maximum
(10^4.98 K for the pinned bank), including internal thermal trials. General
hot/cold live admission remains open. See the
[implementation record](snrt_chimes_molecular_coupling_implementation_2026-09-11.md).

Q-GOAL first: does the concrete closure below advance that goal honestly?
Q-LEAN second: remove excessive machinery/gates. Then assess physics.

Already done: CHIMES157 conservative atomic/Auger node solver, native hot
photo/dark split, per-direction N/E live MPI staging and exact restart.

Implementation design to review once (not a per-helper audit):

- Add fixed dust absorption per spectral node to the SAME survival-fraction
  ODE as gas. Accumulate separate grain energy/photon counters. Transport
  must then have neither gas nor dust absorption for this mode; dust IR and
  material use the accepted grain counter, not another optical-depth debit.
  Four D03 C/silicate bins only initially; Fe/PAH/relative motion unchanged.
- Build molecular coefficients on the existing128 nodes from retained
  Leiden H2/CO HDF5. Integrate narrow lines over node support, not sample
  point values. Distinguish nonionizing absorption from dissociation; their
  ratio produces reacting vs fluorescent events. Preserve actual node photon
  energy in counters; no extra photon for each dissociation. This remains a
  coarse spectral closure, not resolved line radiative transfer.
- Existing30 CHIMES empirical molecular/negative-ion channels retain their
  declared Habing-shape approximation using actual 6–13.6eV energy density
  (CHIMES normalization c*5.29e-14). Do not claim measured cross sections for
  those channels. Use the pinned native reactant/product maps for nuclei/charge.
- Retain local H2/CO shielding approximation (T, cell column, b=1km/s), but
  do not reapply dust or HI attenuation already handled by transport/receiver.
  Flag if multiplying this shielding by resolved absorption would double
  count molecular optical depth and requires a different closure.
- Molecular photon energy not supplied as measured gas heat stays an explicit
  unresolved excitation/dissociation/fluorescence budget. H2 direct heat uses
  native6.4e-13erg/dissociation, bounded by available energy; UV pumping needs
  an energy-limited treatment of NON-dissociating absorbed photons, not an
  unconstrained extra heat source. No invented cooling-ray emission.
- Avoid upstream zero_molecular_abundances at1e5K. Candidate approach:
  a separate nonradiative entry using a LOCAL config Tmol_K=1e9 to retain
  molecular evolution across the boundary. Audit rate-domain validity and
  hot grain/cooling effects before admission; leave old grey/hot APIs unchanged.
- New model/data identities only when the full coupled path is verified.
  Reuse native smoke and one short MPI/restart comparison; no new framework.

Please inspect source and say which choices are safe, which need correction,
and if molecular energy/rate data prevent honest runtime admission. Read-only;
do not edit files or launch jobs. Main files: patch/lagRamses/snrt_chimes_photo.cpp,
snrt_chimes_bridge.c, snrt_chimes.f90, snrt_ramses_driver.f90; tool
simulation/snrt/tools/build_chimes_group_tables.py. Upstream local source:
.medium-pah-charge.wBiAfM/chimes/src/{chimes.c,chimes_cooling.c,update_rates.c}.
