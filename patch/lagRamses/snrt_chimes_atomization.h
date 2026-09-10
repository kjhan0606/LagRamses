/* ATcT v1.130, gas-phase ground-state formation enthalpies at 0 K,
 * kJ/mol. H3+ para; CH2 triplet; C2 singlet ground state. Source/IDs:
 * simulation/snrt/data/chimes_atomization_atct1130.json.
 * Indexed as CHIMES sp_H2..sp_O2p. NOT high-temperature rate coefficients. */
static const double atomization_hf0[20]={
  0.,1488.364,1110.304,37.279,-238.902,819.991,0.,827.754,
  592.832,391.050,1099.347,-113.800,1619.753,1393.18,1293.361,
  978.496,605.88,1238.310,986.01,1164.585
};
static const double atomization_atom_hf0[3]={216.034,711.396,246.844};
static const double atomization_ion_hf0[3]={1528.084,1797.849,1560.786};
static const double atomization_kjmol_to_ev=1000./(6.02214076e23*1.602176634e-19);
