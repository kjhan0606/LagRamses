#ifndef SNRT_CHIMES_SPECTRUM_H
#define SNRT_CHIMES_SPECTRUM_H
// Immutable, explicitly loaded native atomic reaction bank. No global
// replacement of CHIMES tables and no CVODE mutation from these functions.
#ifdef __cplusplus
extern "C" {
#endif
int snrt_chimes_band_load(const char *path, void **handle, int *reactions, int *shells, double *identity);
void snrt_chimes_band_free(void *handle);
int snrt_chimes_molecular_load(const char *path,void **handle,double *identity);
void snrt_chimes_molecular_free(void *handle);
int snrt_chimes_band_reactions(void *handle, int reactions, int *mapping);
// N[band][direction] = photons/cm3, E = eV/cm3, direction-integrated shares.
// Output [band][moment][reaction]: sum N*sigma, N*sigma*E and
// N*sigma*(E-shell binding). Multiply by c_hat and target number density
// for event/energy rates. These are NOT finite-inventory absorbed counts.
// Photoelectron moment excludes unprovided cascade/Auger relaxation energy.
int snrt_chimes_band_moments(void *handle, int directions, int reactions,
    const double *number, const double *energy, double *moments);
// Conservative atomic photo operator, NOT a complete CHIMES network step.
// Requires initialized CHIMES stoichiometry and the existing FS2010 tables.
// Species are n_i/nH; N/E use the layout above. c_hat in cm/s, dt in seconds.
// ledger[0:7]: heat, HI/HeI/HeII secondary ionization, excitation,
// unresolved binding/cascade reservoir, absorbed photon energy (all eV/cm3),
// and primary absorbed photon NUMBER (cm^-3). Heat is an increment, not T.
// No recombination, collisional chemistry/cooling, molecular dissociation,
// grain absorption or gas thermal evolution is included in this operator.
// Success publishes all outputs; any failure leaves every output untouched.
int snrt_chimes_band_photo_step(void *handle, int directions, double nH,
    double dt, double c_hat, const double *abundance, const double *number,
    const double *energy, double *next_abundance, double *next_number,
    double *next_energy, double *ledger);
// Same photo operator with competitive fixed-grain absorption. alpha is
// [group][128 nodes] in cm^-1, aligned to the loaded bank, already multiplied
// by physical grain density. No second dust absorption may precede/follow it.
// ledger[0:7] retains GAS-only semantics above; [8] is accepted grain energy
// (eV/cm3), [9] accepted grain photon number (cm^-3). Not a material/IR solve.
int snrt_chimes_band_photo_dust_step(void *handle, int directions, double nH,
    double dt, double c_hat, const double *dust_alpha, const double *abundance,
    const double *number, const double *energy, double *next_abundance,
    double *next_number, double *next_energy, double *ledger);
// Atomic + molecular + grain finite-photon operator. The molecular handle
// declares the coarse line projection / empirical shape closure. shield[2]
// is an explicit fixed effective H2/CO multiplier in [0,1]; unity is the
// optically thin coefficient reference, NOT an implicit local-column law.
// Molecular non-heating energy enters the unresolved ledger[5], including
// fluorescent absorptions that do not change species. H2 direct heat is
// 6.4e-13 erg per dissociation; pumping uses min(native16.9eV/dissociation,
// 2eV/fluorescent capture) times the explicit critical-density factor in [0,1].
// Other photofragment heat is not supplied.
// No dark chemistry, temperature update, IR or live admission in this API.
int snrt_chimes_band_photo_molecular_step(void *handle,void *molecules,int directions,
    double nH,double dt,double c_hat,const double *dust_alpha,const double *shield,double pumping,
    const double *abundance,const double *number,const double *energy,
    double *next_abundance,double *next_number,double *next_energy,double *ledger);
// Same operator, with nine accepted grain N (cm^-3) / E (eV/cm3) counters
// for the existing group-wise material receiver. These sum to ledger[9]/[8].
// No second attenuation or representative-energy substitution is permitted.
// Both arrays and the original outputs remain untouched on failure.
int snrt_chimes_band_photo_molecular_groups(void *handle,void *molecules,int directions,
    double nH,double dt,double c_hat,const double *dust_alpha,const double *shield,double pumping,
    const double *abundance,const double *number,const double *energy,
    double *next_abundance,double *next_number,double *next_energy,double *ledger,
    double *grain_number,double *grain_energy);
#ifdef __cplusplus
}
#endif
#endif
