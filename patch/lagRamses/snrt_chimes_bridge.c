/* Native CHIMES adapter. External LGPL3+ CHIMES remains a separately linked
 * library. No Python interpreter, equilibrium abundance substitution, or
 * process-global per-cell thermochemical state is used here. */
#include <math.h>
#include <float.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <openssl/evp.h>
#include <mpi.h>
#include "chimes_proto.h"
#include "chimes_vars.h"

#define NS 157
#define NG 9
extern int snrt_chimes_receiver_abi(void);
/* ABI 4 remains usable for the existing neutral-solid receiver. ABI 5 adds
 * this optional hook without changing either CHIMES structure layout. */
typedef ChimesFloat (*solid_charge_fn)(const struct globalVariables *);
static solid_charge_fn *solid_charge_hook;
static int receiver_abi;
static struct globalVariables configuration;
static int initialized;
static unsigned char table_identity[32*(NG+1)];
static double *auger_heat_fuv,*auger_heat_euv;
static double *primary_heat_fuv,*primary_heat_euv;
extern int snrt_chimes_fs_ready(void);
extern int snrt_chimes_fs_fractions(double,double,double *);
extern int snrt_chimes_fs_samples(double,double *,double *);
static const double ev_erg=1.602176634e-12;
struct secondary_source {
    int n,stride;
    const int *reactant;
    const double *sigma,*energy;
    double *samples;
};
static struct secondary_source secondary_sources[4];
static double secondary_xi[14];
static double h2_sigma[NG][2],co_sigma[NG][2];
static const double group_edges[NG+1]={.01,1,5.6,11.2,13.6,24.59,54.42,500,2000,10000};
static const double group_energy[NG]={.1,2.3664319132398464,10.56942670578782,12.2954112566775,
    17.662918001204492,34.381528542610674,106.63925876378897,869.6341490248457,4023.594574013186};
static int nuclei[NS][11],charge[NS];
static const int neutral[11]={sp_HI,sp_HeI,sp_CI,sp_NI,sp_OI,sp_NeI,sp_MgI,sp_SiI,sp_SI,sp_CaI,sp_FeI};
static const int atomic_z[11]={1,2,6,7,8,10,12,14,16,20,26};

struct cell_context {
    double solid_q;
    double molecular_temperature_max; /* zero for unchanged legacy/hot calls */
    int outside_molecular_domain;
};

double snrt_chimes_molecular_temperature_max(void)
{
    if(initialized!=1 || chimes_table_bins.N_mol_cool_Temperatures<1)return 0;
    return fmin(1e5,pow(10.,chimes_table_bins.mol_cool_Temperatures[
        chimes_table_bins.N_mol_cool_Temperatures-1]));
}

static ChimesFloat solid_charge(const struct globalVariables *c)
{
    /* hybrid_data belongs to the private configuration copy, never shared
     * mutable state. Positive solid charge supplies additional electrons. */
    return c->hybrid_data ? ((const struct cell_context *)c->hybrid_data)->solid_q : 0;
}

static int file_digest(const char *path,unsigned char *digest)
{
    unsigned char block[65536];unsigned int length=0;size_t n;
    FILE *file=fopen(path,"rb");if(!file)return 1;
    EVP_MD_CTX *ctx=EVP_MD_CTX_new();
    if(!ctx){fclose(file);return 1;}
    int ok=EVP_DigestInit_ex(ctx,EVP_sha256(),NULL);
    while(ok && (n=fread(block,1,sizeof(block),file))>0)ok=EVP_DigestUpdate(ctx,block,n);
    if(ferror(file))ok=0;
    if(ok)ok=EVP_DigestFinal_ex(ctx,digest,&length);
    EVP_MD_CTX_free(ctx);fclose(file);
    return !ok || length!=32;
}

static int read_vector(hid_t file,const char *name,int count,double *values)
{
    hid_t set=-1,space=-1;int ok=0;
    if(H5Lexists(file,name,H5P_DEFAULT)<=0)return 1;
    set=H5Dopen2(file,name,H5P_DEFAULT);
    if(set<0)return 1;
    space=H5Dget_space(set);
    if(space>=0 && H5Sget_simple_extent_npoints(space)==count)
        ok=H5Dread(set,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values)>=0;
    if(space>=0)H5Sclose(space);
    H5Dclose(set);
    if(!ok)return 1;
    for(int i=0;i<count;i++)if(!isfinite(values[i]) || values[i]<0)return 1;
    return 0;
}

int snrt_chimes_identity(double *values)
{
    if(initialized!=1)return 1;
    for(int i=0;i<32*(NG+1);i++)values[i]=table_identity[i];
    return 0;
}

double snrt_chimes_boltzmann(void) { return BOLTZMANNCGS; }
int snrt_chimes_charge_supported(void) { return initialized==1 && receiver_abi==5; }

int snrt_chimes_group_binding(int n,const double *edges,const double *means)
{
    if(initialized!=1 || n!=NG)return 1;
    for(int i=0;i<=NG;i++)if(edges[i]!=group_edges[i])return 2;
    for(int i=0;i<NG;i++)if(means[i]!=group_energy[i])return 2;
    return 0;
}

int snrt_chimes_nuclear_sums(const double *state,double *elements)
{
    if(initialized!=1)return 1;
    memset(elements,0,11*sizeof(double));
    /* Also used for signed face fluxes. No positivity clipping here. */
    for(int s=0;s<NS;s++){
        if(!isfinite(state[s]))return 2;
        for(int e=0;e<11;e++)elements[e]+=nuclei[s][e]*state[s];
    }
    return 0;
}

static void fatal_library_error(void)
{
    int active=0;
    fflush(NULL);
    MPI_Initialized(&active);
    if(active)MPI_Abort(MPI_COMM_WORLD,31);
    /* Upstream's fatal callback must never turn an invalid table into a
     * successful process exit or strand other ranks in a collective. */
    _Exit(31);
}

static void initialize_stoichiometry(void)
{
    for(int e=0;e<11;e++)for(int z=0;z<=atomic_z[e];z++){
        nuclei[neutral[e]+z][e]=1;charge[neutral[e]+z]=z;
    }
    nuclei[sp_Hm][0]=1;charge[sp_Hm]=-1;
    nuclei[sp_Cm][2]=1;charge[sp_Cm]=-1;
    nuclei[sp_Om][4]=1;charge[sp_Om]=-1;
    charge[sp_elec]=-1;
    const int mol_h[20]={2,2,3,1,2,0,0,1,1,2,3,0,1,2,1,2,3,0,1,0};
    const int mol_c[20]={0,0,0,0,0,2,0,1,1,1,1,1,1,1,0,0,0,1,1,0};
    const int mol_o[20]={0,0,0,1,1,0,2,1,0,0,0,1,0,0,1,1,1,1,1,2};
    const int mol_q[20]={0,1,1,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,1,1};
    for(int j=0;j<20;j++){
        nuclei[sp_H2+j][0]=mol_h[j];nuclei[sp_H2+j][2]=mol_c[j];
        nuclei[sp_H2+j][4]=mol_o[j];charge[sp_H2+j]=mol_q[j];
    }
}

/* Fractions of primary electron kinetic energy: heat, HI/HeI/HeII
 * ionization, escaping HI excitation. FS2010 is a primordial atomic-gas
 * closure, not a molecular electron-degradation calculation. Where actual
 * atomic targets fall below its assumed populations, continuously limit
 * those channels and thermalize the unavailable share. No H2 is destroyed
 * to manufacture HI targets. This explicit extension is not a new FS table. */
static void secondary_target_limit(const double *state,double *f)
{
    double xi=fmin(1.,state[sp_HII]);
    double xr=fmax(1e-4,fmin(.999,xi));
    /* Use the original SNRT primordial mass fractions for the target
     * normalization, NOT division by the actual He abundance. The latter
     * would retain a finite He reaction rate as total He tends to zero. */
    const double he_reference=.24/(4.*.76);
    double available[3]={fmin(1.,state[sp_HI]/(1.-xr)),
        fmin(1.,state[sp_HeI]/(he_reference*(1.-xr))),
        fmin(1.,state[sp_HeII]/(he_reference*xr))};
    for(int j=0;j<3;j++)f[j+1]*=available[j];
    f[4]*=available[0];
    f[0]=fmax(0.,1.-f[1]-f[2]-f[3]-f[4]);
}

int snrt_chimes_secondary_partition(double energy,const double *state,double *f)
{
    if(!isfinite(energy) || energy<0)return 1;
    for(int i=0;i<NS;i++)if(!isfinite(state[i]) || state[i]<0)return 1;
    if(snrt_chimes_fs_fractions(energy,fmin(1.,state[sp_HII]),f))return 2;
    secondary_target_limit(state,f);
    return 0;
}

static int secondary_cache_initialize(int nspectra)
{
    /* The xi axis must exist even for a spectrum with no ionizing channels. */
    double axis_samples[84];
    if(nspectra>0 && snrt_chimes_fs_samples(200.,axis_samples,secondary_xi))return 1;
#define SECONDARY_SOURCE(TABLE,ENERGY) \
    (struct secondary_source){TABLE.N_reactions[0],TABLE.N_reactions[1],TABLE.reactants,TABLE.sigmaPhot,ENERGY,NULL}
    secondary_sources[0]=SECONDARY_SOURCE(chimes_table_photoion_fuv,primary_heat_fuv);
    secondary_sources[1]=SECONDARY_SOURCE(chimes_table_photoion_euv,primary_heat_euv);
    secondary_sources[2]=SECONDARY_SOURCE(chimes_table_photoion_auger_fuv,auger_heat_fuv);
    secondary_sources[3]=SECONDARY_SOURCE(chimes_table_photoion_auger_euv,auger_heat_euv);
#undef SECONDARY_SOURCE
    for(int j=0;j<4;j++){
        struct secondary_source *s=&secondary_sources[j];
        s->samples=calloc(NG*s->stride*84,sizeof(double));
        if(!s->samples)return 1;
        for(int b=0;b<nspectra;b++)for(int i=0;i<s->n;i++){
            int k=b*s->stride+i;
            if(s->energy[k]>=10.*ev_erg && s->sigma[k]>0 &&
               snrt_chimes_fs_samples(s->energy[k]/ev_erg,s->samples+84*k,secondary_xi))return 1;
        }
    }
    return 0;
}

static void secondary_channel(double rate,double energy,const double *samples,int ix,double wx,
                              const double *state,double *rates,double *loss)
{
    if(rate<=0 || energy<10.*ev_erg)return; /* Same sub-table all-heat limit. */
    double f[5]={0},v[6],w[3];
    for(int j=0;j<6;j++)v[j]=fmax(0.,(1.-wx)*samples[ix*6+j]+wx*samples[(ix+1)*6+j]);
    w[0]=v[3]*13.6;w[1]=v[4]*24.59;w[2]=v[5]*54.42;
    double total=w[0]+w[1]+w[2];
    if(total>0)for(int j=0;j<3;j++)f[j+1]=v[0]*w[j]/total;
    total=v[1]+f[1]+f[2]+f[3]+v[2];
    if(total<=0 || !isfinite(total))fatal_library_error();
    f[0]=v[1]/total;f[4]=v[2]/total;
    for(int j=1;j<=3;j++)f[j]/=total;
    secondary_target_limit(state,f);
    double power=rate*energy;
    const double thresholds[3]={13.6,24.59,54.42};
    for(int j=0;j<3;j++)rates[j]+=power*f[j+1]/(thresholds[j]*ev_erg);
    *loss+=power*(f[1]+f[2]+f[3]+f[4]);
}

/* Rates per H nucleus (s^-1), energy per H (erg/s). Identical calculation
 * supplies the chemical RHS and the subtraction from primary heating.
 * Secondary electrons do not absorb any further primary photons. */
static void secondary_budget(struct gasVariables *g,const struct globalVariables *c,double *rates,double *loss)
{
    memset(rates,0,3*sizeof(double));*loss=0;
    double state[NS];
    /* CVODE Newton trials can briefly leave the positive cone. This local
     * rate evaluation uses nonnegative reactants, not published clipping. */
    for(int s=0;s<NS;s++)state[s]=fmax(0.,g->abundances[s]);
    int ix=0;
    double xi=fmax(secondary_xi[0],fmin(secondary_xi[13],state[sp_HII]));
    while(ix<12 && xi>secondary_xi[ix+1])ix++;
    double wx=0;
    if(c->N_spectra>0)wx=(xi-secondary_xi[ix])/(secondary_xi[ix+1]-secondary_xi[ix]);
    for(int b=0;b<c->N_spectra;b++){
        double flux=fmax(0.,g->isotropic_photon_density[b])*LIGHTSPEED*g->light_speed_reduction_factor;
        if(flux==0)continue;
        for(int j=0;j<4;j++){
            const struct secondary_source *s=&secondary_sources[j];
            for(int i=0;i<s->n;i++){
                int k=b*s->stride+i;
                secondary_channel(flux*s->sigma[k]*state[s->reactant[i]],s->energy[k],
                                  s->samples+84*k,ix,wx,state,rates,loss);
            }
        }
    }
}

static void secondary_rates(struct Species_Structure *species,struct gasVariables *g,const struct globalVariables *c)
{
    double rates[3],loss;
    secondary_budget(g,c,rates,&loss);
    const int from[3]={sp_HI,sp_HeI,sp_HeII},to[3]={sp_HII,sp_HeII,sp_HeIII};
    for(int j=0;j<3;j++){
        species[from[j]].destruction_rate+=rates[j];
        species[to[j]].creation_rate+=rates[j];
        species[sp_elec].creation_rate+=rates[j];
    }
}

/* Local molecular shielding only. Ionizing attenuation and dust extinction
 * are owned by RT, so cellSelfShieldingOn stays zero. The unresolved turbulent
 * width is explicitly 1 km/s (matching the admitted CO line data). */
int snrt_chimes_molecular_factors(double temperature,double nh,double length,const double *state,double *out)
{
    if(initialized!=1 || !state || !out || !isfinite(temperature) || temperature<10 ||
        temperature>snrt_chimes_molecular_temperature_max() ||
       !isfinite(nh) || nh<=0 || !isfinite(length) || length<0)return 1;
    for(int s=0;s<NS;s++)if(!isfinite(state[s]) || state[s]<0)return 1;
    double nh2=nh*length*state[sp_H2],nco=nh*length*state[sp_CO];
    if(!isfinite(nh2) || !isfinite(nco))return 1;
    int it,ih,ib,ic;double ft,fh,fb,fc,sh2=1,sco=1;
    chimes_get_table_index(chimes_table_bins.Temperatures,chimes_table_bins.N_Temperatures,log10(temperature),&it,&ft);
    chimes_get_table_index(chimes_table_bins.H2self_column_densities,chimes_table_bins.N_H2self_column_densities,
                          log10(fmax(nh2,1e-100)),&ih,&fh);
    chimes_get_table_index(chimes_table_bins.b_turbulence,chimes_table_bins.N_b_turbulence,0.,&ib,&fb);
    if(nh2>0)sh2=pow(10.,chimes_interpol_4d_fix_x(chimes_table_H2_photodissoc.self_shielding,
        0,it,ih,ib,ft,fh,fb,chimes_table_bins.N_Temperatures,
        chimes_table_bins.N_H2self_column_densities,chimes_table_bins.N_b_turbulence));
    chimes_get_table_index(chimes_table_bins.COself_column_densities,chimes_table_bins.N_COself_column_densities,
                          log10(fmax(nco,1e-100)),&ic,&fc);
    chimes_get_table_index(chimes_table_bins.H2CO_column_densities,chimes_table_bins.N_H2CO_column_densities,
                          log10(fmax(nh2,1e-100)),&ih,&fh);
    if(nh2>0 || nco>0)sco=pow(10.,chimes_interpol_3d_fix_x(chimes_table_CO_photodissoc.self_shielding,
        0,ic,ih,fc,fh,chimes_table_bins.N_COself_column_densities,chimes_table_bins.N_H2CO_column_densities));
    const double ch=pow(10.,chimes_interpol_1d(chimes_table_H2_collis_dissoc.critical_density_H,it,ft));
    const double ch2=pow(10.,chimes_interpol_1d(chimes_table_H2_collis_dissoc.critical_density_H2,it,ft));
    const double h=state[sp_HI],h2=state[sp_H2];
    if(!isfinite(sh2) || !isfinite(sco) || !isfinite(ch) || !isfinite(ch2) || ch<=0 || ch2<=0)return 2;
    const double critical=(h+h2)>0?(h+h2)/(h/ch+h2/ch2):0.;
    const double result[3]={fmin(1.,fmax(0.,sh2)),fmin(1.,fmax(0.,sco)),nh/(nh+critical)};
    memcpy(out,result,sizeof(result));return 0;
}

static void molecular_coefficients(struct gasVariables *g,const struct globalVariables *c,struct UserData d)
{
    struct chimes_current_rates_struct *r=d.chimes_current_rates;
    /* Use the admitted primary-electron energies for both heat and its
     * secondary partition, including the base channels. No 0/0 floor heat. */
    for(int i=0;i<chimes_table_photoion_fuv.N_reactions[d.mol_flag_index];i++){
        r->photoion_fuv_heat_rate[i]=0;
        for(int b=0;b<c->N_spectra;b++){
            int k=b*chimes_table_photoion_fuv.N_reactions[1]+i;
            r->photoion_fuv_heat_rate[i]+=g->isotropic_photon_density[b]*LIGHTSPEED*g->light_speed_reduction_factor
                *chimes_table_photoion_fuv.sigmaPhot[k]*primary_heat_fuv[k];
        }
    }
    for(int i=0;i<chimes_table_photoion_euv.N_reactions[d.mol_flag_index];i++){
        r->photoion_euv_heat_rate[i]=0;
        for(int b=0;b<c->N_spectra;b++){
            int k=b*chimes_table_photoion_euv.N_reactions[1]+i;
            r->photoion_euv_heat_rate[i]+=g->isotropic_photon_density[b]*LIGHTSPEED*g->light_speed_reduction_factor
                *chimes_table_photoion_euv.sigmaPhot[k]*primary_heat_euv[k];
        }
    }
    for(int i=0;i<chimes_table_photodissoc_group1.N_reactions[d.mol_flag_index];i++)
        r->photodissoc_group1_rate_coefficient[i]=0;
    for(int i=0;i<chimes_table_photodissoc_group2.N_reactions[d.mol_flag_index];i++)
        r->photodissoc_group2_rate_coefficient[i]=0;
    double sh2=1,sco=1;
    if(d.mol_flag_index){
        double nh2=g->nH_tot*g->cell_size*fmax(g->abundances[sp_H2],0);
        double nco=g->nH_tot*g->cell_size*fmax(g->abundances[sp_CO],0);
        int it,ih,ib,ic;double ft,fh,fb,fc;
        chimes_get_table_index(chimes_table_bins.Temperatures,chimes_table_bins.N_Temperatures,log10(g->temperature),&it,&ft);
        chimes_get_table_index(chimes_table_bins.H2self_column_densities,chimes_table_bins.N_H2self_column_densities,
                              log10(fmax(nh2,1e-100)),&ih,&fh);
        chimes_get_table_index(chimes_table_bins.b_turbulence,chimes_table_bins.N_b_turbulence,0.,&ib,&fb);
        if(nh2>0)sh2=pow(10.,chimes_interpol_4d_fix_x(chimes_table_H2_photodissoc.self_shielding,
            0,it,ih,ib,ft,fh,fb,chimes_table_bins.N_Temperatures,
            chimes_table_bins.N_H2self_column_densities,chimes_table_bins.N_b_turbulence));
        chimes_get_table_index(chimes_table_bins.COself_column_densities,chimes_table_bins.N_COself_column_densities,
                              log10(fmax(nco,1e-100)),&ic,&fc);
        chimes_get_table_index(chimes_table_bins.H2CO_column_densities,chimes_table_bins.N_H2CO_column_densities,
                              log10(fmax(nh2,1e-100)),&ih,&fh);
        if(nh2>0 || nco>0)sco=pow(10.,chimes_interpol_3d_fix_x(chimes_table_CO_photodissoc.self_shielding,
            0,ic,ih,fc,fh,chimes_table_bins.N_COself_column_densities,chimes_table_bins.N_H2CO_column_densities));
        r->H2_photodissoc_shield_factor[0]=fmin(1.,fmax(0.,sh2));
        r->CO_photodissoc_shield_factor[0]=fmin(1.,fmax(0.,sco));
        r->H2_photodissoc_rate_coefficient[0]=0;r->CO_photodissoc_rate_coefficient[0]=0;
    }
    for(int b=0;b<c->N_spectra;b++){
        double flux=g->isotropic_photon_density[b]*LIGHTSPEED*g->light_speed_reduction_factor;
        double G0=flux*g->G0_parameter[b];
        /* Other CHIMES molecular channels retain their documented FUV-shape
         * approximation, but now track the evolving photons and consume them. */
        for(int i=0;i<chimes_table_photodissoc_group1.N_reactions[d.mol_flag_index];i++)
            r->photodissoc_group1_rate_coefficient[i]+=G0*chimes_table_photodissoc_group1.rates[i];
        for(int i=0;i<chimes_table_photodissoc_group2.N_reactions[d.mol_flag_index];i++)
            r->photodissoc_group2_rate_coefficient[i]+=G0*chimes_table_photodissoc_group2.rates[i];
        if(d.mol_flag_index){
            r->H2_photodissoc_rate_coefficient[0]+=flux*h2_sigma[b][1]*r->H2_photodissoc_shield_factor[0];
            r->CO_photodissoc_rate_coefficient[0]+=flux*co_sigma[b][1]*r->CO_photodissoc_shield_factor[0];
        }
    }
}

static ChimesFloat molecular_opacity(struct gasVariables *g,const struct globalVariables *c,struct UserData d,int b)
{
    (void)c;
    double opacity=0;
    for(int i=0;i<chimes_table_photodissoc_group1.N_reactions[d.mol_flag_index];i++)
        opacity+=g->G0_parameter[b]*chimes_table_photodissoc_group1.rates[i]
            *fmax(0,g->abundances[chimes_table_photodissoc_group1.reactants[i]]);
    for(int i=0;i<chimes_table_photodissoc_group2.N_reactions[d.mol_flag_index];i++)
        opacity+=g->G0_parameter[b]*chimes_table_photodissoc_group2.rates[i]
            *fmax(0,g->abundances[chimes_table_photodissoc_group2.reactants[i]]);
    if(d.mol_flag_index){
        opacity+=h2_sigma[b][0]*d.chimes_current_rates->H2_photodissoc_shield_factor[0]*fmax(0,g->abundances[sp_H2]);
        opacity+=co_sigma[b][0]*d.chimes_current_rates->CO_photodissoc_shield_factor[0]*fmax(0,g->abundances[sp_CO]);
    }
    return opacity*g->nH_tot;
}

/* Gas/grain exchange is owned by the live dust energy receiver. The current
 * dust receiver assigns all primary grain absorption to grain enthalpy:
 * disable CHIMES PE until its energy can be debited from that same ledger,
 * not because a separate PE gas heating implementation already exists.
 * CHIMES' callback units are erg cm^3 s^-1, before nH^2.
 * Grain-surface chemical reactions and their formation heat remain active.
 * With RT, CHIMES cellSelfShieldingOn is disabled: extinction=0 here. */
static ChimesFloat remove_duplicate_dust_energy(struct gasVariables *g,const struct globalVariables *c)
{
    // Record every cooling evaluation, including internal solver trials.
    // A cold step that crosses the table limit and returns below it is still
    // rejected; endpoint-only checks cannot certify the evaluated rates.
    struct cell_context *context=c->hybrid_data;
    if(context && context->molecular_temperature_max>0 &&
       (!isfinite(g->temperature) || g->temperature>context->molecular_temperature_max))
        context->outside_molecular_domain=1;
    /* Upstream only heats from base ionization channels. Include primary
     * photoelectrons of Auger-producing channels from the same shell fits.
     * Unspecified cascade electron/fluorescence energy is not invented. */
    ChimesFloat extra_heat=0;
    for(int b=0;b<c->N_spectra;b++){
        double flux=g->isotropic_photon_density[b]*LIGHTSPEED*g->light_speed_reduction_factor;
        int nf=chimes_table_photoion_auger_fuv.N_reactions[1];
        int ne=chimes_table_photoion_auger_euv.N_reactions[1];
        for(int i=0;i<chimes_table_photoion_auger_fuv.N_reactions[g->temperature<=c->Tmol_K];i++)
            extra_heat+=flux*chimes_table_photoion_auger_fuv.sigmaPhot[b*nf+i]*auger_heat_fuv[b*nf+i]
                *g->abundances[chimes_table_photoion_auger_fuv.reactants[i]]/g->nH_tot;
        for(int i=0;i<chimes_table_photoion_auger_euv.N_reactions[g->temperature<=c->Tmol_K];i++)
            extra_heat+=flux*chimes_table_photoion_auger_euv.sigmaPhot[b*ne+i]*auger_heat_euv[b*ne+i]
                *g->abundances[chimes_table_photoion_auger_euv.reactants[i]]/g->nH_tot;
    }
    double secondary[3],loss;
    secondary_budget(g,c,secondary,&loss);
    extra_heat-=loss/g->nH_tot;
    if(g->temperature>c->Tmol_K || g->dust_ratio==0)return extra_heat;
    int i,j;ChimesFloat f,h;
    chimes_get_table_index(chimes_table_bins.Temperatures,chimes_table_bins.N_Temperatures,
                          chimes_log10(g->temperature),&i,&f);
    ChimesFloat correction=chimes_exp10(chimes_interpol_1d(chimes_table_cooling.gas_grain_transfer,i,f))
        *g->dust_ratio*(g->temperature-c->grain_temperature);
    ChimesFloat G0=0,ne=g->nH_tot*g->abundances[sp_elec];
    for(int b=0;b<c->N_spectra;b++)G0+=g->isotropic_photon_density[b]*LIGHTSPEED
        *g->light_speed_reduction_factor*g->G0_parameter[b];
    if(G0>0 && ne>0){
        ChimesFloat psi=chimes_log10(chimes_max(G0*chimes_sqrt(g->temperature)
                                     /chimes_max(ne*.5,CHIMES_FLT_MIN),CHIMES_FLT_MIN));
        chimes_get_table_index(chimes_table_bins.Psi,chimes_table_bins.N_Psi,psi,&j,&h);
        correction-=chimes_exp10(chimes_interpol_2d(chimes_table_cooling.photoelectric_heating,
                             i,j,f,h,chimes_table_bins.N_Psi))*G0*g->dust_ratio/g->nH_tot;
    }
    return correction+extra_heat;
}

int snrt_chimes_initialize(const char *main_path,int nspectra,const char *paths,int stride)
{
    receiver_abi=snrt_chimes_receiver_abi();
    if(initialized || NS!=CHIMES_TOTSIZE || sizeof(ChimesFloat)!=sizeof(double) ||
       (receiver_abi!=4 && receiver_abi!=5))return 1;
    /* Lookup avoids ELF copy relocations turning an optional data symbol
     * into a mandatory runtime dependency when loading an ABI4 library. */
    if(receiver_abi==5){
        void *handle=dlopen(NULL,RTLD_NOW);
        if(!handle)return 1;
        solid_charge_hook=(solid_charge_fn *)dlsym(handle,"snrt_chimes_solid_charge");
        dlclose(handle); /* CHIMES itself remains a linked dependency. */
        if(!solid_charge_hook)return 1;
    }
    if(nspectra>0 && !snrt_chimes_fs_ready())return 1;
    if(!main_path || strlen(main_path)>=500 || (nspectra!=0 && nspectra!=NG))return 2;
    if(nspectra && (!paths || stride<1 || stride>500))return 2;
    if(file_digest(main_path,table_identity))return 2;
    char digest_hex[65];
    for(int i=0;i<32;i++)sprintf(digest_hex+2*i,"%02x",table_identity[i]);
    if(strcmp(digest_hex,"8bde78faacd59249fad5e810cae43311ed03ef09131c62b0a0a5c39e8407cb0a"))return 2;
    memset(&configuration,0,sizeof(configuration));
    strcpy(configuration.MainDataTablePath,main_path);
    strcpy(configuration.EqAbundanceTablePath,"None");
    configuration.N_spectra=nspectra;
    for(int b=0;b<nspectra;b++){
        size_t len=strnlen(paths+b*stride,(size_t)stride);
        if(len==(size_t)stride || len>=500)return 2;
        memcpy(configuration.PhotoIonTablePath[b],paths+b*stride,len+1);
        if(file_digest(configuration.PhotoIonTablePath[b],table_identity+32*(b+1)))return 2;
        hid_t file=H5Fopen(configuration.PhotoIonTablePath[b],H5F_ACC_RDONLY,H5P_DEFAULT);
        if(file<0)return 2;
        double edges[2],mean;
        int bad=read_vector(file,"Header/snrt_group_edges_ev",2,edges);
        bad|=read_vector(file,"Header/snrt_mean_energy_ev",1,&mean);
        H5Fclose(file);
        if(bad || edges[0]!=group_edges[b] || edges[1]!=group_edges[b+1] || mean!=group_energy[b])return 2;
    }
    configuration.rt_update_flux=nspectra>0;
    configuration.rt_use_on_the_spot_approx=1;
    configuration.redshift_dependent_UVB_index=-1;
    configuration.StaticMolCooling=1;
    configuration.Tmol_K=1e5;
    configuration.grain_temperature=20;
    configuration.cmb_temperature=2.727;
    configuration.relativeTolerance=1e-5;
    configuration.absoluteTolerance=1e-12;
    configuration.explicitTolerance=1e-6;
    configuration.rt_density_absoluteTolerance=1e-30;
    configuration.rt_flux_absoluteTolerance=1e-10;
    configuration.scale_metal_tolerances=1;
    for(int e=0;e<9;e++)configuration.element_included[e]=1;
    configuration.hybrid_cooling_mode=1;
    chimes_exit=fatal_library_error;
    initialized=-1; /* A failed library initialization cannot be retried. */
    init_chimes(&configuration);
    /* init_chimes intentionally resets these callback pointers. Install the
     * receiver only after initialization has finished. */
    configuration.hybrid_cooling_fn=remove_duplicate_dust_energy;
    if(receiver_abi==5)*solid_charge_hook=solid_charge;
    if(configuration.totalNumberOfSpecies!=NS)return 3;
    for(int i=0;i<NS;i++)if(configuration.speciesIndices[i]!=i)return 3;
    int nf=chimes_table_photoion_auger_fuv.N_reactions[1];
    int ne=chimes_table_photoion_auger_euv.N_reactions[1];
    auger_heat_fuv=calloc(NG*nf,sizeof(double));auger_heat_euv=calloc(NG*ne,sizeof(double));
    int pf=chimes_table_photoion_fuv.N_reactions[1],pe=chimes_table_photoion_euv.N_reactions[1];
    primary_heat_fuv=calloc(NG*pf,sizeof(double));primary_heat_euv=calloc(NG*pe,sizeof(double));
    if(!auger_heat_fuv || !auger_heat_euv || !primary_heat_fuv || !primary_heat_euv)return 3;
    for(int b=0;b<nspectra;b++){
        hid_t file=H5Fopen(configuration.PhotoIonTablePath[b],H5F_ACC_RDONLY,H5P_DEFAULT);
        if(file<0)return 3;
        int bad=read_vector(file,"snrt_primary_photoelectron_energy_erg/photoion_auger_fuv",nf,auger_heat_fuv+b*nf);
        bad|=read_vector(file,"snrt_primary_photoelectron_energy_erg/photoion_auger_euv",ne,auger_heat_euv+b*ne);
        bad|=read_vector(file,"snrt_primary_photoelectron_energy_erg/photoion_fuv",pf,primary_heat_fuv+b*pf);
        bad|=read_vector(file,"snrt_primary_photoelectron_energy_erg/photoion_euv",pe,primary_heat_euv+b*pe);
        bad|=read_vector(file,"snrt_molecule_sigma/H2",2,h2_sigma[b]);
        bad|=read_vector(file,"snrt_molecule_sigma/CO",2,co_sigma[b]);
        H5Fclose(file);
        if(bad)return 3;
        if(h2_sigma[b][1]>h2_sigma[b][0] || co_sigma[b][1]>co_sigma[b][0])return 3;
    }
    if(secondary_cache_initialize(nspectra))return 3;
    snrt_chimes_molecular_coefficients=molecular_coefficients;
    snrt_chimes_molecular_opacity=molecular_opacity;
    snrt_chimes_secondary_rates=secondary_rates;
    initialize_stoichiometry();initialized=1;
    return 0;
}

int snrt_chimes_budget(const double *abundance,double *elements,double *net_charge)
{
    if(initialized!=1)return 1;
    for(int e=0;e<11;e++)elements[e]=0;
    *net_charge=0;
    for(int s=0;s<NS;s++){
        if(!isfinite(abundance[s]) || abundance[s]<0)return 2;
        *net_charge+=charge[s]*abundance[s];
        for(int e=0;e<11;e++)elements[e]+=nuclei[s][e]*abundance[s];
    }
    return 0;
}

int snrt_chimes_neutral(const double *elements,double *abundance)
{
    if(initialized!=1)return 1;
    memset(abundance,0,NS*sizeof(double));
    for(int e=0;e<11;e++){
        if(!isfinite(elements[e]) || elements[e]<0)return 2;
        abundance[neutral[e]]=elements[e];
    }
    return 0;
}

/* Reconcile advected species with the gas-phase element ledger after neutral
 * ejecta injection or grain exchange. Molecules are never silently broken to
 * supply grain growth: the grain receiver must reserve their locked atoms.
 * Both inputs can be number densities instead of normalized abundances. */
int snrt_chimes_reconcile_charged(const double *elements,const double *old,double solid_q,double *next)
{
    double measured[11],q,work[NS],locked[11]={0};
    memcpy(next,old,NS*sizeof(double));
    if(!isfinite(solid_q))return 1;
    if(snrt_chimes_budget(old,measured,&q))return 1;
    memcpy(work,old,sizeof(work));
    for(int s=sp_H2;s<NS;s++)for(int e=0;e<11;e++)locked[e]+=nuclei[s][e]*old[s];
    for(int e=0;e<11;e++){
        if(!isfinite(elements[e]) || elements[e]<0 ||
           elements[e]<locked[e]-1e-12*fmax(elements[e],locked[e]))return 2;
        double free_old=measured[e]-locked[e],free_new=fmax(0,elements[e]-locked[e]);
        if(free_new>=free_old)work[neutral[e]]+=free_new-free_old;
        else if(free_old>0){
            double scale=free_new/free_old;
            for(int s=1;s<sp_H2;s++)if(nuclei[s][e])work[s]*=scale;
        }
    }
    work[sp_elec]=solid_q;
    for(int s=1;s<NS;s++)work[sp_elec]+=charge[s]*work[s];
    /* No physical abundance floor. At underflow, the sum of individually
     * nonnegative solver populations can require a negative subnormal
     * electron population (observed -7.9e-315). Keep EVERY ion/molecule and
     * accept only this <DBL_MIN charge residual as zero electrons. Normal
     * negative populations, even 1e-300, still reject unchanged. */
    if(work[sp_elec]<0 && work[sp_elec]>-DBL_MIN)work[sp_elec]=0;
    if(work[sp_elec]<0)return 3;
    if(snrt_chimes_budget(work,measured,&q))return 3;
    for(int e=0;e<11;e++)if(fabs(measured[e]-elements[e])>
        1e-10*fmax(elements[e],1e-30))return 3;
    memcpy(next,work,sizeof(work));
    return 0;
}

int snrt_chimes_reconcile(const double *elements,const double *old,double *next)
{
    return snrt_chimes_reconcile_charged(elements,old,0,next);
}

int snrt_chimes_locked(const double *abundance,double *elements)
{
    if(initialized!=1)return 1;
    memset(elements,0,11*sizeof(double));
    for(int s=sp_H2;s<NS;s++){
        if(!isfinite(abundance[s]) || abundance[s]<0)return 2;
        for(int e=0;e<11;e++)elements[e]+=nuclei[s][e]*abundance[s];
    }
    return 0;
}

/* Controls: nH [cm^-3], Tgas, Tdust [K], dt [s], length [cm], D/D_MW,
 * grain surface correction, CR HI rate [s^-1], chat/c. Abundances n_i/nH.
 * Return only complete trial states. Native caller owns their publication. */
static int cell_charged(const double *controls,const double *elements,const double *old_abundance,
                     double solid_q,
                     const double *old_photons,double *temperature,double *abundance,double *photons,int cold_spectral)
{
    *temperature=controls[1];
    memcpy(abundance,old_abundance,NS*sizeof(double));
    memcpy(photons,old_photons,NG*sizeof(double));
    if(initialized!=1)return 1;
    if(!isfinite(solid_q))return 2;
    if(solid_q!=0 && receiver_abi<5)return 4;
    for(int k=0;k<9;k++)if(!isfinite(controls[k]) || controls[k]<0)return 2;
    if(controls[0]<=0 || controls[1]<10 || controls[1]>1e9 || controls[2]<1 || controls[2]>1e4 ||
       controls[8]<=0 || controls[8]>1)return 2;
    double measured[11],q;
    if(snrt_chimes_budget(old_abundance,measured,&q))return 3;
    for(int e=0;e<11;e++)if(!isfinite(elements[e]) || elements[e]<0 ||
        fabs(measured[e]-elements[e])>1e-8*fmax(elements[e],1e-20))return 3;
    if(fabs(elements[0]-1)>1e-12 || fabs(q+solid_q)>1e-10)return 3;
    for(int b=0;b<NG;b++)if(!isfinite(old_photons[b]) || old_photons[b]<0)return 3;
    if(configuration.N_spectra==0){
        for(int b=0;b<NG;b++)if(old_photons[b]!=0)return 3;
        if(controls[8]!=1)return 3;
    }
    if(controls[3]==0)return 0;
    struct globalVariables c=configuration; /* Tdust is local, thread-safe. */
    if(cold_spectral){
        if(controls[1]>snrt_chimes_molecular_temperature_max())return 2;
        // Resolve the molecular photo/dark split below the default grey
        // network's thermal truncation error; do not change global tables,
        // tolerances or Tmol_K for other receivers.
        c.relativeTolerance=1e-10;c.absoluteTolerance=1e-17;c.explicitTolerance=1e-10;
    }
    struct cell_context context={solid_q,cold_spectral?snrt_chimes_molecular_temperature_max():0,0};
    c.hybrid_data=&context;
    struct gasVariables g={0};
    double work[NS],radiation[NG],reduction[NG],G0[NG],dissoc[NG],source[NG]={0};
    memcpy(work,old_abundance,sizeof(work));memcpy(radiation,old_photons,sizeof(radiation));
    c.grain_temperature=controls[2];
    g.nH_tot=controls[0];g.temperature=controls[1];g.hydro_timestep=controls[3];
    g.cell_size=controls[4];g.dust_ratio=controls[5];g.dust_boost_factor=controls[6];
    g.cr_rate=controls[7];g.light_speed_reduction_factor=controls[8];
    memcpy(g.element_abundances,elements+1,10*sizeof(double));
    g.TempFloor=10;g.ThermEvolOn=1;g.temp_floor_mode=0;
    g.abundances=work;g.isotropic_photon_density=radiation;g.flux_reduction_factor=reduction;
    g.G0_parameter=G0;g.H2_dissocJ=dissoc;g.rt_source_term=source;
    for(int b=0;b<c.N_spectra;b++){
        reduction[b]=1;
        G0[b]=chimes_table_spectra.G0_parameter[b];
        dissoc[b]=chimes_table_spectra.H2_dissocJ[b];
    }
    /* Trace O injected into C-rich gas can exceed the nucleus tolerance
     * even when CVODE satisfies each species' local error estimate. Retry
     * only these failed cells from the untouched input, at tighter solver
     * tolerances. The conservation acceptance criterion never changes. */
    const struct gasVariables initial_g=g;
    int status=0;
    for(int attempt=0;attempt<3;attempt++){
        g=initial_g;
        memcpy(work,old_abundance,sizeof(work));
        memcpy(radiation,old_photons,sizeof(radiation));
        for(int b=0;b<c.N_spectra;b++)reduction[b]=1;
        status=chimes_network(&g,&c);
        if(context.outside_molecular_domain)return 50;
        if(status || !isfinite(g.temperature) || g.temperature<10 || g.temperature>1e9)break;
        if(snrt_chimes_budget(work,measured,&q))break;
        int bad=fabs(q+solid_q)>1e-8;
        for(int e=0;e<11;e++)
            if(fabs(measured[e]-elements[e])>1e-5*fmax(elements[e],1e-20))bad=1;
        if(!bad)break;
        c.relativeTolerance*=.01;
        c.absoluteTolerance*=.01;
        c.explicitTolerance*=.01;
    }
    /* Distinct failure codes make a rejected live cell diagnosable without
     * dumping its 157-species state or weakening conservation checks. */
    if(status)return 41;
    if(!isfinite(g.temperature) || g.temperature<10 || g.temperature>1e9)return 42;
    if(snrt_chimes_budget(work,measured,&q))return 43;
    for(int e=0;e<11;e++)if(fabs(measured[e]-elements[e])>1e-5*fmax(elements[e],1e-20)){
        static int reported=0;
        flockfile(stderr);
        if(!reported){
            fprintf(stderr,"CHIMES nucleus rejection: element=%d expected=%.17g actual=%.17g\n",
                    e,elements[e],measured[e]);
            reported=1;
        }
        funlockfile(stderr);
        return 44;
    }
    if(fabs(q+solid_q)>1e-8)return 45;
    /* Publish a state satisfying the next step's stricter input contract.
     * Preserve thermal energy when this roundoff correction changes n_tot. */
    double corrected[NS],before=0,after=0;
    int reconcile_status=snrt_chimes_reconcile_charged(elements,work,solid_q,corrected);
    if(reconcile_status){
        static int reported=0;
        flockfile(stderr);
        if(!reported){
            double locked[11]={0},electron=solid_q;
            for(int s=1;s<NS;s++)electron+=charge[s]*work[s];
            for(int s=sp_H2;s<NS;s++)for(int e=0;e<11;e++)locked[e]+=nuclei[s][e]*work[s];
            fprintf(stderr,"CHIMES reconciliation rejected: subcode=%d electron=%.17g\n",reconcile_status,electron);
            for(int e=0;e<11;e++)if(locked[e]>elements[e])
                fprintf(stderr,"CHIMES locked excess: element=%d target=%.17g locked=%.17g\n",e,elements[e],locked[e]);
            reported=1;
        }
        funlockfile(stderr);return 46;
    }
    for(int s=0;s<NS;s++){before+=work[s];after+=corrected[s];}
    if(after<=0)return 47;
    g.temperature*=before/after;
    if(g.temperature<10 && g.temperature>=10*(1-1e-10))g.temperature=10;
    if(!isfinite(g.temperature) || g.temperature<10 || g.temperature>1e9)return 48;
    memcpy(work,corrected,sizeof(work));
    for(int b=0;b<c.N_spectra;b++)if(!isfinite(radiation[b]) || radiation[b]<0 ||
        radiation[b]>old_photons[b]*(1+1e-5)+1e-30)return 49;
    /* Case B, zero source term: a CVODE logarithm floor is not emission. */
    for(int b=0;b<c.N_spectra;b++)if(old_photons[b]==0)radiation[b]=0;
    *temperature=g.temperature;
    memcpy(abundance,work,sizeof(work));memcpy(photons,radiation,sizeof(radiation));
    return 0;
}

int snrt_chimes_cell_charged(const double *controls,const double *elements,const double *old,
    double solid_q,const double *number,double *temperature,double *next,double *next_number)
{
    return cell_charged(controls,elements,old,solid_q,number,temperature,next,next_number,0);
}

int snrt_chimes_cell_cold_dark(const double *controls,const double *elements,const double *old,
    double *temperature,double *next)
{
    double zero[NG]={0},out[NG];
    return cell_charged(controls,elements,old,0,zero,temperature,next,out,1);
}

int snrt_chimes_cell(const double *controls,const double *elements,const double *old_abundance,
                     const double *old_photons,double *temperature,double *abundance,double *photons)
{
    return snrt_chimes_cell_charged(controls,elements,old_abundance,0,
                                  old_photons,temperature,abundance,photons);
}
