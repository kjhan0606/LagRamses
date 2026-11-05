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
#include "snrt_group_profile.h"
#include "snrt_chimes_atomization.h"
#ifdef SNRT_CHIMES_CUDA
#include "snrt_chimes_cvode_cuda.h"
#include "snrt_chimes_rhs_cuda.h"
#endif

#define NS 157
#define NG SNRT_GROUP_COUNT
#if SNRT_GROUP_COUNT > CHIMES_MAX_UV_SPECTRA
#error "SNRT group profile exceeds CHIMES_MAX_UV_SPECTRA; increase CHIMES capacity and rebuild all ABI consumers"
#endif
#define group_edges snrt_profile_group_edges
extern int snrt_chimes_receiver_abi(void);
/* ABI 4 remains usable for the existing neutral-solid receiver. ABI 5 adds
 * this optional hook without changing either CHIMES structure layout. */
typedef ChimesFloat (*solid_charge_fn)(const struct globalVariables *);
static solid_charge_fn *solid_charge_hook;
typedef double (*thermal_ceiling_fn)(const struct globalVariables *);
typedef void (*thermal_root_fn)(const struct globalVariables *,double);
static thermal_ceiling_fn *thermal_ceiling_hook;
static thermal_root_fn *thermal_root_hook;
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
static double group_energy[NG];
static int nuclei[NS][11],charge[NS];
static const int neutral[11]={sp_HI,sp_HeI,sp_CI,sp_NI,sp_OI,sp_NeI,sp_MgI,sp_SiI,sp_SI,sp_CaI,sp_FeI};
static const int atomic_z[11]={1,2,6,7,8,10,12,14,16,20,26};
int snrt_chimes_budget(const double *,double *,double *);

/* Production-facing selector for the complete CHIMES cell solve. The two
 * lower-level switches remain available for isolated backend experiments. */
static int resolve_compute_backend(int *selected)
{
    *selected=-1;
    const char *value=getenv("SNRT_CHIMES_COMPUTE_BACKEND");
    if(!value || !*value)return 0;

    const char *cvode_expected=NULL,*rhs_expected=NULL;
    if(!strcmp(value,"cpu")) {
        *selected=0;cvode_expected="cpu";rhs_expected="cpu";
    } else if(!strcmp(value,"hybrid")) {
        *selected=2;cvode_expected="cpu";rhs_expected="cuda_dark_reactions";
        const char *scheduler=getenv("SNRT_CHIMES_CELL_SCHEDULER");
        if(scheduler && *scheduler && strcmp(scheduler,"level_queue")) {
            fprintf(stderr,"CHIMES hybrid requires level_queue scheduling\n");
            return 1;
        }
    } else if(!strcmp(value,"cuda_integrated")) {
        *selected=1;cvode_expected="cuda_batched_lu";rhs_expected="cuda_dark_reactions";
    } else {
        fprintf(stderr,"Unknown SNRT_CHIMES_COMPUTE_BACKEND=%s; expected cpu|hybrid (cuda_integrated is experimental)\n",value);
        return 1;
    }

    const char *old_cvode=getenv("SNRT_CHIMES_CVODE_BACKEND");
    const char *old_rhs=getenv("SNRT_CHIMES_RHS_BACKEND");
    if((old_cvode && *old_cvode && strcmp(old_cvode,cvode_expected)) ||
       (old_rhs && *old_rhs && strcmp(old_rhs,rhs_expected))) {
        fprintf(stderr,"SNRT_CHIMES_COMPUTE_BACKEND conflicts with legacy CHIMES backend selectors\n");
        return 1;
    }
    return 0;
}

/* Serial startup query: use the same selector admission as initialization.
 * It does not configure CHIMES or create a device context. */
int snrt_chimes_cuda_requested(void)
{
    int selected;
    if(resolve_compute_backend(&selected))return -1;
    if(selected>=0)return selected!=0;
    const char *cvode=getenv("SNRT_CHIMES_CVODE_BACKEND");
    const char *rhs=getenv("SNRT_CHIMES_RHS_BACKEND");
    return (cvode && *cvode && strcmp(cvode,"cpu")) ||
           (rhs && *rhs && strcmp(rhs,"cpu"));
}

struct cell_context {
    double solid_q;
    double molecular_temperature_max; /* zero for unchanged legacy/hot calls */
    int outside_molecular_domain;
    double event_ceiling,event_time;
};

#ifdef SNRT_CHIMES_JACOBIAN_TESTING
/* Private copies for the bounded parallel-Jacobian experiment. In particular,
 * the molecular-domain diagnostic must not be written by concurrent workers. */
int snrt_chimes_jac_context_clone(const struct globalVariables *src,
                                 struct globalVariables *dst)
{
    *dst=*src;
    if(!src->hybrid_data)return 0;
    struct cell_context *copy=malloc(sizeof(*copy));
    if(!copy)return -1;
    *copy=*(const struct cell_context *)src->hybrid_data;
    dst->hybrid_data=copy;
    return 0;
}
void snrt_chimes_jac_context_release(const struct globalVariables *src,
                                   struct globalVariables *copy)
{
    /* Trial flags are not published: the experiment performs the reference
     * RHS/serial fallback on the owner for side effects before returning. */
    if(copy->hybrid_data!=src->hybrid_data)free(copy->hybrid_data);
    copy->hybrid_data=NULL;
}
int snrt_chimes_jac_context_flag(const struct globalVariables *copy)
{
    const struct cell_context *p=copy->hybrid_data;
    return p && p->outside_molecular_domain;
}
#endif

static double thermal_ceiling(const struct globalVariables *c)
{
    const struct cell_context *p=c->hybrid_data;
    return p?p->event_ceiling:0;
}
static void thermal_root(const struct globalVariables *c,double time)
{
    struct cell_context *p=c->hybrid_data;
    if(p)p->event_time=time;
}
int snrt_chimes_transition_supported(void) { return initialized==1 && receiver_abi==6; }

/* Energy/nucleus/charge preserving instantaneous atomization. All arrays
 * stay untouched on failure. Chemical energy cost is per gas H nucleus. */
int snrt_chimes_atomize(double temperature,const double *old,double *next,double *new_temperature,double *cost_ev)
{
    if(initialized!=1 || !isfinite(temperature) || temperature<10 || temperature>1e9)return 1;
    double e0[11],e1[11],q0,q1,work[NS],cost=0,before=0,after=0;
    if(snrt_chimes_budget(old,e0,&q0))return 1;
    memcpy(work,old,sizeof(work));
    const int element[3]={0,2,4};
    for(int s=sp_H2;s<NS;s++) {
        double channel=-atomization_hf0[s-sp_H2];
        int ion=-1;
        for(int k=0;k<3;k++) {
            int count=nuclei[s][element[k]];
            work[neutral[element[k]]]+=count*old[s];
            channel+=count*atomization_atom_hf0[k];
            if(count && (ion<0 || atomization_ion_hf0[k]-atomization_atom_hf0[k]
                <atomization_ion_hf0[ion]-atomization_atom_hf0[ion]))ion=k;
        }
        if(charge[s]==1 && ion>=0) {
            work[neutral[element[ion]]]-=old[s];
            work[neutral[element[ion]]+1]+=old[s];
            channel+=atomization_ion_hf0[ion]-atomization_atom_hf0[ion];
        } else if(charge[s]!=0)return 2;
        if(channel<=0 || !isfinite(channel))return 2;
        cost+=old[s]*channel*atomization_kjmol_to_ev;
        work[s]=0;
    }
    if(snrt_chimes_budget(work,e1,&q1))return 2;
    for(int e=0;e<11;e++)if(fabs(e1[e]-e0[e])>1e-12*fmax(e0[e],1e-30))return 2;
    if(fabs(q1-q0)>1e-12 || work[sp_elec]!=old[sp_elec])return 2;
    for(int s=0;s<NS;s++){before+=old[s];after+=work[s];}
    if(before<=0 || after<=0)return 2;
    double t=(1.5*BOLTZMANNCGS*before*temperature-cost*ev_erg)/(1.5*BOLTZMANNCGS*after);
    if(!isfinite(t) || t<10 || t>1e9)return 3;
    memcpy(next,work,sizeof(work));*new_temperature=t;*cost_ev=cost;return 0;
}

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

static int read_string_attribute(hid_t file,const char *name,char *value,size_t capacity)
{
    if(!name || !value || capacity<2 || H5Aexists(file,name)<=0)return 1;
    hid_t attr=H5Aopen(file,name,H5P_DEFAULT);
    if(attr<0)return 1;
    hid_t type=H5Aget_type(attr);
    if(type<0){H5Aclose(attr);return 1;}
    int ok=0;
    if(H5Tget_class(type)==H5T_STRING && H5Tis_variable_str(type)>0){
        char *raw=NULL;
        if(H5Aread(attr,type,&raw)>=0 && raw && strlen(raw)<capacity){
            memcpy(value,raw,strlen(raw)+1);ok=1;
        }
        if(raw)H5free_memory(raw);
    }else if(H5Tget_class(type)==H5T_STRING){
        const size_t n=H5Tget_size(type);
        if(n>0 && n<capacity){
            char raw[512]={0};
            if(n<sizeof(raw) && H5Aread(attr,type,raw)>=0){
                raw[n]='\0';memcpy(value,raw,n+1);ok=1;
            }
        }
    }
    H5Tclose(type);H5Aclose(attr);return !ok;
}

int snrt_chimes_identity(double *values)
{
    if(initialized!=1)return 1;
    for(int i=0;i<32*(NG+1);i++)values[i]=table_identity[i];
    return 0;
}

double snrt_chimes_boltzmann(void) { return BOLTZMANNCGS; }
/* Set once on the serial level boundary, before any OpenMP cell callbacks.
 * Each cell copies configuration privately below. Never call this setter
 * concurrently with a cell solve. a=1 retains the historical noncosmo value. */
int snrt_chimes_set_expansion(double a)
{
    if(initialized!=1 || !isfinite(a) || a<=0)return 1;
    const double temperature=2.727/a;
    if(!isfinite(temperature) || temperature<=0)return 1;
    configuration.cmb_temperature=temperature;
    return 0;
}
double snrt_chimes_cmb_temperature(void)
{
    return initialized==1 ? configuration.cmb_temperature : NAN;
}
int snrt_chimes_charge_supported(void) { return initialized==1 && receiver_abi>=5; }

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
    if(c->N_spectra==0)return; /* No primary photoelectron source. */
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

#ifdef SNRT_CHIMES_JACOBIAN_TESTING
/* The reaction dependency experiment must not silently omit an unknown
 * extension. This registered photoelectron hook is a no-op in the dark. */
int snrt_chimes_jac_secondary_is_dark_noop(const struct globalVariables *c)
{
    return c->N_spectra==0 && snrt_chimes_secondary_rates==secondary_rates;
}
#endif

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
    if(c->N_spectra==0){
        /* Retain exact-zero photorates; skip only photon shielding, never
         * molecular cooling columns or collisional dissociation. */
        if(d.mol_flag_index){
            r->H2_photodissoc_rate_coefficient[0]=0;
            r->CO_photodissoc_rate_coefficient[0]=0;
            r->H2_photodissoc_shield_factor[0]=1;
            r->CO_photodissoc_shield_factor[0]=1;
        }
        return;
    }
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
    if(initialized)return 1;
    int compute_backend=-1;
    if(resolve_compute_backend(&compute_backend))return 5;
    int unified_cuda=compute_backend==1;
    const char *cvode_backend=getenv("SNRT_CHIMES_CVODE_BACKEND");
    int cvode_cuda=compute_backend>=0?unified_cuda:
        (cvode_backend && *cvode_backend && strcmp(cvode_backend,"cpu"));
    const char *rhs_backend=getenv("SNRT_CHIMES_RHS_BACKEND");
    int rhs_cuda=compute_backend>=0?(compute_backend!=0):
        (rhs_backend && *rhs_backend && strcmp(rhs_backend,"cpu"));
    if(cvode_cuda) {
#ifdef SNRT_CHIMES_CUDA
        /* Reject a stale external CHIMES library instead of silently running
         * CPU under a requested GPU mode. No change to the structure ABI. */
        void *factory_handle=dlopen(NULL,RTLD_NOW);
        if(!factory_handle)return 5;
        int (*factory_abi)(void)=(int (*)(void))dlsym(factory_handle,"snrt_chimes_cvode_factory_abi");
        int factory_ok=factory_abi && factory_abi()==1;
        dlclose(factory_handle);
        if(!factory_ok)return 5;
#else
        fprintf(stderr,"CHIMES CUDA backend requested but CHIMES_CUDA=1 is absent\n");
        return 5;
#endif
    }
    if(rhs_cuda) {
#ifdef SNRT_CHIMES_CUDA
        void *handle=dlopen(NULL,RTLD_NOW);
        if(!handle)return 5;
        int (*abi)(void)=(int (*)(void))dlsym(handle,"snrt_chimes_dark_rhs_abi");
        int ok=abi && abi()==1;dlclose(handle);
        if(!ok)return 5;
#else
        fprintf(stderr,"CHIMES device RHS requires CHIMES_CUDA=1\n");return 5;
#endif
    }
#ifdef SNRT_CHIMES_CUDA
    if(compute_backend>=0) {
        if(snrt_chimes_cvode_configure_mode(unified_cuda) ||
           snrt_chimes_rhs_configure_mode(rhs_cuda))return 5;
        fprintf(stderr,"CHIMES integrated compute backend: %s\n",
            compute_backend==2?"hybrid (host CVODE/LU, asynchronous device RHS/cooling)":
            unified_cuda?"cuda_integrated":"cpu");
    } else if(snrt_chimes_cvode_configure() || snrt_chimes_rhs_configure())return 5;
#endif
    receiver_abi=snrt_chimes_receiver_abi();
    if(initialized || NS!=CHIMES_TOTSIZE || sizeof(ChimesFloat)!=sizeof(double) ||
       (receiver_abi!=4 && receiver_abi!=5 && receiver_abi!=6))return 1;
    /* Lookup avoids ELF copy relocations turning an optional data symbol
     * into a mandatory runtime dependency when loading an ABI4 library. */
    if(receiver_abi>=5){
        void *handle=dlopen(NULL,RTLD_NOW);
        if(!handle)return 1;
        solid_charge_hook=(solid_charge_fn *)dlsym(handle,"snrt_chimes_solid_charge");
        if(receiver_abi==6){
            thermal_ceiling_hook=(thermal_ceiling_fn *)dlsym(handle,"snrt_chimes_thermal_ceiling");
            thermal_root_hook=(thermal_root_fn *)dlsym(handle,"snrt_chimes_thermal_root_found");
        }
        dlclose(handle); /* CHIMES itself remains a linked dependency. */
        if(!solid_charge_hook)return 1;
        if(receiver_abi==6 && (!thermal_ceiling_hook || !thermal_root_hook))return 1;
    }
    if(nspectra>0 && !snrt_chimes_fs_ready())return 1;
    if(nspectra<0 || nspectra>CHIMES_MAX_UV_SPECTRA){
        fprintf(stderr,"CHIMES spectrum count %d exceeds supported range [0,%d]\n",
                nspectra,CHIMES_MAX_UV_SPECTRA);
        return 2;
    }
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
        if(SNRT_GROUP_COUNT!=9 || H5Aexists(file,"snrt_group_set_id")>0){
            char profile_id[128]={0},edge_digest[128]={0};
            bad|=read_string_attribute(file,"snrt_group_set_id",profile_id,sizeof(profile_id));
            bad|=read_string_attribute(file,"snrt_group_edges_sha256",edge_digest,sizeof(edge_digest));
            if(strcmp(profile_id,SNRT_GROUP_SET_ID) || strcmp(edge_digest,SNRT_GROUP_EDGES_SHA256))bad=1;
        }
        H5Fclose(file);
        if(bad || edges[0]!=group_edges[b] || edges[1]!=group_edges[b+1] ||
           mean<edges[0] || mean>edges[1])return 2;
        group_energy[b]=mean;
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
    if(configuration.totalNumberOfSpecies!=NS)return 3;
    fprintf(stderr,"CHIMES cross-section initialization complete spectra=%d capacity=%d first_table=%s\n",
            configuration.N_spectra,CHIMES_MAX_UV_SPECTRA,
            configuration.N_spectra>0?configuration.PhotoIonTablePath[0]:"none");
    /* init_chimes intentionally resets these callback pointers. Install the
     * receiver only after initialization has finished. */
    configuration.hybrid_cooling_fn=remove_duplicate_dust_energy;
    if(receiver_abi>=5)*solid_charge_hook=solid_charge;
    if(receiver_abi==6){*thermal_ceiling_hook=thermal_ceiling;*thermal_root_hook=thermal_root;}
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
static void report_reconcile_failure(const char *branch,const double *elements,const double *old,
                                     const double *state,double solid_q)
{
    static int reported=0;
    long double old_pos=0,old_neg=0,old_net=0,post_pos=0,post_neg=0,post_net=0;
    long double derived_e=solid_q,charge_scale=fabsl((long double)solid_q),max_element_rel=0;
    for(int s=0;s<NS;s++){
        long double old_term=(long double)charge[s]*old[s];
        long double post_term=(long double)charge[s]*state[s];
        old_net+=old_term;post_net+=post_term;
        if(old_term>0)old_pos+=old_term;else old_neg+=old_term;
        if(post_term>0)post_pos+=post_term;else post_neg+=post_term;
        if(s!=sp_elec){derived_e+=post_term;charge_scale+=fabsl(post_term);}
    }
    for(int e=0;e<11;e++){
        long double measured=0,target=elements[e];
        for(int s=0;s<NS;s++)measured+=(long double)nuclei[s][e]*state[s];
        long double relative=fabsl(measured-target)/fmaxl(fabsl(target),1e-30L);
        if(relative>max_element_rel)max_element_rel=relative;
    }
    flockfile(stderr);
    if(!reported){
        fprintf(stderr,"CHIMES reconcile failure detail branch=%s solid_q=%.17g old_e=%.17g ",
                branch,solid_q,old[sp_elec]);
        fprintf(stderr,"old_charge(pos/neg/net)=%.21Lg/%.21Lg/%.21Lg ",old_pos,old_neg,old_net);
        fprintf(stderr,"post_charge(pos/neg/net)=%.21Lg/%.21Lg/%.21Lg ",post_pos,post_neg,post_net);
        fprintf(stderr,"derived_e=%.21Lg charge_scale=%.21Lg ",derived_e,charge_scale);
        fprintf(stderr,"old_Hm/Cm/Om=%.17g/%.17g/%.17g post_Hm/Cm/Om=%.17g/%.17g/%.17g ",
                old[sp_Hm],old[sp_Cm],old[sp_Om],state[sp_Hm],state[sp_Cm],state[sp_Om]);
        fprintf(stderr,"max_element_relative_error=%.21Lg\n",max_element_rel);
        reported=1;
    }
    funlockfile(stderr);
}

static int reconcile_charged_impl(const double *elements,const double *negative_scale,const double *old,
                                  double solid_q,double *next,int transport_roundoff)
{
    double measured[11],q,work[NS],locked[11]={0};
    memcpy(next,old,NS*sizeof(double));
    if(!isfinite(solid_q))return 1;
    memcpy(work,old,sizeof(work));
    for(int e=0;e<11;e++)if(!isfinite(elements[e]) || elements[e]<0)return 2;
    if(transport_roundoff){
        if(!negative_scale)return 2;
        for(int e=0;e<11;e++)if(!isfinite(negative_scale[e]) || negative_scale[e]<0)return 2;
    }
    /* Direct reconciliation keeps its gas-phase-reservoir bound. At the
     * transport/preparation boundary, however, gas-phase depletion into
     * grains can make that reservoir zero while the total element inventory
     * is finite. Permit only atomic roundoff residues there, bounded by the
     * summed elemental debit and total (gas+solid) inventory. */
    double negative_debit[11]={0};
    for(int s=0;s<NS;s++){
        if(!isfinite(work[s]))return 1;
        if(work[s]>=0)continue;
        if(s==sp_elec)return 1;
        int has_nucleus=0,n_element=0,n_total=0;
        for(int e=0;e<11;e++)if(nuclei[s][e]>0){
            has_nucleus=1;
            n_element++;
            n_total+=nuclei[s][e];
            double debit=-work[s]*nuclei[s][e];
            if(!isfinite(debit))return 1;
            if(transport_roundoff){
                negative_debit[e]+=debit;
                if(!isfinite(negative_debit[e]))return 1;
            }else{
                double tolerance=512*DBL_EPSILON*fmax(elements[e],DBL_MIN);
                if(debit>tolerance)return 1;
            }
        }
        if(!has_nucleus)return 1;
        if(transport_roundoff){
            /* H2 and all multi-nucleus molecules retain strict rejection. */
            if(n_element!=1 || n_total!=1)return 1;
        }else{
            work[s]=0;
        }
    }
    if(transport_roundoff){
        for(int e=0;e<11;e++){
            double tolerance=512*DBL_EPSILON*fmax(negative_scale[e],DBL_MIN);
            if(negative_debit[e]>tolerance)return 1;
        }
        for(int s=0;s<NS;s++)if(work[s]<0)work[s]=0;
    }
    if(snrt_chimes_budget(work,measured,&q))return 1;
    for(int s=sp_H2;s<NS;s++)for(int e=0;e<11;e++)locked[e]+=nuclei[s][e]*work[s];
    /* A gas element can be exactly depleted in the mass carrier while the
     * transported molecular state retains a sub-ULP positive tail. Such a
     * tail cannot be represented in total-minus-solid mass, but the strict
     * molecular lock below would reject every positive value against a zero
     * reservoir. Remove it only when its entire locked inventory is within
     * four ULPs of a nonzero total (gas+solid) element carrier. Transfer the
     * other nuclei in those trace molecules to neutral atoms, preserving
     * their elemental inventories; the depleted element itself has no
     * representable gas reservoir. Larger molecular excesses still reject.
     */
    if(transport_roundoff){
        int depleted[11]={0},have_depleted=0;
        for(int e=0;e<11;e++){
            if(elements[e]!=0.0 || negative_scale[e]<=0.0 || locked[e]<=0.0)continue;
            const double ulp=nextafter(negative_scale[e],INFINITY)-negative_scale[e];
            const double bound=4.0*ulp;
            if(!isfinite(bound) || locked[e]>bound)return 2;
            depleted[e]=1;have_depleted=1;
        }
        if(have_depleted){
            for(int s=sp_H2;s<NS;s++){
                int touches_depleted=0;
                for(int e=0;e<11;e++)if(depleted[e] && nuclei[s][e]>0){
                    touches_depleted=1;break;
                }
                if(!touches_depleted || work[s]==0.0)continue;
                const double population=work[s];
                for(int e=0;e<11;e++)if(!depleted[e] && nuclei[s][e]>0)
                    work[neutral[e]]+=population*nuclei[s][e];
                work[s]=0.0;
            }
            if(snrt_chimes_budget(work,measured,&q))return 1;
            memset(locked,0,sizeof(locked));
            for(int s=sp_H2;s<NS;s++)for(int e=0;e<11;e++)locked[e]+=nuclei[s][e]*work[s];
        }
    }
    /* Lock only the accepted post-roundoff species state. Using old here
     * would let a clipped negative molecule perturb the amount of atomic
     * material available for grain exchange, leaving a small element-ledger
     * drift when the requested gas-phase reservoir also changes. */
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
    if(work[sp_elec]<0){
        report_reconcile_failure("negative_electron",elements,old,work,solid_q);
        return 3;
    }
    if(snrt_chimes_budget(work,measured,&q)){
        report_reconcile_failure("corrected_budget",elements,old,work,solid_q);
        return 3;
    }
    for(int e=0;e<11;e++)if(fabs(measured[e]-elements[e])>
        1e-10*fmax(elements[e],1e-30)){
        report_reconcile_failure("element_closure",elements,old,work,solid_q);
        return 3;
    }
    memcpy(next,work,sizeof(work));
    return 0;
}

int snrt_chimes_reconcile_charged(const double *elements,const double *old,double solid_q,double *next)
{
    return reconcile_charged_impl(elements,elements,old,solid_q,next,0);
}

/* The mode-3 dark transition can leave solver-scale H- and C- populations
 * after its algebraic electron abundance reaches zero. This private-state
 * projection is deliberately narrower than the public strict reconciler:
 * only an exactly dark, uncharged input and only tolerance-bounded H-/C-
 * tails are eligible. Return 0 only after the unchanged strict reconciler
 * accepts the projected state; otherwise leave `next` untouched. */
int snrt_chimes_reconcile_dark_mode3(const double *elements,const double *state,
    double solid_q,int original_zero_incident,double accepted_absolute_tolerance,
    int scale_metal_tolerances,double *next)
{
    if(initialized!=1 || !elements || !state || !next ||
       original_zero_incident!=1 || !isfinite(solid_q) || solid_q!=0.0 ||
       !isfinite(accepted_absolute_tolerance) || accepted_absolute_tolerance<=0.0 ||
       (scale_metal_tolerances!=0 && scale_metal_tolerances!=1))return 1;
    for(int e=0;e<11;e++)if(!isfinite(elements[e]) || elements[e]<0.0)return 1;
    for(int s=0;s<NS;s++)if(!isfinite(state[s]) || state[s]<0.0)return 1;
    if(state[sp_elec]!=0.0)return 1;

    const double hminus=state[sp_Hm],cminus=state[sp_Cm];
    if((hminus>0.0 && state[sp_HI]<=0.0) ||
       (cminus>0.0 && state[sp_CI]<=0.0))return 1;

    /* Reject every unsupported anion, and sum the remaining positive charge
     * without cancellation. The negative required electron must be fully
     * explained by the two permitted tails. */
    long double positive_charge=0.0L;
    for(int s=1;s<NS;s++){
        if(s==sp_Hm || s==sp_Cm || charge[s]==0 || state[s]==0.0)continue;
        if(charge[s]<0)return 1;
        positive_charge+=(long double)charge[s]*(long double)state[s];
    }
    const long double negative_tail=(long double)hminus+(long double)cminus;
    const long double required_electron=positive_charge-negative_tail;
    if(!isfinite((double)positive_charge) || required_electron>=0.0L ||
       positive_charge>negative_tail)return 1;

    const double h_tolerance=accepted_absolute_tolerance*
        (scale_metal_tolerances==1?elements[0]:1.0);
    const double c_tolerance=accepted_absolute_tolerance*
        (scale_metal_tolerances==1?elements[2]:1.0);
    if(!isfinite(h_tolerance) || !isfinite(c_tolerance))return 1;
    const double h_cap=fmin(.01*h_tolerance,512.0*DBL_EPSILON*elements[0]);
    const double c_cap=fmin(.01*c_tolerance,512.0*DBL_EPSILON*elements[2]);
    const double electron_cap=.01*accepted_absolute_tolerance;
    if(!isfinite(electron_cap) || hminus>h_cap || cminus>c_cap ||
       negative_tail>(long double)electron_cap)return 1;

    double before_elements[11],after_elements[11],state_charge,after_charge;
    if(snrt_chimes_budget(state,before_elements,&state_charge) ||
       !isfinite(state_charge) || state_charge>=0.0)return 1;
    double candidate[NS],corrected[NS];
    memcpy(candidate,state,sizeof(candidate));
    candidate[sp_HI]+=hminus;candidate[sp_Hm]=0.0;
    candidate[sp_CI]+=cminus;candidate[sp_Cm]=0.0;
    for(int s=0;s<NS;s++)if(!isfinite(candidate[s]) || candidate[s]<0.0)return 2;
    if(snrt_chimes_budget(candidate,after_elements,&after_charge))return 2;
    for(int e=0;e<11;e++){
        const double scale=fmax(fmax(before_elements[e],after_elements[e]),DBL_MIN);
        if(fabs(after_elements[e]-before_elements[e])>64.0*DBL_EPSILON*scale)return 2;
    }

    /* Keep the public acceptance rule unchanged. In particular, this helper
     * cannot broaden the element-ledger tolerance or repair the original state. */
    if(snrt_chimes_reconcile_charged(elements,candidate,solid_q,corrected))return 2;
    if(snrt_chimes_budget(corrected,after_elements,&after_charge))return 2;
    long double charge_scale=0.0L,stable_charge=(long double)solid_q-(long double)corrected[sp_elec];
    charge_scale=fabsl((long double)corrected[sp_elec]);
    for(int s=1;s<NS;s++){
        const long double term=(long double)charge[s]*(long double)corrected[s];
        stable_charge+=term;charge_scale+=fabsl(term);
    }
    if(!isfinite(after_charge) || !isfinite((double)stable_charge) ||
       fabs(after_charge)>512.0*DBL_EPSILON*fmax((double)charge_scale,DBL_MIN) ||
       fabsl(stable_charge)>512.0L*(long double)DBL_EPSILON*fmaxl(charge_scale,DBL_MIN))return 2;
    memcpy(next,corrected,sizeof(corrected));
    return 0;
}

int snrt_chimes_reconcile_transported(const double *elements,const double *total_element_scale,
                                      const double *old,double solid_q,double *next)
{
    return reconcile_charged_impl(elements,total_element_scale,old,solid_q,next,1);
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
static int cell_charged_impl(const double *controls,const double *elements,const double *old_abundance,
                     double solid_q,
                     const double *old_photons,double *temperature,double *abundance,double *photons,
                     int mode,int allow_unilluminated_cminus_projection,
                     int original_zero_incident,double *elapsed)
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
#ifdef SNRT_CHIMES_CVODE_DIAGNOSTICS
    const char *cv_diag=getenv("SNRT_CHIMES_CVODE_DIAGNOSTICS");
    if(cv_diag && !strcmp(cv_diag,"1") && fabs(controls[0]-.11)<1e-14 &&
       controls[1]==1e4 && controls[3]>=1e10 && old_abundance[sp_H2]==.05)c.chimes_debug=1;
#endif
    /* The split dark entrypoints pass identically zero radiation. Remove
     * its 2*N_spectra ODE variables and photo-rate work, NOT the species
     * or thermal network. Keep rt_update_flux/on-the-spot unchanged: they
     * select case B even though this local solve has no photon variables.
     * Never apply this to the general (possibly illuminated) mode 0.
     * Without on-the-spot recombination the zero-radiation subspace need
     * not be invariant, so retain the original system in that case. */
    /* Retain the original hot/atomic system (mode 2): neutral hot probes
     * can trigger its strict negative-electron admission under the changed
     * CVODE norm. Do not relax that gate to obtain a performance result. */
    if((mode==1 || mode==3) && c.rt_update_flux==1 && c.rt_use_on_the_spot_approx==1){
        for(int b=0;b<NG;b++)if(old_photons[b]!=0)return 3;
        c.N_spectra=0;
    }
    if(mode==1 || mode==3){
        if(controls[1]>snrt_chimes_molecular_temperature_max())return 2;
        // Resolve the molecular photo/dark split below the default grey
        // network's thermal truncation error; do not change global tables,
        // tolerances or Tmol_K for other receivers.
        c.relativeTolerance=1e-10;c.absoluteTolerance=1e-17;c.explicitTolerance=1e-10;
    }
    if(mode>=2 && receiver_abi!=6)return 4;
    if(mode==2){
        for(int s=sp_H2;s<NS;s++)if(old_abundance[s]!=0)return 3;
        c.Tmol_K=0; /* Atomic remainder may cool below the event temperature. */
    }
    if(mode==3){c.Tmol_K=snrt_chimes_molecular_temperature_max();c.explicitTolerance=0;}
#ifdef SNRT_CHIMES_CVODE_DIAGNOSTICS
    /* Diagnostic executable only: tighten the cold-dark tolerance vector.
     * Environment is set before starting the fixed-worker parallel region. */
    const char *tol_factor=getenv("SNRT_CHIMES_TEST_TOLERANCE_FACTOR");
    if(mode==1 && tol_factor){
        char *end=NULL;
        const double factor=strtod(tol_factor,&end);
        if(end==tol_factor || *end || !isfinite(factor) || factor<=0 || factor>1)return 2;
        c.relativeTolerance*=factor;
        c.absoluteTolerance*=factor;
        c.explicitTolerance*=factor;
    }
#endif
    struct cell_context context={solid_q,mode==1?snrt_chimes_molecular_temperature_max():0,0,
        mode==3?snrt_chimes_molecular_temperature_max():0,-1};
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
        context.event_time=-1;
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
    for(int s=0;s<NS;s++)before+=work[s];
    /* In a hot-atomic solve with exactly zero incident N and E, CVODE can
     * retain a roundoff-sized C- tail after its electron abundance reaches
     * zero. Repair only this bounded C- -> C I residual. The explicit flag
     * comes from the full photo+dark split input, not the dark operator's
     * deliberately zero photon vector. */
    if(mode==0 && allow_unilluminated_cminus_projection==1 && solid_q==0 && work[sp_elec]==0){
        int zero_solver_radiation=1;
        for(int b=0;b<NG;b++)if(old_photons[b]!=0)zero_solver_radiation=0;
        if(zero_solver_radiation){
            double electron_required=solid_q;
            for(int s=1;s<NS;s++)electron_required+=charge[s]*work[s];
            const double carbon=elements[2];
            const double effective_abs_tol=c.absoluteTolerance*
                (c.scale_metal_tolerances==1?carbon:1.0);
            const double carbon_roundoff=512*DBL_EPSILON*carbon;
            const double projection_bound=fmin(carbon_roundoff,.01*effective_abs_tol);
            const double cminus=work[sp_Cm];
            if(electron_required<0 && cminus>0){
                const int eligible=isfinite(cminus) && isfinite(carbon) && carbon>=0 &&
                    isfinite(effective_abs_tol) && isfinite(projection_bound) &&
                    cminus<=projection_bound && cminus<=1e-10 && cminus>=-electron_required;
                static int reported=0;
                flockfile(stderr);
                if(!reported){
                    fprintf(stderr,"CHIMES trace C- projection incident_zero=1 eligible=%d amount=%.17g carbon=%.17g ",
                            eligible,cminus,carbon);
                    fprintf(stderr,"solver_abs_tol=%.17g effective_abs_tol=%.17g bound=%.17g ",
                            c.absoluteTolerance,effective_abs_tol,projection_bound);
                    fprintf(stderr,"electron_required_before=%.17g charge_adjustment=%.17g\n",
                            electron_required,cminus);
                    reported=1;
                }
                funlockfile(stderr);
                if(eligible){
                    work[sp_CI]+=cminus;
                    work[sp_Cm]=0;
                }
            }
        }
    }
    int reconcile_status=1;
    if(mode==3)
        reconcile_status=snrt_chimes_reconcile_dark_mode3(elements,work,solid_q,
            original_zero_incident,c.absoluteTolerance,c.scale_metal_tolerances,corrected);
    if(reconcile_status!=0)
        reconcile_status=snrt_chimes_reconcile_charged(elements,work,solid_q,corrected);
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
    for(int s=0;s<NS;s++)after+=corrected[s];
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
    if(elapsed)*elapsed=context.event_time>=0?context.event_time:controls[3];
    return mode==3 && context.event_time>=0?51:0; /* accepted thermal root */
}

static int cell_charged(const double *controls,const double *elements,const double *old,
    double solid_q,const double *old_photons,double *temperature,double *next,double *photons,
    int mode,int allow_unilluminated_cminus_projection,int original_zero_incident,double *elapsed)
{
#ifdef SNRT_CHIMES_CUDA
    snrt_chimes_rhs_cell_begin();
#endif
    int status=cell_charged_impl(controls,elements,old,solid_q,old_photons,temperature,next,photons,
                                 mode,allow_unilluminated_cminus_projection,original_zero_incident,elapsed);
#ifdef SNRT_CHIMES_CUDA
    if(snrt_chimes_rhs_cell_retry()) {
        /* Fresh private gas/rates/CVODE state and original input/tolerances.
         * Never mix CPU/GPU RHS arithmetic inside a numerical Jacobian. */
        status=cell_charged_impl(controls,elements,old,solid_q,old_photons,temperature,next,photons,
                                 mode,allow_unilluminated_cminus_projection,original_zero_incident,elapsed);
    }
    snrt_chimes_rhs_cell_end();
#endif
    return status;
}

int snrt_chimes_cell_charged(const double *controls,const double *elements,const double *old,
    double solid_q,const double *number,double *temperature,double *next,double *next_number)
{
    return cell_charged(controls,elements,old,solid_q,number,temperature,next,next_number,0,0,0,NULL);
}

int snrt_chimes_cell_charged_unilluminated(const double *controls,const double *elements,const double *old,
    double solid_q,const double *number,int zero_incident,double *temperature,double *next,double *next_number)
{
    return cell_charged(controls,elements,old,solid_q,number,temperature,next,next_number,
                        0,zero_incident==1,zero_incident==1,NULL);
}

int snrt_chimes_cell_cold_dark(const double *controls,const double *elements,const double *old,
    double *temperature,double *next)
{
    double zero[NG]={0},out[NG];
    return cell_charged(controls,elements,old,0,zero,temperature,next,out,1,0,1,NULL);
}

int snrt_chimes_radiation_is_exactly_zero(int nd,const double *number,const double *energy)
{
    if(nd<1 || nd>720 || !number || !energy)return 0;
    const size_t n=(size_t)nd*(size_t)NG;
    for(size_t i=0;i<n;i++)
        if(!isfinite(number[i]) || !isfinite(energy[i]) || number[i]!=0.0 || energy[i]!=0.0)return 0;
    return 1;
}

int snrt_chimes_cell_transition_dark(const double *controls,const double *elements,const double *old,
    int atomic,int nd,const double *original_number,const double *original_energy,double solid_q,
    double *temperature,double *next,double *elapsed)
{
    double zero[NG]={0},out[NG];
    if(elapsed)*elapsed=controls[3];
    return cell_charged(controls,elements,old,solid_q,zero,temperature,next,out,
                        atomic?2:3,0,
                        snrt_chimes_radiation_is_exactly_zero(nd,original_number,original_energy),elapsed);
}

int snrt_chimes_cell(const double *controls,const double *elements,const double *old_abundance,
                     const double *old_photons,double *temperature,double *abundance,double *photons)
{
    return snrt_chimes_cell_charged(controls,elements,old_abundance,0,
                                  old_photons,temperature,abundance,photons);
}
