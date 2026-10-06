/* Bounded semi-analytic dark CHIMES Jacobian, SUNDIALS 5.8 experiment.
 * Mass-action derivatives for atomic metals; all other columns use upstream
 * difference quotients. By default a metal column is admitted only if its DQ
 * perturbation leaves temperature exactly unchanged. The optional experimental
 * thermal-chain mode adds the ideal-gas EOS derivative and a shared numerical
 * temperature response, including energy. Fixed-T energy dependence is still
 * differenced with full cooling/heating and fresh rates, including CR heating.
 * No abundance division, abundance floors or dropped rates.
 */
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <dlfcn.h>
#include "cvode_impl.h"
#include "cvode_ls_impl.h"
#include <sunmatrix/sunmatrix_dense.h>
#include "chimes_proto.h"
#include "chimes_vars.h"
#if SUNDIALS_VERSION_MAJOR != 5 || SUNDIALS_VERSION_MINOR != 8
#error "Requires SUNDIALS 5.8"
#endif
typedef int (*jac_fn)(realtype,N_Vector,N_Vector,SUNMatrix,CVodeMem,N_Vector);
static _Thread_local jac_fn upstream;
static _Thread_local long calls,analytic_columns,numeric_columns,rhs_calls,cooling_calls,checked,fallbacks;
static _Thread_local double max_difference;
static _Thread_local long thermal_columns;
typedef void (*cool_fn)(struct gasVariables *,const struct globalVariables *,struct UserData);
static _Thread_local cool_fn upstream_cool;
static _Thread_local int cache_active,cache_valid;
static _Thread_local struct chimes_current_rates_struct *cache_owner;
static _Thread_local ChimesFloat cache_T,cache_e,cache_H,cache_Hp;
static _Thread_local long cooling_hits,cooling_misses;
extern int snrt_chimes_jac_secondary_is_dark_noop(const struct globalVariables *) __attribute__((weak));

/* Same dependency-exact coefficient reuse as the earlier DQ experiment.
 * Only interpolation is cached; the full cooling/heating sum and all hooks
 * still run for every column. Never reuse across Jacobians or buffer owners. */
void update_cooling_rates(struct gasVariables *g,const struct globalVariables *c,struct UserData d)
{
    if(!upstream_cool){upstream_cool=(cool_fn)dlsym(RTLD_NEXT,"update_cooling_rates");if(!upstream_cool)abort();}
    if(!cache_active){upstream_cool(g,c,d);return;}
    ChimesFloat e=g->abundances[c->speciesIndices[sp_elec]]*g->nH_tot;
    ChimesFloat h=g->abundances[c->speciesIndices[sp_HI]]*g->nH_tot;
    ChimesFloat hp=g->abundances[c->speciesIndices[sp_HII]]*g->nH_tot;
    if(cache_valid && cache_owner==d.chimes_current_rates && cache_T==g->temperature &&
       cache_e==e && cache_H==h && cache_Hp==hp){cooling_hits++;return;}
    upstream_cool(g,c,d);cooling_misses++;
    cache_owner=d.chimes_current_rates;cache_T=g->temperature;
    cache_e=e;cache_H=h;cache_Hp=hp;cache_valid=1;
}

static void reactions(SUNMatrix J,const int *map,const ChimesFloat *x,int count,
                      int arity,int np,const int *reactants,const int *products,
                      const ChimesFloat *coeff,ChimesFloat multiplier)
{
    for(int k=0;k<count;k++) {
        const int *r=reactants+k*arity,*p=products+k*np;
        int actual=arity;for(int a=0;a<arity;a++)if(r[a]<0){actual=a;break;}
        ChimesFloat factor=coeff[k]*multiplier;
        if(factor==0)continue;
        /* Third-body tables carry one additional nH; passed separately below. */
        for(int a=0;a<actual;a++) {
            int species=r[a],col=map[species];
            if(col<0 || species<sp_CI || species>=sp_H2)continue;
            ChimesFloat derivative=factor;
            for(int b=0;b<actual;b++)if(b!=a)derivative*=x[r[b]];
            if(derivative==0)continue;
            for(int b=0;b<actual;b++)if(map[r[b]]>=0)SM_ELEMENT_D(J,map[r[b]],col)-=derivative;
            for(int b=0;b<np && p[b]>=0;b++)if(map[p[b]]>=0)SM_ELEMENT_D(J,map[p[b]],col)+=derivative;
        }
    }
}
static void chemical_matrix(SUNMatrix J,struct UserData *d,const int *map)
{
    SUNMatZero(J);
    const struct gasVariables *g=d->myGasVars;
    const struct chimes_current_rates_struct *c=d->chimes_current_rates;
    int m=d->mol_flag_index;
    for(int k=0;k<chimes_table_T_dependent.N_reactions[m];k++) {
        const int *r=chimes_table_T_dependent.reactants+3*k;
        double density=g->nH_tot;if(r[2]>=0)density*=g->nH_tot;
        reactions(J,map,g->abundances,1,3,3,r,chimes_table_T_dependent.products+3*k,
                  c->T_dependent_rate_coefficient+k,density);
    }
#define MASS(table,a,b,coeff,mult) reactions(J,map,g->abundances,table.N_reactions[m],a,b,table.reactants,table.products,coeff,mult)
    MASS(chimes_table_constant,2,3,chimes_table_constant.rates,g->nH_tot);
    MASS(chimes_table_recombination_AB,2,1,c->recombination_AB_rate_coefficient,g->nH_tot);
    MASS(chimes_table_grain_recombination,2,1,c->grain_recombination_rate_coefficient,g->nH_tot*g->dust_ratio*g->dust_boost_factor);
    if(g->cr_rate>0)MASS(chimes_table_cosmic_ray,1,3,chimes_table_cosmic_ray.rates,g->cr_rate);
    if(m)MASS(chimes_table_H2_collis_dissoc,2,3,c->H2_collis_dissoc_rate_coefficient,g->nH_tot);
#undef MASS
    /* H2 formation-on-dust and CO CR rates have no atomic-metal dependency.
     * Electron/H/He/molecular columns, including coefficient derivatives,
     * are evaluated numerically, not approximated by this matrix. */
}
static ChimesFloat perturbed_temperature(struct UserData *d,N_Vector y)
{
    struct gasVariables *g=d->myGasVars;
    ChimesFloat T=chimes_max(NV_Ith_S(y,d->network_size)/(1.5f*
        calculate_total_number_density(g->abundances,g->nH_tot,d->myGlobalVars)*BOLTZMANNCGS),10.1f);
    if(snrt_chimes_thermal_ceiling) {
        double ceiling=snrt_chimes_thermal_ceiling(d->myGlobalVars);
        if(ceiling>0 && T>ceiling)T=ceiling;
    }
    return T;
}
static int evaluate(CVodeMem mem,realtype t,N_Vector y,N_Vector out,struct UserData *d)
{
    rhs_calls++;((CVLsMem)mem->cv_lmem)->nfeDQ++;
    return mem->cv_f(t,y,out,d);
}
int cvLsDenseDQJac(realtype t,N_Vector y,N_Vector fy,SUNMatrix J,CVodeMem mem,N_Vector tmp)
{
    if(!upstream){upstream=(jac_fn)dlsym(RTLD_NEXT,"cvLsDenseDQJac");if(!upstream)abort();}
    struct UserData *d=mem->cv_user_data;struct gasVariables *g=d->myGasVars;
    if(d->myGlobalVars->N_spectra || g->ThermEvolOn!=1 ||
       d->myGlobalVars->totalNumberOfSpecies!=CHIMES_TOTSIZE ||
       getenv("SNRT_CHIMES_RHS_BACKEND") || getenv("SNRT_CHIMES_CVODE_BACKEND") ||
       !snrt_chimes_jac_secondary_is_dark_noop || !snrt_chimes_jac_secondary_is_dark_noop(d->myGlobalVars)) {
        fallbacks++;return upstream(t,y,fy,J,mem,tmp);
    }
    for(int s=0;s<CHIMES_TOTSIZE;s++)if(d->myGlobalVars->speciesIndices[s]!=s){
        fallbacks++;return upstream(t,y,fy,J,mem,tmp);
    }
    calls++;
    cache_active=1;cache_valid=0;
    int map[CHIMES_TOTSIZE],species[CHIMES_TOTSIZE],index=0;
    for(int s=0;s<CHIMES_TOTSIZE;s++) {
        map[s]=-1;
        if(d->species[s].include_species){map[s]=index;species[index++]=s;}
    }
    if(index!=d->network_size)abort();
    int n=d->network_size+1;
    ChimesFloat original_abundances[CHIMES_TOTSIZE];
    N_Vector base=N_VClone(y),column=N_VCloneEmpty(y);
    if(!base || !column)abort();
    int rc=evaluate(mem,t,y,base,d);
    if(rc)goto cleanup;
    memcpy(original_abundances,g->abundances,sizeof(original_abundances));
    ChimesFloat T0=g->temperature;
    chemical_matrix(J,d,map);
    realtype norm=N_VWrmsNorm(fy,mem->cv_ewt);
    realtype minimum=norm!=0?1000.0*fabs(mem->cv_h)*mem->cv_uround*n*norm:1.0;
    /* One numerical temperature response shared by atomic-metal columns.
     * Keep explicit fixed-T chemistry analytic, and retain the full energy
     * response (including all cooling hooks). This is a chain-rule hybrid,
     * NOT an analytic derivative of the temperature interpolation tables. */
    realtype thermal_response[CHIMES_TOTSIZE+1]={0};
    realtype eos_slope=-T0*g->nH_tot/calculate_total_number_density(original_abundances,g->nH_tot,d->myGlobalVars);
    int thermal_ok=0;
    if(getenv("SNRT_CHIMES_JAC_THERMAL_CHAIN")) {
        int e=d->network_size;
        realtype saved=NV_Ith_S(y,e);
        realtype inc=fmax(sqrt(mem->cv_uround)*fabs(saved),minimum/NV_Ith_S(mem->cv_ewt,e));
        if(mem->cv_constraintsSet){realtype c=NV_Ith_S(mem->cv_constraints,e);
            if((fabs(c)==1 && (saved+inc)*c<0)||(fabs(c)==2 && (saved+inc)*c<=0))inc=-inc;}
        NV_Ith_S(y,e)=saved+inc;
        rc=evaluate(mem,t,y,tmp,d);
        NV_Ith_S(y,e)=saved;
        if(rc)goto cleanup;
        realtype delta=g->temperature-T0;
        if(isfinite(delta) && delta!=0) {
            thermal_ok=1;
            for(int i=0;i<n;i++)thermal_response[i]=(NV_Ith_S(tmp,i)-NV_Ith_S(base,i))/delta;
        }
    }
    for(int j=0;j<n;j++) {
        realtype saved=NV_Ith_S(y,j);
        realtype inc=fmax(sqrt(mem->cv_uround)*fabs(saved),minimum/NV_Ith_S(mem->cv_ewt,j));
        if(getenv("SNRT_CHIMES_JAC_LOCAL_CO") && j<d->network_size && species[j]==sp_CO && saved>0) {
            realtype local=fmin(inc,1e-4*saved);
            if(local>0 && saved+local!=saved)inc=local;
        }
        if(mem->cv_constraintsSet){realtype c=NV_Ith_S(mem->cv_constraints,j);
            if((fabs(c)==1 && (saved+inc)*c<0)||(fabs(c)==2 && (saved+inc)*c<=0))inc=-inc;}
        NV_Ith_S(y,j)=saved+inc;
        memcpy(g->abundances,original_abundances,sizeof(original_abundances));g->temperature=T0;
        int s=j<d->network_size?species[j]:-1;
        int analytic=s>=sp_CI && s<sp_H2 && isfinite(saved+inc) &&
            saved+inc>=-1.0e-8f && saved+inc<=CHIMES_MAX_ABUNDANCE_LIMIT;
        if(getenv("SNRT_CHIMES_JAC_NUMERICAL_ONLY"))analytic=0;
        realtype temperature_slope=0;
        if(analytic) {
            g->abundances[s]=saved+inc;
            realtype trial_temperature=perturbed_temperature(d,y);
            analytic=trial_temperature==T0 || thermal_ok;
            if(thermal_ok) {
                /* Differentiate the ideal-gas EOS itself, not the rounded
                 * change of T under a tiny abundance perturbation. On a
                 * numerical temperature clamp retain upstream DQ instead. */
                double ceiling=snrt_chimes_thermal_ceiling?snrt_chimes_thermal_ceiling(d->myGlobalVars):0;
                if(T0<=10.1f || (ceiling>0 && T0>=ceiling))analytic=trial_temperature==T0;
                else temperature_slope=eos_slope;
            }
            /* Secondary CR coefficients depend on HII. Do not assume the
             * table's special base channels can never refer to a metal. */
            for(int k=0;k<2;k++)if(s==chimes_table_cosmic_ray.reactants[chimes_table_cosmic_ray.secondary_base_reaction[k]])analytic=0;
            if(d->mol_flag_index) {
                if(s==chimes_table_H2_dust_formation.reactants[0])analytic=0;
                for(int k=0;k<chimes_table_CO_cosmic_ray.N_reactions[d->mol_flag_index];k++)
                    if(s==chimes_table_CO_cosmic_ray.reactants[k])analytic=0;
            }
        }
        if(analytic) {
            update_rate_coefficients(g,d->myGlobalVars,*d,1);
            update_rates(g,d->myGlobalVars,*d);
            ChimesFloat floor=g->TempFloor;
            if(floor<=0)floor=-floor/(1.5f*calculate_total_number_density(g->abundances,g->nH_tot,d->myGlobalVars)*BOLTZMANNCGS);
            if(T0<=floor && g->temp_floor_mode!=0)rc=-1;
            else {
                double energy=-calculate_total_cooling_rate(g,d->myGlobalVars,*d,0);cooling_calls++;
                if(T0<=floor)energy=chimes_max(energy,0.0f);
                realtype reciprocal=1.0/inc;
                SM_ELEMENT_D(J,d->network_size,j)=reciprocal*energy+(-reciprocal)*NV_Ith_S(fy,d->network_size);
                if(temperature_slope!=0) {
                    for(int i=0;i<n;i++)SM_ELEMENT_D(J,i,j)+=thermal_response[i]*temperature_slope;
                    thermal_columns++;
                }
                analytic_columns++;
            }
        } else {
            rc=evaluate(mem,t,y,tmp,d);numeric_columns++;
            if(!rc){N_VSetArrayPointer(SUNDenseMatrix_Column(J,j),column);N_VLinearSum(1.0/inc,tmp,-1.0/inc,fy,column);}
        }
        NV_Ith_S(y,j)=saved;
        if(rc)break;
    }
    if(!rc && getenv("SNRT_CHIMES_JAC_VERIFY")) {
        cache_active=0;cache_valid=0;
        SUNMatrix ref=SUNMatClone(J);if(!ref)abort();
        rc=upstream(t,y,fy,ref,mem,tmp);
        if(!rc)for(int j=0;j<n;j++)for(int i=0;i<n;i++) {
            double delta=fabs(SM_ELEMENT_D(J,i,j)-SM_ELEMENT_D(ref,i,j));
            if(!isfinite(delta))abort();if(delta>max_difference)max_difference=delta;
        }
        SUNMatDestroy(ref);checked++;
    }
cleanup:
    cache_active=0;cache_valid=0;
    N_VSetArrayPointer(NULL,column);N_VDestroy(column);N_VDestroy(base);
    if(rc){fallbacks++;return upstream(t,y,fy,J,mem,tmp);}
    return 0;
}
__attribute__((destructor)) static void report(void)
{
    if(calls || fallbacks)fprintf(stderr,"SEMIANALYTIC_JAC calls=%ld analytic_columns=%ld numeric_columns=%ld rhs_calls=%ld cooling_calls=%ld matrices_checked=%ld max_abs_J_difference=%.17g fallback=%ld\n",
        calls,analytic_columns,numeric_columns,rhs_calls,cooling_calls,checked,max_difference,fallbacks);
    if(calls)fprintf(stderr,"SEMIANALYTIC_COOLING hits=%ld misses=%ld\n",cooling_hits,cooling_misses);
    if(calls)fprintf(stderr,"SEMIANALYTIC_THERMAL columns=%ld\n",thermal_columns);
}
