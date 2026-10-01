// Bounded real-CHIMES dark test + explicit LU correctness/performance probe.
// FS photochemistry is deliberately unavailable: any accidental call aborts.
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "snrt_chimes_cvode_cuda.h"
#include "snrt_chimes_rhs_cuda.h"
extern "C" {
#include "chimes_proto.h"
#include "chimes_vars.h"
SUNLinearSolver snrt_chimes_cvode_linear_solver(N_Vector,SUNMatrix);
void cuda_pool_init(int,int);
int cuda_acquire_stream(void);
void cuda_release_stream(int);
void cuda_pool_finalize(void);
// This executable never initializes MG/scalar/particle device storage.
// The real stream pool's finalizer still references these cleanup hooks.
void cuda_mg_finalize(void){}
void cuda_scal_finalize(void){}
void cuda_pm_finalize(void){}
int snrt_chimes_initialize(const char*,int,const char*,int);
int snrt_chimes_neutral(const double*,double*);
int snrt_chimes_budget(const double*,double*,double*);
int snrt_chimes_cell_cold_dark(const double*,const double*,const double*,double*,double*);
int snrt_chimes_fs_ready(void){return 0;}
int snrt_chimes_fs_fractions(double,double,double*){std::abort();}
int snrt_chimes_fs_samples(double,double*,double*){std::abort();}
}
static void require(bool ok,const char *msg){if(!ok){std::fprintf(stderr,"FAIL %s\n",msg);std::exit(1);}}
static void mode(const char *value){setenv("SNRT_CHIMES_CVODE_BACKEND",value,1);require(!snrt_chimes_cvode_configure(),"configure");}
static void counts(const char *label){unsigned long long c[6];snrt_chimes_cvode_counts(c);
 std::printf("COUNTS %s cpu=%llu gpu=%llu batches=%llu busy=%llu errors=%llu largest=%llu\n",label,c[0],c[1],c[2],c[3],c[4],c[5]);}

static double lu_test(int threads,int repeats,bool singular) {
    int bad=0;double maxerr=0,start=omp_get_wtime();
#pragma omp parallel num_threads(threads) reduction(+:bad) reduction(max:maxerr)
    {
        const int n=158,tid=omp_get_thread_num();
        N_Vector y=N_VNew_Serial(n),b=N_VNew_Serial(n);
        SUNMatrix a=SUNDenseMatrix(n,n);
        SUNLinearSolver s=snrt_chimes_cvode_linear_solver(y,a);
        if(!s)std::abort();
        SUNLinSolInitialize(s);
        for(int rep=0;rep<repeats;rep++) {
            // Permuted, diagonally dominant matrix requires partial pivoting.
            for(int j=0;j<n;j++)for(int i=0;i<n;i++)
                SM_ELEMENT_D(a,(i+13)%n,j)=i==j ? n+1.+tid : .01*std::sin(i+3*j+tid);
            for(int i=0;i<n;i++) {
                double sum=0;for(int j=0;j<n;j++)sum+=SM_ELEMENT_D(a,i,j)*(1.+j/(double)n);
                NV_Ith_S(b,i)=sum;
            }
            if(singular && tid==0)for(int j=0;j<n;j++)SM_ELEMENT_D(a,0,j)=0;
#pragma omp barrier
            int status=SUNLinSolSetup(s,a);
            if(singular && tid==0){if(status!=SUNLS_LUFACT_FAIL)bad++;continue;}
            if(status || SUNLinSolSolve(s,a,y,b,0)){bad++;continue;}
            for(int i=0;i<n;i++)maxerr=std::max(maxerr,std::abs(NV_Ith_S(y,i)-(1.+i/(double)n)));
        }
        SUNLinSolFree(s);SUNMatDestroy(a);N_VDestroy(y);N_VDestroy(b);
    }
    require(!bad && maxerr<1e-11,"LU pivot/status parity");
    double wall=omp_get_wtime()-start;
    std::printf("LU threads=%d repeats=%d singular=%d max_abs_error=%.17g wall_s=%.9g\n",threads,repeats,singular,maxerr,wall);
    return wall;
}
constexpr int nc=24,ns=157;
struct Cell {double controls[9],elements[11],old[ns],next[ns],temperature;int status;};
static Cell input[nc];
static void make_cells() {
    double dt_factor=1;
    if(const char *v=std::getenv("SNRT_CHIMES_TEST_DT_FACTOR")){
        char *end=nullptr;dt_factor=std::strtod(v,&end);
        require(end && !*end && std::isfinite(dt_factor) && dt_factor>0 && dt_factor<=1,"fixture dt factor");
    }
    std::printf("FIXTURE_DT_S=%.17g\n",1e10*dt_factor);
    for(int i=0;i<nc;i++) {
        auto &c=input[i];
        c.elements[0]=1;c.elements[1]=.079;c.elements[2]=1e-5;c.elements[4]=2e-5;
        require(!snrt_chimes_neutral(c.elements,c.old),"initial species");
        c.old[sp_HI]-=1e-4;c.old[sp_HII]+=1e-4;c.old[sp_elec]+=1e-4;
        if(i%2){c.old[sp_HI]-=.1;c.old[sp_H2]=.05;}
        c.controls[0]=.1*(1.+.1*(i%4));c.controls[1]=i%3==0?160:(i%3==1?1000:1e4);
        c.controls[2]=20;c.controls[3]=1e10*dt_factor;c.controls[4]=3e18;c.controls[5]=.01;
        c.controls[6]=1;c.controls[7]=1e-17;c.controls[8]=1;
    }
}
static double evolve(Cell *out,int threads) {
    std::memcpy(out,input,sizeof(input));double t=omp_get_wtime();
#pragma omp parallel for num_threads(threads) schedule(static)
    for(int i=0;i<nc;i++) {
        auto &c=out[i];
        c.status=snrt_chimes_cell_cold_dark(c.controls,c.elements,c.old,&c.temperature,c.next);
    }
    return omp_get_wtime()-t;
}
static int compare(const Cell *ref,const Cell *trial,const char *label,bool strict_status=true) {
    double species=0,temp=0,energy=0,nuclei=0,charge=0;int failed=0,mismatch=0,wi=-1,ws=-1;
    for(int i=0;i<nc;i++) {
        if(ref[i].status || trial[i].status){
            std::printf("CELL_FAILURE %s i=%d cpu=%d trial=%d\n",label,i,ref[i].status,trial[i].status);
            if(ref[i].status!=trial[i].status)mismatch++;
            if(strict_status)require(ref[i].status==trial[i].status,"matching rejection status");
            if(ref[i].status)require(!std::memcmp(ref[i].next,input[i].old,sizeof(input[i].old)),"CPU rejected state unchanged");
            if(trial[i].status)require(!std::memcmp(trial[i].next,input[i].old,sizeof(input[i].old)),"GPU rejected state unchanged");
            failed++;continue;
        }
        double el[11],q=0,a=0,b=0;
        require(!snrt_chimes_budget(trial[i].next,el,&q),"positive finite species");
        charge=std::max(charge,std::abs(q));
        for(int j=0;j<11;j++)nuclei=std::max(nuclei,std::abs(el[j]-input[i].elements[j])/std::max(input[i].elements[j],1e-20));
        for(int j=0;j<ns;j++) {
            double discrepancy=std::abs(ref[i].next[j]-trial[i].next[j])/(1e-17+1e-6*std::abs(ref[i].next[j]));
            if(discrepancy>species){species=discrepancy;wi=i;ws=j;}
            a+=ref[i].next[j];b+=trial[i].next[j];
        }
        temp=std::max(temp,std::abs(trial[i].temperature/ref[i].temperature-1));
        energy=std::max(energy,std::abs(b*trial[i].temperature/(a*ref[i].temperature)-1));
    }
    std::printf("PHYSICS %s failed=%d species_scaled_error=%.17g temperature_rel=%.17g energy_rel=%.17g nuclei_rel=%.17g charge_abs=%.17g\n",label,failed,species,temp,energy,nuclei,charge);
    std::printf("STATUSES %s mismatch=%d\n",label,mismatch);
    if(wi>=0)std::printf("WORST_SPECIES %s cell=%d species=%d cpu=%.17g gpu=%.17g abs_diff=%.17g\n",label,wi,ws,ref[wi].next[ws],trial[wi].next[ws],std::abs(ref[wi].next[ws]-trial[wi].next[ws]));
    require(failed<nc && species<=1 && temp<1e-6 && energy<1e-6 && nuclei<1e-8 && charge<1e-10,"dark CPU/GPU physics parity");
    return mismatch;
}
static void rhs_test() {
    std::printf("RHS_EXECUTION_MODE=%s\n",std::getenv("SNRT_CHIMES_RHS_VERIFY")?"paired_cpu_reference":"device_results");
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");setenv("SNRT_CHIMES_RHS_BACKEND","cpu",1);
    require(!snrt_chimes_rhs_configure(),"CPU RHS config");
    double cpu_s=evolve(reference,8);
    double start=omp_get_wtime();cuda_pool_init(0,3);
    std::printf("RHS_POOL_INIT_S=%.9g\n",omp_get_wtime()-start);
    setenv("SNRT_CHIMES_RHS_BACKEND","cuda_dark_reactions",1);
    require(!snrt_chimes_rhs_configure(),"GPU RHS config");
    double cold=evolve(trial,8);int mismatch=compare(reference,trial,"rhs-cold",false);
    unsigned long long a[5],b[5];snrt_chimes_rhs_counts(a);
    double warm=evolve(trial,8);mismatch+=compare(reference,trial,"rhs-warm",false);
    snrt_chimes_rhs_counts(b);
    std::printf("RHS_CPU_S=%.9g RHS_GPU_COLD_S=%.9g RHS_GPU_WARM_S=%.9g SPEEDUP=%.9g GPU_CALLS=%llu CPU_FALLBACK=%llu BATCHES=%llu ERRORS=%llu TABLE_UPLOADS=%llu WARM_GPU_CALLS=%llu WARM_UPLOADS=%llu\n",
        cpu_s,cold,warm,cpu_s/warm,b[0],b[1],b[2],b[3],b[4],b[0]-a[0],b[4]-a[4]);
    require(b[0]>a[0] && b[3]==0 && b[4]==a[4],"real resident reaction callbacks");
    int held[3];for(int i=0;i<3;i++){held[i]=cuda_acquire_stream();require(held[i]>=0,"busy RHS lease");}
    evolve(trial,8);compare(reference,trial,"rhs-busy-fallback");
    for(int i=0;i<3;i++)cuda_release_stream(held[i]);
    snrt_chimes_rhs_counts(a);require(a[0]==b[0] && a[1]>b[1],"busy uses CPU RHS");
    evolve(trial,1);compare(reference,trial,"rhs-serial-fallback");
#ifdef SNRT_CHIMES_RHS_TESTING
    snrt_chimes_rhs_test_drop_after(3);
    evolve(trial,8);compare(reference,trial,"rhs-midcell-retry");
    unsigned long long drops=snrt_chimes_rhs_test_drops();
    std::printf("RHS_MIDCELL_RETRIES=%llu\n",drops);require(drops>0,"midcell retry exercised");
    snrt_chimes_rhs_test_drop_after(0);
#endif
    double paired[3];snrt_chimes_rhs_verification(paired);
    std::printf("RHS_PAIRED_MAX rates_rel=%.17g species_flux_rel=%.17g net_cooling_rel=%.17g\n",paired[0],paired[1],paired[2]);
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    require(!mismatch,"RHS solver status parity");
    puts("CHIMES_DEVICE_DARK_REACTIONS_PARITY_PASS");
}
int main(int argc,char **argv) {
    require(argc==2 || argc==3,"main table argument");omp_set_dynamic(0);
    double start=omp_get_wtime();
    int init_status=snrt_chimes_initialize(argv[1],0,nullptr,0);
    if(argc==3 && !std::strcmp(argv[2],"rhs")){require(!init_status,"rhs initialization");rhs_test();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"profile")) {
        require(!init_status,"profile initialization");make_cells();Cell reference[nc];
        double wall=evolve(reference,1);int failed=0;
        for(auto &c:reference)if(c.status)failed++;
        std::printf("SERIAL_RHS_PROFILE cells=%d rejected=%d wall_s=%.9g\n",nc,failed,wall);return 0;
    }
    if(argc==3) {std::printf("BRIDGE_CONFIG_STATUS=%d\n",init_status);return init_status;}
    require(!init_status,"bridge initialize");
    std::printf("TABLE_INIT_S=%.9g\n",omp_get_wtime()-start);
    start=omp_get_wtime();cuda_pool_init(0,3);
    std::printf("CUDA_POOL_INIT_S=%.9g\n",omp_get_wtime()-start);
    mode("cpu");double cpu_lu=lu_test(8,100,false);
    mode("cuda_batched_lu");
    double cold_lu=lu_test(8,1,false);double gpu_lu=lu_test(8,100,false);
    lu_test(8,1,true);lu_test(1,2,false);counts("lu");
    std::printf("LU_CPU_GPU_SPEEDUP=%.9g LU_COLD_S=%.9g\n",cpu_lu/gpu_lu,cold_lu);
    unsigned long long before[6],after[6];snrt_chimes_cvode_counts(before);
    require(before[1]>0,"GPU LU actually used");
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");double cpu_s=evolve(reference,8);
    mode("cuda_batched_lu");double gpu_s=evolve(trial,8);
    snrt_chimes_cvode_counts(after);counts("physics");
    std::printf("CHIMES_CPU_S=%.9g GPU_S=%.9g SPEEDUP=%.9g GPU_FACTORS=%llu\n",cpu_s,gpu_s,cpu_s/gpu_s,after[1]-before[1]);
    compare(reference,trial,"gpu-omp8");
    int held[3];for(int i=0;i<3;i++){held[i]=cuda_acquire_stream();require(held[i]>=0,"reserve test lease");}
    double busy_s=evolve(trial,8);compare(reference,trial,"busy-fallback");
    for(int i=0;i<3;i++)cuda_release_stream(held[i]);
    counts("busy");std::printf("BUSY_FALLBACK_S=%.9g\n",busy_s);
    Cell invalid=input[0];invalid.controls[0]=-1;
    require(snrt_chimes_cell_cold_dark(invalid.controls,invalid.elements,invalid.old,&invalid.temperature,invalid.next)!=0,"reject invalid input");
    require(!std::memcmp(invalid.next,invalid.old,sizeof(invalid.old)),"failure leaves input state");
    snrt_chimes_cvode_cuda_finalize();cuda_pool_finalize();
    require(after[1]>before[1],"real CHIMES callback batching exercised");
    puts("CHIMES_CVODE_BATCH_LU_DARK_PASS");
}
