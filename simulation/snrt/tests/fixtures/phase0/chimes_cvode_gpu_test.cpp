// Bounded real-CHIMES dark test + explicit LU correctness/performance probe.
// FS photochemistry is deliberately unavailable: any accidental call aborts.
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <algorithm>
#include <dlfcn.h>
#include <openssl/evp.h>
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
#ifdef SNRT_CHIMES_CPU_ONLY_TEST_HOOKS
// The isolated CPU-SPGMR profile links an older CHIMES updater that has no
// device-rate cache. Its GPU test object still references this hook, but the
// profile runs with every GPU RHS selector unset, so the hook is unreachable.
void snrt_chimes_rate_cache_invalidate(void){std::abort();}
#endif
}
static void require(bool ok,const char *msg){if(!ok){std::fprintf(stderr,"FAIL %s\n",msg);std::exit(1);}}
static void mode(const char *value){setenv("SNRT_CHIMES_CVODE_BACKEND",value,1);require(!snrt_chimes_cvode_configure(),"configure");}
static void counts(const char *label){unsigned long long c[10];snrt_chimes_cvode_counts(c);
 std::printf("COUNTS %s cpu=%llu gpu=%llu batches=%llu busy=%llu errors=%llu largest=%llu peak_concurrent=%llu streams_used=%llu gpu_solves=%llu solve_batches=%llu\n",
     label,c[0],c[1],c[2],c[3],c[4],c[5],c[6],c[7],c[8],c[9]);}

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
#ifndef SNRT_CHIMES_TEST_NCELL
#define SNRT_CHIMES_TEST_NCELL 24
#endif
constexpr int nc=SNRT_CHIMES_TEST_NCELL,ns=157;
struct Cell {double controls[9],elements[11],old[ns],next[ns],temperature;int status,tile_worker,tile_gpu_lease;};
static Cell input[nc];
static int fixture_cell_index[nc];
static const bool diagnostic_mode=std::getenv("SNRT_CHIMES_RHS_DIAGNOSTIC")!=nullptr;
static int gpu_local_rank() {
    const char *value=std::getenv("SLURM_LOCALID");
    if(!value || !*value)return 0;
    char *end=nullptr;long rank=std::strtol(value,&end,10);
    require(end && !*end && rank>=0 && rank<=INT_MAX,"Slurm local GPU rank");
    return (int)rank;
}
static void make_cells() {
    const char *composition=std::getenv("SNRT_CHIMES_TEST_COMPOSITION");
    require(!composition || !std::strcmp(composition,"mixed"),"fixture composition");
    long cell_offset=0,cell_stride=1;
    if(const char *v=std::getenv("SNRT_CHIMES_TEST_CELL_OFFSET")) {
        char *end=nullptr;cell_offset=std::strtol(v,&end,10);
        require(end && !*end && cell_offset>=0,"fixture cell offset");
    }
    if(const char *v=std::getenv("SNRT_CHIMES_TEST_CELL_STRIDE")) {
        char *end=nullptr;cell_stride=std::strtol(v,&end,10);
        require(end && !*end && cell_stride>0,"fixture cell stride");
    }
    double dt_factor=1;
    if(const char *v=std::getenv("SNRT_CHIMES_TEST_DT_FACTOR")){
        char *end=nullptr;dt_factor=std::strtod(v,&end);
        require(end && !*end && std::isfinite(dt_factor) && dt_factor>0 && dt_factor<=1,"fixture dt factor");
    }
    std::printf("FIXTURE_DT_S=%.17g\n",1e10*dt_factor);
    for(int i=0;i<nc;i++) {
        const long logical=cell_offset+i*cell_stride;
        require(logical<=INT_MAX,"fixture logical cell range");
        fixture_cell_index[i]=(int)logical;
        const int k=fixture_cell_index[i];
        auto &c=input[i];
        c.elements[0]=1;c.elements[1]=.079;c.elements[2]=1e-5;c.elements[4]=2e-5;
        if(composition) {
            // Exercise buffer reuse across full, absent and tiny nonzero
            // metal budgets. Tiny nonzero budgets must never be skipped.
            int kind=k%4;
            if(kind)for(int k=2;k<11;k++)c.elements[k]=kind==1?1e-6:(kind==2?0:1e-40);
        }
        require(!snrt_chimes_neutral(c.elements,c.old),"initial species");
        c.old[sp_HI]-=1e-4;c.old[sp_HII]+=1e-4;c.old[sp_elec]+=1e-4;
        if(k%2){c.old[sp_HI]-=.1;c.old[sp_H2]=.05;}
        c.controls[0]=.1*(1.+.1*(k%4));c.controls[1]=k%3==0?160:(k%3==1?1000:1e4);
        c.controls[2]=20;c.controls[3]=1e10*dt_factor;c.controls[4]=3e18;c.controls[5]=.01;
        c.controls[6]=1;c.controls[7]=1e-17;c.controls[8]=1;
        if(std::getenv("SNRT_CHIMES_TEST_CELL_OFFSET") || std::getenv("SNRT_CHIMES_TEST_CELL_STRIDE"))
            std::printf("FIXTURE_CELL local=%d logical=%d Tinit=%.17g nH=%.17g H2=%.17g\n",
                i,k,c.controls[1],c.controls[0],c.old[sp_H2]);
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
static void tile_cell(int i,void *context){
    auto &c=static_cast<Cell*>(context)[i];
    c.tile_worker=omp_get_thread_num();
#ifdef SNRT_CHIMES_RHS_TESTING
    c.tile_gpu_lease=snrt_chimes_rhs_tile_gpu_selected();
#endif
    c.status=snrt_chimes_cell_cold_dark(c.controls,c.elements,c.old,&c.temperature,c.next);
}
static double evolve_tile(Cell *out){
    std::memcpy(out,input,sizeof(input));double start=omp_get_wtime();
    require(!snrt_chimes_rhs_tile_run(nc,tile_cell,out),"tile dispatcher");
    return omp_get_wtime()-start;
}
static void evolve_parallel_tiles(Cell *out){
    std::memcpy(out,input,sizeof(input));
#pragma omp parallel for num_threads(8) schedule(dynamic,1)
    for(int first=0;first<nc;first+=6)
        require(!snrt_chimes_rhs_tile_run(std::min(6,nc-first),tile_cell,out+first),"parallel tile dispatcher");
}
static int compare(const Cell *ref,const Cell *trial,const char *label,bool strict_status=true) {
    double species=0,temp=0,energy=0,nuclei=0,charge=0;int failed=0,mismatch=0,cpu_rejected=0,trial_rejected=0,wi=-1,ws=-1;
    for(int i=0;i<nc;i++) {
        if(ref[i].status || trial[i].status){
            std::printf("CELL_FAILURE %s i=%d logical=%d cpu=%d trial=%d\n",label,i,fixture_cell_index[i],ref[i].status,trial[i].status);
            if(ref[i].status!=trial[i].status)mismatch++;
            if(ref[i].status)cpu_rejected++;
            if(trial[i].status)trial_rejected++;
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
    std::printf("STATUSES %s mismatch=%d cpu_rejected=%d trial_rejected=%d\n",label,mismatch,cpu_rejected,trial_rejected);
    if(wi>=0)std::printf("WORST_SPECIES %s cell=%d species=%d cpu=%.17g gpu=%.17g abs_diff=%.17g\n",label,wi,ws,ref[wi].next[ws],trial[wi].next[ws],std::abs(ref[wi].next[ws]-trial[wi].next[ws]));
    require(failed<nc || diagnostic_mode,"all fixture cells rejected");
    if(diagnostic_mode)require(nuclei<1e-8 && charge<1e-10,"dark CPU/GPU conservation diagnostic");
    else require(species<=1 && temp<1e-6 && energy<1e-6 && nuclei<1e-8 && charge<1e-10,"dark CPU/GPU physics parity");
    return mismatch;
}
static void rhs_test() {
    const char *verify_mode=std::getenv("SNRT_CHIMES_RHS_VERIFY");
    std::printf("RHS_EXECUTION_MODE=%s paired_return=%s diagnostic_only=%d\n",
        verify_mode?"paired_rhs":"device_results",verify_mode?verify_mode:"none",diagnostic_mode?1:0);
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");setenv("SNRT_CHIMES_RHS_BACKEND","cpu",1);
    require(!snrt_chimes_rhs_configure(),"CPU RHS config");
    double cpu_s=evolve(reference,8);
    double start=omp_get_wtime();cuda_pool_init(gpu_local_rank(),3);
    std::printf("RHS_POOL_INIT_S=%.9g\n",omp_get_wtime()-start);
    setenv("SNRT_CHIMES_RHS_BACKEND","cuda_dark_reactions",1);
    require(!snrt_chimes_rhs_configure(),"GPU RHS config");
#ifdef SNRT_CHIMES_RHS_TESTING
    require(snrt_chimes_rhs_test_retry_sticky(),"pending whole-cell retry remains sticky");
    puts("RHS_PENDING_RETRY_STICKY_PASS");
#endif
    double cold=evolve(trial,8);int mismatch=compare(reference,trial,"rhs-cold",false);
    unsigned long long a[7],b[7];snrt_chimes_rhs_counts(a);
    double warm=evolve(trial,8);mismatch+=compare(reference,trial,"rhs-warm",false);
    snrt_chimes_rhs_counts(b);
    std::printf("RHS_CPU_S=%.9g RHS_GPU_COLD_S=%.9g RHS_GPU_WARM_S=%.9g SPEEDUP=%.9g GPU_CALLS=%llu CPU_FALLBACK=%llu BATCHES=%llu ERRORS=%llu TABLE_UPLOADS=%llu WARM_GPU_CALLS=%llu WARM_UPLOADS=%llu PEAK_CONCURRENT_BATCHES=%llu STREAMS_USED=%llu\n",
        cpu_s,cold,warm,cpu_s/warm,b[0],b[1],b[2],b[3],b[4],b[0]-a[0],b[4]-a[4],b[5],b[6]);
    require(b[0]>a[0] && b[3]==0 && b[4]>=a[4],"real resident reaction callbacks");
    require(b[5]>1 && b[6]>1,"multiple RHS streams execute concurrently");
    if(!diagnostic_mode){
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
    }
    double paired[9];snrt_chimes_rhs_verification(paired);
    std::printf("RHS_PAIRED_MAX rates_rel=%.17g gross_species_rel=%.17g cooling_rel=%.17g net_species_rel=%.17g CO_net_abs=%.17g CO_net_rel=%.17g nonfinite_gpu=%llu nonfinite_cpu=%llu invalid_gpu_outputs=%llu\n",
        paired[0],paired[1],paired[2],paired[3],paired[4],paired[5],
        (unsigned long long)paired[6],(unsigned long long)paired[7],(unsigned long long)paired[8]);
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    if(!diagnostic_mode)require(!mismatch,"RHS solver status parity");
    puts(diagnostic_mode?"CHIMES_DEVICE_DARK_REACTIONS_DIAGNOSTIC_COMPLETE":"CHIMES_DEVICE_DARK_REACTIONS_PARITY_PASS");
}
static void tile_test(){
    require(!diagnostic_mode,"tile correctness cannot bypass acceptance");
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");require(!snrt_chimes_rhs_configure_mode(0),"tile CPU configure");
    const double cpu=evolve(reference,8);
    cuda_pool_init(gpu_local_rank(),3);
    require(!snrt_chimes_rhs_configure_mode(1),"tile GPU configure");
    unsigned long long before[7],after[7];snrt_chimes_rhs_counts(before);
    const double gpu=evolve_tile(trial);compare(reference,trial,"tile-staged");
    snrt_chimes_rhs_counts(after);
    require(after[0]>before[0] && after[3]==before[3],"tile real GPU work without device errors");
    int held[3];for(int i=0;i<3;i++){held[i]=cuda_acquire_stream();require(held[i]>=0,"tile busy lease");}
    evolve_tile(trial);compare(reference,trial,"tile-whole-CPU-fallback");
    for(int i=0;i<3;i++)cuda_release_stream(held[i]);
    unsigned long long busy[7];snrt_chimes_rhs_counts(busy);
    require(busy[0]==after[0],"busy tile never dispatches GPU callbacks");
    require(!snrt_chimes_rhs_tile_run(0,tile_cell,trial),"empty tile");
    require(snrt_chimes_rhs_tile_run(129,tile_cell,trial)!=0,"bounded tile capacity");
    // Match production's OMP cohort ownership: one available lease and
    // independent CPU cohorts; active GPU cohorts keep host LU on CPU.
    mode("cuda_batched_lu");
    unsigned long long lu_before[10],lu_after[10];snrt_chimes_cvode_counts(lu_before);
    for(int i=0;i<2;i++){held[i]=cuda_acquire_stream();require(held[i]>=0,"parallel tile held lease");}
    evolve_parallel_tiles(trial);compare(reference,trial,"tile-parallel-hybrid");
    for(int i=0;i<2;i++)cuda_release_stream(held[i]);
    snrt_chimes_cvode_counts(lu_after);
    require(lu_after[1]==lu_before[1],"tile-owned lease cannot be reacquired by host LU");
#ifdef SNRT_CHIMES_RHS_TESTING
    const auto drops_before=snrt_chimes_rhs_test_drops();
    snrt_chimes_rhs_test_drop_after(3);
    evolve_tile(trial);compare(reference,trial,"tile-midcell-retry");
    snrt_chimes_rhs_test_drop_after(0);
    require(snrt_chimes_rhs_test_drops()>drops_before,"tile forced midcell retry executed");
#endif
    // Same queue in CPU-only and broker arms; validate ownership, completion,
    // busy-pool CPU fallback and fresh-state retry without changing tolerances.
    setenv("SNRT_CHIMES_CELL_SCHEDULER","level_queue",1);
    setenv("SNRT_CHIMES_GPU_WORKER","thread0",1);
    require(!snrt_chimes_rhs_configure_mode(0),"level CPU configure");
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"level CPU queue");
    compare(reference,trial,"level-native-CPU");
    require(!snrt_chimes_rhs_configure_mode(1),"level broker configure");
    snrt_chimes_rhs_counts(before);
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"level broker queue");
    compare(reference,trial,"level-broker");
    snrt_chimes_rhs_counts(after);
    require(after[0]>before[0] && after[3]==before[3],"level real GPU work without errors");
    for(int i=0;i<3;i++){held[i]=cuda_acquire_stream();require(held[i]>=0,"level busy lease");}
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"level busy CPU queue");
    compare(reference,trial,"level-busy-CPU");
    for(int i=0;i<3;i++)cuda_release_stream(held[i]);
    snrt_chimes_rhs_counts(busy);
    require(busy[0]==after[0],"level busy queue does no GPU work");
#ifdef SNRT_CHIMES_RHS_TESTING
    const auto level_drops=snrt_chimes_rhs_test_drops();
    snrt_chimes_rhs_test_drop_after(3);
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"level retry queue");
    compare(reference,trial,"level-midcell-retry");
    snrt_chimes_rhs_test_drop_after(0);
    require(snrt_chimes_rhs_test_drops()>level_drops,"level forced midcell retry");
#endif
    // Four independent host brokers, each owning a distinct stream/buffers.
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    cuda_pool_init(gpu_local_rank(),4);
    setenv("SNRT_CHIMES_GPU_WORKER","any",1);
    setenv("SNRT_CHIMES_GPU_BROKERS","4",1);
    require(!snrt_chimes_rhs_configure_mode(1),"four broker configure");
    snrt_chimes_rhs_counts(before);
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"four broker queue");
    compare(reference,trial,"level-four-brokers");
    snrt_chimes_rhs_counts(after);
    require(after[0]>before[0] && after[3]==before[3],"four brokers real GPU work without errors");
    require((snrt_chimes_rhs_worker_mask()&15)==15,"four GPU workers executed");
    unsigned long long transfer[3];snrt_chimes_rhs_result_counts(transfer);
    require(transfer[0]>0,"compact result batches actually executed");
    puts("CHIMES_LEVEL_FOUR_BROKERS_PARITY_PASS");
    setenv("SNRT_CHIMES_GPU_RESULT","full",1);
    require(!snrt_chimes_rhs_configure_mode(1),"full result configure");
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"full result level queue");
    compare(reference,trial,"level-full-result-reference");
    unsetenv("SNRT_CHIMES_GPU_RESULT");
    require(!snrt_chimes_rhs_configure_mode(1),"compact result restore");
    puts("CHIMES_LEVEL_COMPACT_FULL_PARITY_PASS");
#ifdef SNRT_CHIMES_RHS_TESTING
    // Fail a submission after earlier batches have recorded reusable events.
    // Its stale event must not cause consumption of an unrecorded candidate.
    snrt_chimes_rhs_test_submit_fail_after(2);
    snrt_chimes_rhs_counts(before);
    std::memcpy(trial,input,sizeof(input));
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"level submission failure queue");
    compare(reference,trial,"level-submission-failure-CPU-retry");
    snrt_chimes_rhs_counts(after);
    require(after[3]>before[3],"submission failure path executed");
    snrt_chimes_rhs_test_submit_fail_after(0);
    puts("CHIMES_LEVEL_SUBMISSION_FAILURE_PARITY_PASS");
#endif
    unsetenv("SNRT_CHIMES_GPU_BROKERS");
    unsetenv("SNRT_CHIMES_CELL_SCHEDULER");unsetenv("SNRT_CHIMES_GPU_WORKER");
    puts("CHIMES_LEVEL_BROKER_PARITY_PASS");
    std::printf("TILE_STAGED_TIME cpu_omp8_s=%.9g tile_gpu_s=%.9g gpu_calls=%llu batches=%llu errors=%llu\n",
        cpu,gpu,after[0]-before[0],after[2]-before[2],after[3]-before[3]);
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_TILE_STAGED_PARITY_PASS");
}
static void tile_rhs_diagnostic(){
    require(diagnostic_mode,"tile RHS diagnostic requires diagnostic mode");
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");
    setenv("SNRT_CHIMES_RHS_BACKEND","cpu",1);
    require(!snrt_chimes_rhs_configure_mode(0),"tile RHS diagnostic CPU configure");
    const double cpu=evolve(reference,8);
    cuda_pool_init(gpu_local_rank(),3);
    setenv("SNRT_CHIMES_RHS_BACKEND","cuda_dark_reactions",1);
    require(!snrt_chimes_rhs_configure_mode(1),"tile RHS diagnostic GPU configure");
    const double gpu=evolve_tile(trial);
    const int mismatch=compare(reference,trial,"tile-rhs-paired",false);
    int cpu_rejected=0,gpu_rejected=0;
    for(int i=0;i<nc;i++){cpu_rejected+=reference[i].status!=0;gpu_rejected+=trial[i].status!=0;}
    double paired[9];snrt_chimes_rhs_verification(paired);
    unsigned long long work[7];snrt_chimes_rhs_counts(work);
    std::printf("TILE_RHS_DIAGNOSTIC cells=%d cpu_s=%.9g gpu_tile_s=%.9g cpu_rejected=%d gpu_rejected=%d status_mismatch=%d gpu_calls=%llu cpu_fallback=%llu batches=%llu errors=%llu\n",
        nc,cpu,gpu,cpu_rejected,gpu_rejected,mismatch,work[0],work[1],work[2],work[3]);
    std::printf("TILE_RHS_PAIRED_MAX rates_rel=%.17g gross_species_rel=%.17g cooling_rel=%.17g net_species_rel=%.17g CO_net_abs=%.17g CO_net_rel=%.17g nonfinite_gpu=%.0f nonfinite_cpu=%.0f invalid_gpu_outputs=%.0f\n",
        paired[0],paired[1],paired[2],paired[3],paired[4],paired[5],paired[6],paired[7],paired[8]);
    require(work[0]>0,"tile diagnostic executed GPU RHS");
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_TILE_PAIRED_RHS_DIAGNOSTIC_COMPLETE");
}
static void tile_route_diagnostic(){
    require(diagnostic_mode,"tile route diagnostic requires diagnostic mode");
    make_cells();Cell reference[nc],trial[nc];
    mode("cpu");
    require(!snrt_chimes_rhs_configure_mode(0),"tile route CPU configure");
    const double cpu=evolve(reference,8);
    cuda_pool_init(gpu_local_rank(),3);
    setenv("SNRT_CHIMES_GPU_WORKER","any",1);
    require(!snrt_chimes_rhs_configure_mode(1),"tile route GPU configure");
    unsigned long long before[7],after[7];snrt_chimes_rhs_counts(before);
    const double start=omp_get_wtime();
    evolve_parallel_tiles(trial);
    const double hybrid=omp_get_wtime()-start;
    const int mismatch=compare(reference,trial,"tile-route-status",false);
    snrt_chimes_rhs_counts(after);
    int thread_tiles[8]={0},thread_gpu_tiles[8]={0},thread_mismatch[8]={0};
    int gpu_tiles=0,cpu_tiles=0,mismatch_tiles=0,mixed_tiles=0;
    for(int first=0,tile=0;first<nc;first+=6,tile++){
        const int last=std::min(nc,first+6),worker=trial[first].tile_worker;
        const int gpu_lease=trial[first].tile_gpu_lease;
        int tile_mismatch=0,tile_bad_route=0;
        for(int i=first;i<last;i++){
            if(trial[i].tile_worker!=worker || trial[i].tile_gpu_lease!=gpu_lease)tile_bad_route++;
            if(reference[i].status!=trial[i].status){
                tile_mismatch++;
                std::printf("TILE_ROUTE_MISMATCH tile=%d local=%d logical=%d worker=%d gpu_lease=%d cpu_status=%d trial_status=%d\n",
                    tile,i,fixture_cell_index[i],trial[i].tile_worker,trial[i].tile_gpu_lease,
                    reference[i].status,trial[i].status);
            }
        }
        if(worker>=0 && worker<8){thread_tiles[worker]++;thread_gpu_tiles[worker]+=gpu_lease;thread_mismatch[worker]+=tile_mismatch;}
        gpu_tiles+=gpu_lease;cpu_tiles+=!gpu_lease;mismatch_tiles+=tile_mismatch>0;mixed_tiles+=tile_bad_route>0;
    }
    std::printf("TILE_ROUTE_SUMMARY cells=%d tiles=%d gpu_lease_tiles=%d cpu_fallback_tiles=%d tiles_with_mismatch=%d mixed_route_tiles=%d status_mismatch_cells=%d gpu_rhs_calls=%llu cpu_fallback_calls=%llu gpu_batches=%llu gpu_errors=%llu worker_mask=0x%llx cpu_s=%.9g mixed_s=%.9g\n",
        nc,(nc+5)/6,gpu_tiles,cpu_tiles,mismatch_tiles,mixed_tiles,mismatch,
        after[0]-before[0],after[1]-before[1],after[2]-before[2],after[3]-before[3],
        snrt_chimes_rhs_worker_mask(),cpu,hybrid);
    for(int thread=0;thread<8;thread++)if(thread_tiles[thread])
        std::printf("TILE_THREAD_ROUTE thread=%d tiles=%d gpu_lease_tiles=%d cpu_fallback_tiles=%d mismatch_cells=%d\n",
            thread,thread_tiles[thread],thread_gpu_tiles[thread],thread_tiles[thread]-thread_gpu_tiles[thread],thread_mismatch[thread]);
    require(after[0]>before[0],"tile route diagnostic exercised GPU RHS");
    snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_TILE_ROUTE_DIAGNOSTIC_COMPLETE");
}
#ifndef SNRT_CHIMES_CPU_BUILD_ONLY
static void bundle_test(bool rhs,bool benchmark,bool integrated=false) {
    require(!diagnostic_mode,"qualification cannot use diagnostic bypass");
    make_cells();Cell reference[nc],trial[nc];
    auto select=[&](bool gpu){
        if(integrated) {
            require(!snrt_chimes_cvode_configure_mode(gpu?1:0),"integrated CVODE configure");
            require(!snrt_chimes_rhs_configure_mode(gpu?1:0),"integrated RHS configure");
        } else {
            mode(gpu && !rhs?"cuda_batched_lu":"cpu");
            setenv("SNRT_CHIMES_RHS_BACKEND",gpu && rhs?"cuda_dark_reactions":"cpu",1);
            require(!snrt_chimes_rhs_configure(),"bundle RHS configure");
        }
    };
    auto accepted=[](const Cell *cells){for(int i=0;i<nc;i++)require(cells[i].status==0,"bundle requires every cell accepted");};
    double start=omp_get_wtime();cuda_pool_init(gpu_local_rank(),3);
    std::printf("BUNDLE_POOL_INIT_S=%.9g\n",omp_get_wtime()-start);
    double cold=0;
    if(integrated) {
        const char *selector=std::getenv("SNRT_CHIMES_COMPUTE_BACKEND");
        require(selector && !std::strcmp(selector,"cuda_integrated"),"integrated selector initialized the full GPU path");
        unsigned long long lu_before[10],lu_after[10],rhs_before[7],rhs_after[7];
        snrt_chimes_cvode_counts(lu_before);snrt_chimes_rhs_counts(rhs_before);
        cold=evolve(trial,8);accepted(trial);
        snrt_chimes_cvode_counts(lu_after);snrt_chimes_rhs_counts(rhs_after);
        require(rhs_after[0]>rhs_before[0],"integrated solve exercised GPU RHS");
        require(lu_after[4]==lu_before[4] && rhs_after[3]==rhs_before[3],"integrated cold solve has no device errors");
        select(false);evolve(reference,8);accepted(reference);
        compare(reference,trial,"integrated-cold");
        std::printf("INTEGRATED_COLD_S=%.9g gpu_lu=%llu gpu_rhs=%llu\n",cold,
            lu_after[1]-lu_before[1],rhs_after[0]-rhs_before[0]);
    } else {
        select(false);evolve(reference,8);accepted(reference);
        select(true);cold=evolve(trial,8);accepted(trial);compare(reference,trial,"bundle-cold");
        std::printf("BUNDLE_COLD_S=%.9g\n",cold);
    }
    unsigned long long lu0[10],lu1[10],r0[7],r1[7];
    snrt_chimes_cvode_counts(lu0);snrt_chimes_rhs_counts(r0);
    int samples=benchmark?3:1;
    for(int sample=0;sample<samples;sample++)for(int k=0;k<2;k++) {
        bool gpu=(k==(sample%2==0?1:0));select(gpu);
        bool capture=gpu && sample==0 && std::getenv("SNRT_CHIMES_BUNDLE_PROFILE");
        auto profiler=[&](const char *symbol){auto fn=(int(*)())dlsym(RTLD_DEFAULT,symbol);require(fn && fn()==0,"CUDA profiler range");};
        if(capture)profiler("cudaProfilerStart");
        double seconds=evolve(trial,8);
        if(capture)profiler("cudaProfilerStop");
        accepted(trial);compare(reference,trial,gpu?"bundle-warm-gpu":"bundle-warm-cpu");
        std::printf("BUNDLE_TIME backend=%s arm=%s sample=%d cells=%d seconds=%.9g\n",
            integrated?"integrated":(rhs?"rhs":"lu"),gpu?"gpu":"cpu",sample+1,nc,seconds);
    }
    snrt_chimes_cvode_counts(lu1);snrt_chimes_rhs_counts(r1);
    std::printf("BUNDLE_GPU_WORK lu_factors=%llu lu_batches=%llu lu_errors=%llu lu_peak_concurrent=%llu lu_streams_used=%llu rhs_calls=%llu rhs_batches=%llu rhs_errors=%llu warm_table_uploads=%llu rhs_peak_concurrent=%llu rhs_streams_used=%llu\n",
        lu1[1]-lu0[1],lu1[2]-lu0[2],lu1[4]-lu0[4],lu1[6],lu1[7],r1[0]-r0[0],r1[2]-r0[2],r1[3]-r0[3],r1[4]-r0[4],r1[5],r1[6]);
    if(integrated)require(r1[0]>r0[0],"integrated warm GPU RHS exercised");
    else require(rhs?r1[0]>r0[0]:lu1[1]>lu0[1],"GPU work exercised");
    if(rhs||integrated)require(r1[5]>1 && r1[6]>1,"multiple RHS streams execute concurrently");
    if(!rhs && !integrated)require(lu1[6]>1 && lu1[7]>1,"multiple LU streams execute concurrently");
    require(lu1[4]==lu0[4] && r1[3]==r0[3],"warm device errors");
    snrt_chimes_rhs_finalize();snrt_chimes_cvode_cuda_finalize();cuda_pool_finalize();
    puts("CHIMES_GPU_BUNDLE_PASS");
}
#endif
static void cpu_build_test(bool benchmark) {
    make_cells();Cell reference[nc],trial[nc];
    evolve(reference,1);
    if(const char *path=std::getenv("SNRT_CHIMES_TEST_REFERENCE")) {
        FILE *file=std::fopen(path,"r");require(file,"open reference state");
        int cells=0,species=0;require(std::fscanf(file,"%d%d",&cells,&species)==2 && cells==nc && species==ns,"reference dimensions");
        Cell saved[nc];std::memcpy(saved,input,sizeof(input));
        for(auto &c:saved) {
            require(std::fscanf(file,"%d%lf",&c.status,&c.temperature)==2,"reference cell");
            for(double &x:c.next)require(std::fscanf(file,"%lf",&x)==1,"reference abundance");
        }
        require(!std::fclose(file),"reference close");compare(saved,reference,"cpu-build-library");
    }
    if(const char *path=std::getenv("SNRT_CHIMES_TEST_STATE_FILE")) {
        FILE *file=std::fopen(path,"wx");require(file,"create new state file");
        std::fprintf(file,"%d %d\n",nc,ns);
        for(auto &c:reference) {
            std::fprintf(file,"%d %.17g",c.status,c.temperature);
            for(double x:c.next)std::fprintf(file," %.17g",x);
            std::fprintf(file,"\n");
        }
        require(!std::fclose(file),"state close");
    }
    for(int threads:{1,8}) {
        evolve(trial,threads);compare(reference,trial,"cpu-build-threads");
        EVP_MD_CTX *hash=EVP_MD_CTX_new();require(hash,"state hash allocation");
        require(EVP_DigestInit_ex(hash,EVP_sha256(),nullptr)==1,"state hash init");
        int rejected=0;
        for(auto &c:trial) {
            if(c.status)rejected++;
            require(EVP_DigestUpdate(hash,&c.status,sizeof(c.status))==1 &&
                    EVP_DigestUpdate(hash,&c.temperature,sizeof(c.temperature))==1 &&
                    EVP_DigestUpdate(hash,c.next,sizeof(c.next))==1,"state hash update");
        }
        unsigned char bytes[EVP_MAX_MD_SIZE];unsigned int nbytes=0;
        require(EVP_DigestFinal_ex(hash,bytes,&nbytes)==1,"state hash finish");EVP_MD_CTX_free(hash);
        std::printf("CPU_STATE threads=%d rejected=%d sha256=",threads,rejected);
        for(unsigned int i=0;i<nbytes;i++)std::printf("%02x",bytes[i]);std::puts("");
        if(benchmark) {
            constexpr int repeats=20;double start=omp_get_wtime();
            for(int rep=0;rep<repeats;rep++)evolve(trial,threads);
            double wall=omp_get_wtime()-start;
            compare(reference,trial,"cpu-build-repeat");
            std::printf("CPU_BUILD_TIME threads=%d repeats=%d wall_s=%.9g\n",threads,repeats,wall);
        }
    }
}
static void cuda_nvector_trial() {
    require(!diagnostic_mode,"CUDA N_Vector trial requires production tolerances");
    make_cells();Cell reference[nc],trial[nc];
    std::fprintf(stderr,"NVECTOR_TRIAL_PHASE cpu-baseline-start\n");
    setenv("SNRT_CHIMES_CVODE_BACKEND","cpu",1);
    setenv("SNRT_CHIMES_CVODE_NVECTOR","cpu",1);
    setenv("SNRT_CHIMES_RHS_BACKEND","cpu",1);
    mode("cpu");
    require(!snrt_chimes_rhs_configure_mode(0),"CPU RHS configure");
    const double cpu_s=evolve(reference,8);
    std::fprintf(stderr,"NVECTOR_TRIAL_PHASE cpu-baseline-complete\n");

    cuda_pool_init(gpu_local_rank(),3);
    setenv("SNRT_CHIMES_CVODE_NVECTOR","cuda_managed",1);
    mode("cpu"); // disable batched LU; test CUDA N_Vector with CPU RHS/control
    require(!snrt_chimes_rhs_configure_mode(0),"CPU RHS configure for managed vector");
    std::fprintf(stderr,"NVECTOR_TRIAL_PHASE cuda-vector-start\n");
    const double cuda_s=evolve(trial,8);
    std::fprintf(stderr,"NVECTOR_TRIAL_PHASE cuda-vector-complete\n");
    const int status_mismatch=compare(reference,trial,"cuda-managed-nvector",false);
    require(status_mismatch==0,"matching rejection status");
    std::printf("NVECTOR_CPU_S=%.9g NVECTOR_CUDA_MANAGED_S=%.9g SPEEDUP=%.9g\n",
        cpu_s,cuda_s,cpu_s/cuda_s);
    unsetenv("SNRT_CHIMES_CVODE_NVECTOR");
    snrt_chimes_cvode_cuda_finalize();snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_CUDA_MANAGED_NVECTOR_PARITY_PASS");
}
static void cuda_dense_solve_trial() {
    require(!diagnostic_mode,"CUDA dense-solve trial requires production tolerances");
    make_cells();Cell reference[nc],trial[nc];
    setenv("SNRT_CHIMES_CVODE_NVECTOR","cpu",1);
    mode("cpu");
    require(!snrt_chimes_rhs_configure_mode(0),"CPU RHS configure");
    const double cpu_s=evolve(reference,8);
    cuda_pool_init(gpu_local_rank(),3);
    mode("cuda_batched_lu_solve");
    require(!snrt_chimes_rhs_configure_mode(0),"CPU RHS configure for GPU dense solve");
    unsigned long long before[10],after[10];snrt_chimes_cvode_counts(before);
    const double gpu_s=evolve(trial,8);
    snrt_chimes_cvode_counts(after);
    const int status_mismatch=compare(reference,trial,"cuda-batched-dense-solve",false);
    std::printf("DENSE_SOLVE_CPU_S=%.9g DENSE_SOLVE_GPU_S=%.9g gpu_factors=%llu gpu_solves=%llu solve_batches=%llu\n",
        cpu_s,gpu_s,after[1]-before[1],after[8]-before[8],after[9]-before[9]);
    require(status_mismatch==0,"matching rejection status");
    require(after[1]>before[1] && after[8]>before[8] && after[9]>before[9],"GPU factor and solve work executed");
    require(after[4]==before[4],"no GPU dense-solve device errors");

    // Exercise both RAMSES dispatch contracts with the experimental solver
    // enabled: fiber/tile cohorts retain host dense solves, and native CPU
    // level-queue workers must never enter the CUDA factor/solve queues.
    setenv("SNRT_CHIMES_GPU_WORKER","any",1);
    setenv("SNRT_CHIMES_GPU_BROKERS","3",1);
    require(!snrt_chimes_rhs_configure_mode(1),"tile RHS configure for solver isolation");
    unsigned long long rhs_before[7],rhs_after[7];
    snrt_chimes_rhs_counts(rhs_before);
    snrt_chimes_cvode_counts(before);
    evolve_parallel_tiles(trial);
    snrt_chimes_cvode_counts(after);
    snrt_chimes_rhs_counts(rhs_after);
    require(rhs_after[0]>rhs_before[0],"tile dispatch exercised real GPU RHS work");
    require(after[0]==before[0] && after[1]==before[1] && after[8]==before[8] && after[9]==before[9],
        "tile fibers retain host dense linear algebra");
    puts("CHIMES_CVODE_TILE_SOLVER_ROUTE_PASS");

    require(!snrt_chimes_rhs_configure_mode(0),"disable RHS GPU for CPU-direct lane test");
    setenv("SNRT_CHIMES_CELL_SCHEDULER","level_queue",1);
    std::memcpy(trial,input,sizeof(input));
    snrt_chimes_cvode_counts(before);
    require(!snrt_chimes_rhs_level_run(nc,tile_cell,trial),"native CPU level queue with experimental solver configured");
    compare(reference,trial,"solver-cpu-direct-lane");
    snrt_chimes_cvode_counts(after);
    require(after[0]==before[0] && after[1]==before[1] && after[8]==before[8] && after[9]==before[9],
        "CPU-direct workers never enter CUDA dense factor/solve queues");
    puts("CHIMES_CVODE_DISPATCH_ISOLATION_PASS");

    snrt_chimes_cvode_cuda_finalize();snrt_chimes_rhs_finalize();cuda_pool_finalize();
    puts("CHIMES_CUDA_BATCHED_DENSE_SOLVE_PARITY_PASS");
}
int main(int argc,char **argv) {
    require(argc==2 || argc==3,"main table argument");omp_set_dynamic(0);
    double start=omp_get_wtime();
    int init_status=snrt_chimes_initialize(argv[1],0,nullptr,0);
    if(argc==3 && (!std::strcmp(argv[2],"cpu-bench") || !std::strcmp(argv[2],"cpu-state"))) {
        require(!init_status,"CPU build initialization");cpu_build_test(!std::strcmp(argv[2],"cpu-bench"));return 0;
    }
#ifndef SNRT_CHIMES_CPU_BUILD_ONLY
    if(argc==3 && (!std::strcmp(argv[2],"check-lu") || !std::strcmp(argv[2],"check-rhs") ||
                   !std::strcmp(argv[2],"bench-lu") || !std::strcmp(argv[2],"bench-rhs") ||
                   !std::strcmp(argv[2],"check-integrated") || !std::strcmp(argv[2],"bench-integrated"))) {
        require(!init_status,"bundle initialization");
        bool integrated=std::strstr(argv[2],"integrated")!=nullptr;
        bundle_test(std::strstr(argv[2],"rhs")!=nullptr,std::strncmp(argv[2],"bench",5)==0,integrated);return 0;
    }
    if(argc==3 && !std::strcmp(argv[2],"cuda-nvector")){require(!init_status,"CUDA N_Vector initialization");cuda_nvector_trial();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"cuda-direct-solve")){require(!init_status,"CUDA direct-solve initialization");cuda_dense_solve_trial();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"rhs")){require(!init_status,"rhs initialization");rhs_test();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"tile")){require(!init_status,"tile initialization");tile_test();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"tile-rhs-diagnostic")){require(!init_status,"tile RHS diagnostic initialization");tile_rhs_diagnostic();return 0;}
    if(argc==3 && !std::strcmp(argv[2],"tile-route-diagnostic")){require(!init_status,"tile route diagnostic initialization");tile_route_diagnostic();return 0;}
#endif
    if(argc==3 && !std::strcmp(argv[2],"profile")) {
        require(!init_status,"profile initialization");make_cells();Cell reference[nc];
        double wall=evolve(reference,1);int failed=0;
        for(auto &c:reference)if(c.status)failed++;
        std::printf("SERIAL_RHS_PROFILE cells=%d rejected=%d wall_s=%.9g\n",nc,failed,wall);
        for(int i:{5,17})std::printf("PROFILE_CELL i=%d status=%d T=%.17g H2=%.17g CO=%.17g\n",
            i,reference[i].status,reference[i].temperature,reference[i].next[sp_H2],reference[i].next[sp_CO]);
        return 0;
    }
#ifdef SNRT_CHIMES_CPU_BUILD_ONLY
    return 2;
#else
    if(argc==3) {std::printf("BRIDGE_CONFIG_STATUS=%d\n",init_status);return init_status;}
    require(!init_status,"bridge initialize");
    std::printf("TABLE_INIT_S=%.9g\n",omp_get_wtime()-start);
    start=omp_get_wtime();cuda_pool_init(gpu_local_rank(),3);
    std::printf("CUDA_POOL_INIT_S=%.9g\n",omp_get_wtime()-start);
    mode("cpu");double cpu_lu=lu_test(8,100,false);
    mode("cuda_batched_lu");
    double cold_lu=lu_test(8,1,false);double gpu_lu=lu_test(8,100,false);
    lu_test(8,1,true);lu_test(1,2,false);counts("lu");
    std::printf("LU_CPU_GPU_SPEEDUP=%.9g LU_COLD_S=%.9g\n",cpu_lu/gpu_lu,cold_lu);
    unsigned long long before[10],after[10];snrt_chimes_cvode_counts(before);
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
#endif
}
