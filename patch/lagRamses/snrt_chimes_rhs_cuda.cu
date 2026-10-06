/* SPDX-License-Identifier: LGPL-3.0-or-later
 * CHIMES reaction equations: Copyright (c) 2020 Alexander Richings.
 * CHIMES dark reaction stage, FP64. Equations/operation order follow the
 * LGPL3+ CHIMES update_rate_coefficients/update_rates/update_rate_vector
 * (Alexander Richings 2020), inherited source a58e5c0 plus local patches.
 * Tables/topology stay on device; each independent CVODE callback retains
 * ownership of its abundance vector, rates and error control on the host.
 */
#include "snrt_chimes_rhs_cuda.h"
#include "../cuRamses/cuda_stream_pool.h"
extern "C" {
#include "chimes_proto.h"
#include "chimes_vars.h"
void snrt_chimes_rate_cache_invalidate(void);
}
#include <omp.h>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <mutex>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstdint>
#include <cstddef>
#include <thread>
#include <cuda_profiler_api.h>
#include <nvtx3/nvToolsExt.h>
static_assert(sizeof(ChimesFloat)==8 && sizeof(sunindextype)==8,
              "Device reaction slice requires FP64 / 64-bit SUNDIALS indices");
/* SysV AMD64: only callee-saved registers, plus the FP control words.
 * The saved stack pointer addresses the return slot (RSP % 16 == 8).
 * A fresh fiber stores its trampoline at a 16-byte-aligned slot so the
 * first ret enters that function with the ABI's entry alignment. */
struct alignas(64) FastCtx {
  std::uint64_t sp,rbx,rbp,r12,r13,r14,r15;
  std::uint32_t mxcsr;
  std::uint16_t fcw,pad;
};
static_assert(sizeof(FastCtx)==64,"fast context");
static_assert(offsetof(FastCtx,mxcsr)==56 && offsetof(FastCtx,fcw)==60,"fast context layout");
extern "C" __attribute__((naked)) void snrt_chimes_fast_swap(FastCtx *,FastCtx *){
  __asm__(
    "movq %rsp, (%rdi)\n\t"
    "movq %rbx, 8(%rdi)\n\t"
    "movq %rbp, 16(%rdi)\n\t"
    "movq %r12, 24(%rdi)\n\t"
    "movq %r13, 32(%rdi)\n\t"
    "movq %r14, 40(%rdi)\n\t"
    "movq %r15, 48(%rdi)\n\t"
    "stmxcsr 56(%rdi)\n\t"
    "fnstcw 60(%rdi)\n\t"
    "movq (%rsi), %rsp\n\t"
    "movq 8(%rsi), %rbx\n\t"
    "movq 16(%rsi), %rbp\n\t"
    "movq 24(%rsi), %r12\n\t"
    "movq 32(%rsi), %r13\n\t"
    "movq 40(%rsi), %r14\n\t"
    "movq 48(%rsi), %r15\n\t"
    "ldmxcsr 56(%rsi)\n\t"
    "fldcw 60(%rsi)\n\t"
    "ret\n\t");
}
namespace {
void capture_fp(std::uint32_t &mxcsr,std::uint16_t &fcw){
  asm volatile("stmxcsr %0" : "=m"(mxcsr));
  asm volatile("fnstcw %0" : "=m"(fcw));
}
constexpr int ns=157, nb_max=128, nrmax=1024, groups=8;
int nb=64;
enum {TD,CON,REC,GRAIN,CR,H2D,H2C,CO};
enum {AXT,AXD,AXX,TDT,CONT,RECT,GRAINT,CRT,CRSEC,H2DT,KZERO,KLTE,CRITH,CRITH2,CRITHE,COT,NTAB};
struct Reaction {int kind,local,a[3];};
struct Tables {
 int nt,nd,nx,np,nr,count[groups][2],base[groups],offset[NTAB],secondary[2];
 int cooling_index[4];bool compact_safe=false;
 int begin[2*ns+1]; // CSR: destruction first, creation second; duplicate nuclei retained
};
struct Input {
 double x[ns],t,td,nh,dust,boost,cr;
 double cell_size,div_vel,extinction,constant_heating,cmb,tmol;
 double co_column,h2o_column,oh_column,elem[10];
 int mol,ab[2],static_mol,hybrid_mode,elem_mask,rewrite_columns;
};
struct Output {
 double rate[nrmax],coefficient[nrmax],crit[3],destroy[ns],create[ns];
 double net_cooling;int cooling_valid;
};
// The admitted CHIMES 2.4 main table has 136 molecular cosmic-ray rates.
// compact_safe keeps larger tables on the full-result path rather than
// truncating this payload or overrunning the fixed compact record.
constexpr int compact_cr_max=136;
// Keep full reaction scratch on device. Host cooling in the admitted dark
// network reads only these five molecular values and cosmic-ray rates.
struct CompactOutput {
 double destroy[ns],create[ns],crit[3],molecular[5],cosmic[compact_cr_max];
 double net_cooling;
 int valid,cooling_valid;
};
static_assert(sizeof(CompactOutput)<sizeof(Output)/5,"Compact result must reduce D2H payload");
struct HostTables {Tables t{};std::vector<double> v;std::vector<Reaction> r;std::vector<int> entries;bool valid=false;};
HostTables host;
std::once_flag tables_once;
bool tables_ready=false;
bool enabled=false;
bool gpu_thread0_only=false;
bool level_queue_enabled=false;
bool compact_results=true;
bool gpu_cooling=true;
int level_gpu_brokers=1;
bool brokers_explicit=false;
thread_local bool cpu_direct=false;
thread_local double gpu_net_cooling=0;
thread_local int gpu_net_cooling_ready=0;
std::atomic<bool> timeline_claimed{false};
bool verify=false;
bool verify_return_gpu=false;
thread_local int cell_backend=-2; // -2 outside bridge, -1 undecided, 0 CPU, 1 GPU
thread_local bool retry_cell=false;
#ifdef SNRT_CHIMES_RHS_TESTING
int drop_after=0;thread_local int cell_calls=0;
int submit_fail_after=0;
std::atomic<int> submit_attempts{0};
std::atomic<unsigned long long> drops{0};
#endif
std::mutex verification_mutex;
double verify_rates=0,verify_species=0,verify_cooling=0;
double verify_net=0,verify_co_abs=0,verify_co_rel=0;
unsigned long long verify_nonfinite_device=0,verify_nonfinite_cpu=0;
std::atomic<unsigned long long> used{0},fallback{0},batches{0},errors{0},uploads{0};
std::atomic<unsigned long long> invalid_device_outputs{0};
std::atomic<unsigned long long> active_batches{0},peak_active_batches{0},used_stream_mask{0};
std::atomic<unsigned long long> used_worker_mask{0};
std::atomic<unsigned long long> compact_batches{0},full_batches{0},result_bytes{0};

struct CoolTables {
 int enabled=0,n1=0,n2=0,n4=0;
 int nT2=0,nNe2=0,nThi2=0,nT4=0,nHI4=0,nNe4=0,nHII4=0,nThi4=0;
 int nTmol=0,nCOrot=0,nCOvib=0,nH2Ohi=0,nH2Olo=0,nH2Orot=0,nH2Ovib=0;
 int have_co=0,have_h2o=0,base1=0,base2=0,base4=0;
 int o_rates=0,o_rates2=0,o_rates4=0,o_hi2=0,o_hi4=0,o_gg=0;
 int o_T2=0,o_ne2=0,o_Thi2=0,o_T4=0,o_nHI=0,o_ne4=0,o_nHII=0,o_Thi4=0;
 int o_tmol=0,o_h2h2=0,o_h2hi=0,o_h2hii=0,o_h2hei=0,o_h2e=0,o_h2lte=0;
 int o_corN=0,o_covN=0,o_corL0=0,o_corLlte=0,o_corNh=0,o_corA=0,o_covL0=0,o_covLlte=0;
 int o_hiT=0,o_loT=0,o_rotN=0,o_vibN=0;
 int o_hiL0=0,o_hiLlte=0,o_hiNh=0,o_hiA=0;
 int o_orL0=0,o_orLlte=0,o_orNh=0,o_orA=0;
 int o_paL0=0,o_paLlte=0,o_paNh=0,o_paA=0;
 int o_vibL0=0,o_vibLlte=0;
};
CoolTables cool;
bool put_doubles(std::vector<double> &dst,const double *p,size_t n,int &off){
 if(n && !p)return false;
 if(dst.size()+n>(size_t)INT_MAX)return false;
 off=(int)dst.size();if(n)dst.insert(dst.end(),p,p+n);return true;
}
bool put_species(std::vector<int> &dst,const int *s,int n,int &base){
 if(n<0 || (n && !s))return false;
 for(int i=0;i<n;i++)if(s[i]<0 || s[i]>=ns)return false;
 base=(int)dst.size();if(n)dst.insert(dst.end(),s,s+n);return true;
}
void append_cooling(HostTables &h){
 cool={};
 if(!gpu_cooling)return;
 auto &bins=chimes_table_bins;auto &tab=chimes_table_cooling;
 const int nT=bins.N_Temperatures,n1=tab.N_coolants,n2=tab.N_coolants_2d,n4=tab.N_coolants_4d;
 if(nT<2 || n1<0 || n2<0 || n4<0 || (n1 && !tab.rates) || !tab.gas_grain_transfer)return;
 if(bins.N_mol_cool_Temperatures<2 || !bins.mol_cool_Temperatures)return;
 if(!tab.H2_cool_lowDens_H2 || !tab.H2_cool_lowDens_HI || !tab.H2_cool_lowDens_HII ||
    !tab.H2_cool_lowDens_HeI || !tab.H2_cool_lowDens_elec || !tab.H2_cool_LTE)return;
 CoolTables next{};
 next.n1=n1;next.n2=n2;next.n4=n4;next.nTmol=bins.N_mol_cool_Temperatures;
 if(!put_species(h.entries,tab.coolants,n1,next.base1))return;
 if(!put_species(h.entries,tab.coolants_2d,n2,next.base2))return;
 if(!put_species(h.entries,tab.coolants_4d,n4,next.base4))return;
 if(!put_doubles(h.v,tab.rates,(size_t)n1*nT,next.o_rates))return;
 if(!put_doubles(h.v,tab.gas_grain_transfer,(size_t)nT,next.o_gg))return;
 if(!put_doubles(h.v,bins.mol_cool_Temperatures,(size_t)next.nTmol,next.o_tmol))return;
 const size_t nm=(size_t)next.nTmol;
 if(!put_doubles(h.v,tab.H2_cool_lowDens_H2,nm,next.o_h2h2))return;
 if(!put_doubles(h.v,tab.H2_cool_lowDens_HI,nm,next.o_h2hi))return;
 if(!put_doubles(h.v,tab.H2_cool_lowDens_HII,nm,next.o_h2hii))return;
 if(!put_doubles(h.v,tab.H2_cool_lowDens_HeI,nm,next.o_h2hei))return;
 if(!put_doubles(h.v,tab.H2_cool_lowDens_elec,nm,next.o_h2e))return;
 if(!put_doubles(h.v,tab.H2_cool_LTE,nm,next.o_h2lte))return;
 if(n2){
   if(bins.N_cool_2d_Temperatures<2 || bins.N_cool_2d_ElectronDensities<2 || bins.N_cool_hiT_2d_Temperatures<2 ||
      !bins.cool_2d_Temperatures || !bins.cool_2d_ElectronDensities || !bins.cool_hiT_2d_Temperatures ||
      !tab.rates_2d || !tab.rates_hiT_2d)return;
   next.nT2=bins.N_cool_2d_Temperatures;next.nNe2=bins.N_cool_2d_ElectronDensities;next.nThi2=bins.N_cool_hiT_2d_Temperatures;
   if(!put_doubles(h.v,bins.cool_2d_Temperatures,(size_t)next.nT2,next.o_T2))return;
   if(!put_doubles(h.v,bins.cool_2d_ElectronDensities,(size_t)next.nNe2,next.o_ne2))return;
   if(!put_doubles(h.v,bins.cool_hiT_2d_Temperatures,(size_t)next.nThi2,next.o_Thi2))return;
   if(!put_doubles(h.v,tab.rates_2d,(size_t)n2*next.nT2*next.nNe2,next.o_rates2))return;
   if(!put_doubles(h.v,tab.rates_hiT_2d,(size_t)n2*next.nThi2,next.o_hi2))return;
 }
 if(n4){
   if(bins.N_cool_4d_Temperatures<2 || bins.N_cool_4d_HIDensities<2 || bins.N_cool_4d_ElectronDensities<2 ||
      bins.N_cool_4d_HIIDensities<2 || bins.N_cool_hiT_4d_Temperatures<2 || !bins.cool_4d_Temperatures ||
      !bins.cool_4d_HIDensities || !bins.cool_4d_ElectronDensities || !bins.cool_4d_HIIDensities ||
      !bins.cool_hiT_4d_Temperatures || !tab.rates_4d || !tab.rates_hiT_4d)return;
   next.nT4=bins.N_cool_4d_Temperatures;next.nHI4=bins.N_cool_4d_HIDensities;
   next.nNe4=bins.N_cool_4d_ElectronDensities;next.nHII4=bins.N_cool_4d_HIIDensities;
   next.nThi4=bins.N_cool_hiT_4d_Temperatures;
   if(!put_doubles(h.v,bins.cool_4d_Temperatures,(size_t)next.nT4,next.o_T4))return;
   if(!put_doubles(h.v,bins.cool_4d_HIDensities,(size_t)next.nHI4,next.o_nHI))return;
   if(!put_doubles(h.v,bins.cool_4d_ElectronDensities,(size_t)next.nNe4,next.o_ne4))return;
   if(!put_doubles(h.v,bins.cool_4d_HIIDensities,(size_t)next.nHII4,next.o_nHII))return;
   if(!put_doubles(h.v,bins.cool_hiT_4d_Temperatures,(size_t)next.nThi4,next.o_Thi4))return;
   if(!put_doubles(h.v,tab.rates_4d,(size_t)n4*next.nT4*next.nHI4*next.nNe4*next.nHII4,next.o_rates4))return;
   if(!put_doubles(h.v,tab.rates_hiT_4d,(size_t)n4*next.nThi4,next.o_hi4))return;
 }
 if(tab.CO_cool_rot_L0 && tab.CO_cool_rot_Llte && tab.CO_cool_rot_nhalf && tab.CO_cool_rot_a &&
    tab.CO_cool_vib_L0 && tab.CO_cool_vib_Llte && bins.CO_cool_rot_ColumnDensities && bins.CO_cool_vib_ColumnDensities &&
    bins.N_CO_cool_rot_ColumnDensities>=2 && bins.N_CO_cool_vib_ColumnDensities>=2){
   next.have_co=1;next.nCOrot=bins.N_CO_cool_rot_ColumnDensities;next.nCOvib=bins.N_CO_cool_vib_ColumnDensities;
   if(!put_doubles(h.v,bins.CO_cool_rot_ColumnDensities,(size_t)next.nCOrot,next.o_corN))return;
   if(!put_doubles(h.v,bins.CO_cool_vib_ColumnDensities,(size_t)next.nCOvib,next.o_covN))return;
   if(!put_doubles(h.v,tab.CO_cool_rot_L0,nm,next.o_corL0))return;
   if(!put_doubles(h.v,tab.CO_cool_rot_Llte,nm*next.nCOrot,next.o_corLlte))return;
   if(!put_doubles(h.v,tab.CO_cool_rot_nhalf,nm*next.nCOrot,next.o_corNh))return;
   if(!put_doubles(h.v,tab.CO_cool_rot_a,nm*next.nCOrot,next.o_corA))return;
   if(!put_doubles(h.v,tab.CO_cool_vib_L0,nm,next.o_covL0))return;
   if(!put_doubles(h.v,tab.CO_cool_vib_Llte,nm*next.nCOvib,next.o_covLlte))return;
 }
 if(tab.H2O_cool_rot_hiT_L0 && tab.H2O_cool_rot_hiT_Llte && tab.H2O_cool_rot_hiT_nhalf && tab.H2O_cool_rot_hiT_a &&
    tab.H2Oortho_cool_rot_lowT_L0 && tab.H2Oortho_cool_rot_lowT_Llte && tab.H2Oortho_cool_rot_lowT_nhalf && tab.H2Oortho_cool_rot_lowT_a &&
    tab.H2Opara_cool_rot_lowT_L0 && tab.H2Opara_cool_rot_lowT_Llte && tab.H2Opara_cool_rot_lowT_nhalf && tab.H2Opara_cool_rot_lowT_a &&
    tab.H2O_cool_vib_L0 && tab.H2O_cool_vib_Llte && bins.H2O_cool_hiT_Temperatures && bins.H2O_cool_lowT_Temperatures &&
    bins.H2O_cool_rot_ColumnDensities && bins.H2O_cool_vib_ColumnDensities &&
    bins.N_H2O_cool_hiT_Temperatures>=2 && bins.N_H2O_cool_lowT_Temperatures>=2 &&
    bins.N_H2O_cool_rot_ColumnDensities>=2 && bins.N_H2O_cool_vib_ColumnDensities>=2){
   next.have_h2o=1;
   next.nH2Ohi=bins.N_H2O_cool_hiT_Temperatures;next.nH2Olo=bins.N_H2O_cool_lowT_Temperatures;
   next.nH2Orot=bins.N_H2O_cool_rot_ColumnDensities;next.nH2Ovib=bins.N_H2O_cool_vib_ColumnDensities;
   if(!put_doubles(h.v,bins.H2O_cool_hiT_Temperatures,(size_t)next.nH2Ohi,next.o_hiT))return;
   if(!put_doubles(h.v,bins.H2O_cool_lowT_Temperatures,(size_t)next.nH2Olo,next.o_loT))return;
   if(!put_doubles(h.v,bins.H2O_cool_rot_ColumnDensities,(size_t)next.nH2Orot,next.o_rotN))return;
   if(!put_doubles(h.v,bins.H2O_cool_vib_ColumnDensities,(size_t)next.nH2Ovib,next.o_vibN))return;
   if(!put_doubles(h.v,tab.H2O_cool_rot_hiT_L0,(size_t)next.nH2Ohi,next.o_hiL0))return;
   if(!put_doubles(h.v,tab.H2O_cool_rot_hiT_Llte,(size_t)next.nH2Ohi*next.nH2Orot,next.o_hiLlte))return;
   if(!put_doubles(h.v,tab.H2O_cool_rot_hiT_nhalf,(size_t)next.nH2Ohi*next.nH2Orot,next.o_hiNh))return;
   if(!put_doubles(h.v,tab.H2O_cool_rot_hiT_a,(size_t)next.nH2Ohi*next.nH2Orot,next.o_hiA))return;
   if(!put_doubles(h.v,tab.H2Oortho_cool_rot_lowT_L0,(size_t)next.nH2Olo,next.o_orL0))return;
   if(!put_doubles(h.v,tab.H2Oortho_cool_rot_lowT_Llte,(size_t)next.nH2Olo*next.nH2Orot,next.o_orLlte))return;
   if(!put_doubles(h.v,tab.H2Oortho_cool_rot_lowT_nhalf,(size_t)next.nH2Olo*next.nH2Orot,next.o_orNh))return;
   if(!put_doubles(h.v,tab.H2Oortho_cool_rot_lowT_a,(size_t)next.nH2Olo*next.nH2Orot,next.o_orA))return;
   if(!put_doubles(h.v,tab.H2Opara_cool_rot_lowT_L0,(size_t)next.nH2Olo,next.o_paL0))return;
   if(!put_doubles(h.v,tab.H2Opara_cool_rot_lowT_Llte,(size_t)next.nH2Olo*next.nH2Orot,next.o_paLlte))return;
   if(!put_doubles(h.v,tab.H2Opara_cool_rot_lowT_nhalf,(size_t)next.nH2Olo*next.nH2Orot,next.o_paNh))return;
   if(!put_doubles(h.v,tab.H2Opara_cool_rot_lowT_a,(size_t)next.nH2Olo*next.nH2Orot,next.o_paA))return;
   if(!put_doubles(h.v,tab.H2O_cool_vib_L0,nm,next.o_vibL0))return;
   if(!put_doubles(h.v,tab.H2O_cool_vib_Llte,nm*next.nH2Ovib,next.o_vibLlte))return;
 }
 next.enabled=1;cool=next;
}

bool build_tables_impl() {
 if(host.valid)return true;
 host=HostTables{};auto &h=host;auto &t=h.t;
 t.nt=chimes_table_bins.N_Temperatures;t.nd=chimes_table_bins.N_Dust_Temperatures;
 t.nx=chimes_table_bins.N_secondary_cosmic_ray_xHII;t.np=chimes_table_bins.N_Psi;
 if(t.nt<2||t.nd<2||t.nx<2||t.np<1)return false;
 auto add=[&](int key,const double *p,size_t n){t.offset[key]=(int)h.v.size();if(n && !p)return false;if(n)h.v.insert(h.v.end(),p,p+n);return true;};
 std::vector<int> dst[2*ns];
 auto group=[&](int kind,int na,int nm,const int *a,int ar,const int *p,int pr){
   if(na<0||nm<na||h.r.size()+nm>nrmax|| (nm && (!a||!p)))return false;
   t.count[kind][0]=na;t.count[kind][1]=nm;t.base[kind]=(int)h.r.size();
   for(int i=0;i<nm;i++){
     Reaction r{kind,i,{-1,-1,-1}};int idx=(int)h.r.size();
     for(int j=0;j<ar;j++){int s=a[i*ar+j];if(s<0)break;if(s>=ns)return false;r.a[j]=s;dst[s].push_back(idx);}
     for(int j=0;j<pr;j++){int s=p[i*pr+j];if(s<0)break;if(s>=ns)return false;dst[ns+s].push_back(idx);}
     h.r.push_back(r);
   }return true;
 };
#define GROUP(k,table,ar,pr) if(!group(k,table.N_reactions[0],table.N_reactions[1],table.reactants,ar,table.products,pr))return false
 GROUP(TD,chimes_table_T_dependent,3,3);GROUP(CON,chimes_table_constant,2,3);
 GROUP(REC,chimes_table_recombination_AB,2,1);GROUP(GRAIN,chimes_table_grain_recombination,2,1);
 GROUP(CR,chimes_table_cosmic_ray,1,3);
 if(!group(H2D,0,1,chimes_table_H2_dust_formation.reactants,2,chimes_table_H2_dust_formation.products,1))return false;
 GROUP(H2C,chimes_table_H2_collis_dissoc,2,3);GROUP(CO,chimes_table_CO_cosmic_ray,1,2);
#undef GROUP
 if(t.count[REC][1]!=2 || !chimes_table_cosmic_ray.secondary_base_reaction)return false;
 for(int j=0;j<2;j++){t.secondary[j]=chimes_table_cosmic_ray.secondary_base_reaction[j];if(t.secondary[j]<0||t.secondary[j]>=t.count[CR][1])return false;}
 t.nr=(int)h.r.size();
 t.cooling_index[0]=chimes_table_T_dependent.H2_collis_dissoc_heating_reaction_index;
 t.cooling_index[1]=chimes_table_H2_collis_dissoc.Heating_reaction_index;
 t.cooling_index[2]=chimes_table_T_dependent.H2_form_heating_reaction_index;
 t.cooling_index[3]=chimes_table_constant.H2_form_heating_reaction_index;
 t.compact_safe=t.count[CR][1]<=compact_cr_max;
 const int cooling_group[4]={TD,H2C,TD,CON};
 for(int j=0;j<4;j++)t.compact_safe=t.compact_safe &&
   t.cooling_index[j]>=0 && t.cooling_index[j]<t.count[cooling_group[j]][1];
 for(int i=0;i<2*ns;i++){t.begin[i]=(int)h.entries.size();h.entries.insert(h.entries.end(),dst[i].begin(),dst[i].end());}
 t.begin[2*ns]=(int)h.entries.size();
#define ADD(k,p,n) if(!add(k,p,(size_t)(n)))return false
 ADD(AXT,chimes_table_bins.Temperatures,t.nt);ADD(AXD,chimes_table_bins.Dust_Temperatures,t.nd);
 ADD(AXX,chimes_table_bins.secondary_cosmic_ray_xHII,t.nx);
 ADD(TDT,chimes_table_T_dependent.rates,t.count[TD][1]*t.nt);
 ADD(CONT,chimes_table_constant.rates,t.count[CON][1]);
 ADD(RECT,chimes_table_recombination_AB.rates,t.count[REC][1]*2*t.nt);
 ADD(GRAINT,chimes_table_grain_recombination.rates,t.count[GRAIN][1]*t.nt*t.np);
 ADD(CRT,chimes_table_cosmic_ray.rates,t.count[CR][1]);ADD(CRSEC,chimes_table_cosmic_ray.secondary_ratio,2*t.nx);
 ADD(H2DT,chimes_table_H2_dust_formation.rates,t.nt*t.nd);
 ADD(KZERO,chimes_table_H2_collis_dissoc.k0,t.count[H2C][1]*t.nt);
 ADD(KLTE,chimes_table_H2_collis_dissoc.kLTE,t.count[H2C][1]*t.nt);
 ADD(CRITH,chimes_table_H2_collis_dissoc.critical_density_H,t.nt);
 ADD(CRITH2,chimes_table_H2_collis_dissoc.critical_density_H2,t.nt);
 ADD(CRITHE,chimes_table_H2_collis_dissoc.critical_density_He,t.nt);
 ADD(COT,chimes_table_CO_cosmic_ray.rates,t.count[CO][1]*t.nt);
#undef ADD
 append_cooling(h);
 h.valid=true;return true;
}
bool build_tables() {
 std::call_once(tables_once,[]{tables_ready=build_tables_impl();});
 return tables_ready;
}
__device__ void index(const double *a,int n,double x,int &i,double &w){
 double den=(a[n-1]-a[0])/(n-1.0);
 if(x<=a[0]){i=0;w=0;}else if(x>=a[n-1]){i=n-2;w=1;}
 else{i=(int)floor((x-a[0])/den);w=(x-a[i])/den;}
}
__device__ double interp(const double *v,int i,double w){return (1.-w)*v[i]+w*v[i+1];}
__device__ double powten(double x){return exp(x*2.30258509299404568402);}
__global__ void reactions(Tables t,const double *v,const Reaction *r,const int *entry,const Input *in,Output *out,CompactOutput *packed){
 const Input &x=in[blockIdx.x];Output &y=out[blockIdx.x];
 __shared__ int ti,di,xi;__shared__ double tw,dw,xw,ncrit;
 __shared__ int valid;
 if(threadIdx.x==0){
   valid=1;
   index(v+t.offset[AXT],t.nt,log10(x.t),ti,tw);
   index(v+t.offset[AXD],t.nd,log10(x.td),di,dw);
   index(v+t.offset[AXX],t.nx,log10(fmax(x.x[sp_HII],DBL_MIN)),xi,xw);
   y.crit[0]=powten(interp(v+t.offset[CRITH],ti,tw));
   y.crit[1]=powten(interp(v+t.offset[CRITH2],ti,tw));
   y.crit[2]=powten(interp(v+t.offset[CRITHE],ti,tw));
   ncrit=x.x[sp_HI]/y.crit[0];ncrit+=2.*x.x[sp_H2]/y.crit[1];
   ncrit+=x.x[sp_HeI]/y.crit[2];ncrit*=x.nh;
   if(packed)for(double z:y.crit)if(!isfinite(z))valid=0;
 }
 __syncthreads();
 for(int k=threadIdx.x;k<t.nr;k+=blockDim.x){
   const Reaction q=r[k];double z=0;int i=q.local;y.coefficient[k]=0;
   if(i<t.count[q.kind][x.mol])switch(q.kind){
     case TD:
       y.coefficient[k]=powten(interp(v+t.offset[TDT]+i*t.nt,ti,tw));
       z=y.coefficient[k]*x.x[q.a[0]]*x.x[q.a[1]]*x.nh;
       if(q.a[2]>=0)z*=x.nh*x.x[q.a[2]];break;
     case CON:z=v[t.offset[CONT]+i]*x.x[q.a[0]]*x.x[q.a[1]]*x.nh;break;
     case REC:z=powten(interp(v+t.offset[RECT]+(i*2+x.ab[i])*t.nt,ti,tw))*x.x[q.a[0]]*x.x[q.a[1]]*x.nh;break;
     case GRAIN:
       // N_spectra=0 fixes Psi to its first axis value; 2D interpolation
       // retains the same nonzero terms as CPU (dPsi=0).
       z=powten((1.-tw)*v[t.offset[GRAINT]+(i*t.nt+ti)*t.np]+tw*v[t.offset[GRAINT]+(i*t.nt+ti+1)*t.np]);
       z=z*x.x[q.a[0]]*x.x[q.a[1]]*x.nh*x.dust*x.boost;break;
     case CR:
       if(x.cr>0){z=x.cr*v[t.offset[CRT]+i]*x.x[q.a[0]];
         for(int j=0;j<2;j++)if(i==t.secondary[j])z*=1.+powten(interp(v+t.offset[CRSEC]+j*t.nx,xi,xw));}break;
     case H2D:{const double *p=v+t.offset[H2DT];
       double a=(1.-tw)*(1.-dw)*p[ti*t.nd+di];a+=(1.-tw)*dw*p[ti*t.nd+di+1];
       a+=tw*(1.-dw)*p[(ti+1)*t.nd+di];a+=tw*dw*p[(ti+1)*t.nd+di+1];
       z=powten(a)*x.x[q.a[0]]*x.dust*x.boost*x.nh;break;}
     case H2C:
       z=powten((ncrit/(1.+ncrit))*interp(v+t.offset[KLTE]+i*t.nt,ti,tw)+(1./(1.+ncrit))*interp(v+t.offset[KZERO]+i*t.nt,ti,tw));
       y.coefficient[k]=z;
       z=z*x.x[q.a[0]]*x.x[q.a[1]]*x.nh;break;
     case CO:if(x.cr>0)z=powten(interp(v+t.offset[COT]+i*t.nt,ti,tw))*x.x[sp_H2]*x.cr*sqrt(fmax(x.x[q.a[0]],0.));break;
   }
   y.rate[k]=z;
   if(packed && (!isfinite(z)||!isfinite(y.coefficient[k])))atomicExch(&valid,0);
 }
 __syncthreads();
 for(int s=threadIdx.x;s<ns;s+=blockDim.x){
   double d=0,c=0;
   for(int j=t.begin[s];j<t.begin[s+1];j++)d+=y.rate[entry[j]];
   for(int j=t.begin[ns+s];j<t.begin[ns+s+1];j++)c+=y.rate[entry[j]];
   y.destroy[s]=d;y.create[s]=c;
   if(packed){
     packed[blockIdx.x].destroy[s]=d;packed[blockIdx.x].create[s]=c;
     if(!isfinite(d)||!isfinite(c))atomicExch(&valid,0);
   }
 }
 if(packed){
   __syncthreads();
   if(threadIdx.x==0){
     auto &p=packed[blockIdx.x];p.valid=valid;
     for(int j=0;j<3;j++)p.crit[j]=y.crit[j];
     for(int j=0;j<compact_cr_max;j++)p.cosmic[j]=j<t.count[CR][x.mol]?y.rate[t.base[CR]+j]:0.;
     for(int j=0;j<5;j++)p.molecular[j]=0.;
     if(x.mol){
       p.molecular[0]=y.coefficient[t.base[TD]+t.cooling_index[0]];
       p.molecular[1]=y.coefficient[t.base[H2C]+t.cooling_index[1]];
       p.molecular[2]=y.rate[t.base[TD]+t.cooling_index[2]];
       p.molecular[3]=y.rate[t.base[CON]+t.cooling_index[3]];
       p.molecular[4]=y.rate[t.base[H2D]];
     }
   }
 }
}
#include "snrt_chimes_dark_cooling.cuh"
struct Buffers {
 int device=-1;bool ready=false;
 double *v=nullptr;Reaction *r=nullptr;int *entry=nullptr;
 Input *hi=nullptr,*di=nullptr;Output *ho=nullptr,*dout=nullptr;
 CompactOutput *hp=nullptr,*dp=nullptr;double *cr_heat=nullptr;bool compact=false;
 void clear(){if(device>=0)cudaSetDevice(device);cudaFree(v);cudaFree(r);cudaFree(entry);cudaFree(di);cudaFree(dout);cudaFree(dp);cudaFree(cr_heat);cudaFreeHost(hi);cudaFreeHost(ho);cudaFreeHost(hp);*this=Buffers{};}
 bool init(cudaStream_t stream){
   if(ready)return true;
   if(cudaGetDevice(&device)!=cudaSuccess)return false;
   auto &h=host;
   bool ok=cudaMalloc(&v,h.v.size()*sizeof(double))==cudaSuccess && cudaMalloc(&r,h.r.size()*sizeof(Reaction))==cudaSuccess &&
     cudaMalloc(&entry,h.entries.size()*sizeof(int))==cudaSuccess && cudaMallocHost(&hi,nb*sizeof(Input))==cudaSuccess &&
     cudaMallocHost(&ho,nb*sizeof(Output))==cudaSuccess && cudaMalloc(&di,nb*sizeof(Input))==cudaSuccess && cudaMalloc(&dout,nb*sizeof(Output))==cudaSuccess &&
     cudaMallocHost(&hp,nb*sizeof(CompactOutput))==cudaSuccess && cudaMalloc(&dp,nb*sizeof(CompactOutput))==cudaSuccess &&
     cudaMalloc(&cr_heat,nb*sizeof(double))==cudaSuccess;
   if(ok)ok=cudaMemcpyAsync(v,h.v.data(),h.v.size()*sizeof(double),cudaMemcpyHostToDevice,stream)==cudaSuccess &&
     cudaMemcpyAsync(r,h.r.data(),h.r.size()*sizeof(Reaction),cudaMemcpyHostToDevice,stream)==cudaSuccess &&
     cudaMemcpyAsync(entry,h.entries.data(),h.entries.size()*sizeof(int),cudaMemcpyHostToDevice,stream)==cudaSuccess;
   if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
   if(!ok){clear();return false;}ready=true;uploads++;return true;
 }
};
Buffers buffers[MAX_CUDA_STREAMS];
Buffers broker_extra[MAX_CUDA_STREAMS];
struct Request {
 UserData *d;Input input;bool claimed=false,done=false,ready=false;int ok=0;bool required=false;
 int cooling_valid=0;double net_cooling=0;
 const CompactOutput *crow=nullptr;const Output *frow=nullptr;std::atomic<int> *consume_left=nullptr;
};
// Each host continuation retains its own CVODE stack and candidate state.
// There are no additional OS threads and no cross-cell physical dependency.
struct TileFiber {
 FastCtx ctx{};std::vector<char> stack;Request *waiting=nullptr;
 bool done=false,retry=false,queued=false;int backend=-2,global_cell=-1;
#ifdef SNRT_CHIMES_RHS_TESTING
 int calls=0;
#endif
};
struct TileRun {
 FastCtx sched{};int slot=-1,active=-1;bool shared=false;
 snrt_chimes_tile_cell_fn callback=nullptr;void *context=nullptr;
};
thread_local TileRun *tile_run=nullptr;
thread_local std::vector<TileFiber> tile_fibers;
std::mutex mutex;std::condition_variable changed;std::vector<Request*> pending;
bool select_compact(const std::vector<Request*> &requests){
 if(!compact_results || verify || !host.t.compact_safe)return false;
 // rt_update_flux may remain enabled to preserve case-B recombination even
 // when this cold dark solve has no photon variables. Its registered host
 // callbacks then write exact-zero photo rates; compact chemistry fields
 // must still be restored before those callbacks and host cooling run.
 for(auto *q:requests)if(q->d->myGlobalVars->N_spectra)return false;
 return true;
}
bool copy_result_async(Buffers &b,int n,cudaStream_t stream){
 return (b.compact?cudaMemcpyAsync(b.hp,b.dp,n*sizeof(CompactOutput),cudaMemcpyDeviceToHost,stream):
   cudaMemcpyAsync(b.ho,b.dout,n*sizeof(Output),cudaMemcpyDeviceToHost,stream))==cudaSuccess;
}
void apply_compact_result(Request &q,const CompactOutput &o){
 if(!o.valid){invalid_device_outputs++;q.ok=0;return;}
 auto &d=*q.d;const auto &t=host.t;auto *c=d.chimes_current_rates;
 std::memcpy(c->cosmic_ray_rate,o.cosmic,t.count[CR][d.mol_flag_index]*sizeof(double));
 if(d.mol_flag_index){
   c->T_dependent_rate_coefficient[t.cooling_index[0]]=o.molecular[0];
   c->H2_collis_dissoc_rate_coefficient[t.cooling_index[1]]=o.molecular[1];
   c->T_dependent_rate[t.cooling_index[2]]=o.molecular[2];
   c->constant_rate[t.cooling_index[3]]=o.molecular[3];
   c->H2_dust_formation_rate=o.molecular[4];
   c->H2_collis_dissoc_crit_H=o.crit[0];c->H2_collis_dissoc_crit_H2=o.crit[1];c->H2_collis_dissoc_crit_He=o.crit[2];
 }
 for(int s=0;s<ns;s++){d.species[s].destruction_rate=o.destroy[s];d.species[s].creation_rate=o.create[s];}
 q.ok=1;q.net_cooling=o.net_cooling;
 q.cooling_valid=o.cooling_valid && std::isfinite(o.net_cooling);
 used++;
}
void apply_full_result(Request &q,const Output &o){
 auto &d=*q.d;auto &t=host.t;auto *c=d.chimes_current_rates;
 bool valid=true;for(int k=0;k<t.nr;k++)if(!std::isfinite(o.rate[k])||!std::isfinite(o.coefficient[k]))valid=false;
 for(double v:o.crit)if(!std::isfinite(v))valid=false;
 for(int s=0;s<ns;s++)if(!std::isfinite(o.destroy[s])||!std::isfinite(o.create[s]))valid=false;
 if(!valid){invalid_device_outputs++;q.ok=0;return;}
#define COPY(field,group) std::memcpy(c->field,o.rate+t.base[group],t.count[group][d.mol_flag_index]*sizeof(double))
 COPY(T_dependent_rate,TD);COPY(constant_rate,CON);COPY(recombination_AB_rate,REC);
 COPY(grain_recombination_rate,GRAIN);COPY(cosmic_ray_rate,CR);
 std::memcpy(c->T_dependent_rate_coefficient,o.coefficient+t.base[TD],t.count[TD][d.mol_flag_index]*sizeof(double));
 if(d.mol_flag_index){
   std::memcpy(c->H2_collis_dissoc_rate_coefficient,o.coefficient+t.base[H2C],t.count[H2C][1]*sizeof(double));
   c->H2_collis_dissoc_crit_H=o.crit[0];c->H2_collis_dissoc_crit_H2=o.crit[1];c->H2_collis_dissoc_crit_He=o.crit[2];
 }
 if(d.mol_flag_index){c->H2_dust_formation_rate=o.rate[t.base[H2D]];COPY(H2_collis_dissoc_rate,H2C);COPY(CO_cosmic_ray_rate,CO);}
#undef COPY
 for(int s=0;s<ns;s++){d.species[s].destruction_rate=o.destroy[s];d.species[s].creation_rate=o.create[s];}
 q.ok=1;q.net_cooling=o.net_cooling;
 q.cooling_valid=o.cooling_valid && std::isfinite(o.net_cooling);
 used++;
}
void consume_outputs(std::vector<Request*> &requests,Buffers &b,bool ok){
 const int n=(int)requests.size();
 if(ok){
   batches++;
   if(b.compact)compact_batches++;else full_batches++;
   result_bytes+=(unsigned long long)n*(b.compact?sizeof(CompactOutput):sizeof(Output));
   for(int i=0;i<n;i++){
     if(b.compact)apply_compact_result(*requests[i],b.hp[i]);
     else apply_full_result(*requests[i],b.ho[i]);
   }
 }else{errors+=n;}
}
void execute(std::vector<Request*> &requests,int slot,bool release=true){
 const int worker=omp_get_thread_num();
 if(gpu_thread0_only && worker!=0)std::abort();
 if(worker<64)used_worker_mask.fetch_or(1ull<<worker);
 auto stream=cuda_get_stream_internal(slot);auto &b=buffers[slot];
 bool ok=build_tables() && b.init(stream);int n=(int)requests.size();
 bool tracking=false;
 if(ok){
   const auto active=active_batches.fetch_add(1)+1;
   auto peak=peak_active_batches.load();
   while(peak<active && !peak_active_batches.compare_exchange_weak(peak,active)){}
   if(slot>=0 && slot<64)used_stream_mask.fetch_or(1ull<<slot);
   tracking=true;
   for(int i=0;i<n;i++)b.hi[i]=requests[i]->input;
   ok=cudaMemcpyAsync(b.di,b.hi,n*sizeof(Input),cudaMemcpyHostToDevice,stream)==cudaSuccess;
   b.compact=select_compact(requests);
   if(ok){reactions<<<n,256,0,stream>>>(host.t,b.v,b.r,b.entry,b.di,b.dout,b.compact?b.dp:nullptr);ok=cudaGetLastError()==cudaSuccess;}
   if(ok)ok=launch_dark_cooling(host.t,cool,b.v,b.entry,b.di,b.dout,b.compact?b.dp:nullptr,b.cr_heat,n,stream);
   if(ok)ok=copy_result_async(b,n,stream);
 }
 if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
 if(tracking)active_batches.fetch_sub(1);
 consume_outputs(requests,b,ok);
 if(!ok)b.clear();
 if(release)cuda_release_stream(slot);
}

void tile_entry(int cell){
 const int global=tile_fibers[cell].global_cell;
 tile_run->callback(global>=0?global:cell,tile_run->context);
 tile_fibers[cell].done=true;
}

__attribute__((noinline,noclone)) static void fiber_trampoline(){
 const int cell=tile_run->active;
 tile_entry(cell);
 snrt_chimes_fast_swap(&tile_fibers[cell].ctx,&tile_run->sched);
}
void fast_reset(FastCtx &ctx,void *stack,size_t bytes,std::uint32_t mxcsr,std::uint16_t fcw){
 uintptr_t base=reinterpret_cast<uintptr_t>(stack);
 uintptr_t top=(base+bytes)&~uintptr_t(15);
 if(top<base+32)std::abort();
 uintptr_t sp=top-16;
 *reinterpret_cast<void(**)()>(sp)=&fiber_trampoline;
 ctx={};ctx.sp=sp;ctx.mxcsr=mxcsr;ctx.fcw=fcw;
}
void fill_input(Input &in,const UserData *d){
 auto *c=d->myGlobalVars;auto *g=d->myGasVars;
 std::memcpy(in.x,g->abundances,sizeof(in.x));
 in.t=g->temperature;in.td=c->grain_temperature;
 in.nh=g->nH_tot;in.dust=g->dust_ratio;in.boost=g->dust_boost_factor;in.cr=g->cr_rate;
 in.mol=d->mol_flag_index;in.ab[0]=d->case_AB_index[0];in.ab[1]=d->case_AB_index[1];
 in.cell_size=g->cell_size;in.div_vel=g->divVel;in.extinction=d->extinction;
 in.constant_heating=g->constant_heating_rate;in.cmb=c->cmb_temperature;in.tmol=c->Tmol_K;
 in.co_column=d->CO_column;in.h2o_column=d->H2O_column;in.oh_column=d->OH_column;
 std::memcpy(in.elem,g->element_abundances,sizeof(in.elem));
 in.static_mol=c->StaticMolCooling;in.hybrid_mode=c->hybrid_cooling_mode;in.elem_mask=0;
 for(int e=0;e<9;e++)if(c->element_included[e])in.elem_mask|=1<<e;
 in.rewrite_columns=snrt_chimes_molecular_coefficients && c->StaticMolCooling ? 1 : 0;
}

bool serve_tile_request(Request &q){
 auto &fiber=tile_fibers[tile_run->active];fiber.waiting=&q;
 // The scheduler restores these thread-local adapter controls before
 // resuming this same continuation. The CHIMES exact-rate cache is invalidated
 // at each switch, never reused across cells or previous candidate states.
 snrt_chimes_fast_swap(&fiber.ctx,&tile_run->sched);
 return q.ok!=0;
}

bool serve_request(Request &q,int team){
 const int nstreams=std::max(1,cuda_get_n_streams());
 const int target=std::max(4,std::min(nb,(team+nstreams-1)/nstreams));
 const auto deadline=std::chrono::steady_clock::now()+std::chrono::microseconds(50);
 std::unique_lock<std::mutex> lock(mutex);
 pending.push_back(&q);
 changed.notify_all();

 for(;;){
   if(q.done)return q.ok!=0;
   if(q.claimed){changed.wait(lock,[&]{return q.done;});continue;}
   if(pending.size()<(size_t)target && std::chrono::steady_clock::now()<deadline){
     changed.wait_until(lock,deadline,[&]{return q.done||q.claimed||pending.size()>=(size_t)target;});
     continue;
   }

   // Reserve one disjoint request group before competing for a stream. Without
   // this reservation, every waiter can race for leases and fall back before
   // any leader has claimed the batch it could have served.
   std::vector<Request*> batch;
   batch.reserve((size_t)target);
   auto own=std::find(pending.begin(),pending.end(),&q);
   if(own!=pending.end())pending.erase(own);
   q.claimed=true;batch.push_back(&q);
   while(batch.size()<(size_t)target && !pending.empty()){
     Request *next=pending.front();pending.erase(pending.begin());
     if(next->claimed||next->done)continue;
     next->claimed=true;batch.push_back(next);
   }
   const bool required=std::any_of(batch.begin(),batch.end(),[](Request *r){return r->required;});
   if(batch.size()<4 && !required){
     // A short tail is cheaper and more predictable on CPU. Complete every
     // queued caller in this tail together; no stream launch is worthwhile.
     for(auto *r:batch)r->done=true;
     while(!pending.empty() && batch.size()<4){
       Request *next=pending.front();pending.erase(pending.begin());
       if(next->claimed||next->done)continue;
       next->claimed=true;next->done=true;batch.push_back(next);
     }
     changed.notify_all();
     return q.ok!=0;
   }

   lock.unlock();
   int slot=cuda_acquire_stream();
   // Legacy per-callback callers retain bounded CPU retry on contention.
   // The staged path instead reserves a lease for the whole cohort, so it
   // never waits here and cannot discard work due to per-cell lease churn.
   if(slot<0){
     lock.lock();
     for(auto *r:batch)r->done=true;
     changed.notify_all();
     return q.ok!=0; // CPU fallback, or full retry if the pool disappeared
   }
   execute(batch,slot);
   lock.lock();
   for(auto *r:batch)r->done=true;
   changed.notify_all();
 }
}
}
extern "C" int snrt_chimes_rhs_configure_mode(int use_cuda){
 if(use_cuda!=0 && use_cuda!=1)return 1;
 nb=64;
 const char *batch_env=std::getenv("SNRT_CHIMES_GPU_BATCH_CELLS");
 if(batch_env && *batch_env){
   char *end=nullptr;const long value=std::strtol(batch_env,&end,10);
   if(*end || value<16 || value>nb_max || value%16){
     std::fprintf(stderr,"Invalid SNRT_CHIMES_GPU_BATCH_CELLS=%s (16..%d, multiple of 16)\n",batch_env,nb_max);return 1;
   }
   nb=(int)value;
 }
 const char *result_env=std::getenv("SNRT_CHIMES_GPU_RESULT");
 if(result_env && *result_env && std::strcmp(result_env,"compact") && std::strcmp(result_env,"full")){
   std::fprintf(stderr,"Unknown SNRT_CHIMES_GPU_RESULT=%s; expected compact|full\n",result_env);return 1;
 }
 compact_results=!result_env || std::strcmp(result_env,"full");
 const char *worker_policy=std::getenv("SNRT_CHIMES_GPU_WORKER");
 if(worker_policy && *worker_policy && std::strcmp(worker_policy,"any") &&
    std::strcmp(worker_policy,"thread0")){
   std::fprintf(stderr,"Unknown SNRT_CHIMES_GPU_WORKER=%s; expected any|thread0\n",worker_policy);
   return 1;
 }
 gpu_thread0_only=worker_policy && !std::strcmp(worker_policy,"thread0");
 const char *scheduler=std::getenv("SNRT_CHIMES_CELL_SCHEDULER");
 if(scheduler && *scheduler && std::strcmp(scheduler,"cohort") && std::strcmp(scheduler,"level_queue")){
   std::fprintf(stderr,"Unknown SNRT_CHIMES_CELL_SCHEDULER=%s; expected cohort|level_queue\n",scheduler);return 1;
 }
 level_queue_enabled=scheduler && !std::strcmp(scheduler,"level_queue");
 const char *broker_env=std::getenv("SNRT_CHIMES_GPU_BROKERS");
 level_gpu_brokers=1;
 brokers_explicit=false;
 if(broker_env && *broker_env){
   char *end=nullptr;const long value=std::strtol(broker_env,&end,10);
   if(*end || value<1 || value>64 || (gpu_thread0_only && value!=1)){
     std::fprintf(stderr,"Invalid SNRT_CHIMES_GPU_BROKERS=%s (1..64; thread0 requires 1)\n",broker_env);return 1;
   }
   level_gpu_brokers=(int)value;
   brokers_explicit=true;
 }
 const char *verify_env=std::getenv("SNRT_CHIMES_RHS_VERIFY");
 verify=verify_env && *verify_env && std::strcmp(verify_env,"0");
 verify_return_gpu=verify && !std::strcmp(verify_env,"gpu");
 const char *cooling_env=std::getenv("SNRT_CHIMES_GPU_COOLING");
 if(cooling_env && *cooling_env && std::strcmp(cooling_env,"0") && std::strcmp(cooling_env,"1")){
   std::fprintf(stderr,"Invalid SNRT_CHIMES_GPU_COOLING=%s; expected 0|1\n",cooling_env);return 1;
 }
 gpu_cooling=!cooling_env || !std::strcmp(cooling_env,"1") || !*cooling_env;
 if(!use_cuda){enabled=false;return 0;}
 enabled=true;
 std::fprintf(stderr,"CHIMES GPU worker policy: %s\n",gpu_thread0_only?"thread0":"any");
 std::fprintf(stderr,"CHIMES GPU batch policy: cells=%d maximum=%d\n",nb,nb_max);
 std::fprintf(stderr,"CHIMES GPU result policy: %s full_bytes=%zu compact_bytes=%zu (paired verification uses full)\n",
   compact_results?"compact":"full",sizeof(Output),sizeof(CompactOutput));
 std::fprintf(stderr,"CHIMES GPU cooling policy: %s\n",gpu_cooling?"device-dark":"host");
 return 0;
}
extern "C" unsigned long long snrt_chimes_rhs_worker_mask(void){return used_worker_mask.load();}
extern "C" void snrt_chimes_rhs_result_counts(unsigned long long out[3]){
 out[0]=compact_batches.load();out[1]=full_batches.load();out[2]=result_bytes.load();
}
extern "C" int snrt_chimes_rhs_cpu_direct_active(void){return cpu_direct;}
extern "C" int snrt_chimes_rhs_level_requested(void){return level_queue_enabled;}
extern "C" int snrt_chimes_rhs_configure(void){
 const char *s=std::getenv("SNRT_CHIMES_RHS_BACKEND");
 if(!s||!*s||!std::strcmp(s,"cpu"))return snrt_chimes_rhs_configure_mode(0);
 if(std::strcmp(s,"cuda_dark_reactions")){std::fprintf(stderr,"Unknown SNRT_CHIMES_RHS_BACKEND=%s\n",s);return 1;}
 return snrt_chimes_rhs_configure_mode(1);
}
extern "C" int snrt_chimes_dark_rhs_gpu(UserData *d){
 gpu_net_cooling_ready=0;
 if(!enabled)return 0;
 if(cpu_direct)return 0; // native CPU lane: no per-RHS shared atomic
 if(gpu_thread0_only && omp_get_thread_num()!=0){cell_backend=0;fallback++;return 0;}
 if(tile_run && tile_run->slot<0 && !tile_run->shared){cell_backend=0;fallback++;return 0;}
 if(cell_backend<=-2 || cell_backend==0){fallback++;return 0;}
 // A discarded GPU trajectory must stay discarded until the bridge starts
 // its whole-cell CPU retry. CVODE may call f again while unwinding an
 // unsuccessful Jacobian/Newton solve; never resume device arithmetic there.
 if(retry_cell)return -1;
#ifdef SNRT_CHIMES_RHS_TESTING
 if(cell_backend==1 && drop_after && ++cell_calls==drop_after){retry_cell=true;drops++;return -1;}
#endif
 auto decline=[&](){if(cell_backend==1){retry_cell=true;return -1;}cell_backend=0;fallback++;return 0;};
 auto *c=d->myGlobalVars;auto *g=d->myGasVars;
 int team=omp_in_parallel()?omp_get_num_threads():1;
 if(c->N_spectra!=0 || c->totalNumberOfSpecies!=ns || g->ThermEvolOn!=1 ||
    (!tile_run && team<4) || (d->mol_flag_index!=0 && d->mol_flag_index!=1))return decline();
 for(int s=0;s<ns;s++)if(c->speciesIndices[s]!=s)return decline();
 if(d->case_AB_index[0]<0||d->case_AB_index[0]>1||d->case_AB_index[1]<0||d->case_AB_index[1]>1)return decline();
 Request q{};q.d=d;q.required=cell_backend==1;fill_input(q.input,d);
 if(!(q.input.t>0) || !(q.input.td>0) || !(q.input.nh>0) || !std::isfinite(q.input.t))return decline();
 if(!(tile_run?serve_tile_request(q):serve_request(q,team))){
   if(cell_backend==1){retry_cell=true;return -1;}
   cell_backend=0;fallback++;
   return 0;
 }else cell_backend=1;
 if(q.ok && q.cooling_valid && (!verify || verify_return_gpu)){
   if(!(c->hybrid_cooling_mode==1 && !c->hybrid_cooling_fn)){
     gpu_net_cooling=q.net_cooling;gpu_net_cooling_ready=1;
     if(c->hybrid_cooling_mode==1)(void)c->hybrid_cooling_fn(g,c);
   }
 }
 if(q.ok && c->rt_update_flux){
   // Preserve host extension callbacks on their owning thread. In dark
   // mode the current receiver resets photo coefficients/adds zero sources.
   if(snrt_chimes_molecular_coefficients)snrt_chimes_molecular_coefficients(g,c,*d);
   if(snrt_chimes_secondary_rates)snrt_chimes_secondary_rates(d->species,g,c);
 }
 if(q.ok && verify){
   // Compare both arithmetic paths at this trajectory's exact CVODE state.
   // verify=cpu retains the CPU RHS; verify=gpu restores the device RHS after
   // measuring the paired CPU result, so the latter follows a GPU trajectory.
   double device_rates[nrmax]={0},device_coeff[nrmax]={0},device_create[ns],device_destroy[ns],device_crit[3]={0};
   auto &t=host.t;
   auto *r=d->chimes_current_rates;
   double *device_fields[]={r->T_dependent_rate,r->constant_rate,r->recombination_AB_rate,
       r->grain_recombination_rate,r->cosmic_ray_rate};
   const int device_groups[]={TD,CON,REC,GRAIN,CR};
   for(int j=0;j<5;j++)std::memcpy(device_rates+t.base[device_groups[j]],device_fields[j],
       t.count[device_groups[j]][d->mol_flag_index]*sizeof(double));
   if(d->mol_flag_index){
     device_rates[t.base[H2D]]=r->H2_dust_formation_rate;
     std::memcpy(device_rates+t.base[H2C],r->H2_collis_dissoc_rate,t.count[H2C][1]*sizeof(double));
     std::memcpy(device_rates+t.base[CO],r->CO_cosmic_ray_rate,t.count[CO][1]*sizeof(double));
   }
   std::memcpy(device_coeff+t.base[TD],r->T_dependent_rate_coefficient,
       t.count[TD][d->mol_flag_index]*sizeof(double));
   if(d->mol_flag_index){
     std::memcpy(device_coeff+t.base[H2C],r->H2_collis_dissoc_rate_coefficient,
         t.count[H2C][1]*sizeof(double));
     device_crit[0]=r->H2_collis_dissoc_crit_H;device_crit[1]=r->H2_collis_dissoc_crit_H2;
     device_crit[2]=r->H2_collis_dissoc_crit_He;
   }
   for(int s=0;s<ns;s++){device_create[s]=d->species[s].creation_rate;device_destroy[s]=d->species[s].destruction_rate;}
   double device_cool=calculate_total_cooling_rate(g,c,*d,0);
   // Cache invalidation is a C function in the patched private CHIMES.
   snrt_chimes_rate_cache_invalidate();
   update_rate_coefficients(g,c,*d,g->ThermEvolOn);
   if(c->rt_update_flux)update_rt_photochemistry_coefficients(g,c,*d);
   update_rates(g,c,*d);update_rate_vector(d->species,g,c,*d);
   double dr=0,ds=0,dn=0,dco_abs=0,dco_rel=0,dc=0;
   unsigned long long nf_device=0,nf_cpu=0;
   auto observe=[&](double a,double b,double &relative,double *absolute=nullptr){
     bool finite_device=std::isfinite(a),finite_cpu=std::isfinite(b);
     if(!finite_device)nf_device++;
     if(!finite_cpu)nf_cpu++;
     if(!finite_device||!finite_cpu)return;
     double delta=std::abs(a-b);
     relative=std::max(relative,delta/std::max({std::abs(a),std::abs(b),1e-60}));
     if(absolute)*absolute=std::max(*absolute,delta);
   };
   auto compare_rates=[&](int kind,const double *cpu,int molecule,double &maximum){
     for(int j=0;j<t.count[kind][molecule];j++)observe(device_rates[t.base[kind]+j],cpu[j],maximum);
   };
   compare_rates(TD,r->T_dependent_rate,d->mol_flag_index,dr);
   compare_rates(CON,r->constant_rate,d->mol_flag_index,dr);
   compare_rates(REC,r->recombination_AB_rate,d->mol_flag_index,dr);
   compare_rates(GRAIN,r->grain_recombination_rate,d->mol_flag_index,dr);
   compare_rates(CR,r->cosmic_ray_rate,d->mol_flag_index,dr);
   for(int j=0;j<t.count[TD][d->mol_flag_index];j++)
     observe(device_coeff[t.base[TD]+j],r->T_dependent_rate_coefficient[j],dr);
   if(d->mol_flag_index){
     observe(device_rates[t.base[H2D]],r->H2_dust_formation_rate,dr);
     compare_rates(H2C,r->H2_collis_dissoc_rate,1,dr);
     compare_rates(CO,r->CO_cosmic_ray_rate,1,dr);
     for(int j=0;j<t.count[H2C][1];j++)
       observe(device_coeff[t.base[H2C]+j],r->H2_collis_dissoc_rate_coefficient[j],dr);
     observe(device_crit[0],r->H2_collis_dissoc_crit_H,dr);
     observe(device_crit[1],r->H2_collis_dissoc_crit_H2,dr);
     observe(device_crit[2],r->H2_collis_dissoc_crit_He,dr);
   }
   for(int s=0;s<ns;s++){
     double cpu_create=d->species[s].creation_rate,cpu_destroy=d->species[s].destruction_rate;
     observe(device_create[s],cpu_create,ds);observe(device_destroy[s],cpu_destroy,ds);
     double device_net=device_create[s]-device_destroy[s],cpu_net=cpu_create-cpu_destroy;
     observe(device_net,cpu_net,dn);
     if(s==sp_CO)observe(device_net,cpu_net,dco_rel,&dco_abs);
   }
   if(q.cooling_valid)observe(q.net_cooling,device_cool,dc);
   observe(device_cool,calculate_total_cooling_rate(g,c,*d,0),dc);
   std::lock_guard<std::mutex> guard(verification_mutex);
   verify_rates=std::max(verify_rates,dr);verify_species=std::max(verify_species,ds);
   verify_net=std::max(verify_net,dn);verify_co_abs=std::max(verify_co_abs,dco_abs);
   verify_co_rel=std::max(verify_co_rel,dco_rel);verify_cooling=std::max(verify_cooling,dc);
   verify_nonfinite_device+=nf_device;verify_nonfinite_cpu+=nf_cpu;
   if(verify_return_gpu){
#define RESTORE(field,k) std::memcpy(r->field,device_rates+t.base[k],t.count[k][d->mol_flag_index]*sizeof(double))
     RESTORE(T_dependent_rate,TD);RESTORE(constant_rate,CON);RESTORE(recombination_AB_rate,REC);
     RESTORE(grain_recombination_rate,GRAIN);RESTORE(cosmic_ray_rate,CR);
     std::memcpy(r->T_dependent_rate_coefficient,device_coeff+t.base[TD],
         t.count[TD][d->mol_flag_index]*sizeof(double));
     if(d->mol_flag_index){
       r->H2_dust_formation_rate=device_rates[t.base[H2D]];
       RESTORE(H2_collis_dissoc_rate,H2C);RESTORE(CO_cosmic_ray_rate,CO);
       std::memcpy(r->H2_collis_dissoc_rate_coefficient,device_coeff+t.base[H2C],
           t.count[H2C][1]*sizeof(double));
       r->H2_collis_dissoc_crit_H=device_crit[0];r->H2_collis_dissoc_crit_H2=device_crit[1];
       r->H2_collis_dissoc_crit_He=device_crit[2];
     }
#undef RESTORE
     for(int s=0;s<ns;s++){d->species[s].creation_rate=device_create[s];d->species[s].destruction_rate=device_destroy[s];}
   }
 }
 return q.ok;
}
extern "C" void snrt_chimes_rhs_verification(double out[9]){
 out[0]=verify_rates;out[1]=verify_species;out[2]=verify_cooling;out[3]=verify_net;
 out[4]=verify_co_abs;out[5]=verify_co_rel;out[6]=(double)verify_nonfinite_device;out[7]=(double)verify_nonfinite_cpu;
 out[8]=(double)invalid_device_outputs.load();
}
extern "C" void snrt_chimes_rhs_cell_begin(void){cell_backend=-1;retry_cell=false;
#ifdef SNRT_CHIMES_RHS_TESTING
 cell_calls=0;
#endif
}
extern "C" int snrt_chimes_rhs_cell_retry(void){int retry=retry_cell;retry_cell=false;cell_backend=0;return retry;}
extern "C" void snrt_chimes_rhs_cell_end(void){cell_backend=-2;retry_cell=false;}
extern "C" int snrt_chimes_rhs_tile_active(void){return tile_run!=nullptr;}
#ifdef SNRT_CHIMES_RHS_TESTING
extern "C" int snrt_chimes_rhs_tile_gpu_selected(void){return tile_run && tile_run->slot>=0;}
#endif
extern "C" int snrt_chimes_rhs_tile_run(int cells,snrt_chimes_tile_cell_fn callback,void *context){
 if(cells<0 || cells>128 || !callback || tile_run)return 1;
 if(!cells)return 0;
 TileRun run;run.callback=callback;run.context=context;
 const int saved_backend=cell_backend;const bool saved_retry=retry_cell;
#ifdef SNRT_CHIMES_RHS_TESTING
 const int saved_calls=cell_calls;
#endif
 std::uint32_t worker_mxcsr=0;std::uint16_t worker_fcw=0;
 capture_fp(worker_mxcsr,worker_fcw);
 run.slot=enabled && (!gpu_thread0_only || omp_get_thread_num()==0)?cuda_acquire_stream():-1;
 if(run.slot>=0){
   // Allocate/reuse bounded stacks before running any callbacks. Failure to
   // prepare a GPU cohort is a whole-tile CPU selection, not partial commit.
   try{
     tile_fibers.resize(cells);
     for(int i=0;i<cells;i++){
       auto &fiber=tile_fibers[i];fiber.stack.resize(256*1024);
       fiber.done=false;fiber.waiting=nullptr;fiber.backend=-2;fiber.retry=false;
       fiber.global_cell=-1;fiber.queued=false;
#ifdef SNRT_CHIMES_RHS_TESTING
       fiber.calls=0;
#endif
       fast_reset(fiber.ctx,fiber.stack.data(),fiber.stack.size(),worker_mxcsr,worker_fcw);
     }
   }catch(...){cuda_release_stream(run.slot);run.slot=-1;}
 }
 tile_run=&run;
 if(run.slot<0){
   for(int i=0;i<cells;i++)callback(i,context);
 }else{
   for(;;){
     std::vector<Request*> requests;requests.reserve(cells);int unfinished=0;
     // CPU phase: advance independent cells to their next CURRENT RHS.
     for(int i=0;i<cells;i++){
       auto &fiber=tile_fibers[i];if(fiber.done)continue;
       run.active=i;cell_backend=fiber.backend;retry_cell=fiber.retry;
#ifdef SNRT_CHIMES_RHS_TESTING
       cell_calls=fiber.calls;
#endif
       snrt_chimes_rate_cache_invalidate();
       snrt_chimes_fast_swap(&run.sched,&fiber.ctx);
       fiber.backend=cell_backend;fiber.retry=retry_cell;
#ifdef SNRT_CHIMES_RHS_TESTING
       fiber.calls=cell_calls;
#endif
       if(!fiber.done){unfinished++;requests.push_back(fiber.waiting);}
     }
     if(!unfinished)break;
     // GPU phase: consume staged candidate states into temporary rate arrays.
     // One lease belongs to this entire tile; no per-cell lease wait/restart.
     for(size_t first=0;first<requests.size();first+=nb){
       const size_t last=std::min(requests.size(),first+nb);
       std::vector<Request*> batch(requests.begin()+first,requests.begin()+last);
       execute(batch,run.slot,false);
     }
     for(auto &fiber:tile_fibers)fiber.waiting=nullptr;
     // Next CPU phase resumes CVODE only with this phase's completed results.
   }
   cuda_release_stream(run.slot);
 }
 tile_run=nullptr;cell_backend=saved_backend;retry_cell=saved_retry;
#ifdef SNRT_CHIMES_RHS_TESTING
 cell_calls=saved_calls;
#endif
 snrt_chimes_rate_cache_invalidate();
 return 0;
}
extern "C" int snrt_chimes_rhs_level_run(int cells,snrt_chimes_tile_cell_fn callback,void *context){
 if(cells<0 || !callback || tile_run)return 1;
 if(!cells)return 0;
 const char *trace_env=std::getenv("SNRT_CHIMES_TIMELINE");
 const bool trace=trace_env && !std::strcmp(trace_env,"1") && !timeline_claimed.exchange(true);
 std::chrono::steady_clock::time_point trace_start;
 std::atomic<int> next{0};
 const int workers=std::max(1,omp_get_max_threads());
 /* One broker thread serializes every CVODE continuation it owns. That loses
  * to a full CPU team unless every worker runs continuations and each has
  * enough cells for GPU latency to hide behind the other fibers' host work.
  * Smaller levels stay on the CPU. */
 int active_brokers=0;
 if(enabled){
   active_brokers=gpu_thread0_only?1:(brokers_explicit?std::min(level_gpu_brokers,workers):workers);
   if(workers>1 && cells<16*workers) active_brokers=0;
 }
 pending.clear();
 std::atomic<int> shared_launcher{0};
 std::atomic<int> shared_healthy{1};
 // Counted up front so a broker cannot observe zero before its siblings enter.
 std::atomic<int> brokers_running{0};
 std::vector<unsigned long long> completed(workers,0);
#pragma omp parallel num_threads(workers)
 {
   const int worker=omp_get_thread_num();
   const int saved_backend=cell_backend;const bool saved_retry=retry_cell;
   const bool saved_direct=cpu_direct;
   std::uint32_t worker_mxcsr=0;std::uint16_t worker_fcw=0;
   capture_fp(worker_mxcsr,worker_fcw);
   /* Job 414163 filled one shared batch to ~64 cells and matched the CPU
    * result, but that single stream ran the 512-cell hybrid at 0.34x the
    * 8-thread CPU wall. Each broker keeps its own stream instead. */
   const bool shared_team=false;
   int slot=-1;
   if(shared_team){
     if(worker==0) slot=cuda_acquire_stream();
   }else if(worker<active_brokers){
     slot=cuda_acquire_stream();
   }
   bool prepared=false;
   int initial[128];std::fill(initial,initial+128,-1);
   cudaEvent_t events[2]={nullptr,nullptr};
   if(slot>=0){
     const auto setup_stream=cuda_get_stream_internal(slot);
     prepared=build_tables() && buffers[slot].init(setup_stream) && broker_extra[slot].init(setup_stream);
     try{
       tile_fibers.resize(128);
       for(auto &fiber:tile_fibers){
         fiber.stack.resize(256*1024);fiber.done=true;fiber.global_cell=-1;
         fiber.waiting=nullptr;fiber.queued=false;
       }
     }catch(...){prepared=false;}
     for(auto &event:events)if(prepared)
       prepared=cudaEventCreateWithFlags(&event,cudaEventDisableTiming)==cudaSuccess;
     if(!prepared){
       cudaStreamSynchronize(setup_stream);
       for(auto &event:events)if(event){cudaEventDestroy(event);event=nullptr;}
       cuda_release_stream(slot);slot=-1;
     }
   }
   if(shared_team && worker<active_brokers && slot<0){
     try{
       tile_fibers.resize(128);
       for(auto &fiber:tile_fibers){
         fiber.stack.resize(256*1024);fiber.done=true;fiber.global_cell=-1;
         fiber.waiting=nullptr;fiber.queued=false;
       }
       prepared=true;
     }catch(...){prepared=false;}
   }
   const int reserve_limit=active_brokers?std::min(128,std::max(1,cells/active_brokers)):0;
   const int submit_at=std::max(1,std::min(nb,std::max(8,reserve_limit/2)));
   if(trace && worker==0){
     if(cudaProfilerStart()!=cudaSuccess)std::abort();
     trace_start=std::chrono::steady_clock::now();
     std::fprintf(stderr,"CHIMES_TIMELINE_BEGIN window_s=20 cells=%d\n",cells);
   }
#pragma omp barrier
   if(worker==0) shared_launcher.store(shared_team && slot>=0 && prepared?1:0);
#pragma omp barrier
   const bool shared=shared_launcher.load()!=0;
   const bool fiber=shared?(worker<active_brokers && prepared):slot>=0;
   if(shared && fiber) brokers_running.fetch_add(1);
#pragma omp barrier
   if(fiber)for(int i=0;i<reserve_limit;i++){
     const int cell=next.fetch_add(1);
     if(cell>=cells)break;
     initial[i]=cell;
   }
   cudaStream_t stream=nullptr;
   if(!fiber){
     cpu_direct=true;
     for(int cell=next.fetch_add(1);cell<cells;cell=next.fetch_add(1)){
       callback(cell,context);++completed[worker];
     }
   }else if(!shared){
     struct Flight {
       Buffers *buffer=nullptr;cudaEvent_t event=nullptr;bool busy=false;
       unsigned long long sequence=0;
       std::vector<Request*> requests;std::vector<int> owners;
     };
     Flight flights[2];
     flights[0].buffer=&buffers[slot];flights[1].buffer=&broker_extra[slot];
     for(int i=0;i<2;i++){
       flights[i].event=events[i];flights[i].requests.reserve(nb);flights[i].owners.reserve(nb);
     }
     TileRun run;run.slot=slot;run.callback=callback;run.context=context;
     tile_run=&run;cpu_direct=false;
     stream=cuda_get_stream_internal(slot);
     used_worker_mask.fetch_or(1ull<<worker);used_stream_mask.fetch_or(1ull<<slot);
     bool healthy=true,exhausted=false;
     bool trace_stopped=false;unsigned long long sequence=0;
     auto finish=[&](Flight &flight,bool ok){
       consume_outputs(flight.requests,*flight.buffer,ok);
       for(size_t i=0;i<flight.owners.size();i++){
         auto &fiber=tile_fibers[flight.owners[i]];
         // The request lives on this suspended fiber's stack. It is never
         // resumed, replaced, or recycled before its own completion event.
         if(fiber.waiting!=flight.requests[i] || !fiber.queued)std::abort();
         fiber.waiting=nullptr;fiber.queued=false;
       }
       flight.busy=false;flight.requests.clear();flight.owners.clear();
     };
     // Launch one idle flight. A full group starts as soon as submit_at
     // fibers are blocked so the remaining fibers' host CVODE overlaps it.
     // Leftovers smaller than that wait until every owned fiber is swept.
     auto submit_flight=[&](bool allow_partial)->bool{
       if(!healthy) return false;
       Flight *flight=nullptr;
       for(auto &candidate:flights) if(!candidate.busy){flight=&candidate;break;}
       if(!flight) return false;
       for(int i=0;i<128 && (int)flight->requests.size()<nb;i++){
         auto &fiber=tile_fibers[i];
         if(fiber.waiting && !fiber.queued){
           flight->requests.push_back(fiber.waiting);flight->owners.push_back(i);fiber.queued=true;
         }
       }
       const int n=(int)flight->requests.size();
       if(n==0 || (!allow_partial && n<submit_at)){
         for(int owner:flight->owners) tile_fibers[owner].queued=false;
         flight->requests.clear();flight->owners.clear();
         return false;
       }
       auto &b=*flight->buffer;
       flight->sequence=++sequence;
       if(trace){
         char label[128];std::snprintf(label,sizeof(label),"CHIMES submit worker=%d flight=%ld seq=%llu n=%d",worker,flight-flights,flight->sequence,n);
         nvtxRangePushA(label);
       }
       for(int i=0;i<n;i++) b.hi[i]=flight->requests[i]->input;
       bool ok=cudaMemcpyAsync(b.di,b.hi,n*sizeof(Input),cudaMemcpyHostToDevice,stream)==cudaSuccess;
       b.compact=select_compact(flight->requests);
       if(ok){reactions<<<n,256,0,stream>>>(host.t,b.v,b.r,b.entry,b.di,b.dout,b.compact?b.dp:nullptr);ok=cudaGetLastError()==cudaSuccess;}
       if(ok) ok=launch_dark_cooling(host.t,cool,b.v,b.entry,b.di,b.dout,b.compact?b.dp:nullptr,b.cr_heat,n,stream);
       if(ok) ok=copy_result_async(b,n,stream);
#ifdef SNRT_CHIMES_RHS_TESTING
       if(ok && submit_fail_after && submit_attempts.fetch_add(1)>=submit_fail_after) ok=false;
#endif
       if(ok) ok=cudaEventRecord(flight->event,stream)==cudaSuccess;
       if(trace) nvtxRangePop();
       flight->busy=true;
       if(!ok) healthy=false;
       else peak_active_batches.store(std::max(peak_active_batches.load(),1ull));
       return ok;
     };
     for(;;){
       bool progressed=false,alive=false;
       // Never query a reused event after a failed submission: it may still
       // describe the previous batch. Drain failures before consuming results.
       if(healthy)for(auto &flight:flights)if(flight.busy){
         const auto status=cudaEventQuery(flight.event);
         if(status==cudaSuccess){
           if(trace)nvtxRangePushA("CHIMES consume completed batch");
           finish(flight,true);progressed=true;
           if(trace)nvtxRangePop();
         }
         else if(status!=cudaErrorNotReady)healthy=false;
       }
       if(!healthy && run.slot>=0){
         cudaStreamSynchronize(stream);
         for(auto &flight:flights)if(flight.busy)finish(flight,false);
         // Fail pending CURRENT candidates; the bridge restarts any already
         // GPU-owned cell from its original input with fresh CPU CVODE state.
         for(auto &fiber:tile_fibers)if(fiber.waiting){
           fiber.waiting->ok=0;fiber.waiting=nullptr;fiber.queued=false;
         }
         run.slot=-1;
       }
       if(trace)nvtxRangePushA("CHIMES host continuation sweep");
       for(int i=0;i<128;i++){
         auto &fiber=tile_fibers[i];
         if(fiber.done){
           if(fiber.global_cell>=0){++completed[worker];fiber.global_cell=-1;}
           if(exhausted && initial[i]<0)continue;
           const int cell=initial[i]>=0?initial[i]:next.fetch_add(1);
           initial[i]=-1;
           if(cell>=cells){exhausted=true;continue;}
           fiber.global_cell=cell;fiber.done=false;fiber.retry=false;fiber.backend=-2;
           fiber.waiting=nullptr;fiber.queued=false;
#ifdef SNRT_CHIMES_RHS_TESTING
           fiber.calls=0;
#endif
           fast_reset(fiber.ctx,fiber.stack.data(),fiber.stack.size(),worker_mxcsr,worker_fcw);
         }
         alive=true;
         if(fiber.waiting)continue;
         run.active=i;cell_backend=fiber.backend;retry_cell=fiber.retry;
#ifdef SNRT_CHIMES_RHS_TESTING
         cell_calls=fiber.calls;
#endif
         snrt_chimes_rate_cache_invalidate();
         snrt_chimes_fast_swap(&run.sched,&fiber.ctx);
         fiber.backend=cell_backend;fiber.retry=retry_cell;
#ifdef SNRT_CHIMES_RHS_TESTING
         fiber.calls=cell_calls;
#endif
         progressed=true;
         if(healthy && fiber.waiting && !fiber.queued){
           int unqueued=0;
           for(int j=0;j<128;j++) if(tile_fibers[j].waiting && !tile_fibers[j].queued) unqueued++;
           if(unqueued>=submit_at && submit_flight(false)) progressed=true;
         }
       }
       if(trace)nvtxRangePop();
       if(!alive)break;
       // Cells still running stay off this launch. A short tail launches only
       // after the sweep, when every owned fiber is blocked or finished.
       if(submit_flight(true)) progressed=true;
       // Wait only when all owned continuations are suspended. CPU workers
       // still draw independent cells from the level queue during this wait.
       if(!progressed){
         Flight *oldest=nullptr;
         for(auto &flight:flights)if(flight.busy && (!oldest || flight.sequence<oldest->sequence))oldest=&flight;
         if(oldest){
           if(trace){
             char label[128];std::snprintf(label,sizeof(label),"CHIMES wait worker=%d flight=%ld seq=%llu",worker,oldest-flights,oldest->sequence);
             nvtxRangePushA(label);
           }
           if(cudaEventSynchronize(oldest->event)!=cudaSuccess)healthy=false;
           if(trace)nvtxRangePop();
         }
       }
       if(trace && worker==0 && !trace_stopped &&
          std::chrono::steady_clock::now()-trace_start>=std::chrono::seconds(20)){
         std::fprintf(stderr,"CHIMES_TIMELINE_END window_s=20\n");
         cudaProfilerStop();trace_stopped=true;
       }
     }
     tile_run=nullptr;
   }else{
     if(worker==0) std::fprintf(stderr,"CHIMES shared batch launcher=0 batch=%d\n",nb);
     struct Flight {
       Buffers *buffer=nullptr;cudaEvent_t event=nullptr;bool busy=false,consuming=false;
       unsigned long long sequence=0;std::atomic<int> consume_left{0};
       std::vector<Request*> requests;
     };
     Flight flights[2];
     if(worker==0){
       stream=cuda_get_stream_internal(slot);
       flights[0].buffer=&buffers[slot];flights[1].buffer=&broker_extra[slot];
       for(int i=0;i<2;i++){flights[i].event=events[i];flights[i].requests.reserve(nb);}
       used_worker_mask.fetch_or(1ull<<worker);used_stream_mask.fetch_or(1ull<<slot);
     }
     TileRun run;run.slot=worker==0?slot:-1;run.shared=true;run.callback=callback;run.context=context;
     tile_run=&run;cpu_direct=false;
     bool healthy=worker!=0 || slot>=0,exhausted=false;
     unsigned long long sequence=0;int stall=0,last_have=-1;
     auto fail_request=[&](Request *r){
       r->ok=0;r->crow=nullptr;r->frow=nullptr;r->consume_left=nullptr;r->ready=true;
     };
     auto pending_have=[&](){
       std::lock_guard<std::mutex> lock(mutex);return (int)pending.size();
     };
     auto publish_fiber=[&](TileFiber &one){
       if(!one.waiting || one.queued) return;
       std::lock_guard<std::mutex> lock(mutex);
       if(!one.waiting || one.queued) return;
       auto *r=one.waiting;
       r->done=false;r->claimed=true;r->ok=0;r->ready=false;r->crow=nullptr;r->frow=nullptr;r->consume_left=nullptr;
       pending.push_back(r);one.queued=true;
     };
     auto reclaim_index=[&](int index)->bool{
       auto &one=tile_fibers[index];
       if(!one.waiting || !one.queued) return false;
       const CompactOutput *crow=nullptr;const Output *frow=nullptr;std::atomic<int> *left=nullptr;Request *req=nullptr;
       {
         std::lock_guard<std::mutex> lock(mutex);
         if(!one.waiting || !one.waiting->ready) return false;
         req=one.waiting;crow=req->crow;frow=req->frow;left=req->consume_left;req->ready=false;
       }
       if(crow) apply_compact_result(*req,*crow);
       else if(frow) apply_full_result(*req,*frow);
       else req->ok=0;
       if(left) left->fetch_sub(1,std::memory_order_acq_rel);
       one.waiting=nullptr;one.queued=false;return true;
     };
     auto release_consumed=[&](){
       if(worker!=0) return;
       for(auto &flight:flights){
         if(flight.consuming && flight.consume_left.load(std::memory_order_acquire)==0){
           flight.consuming=false;flight.busy=false;flight.requests.clear();
         }
       }
     };
     auto mark_ready=[&](Flight &flight){
       auto &b=*flight.buffer;const int n=(int)flight.requests.size();
       flight.consume_left.store(n,std::memory_order_relaxed);flight.consuming=true;
       {
         std::lock_guard<std::mutex> lock(mutex);
         for(int i=0;i<n;i++){
           Request *r=flight.requests[i];
           r->crow=b.compact?&b.hp[i]:nullptr;r->frow=b.compact?nullptr:&b.ho[i];
           r->consume_left=&flight.consume_left;r->ready=true;
         }
       }
       batches++;if(b.compact)compact_batches++;else full_batches++;
       result_bytes+=(unsigned long long)n*(b.compact?sizeof(CompactOutput):sizeof(Output));
     };
     auto fail_flight=[&](Flight &flight){
       const int n=(int)flight.requests.size();
       {
         std::lock_guard<std::mutex> lock(mutex);
         for(auto *r:flight.requests) fail_request(r);
       }
       if(n) errors+=n;
       flight.requests.clear();flight.busy=false;flight.consuming=false;flight.consume_left.store(0);
     };
     auto fail_pending=[&](){
       std::lock_guard<std::mutex> lock(mutex);
       for(auto *r:pending) fail_request(r);
       pending.clear();
     };
     auto shared_poll=[&]()->bool{
       if(worker!=0 || !healthy) return false;
       bool any=false;
       for(auto &flight:flights){
         if(!flight.busy || flight.consuming) continue;
         const auto status=cudaEventQuery(flight.event);
         if(status==cudaSuccess){mark_ready(flight);any=true;}
         else if(status!=cudaErrorNotReady){healthy=false;shared_healthy.store(0);fail_flight(flight);any=true;}
       }
       if(!healthy){
         cudaStreamSynchronize(stream);
         for(auto &flight:flights) if(flight.busy && !flight.consuming) fail_flight(flight);
         fail_pending();
       }
       return any;
     };
     auto shared_busy=[&](){
       for(auto &flight:flights) if(flight.busy) return true;
       return false;
     };
     auto submit_shared=[&](bool allow_partial)->bool{
       if(worker!=0 || !healthy) return false;
       Flight *flight=nullptr;
       for(auto &candidate:flights) if(!candidate.busy){flight=&candidate;break;}
       if(!flight) return false;
       {
         std::lock_guard<std::mutex> lock(mutex);
         const int have=(int)pending.size();
         if(have==0 || (!allow_partial && have<32)) return false;
         const int n=std::min(have,nb);
         flight->requests.assign(pending.begin(),pending.begin()+n);
         pending.erase(pending.begin(),pending.begin()+n);
       }
       const int n=(int)flight->requests.size();auto &b=*flight->buffer;
       flight->sequence=++sequence;
       for(int i=0;i<n;i++) b.hi[i]=flight->requests[i]->input;
       bool ok=cudaMemcpyAsync(b.di,b.hi,n*sizeof(Input),cudaMemcpyHostToDevice,stream)==cudaSuccess;
       b.compact=select_compact(flight->requests);
       if(ok){reactions<<<n,256,0,stream>>>(host.t,b.v,b.r,b.entry,b.di,b.dout,b.compact?b.dp:nullptr);ok=cudaGetLastError()==cudaSuccess;}
       if(ok) ok=launch_dark_cooling(host.t,cool,b.v,b.entry,b.di,b.dout,b.compact?b.dp:nullptr,b.cr_heat,n,stream);
       if(ok) ok=copy_result_async(b,n,stream);
       if(ok) ok=cudaEventRecord(flight->event,stream)==cudaSuccess;
       if(!ok){healthy=false;shared_healthy.store(0);fail_flight(*flight);fail_pending();return false;}
       flight->busy=true;flight->consuming=false;
       peak_active_batches.store(std::max(peak_active_batches.load(),1ull));
       return true;
     };
     auto wait_oldest=[&]()->bool{
       if(worker!=0 || !healthy) return false;
       Flight *oldest=nullptr;
       for(auto &flight:flights) if(flight.busy && !flight.consuming && (!oldest || flight.sequence<oldest->sequence)) oldest=&flight;
       if(!oldest) return false;
       if(cudaEventSynchronize(oldest->event)!=cudaSuccess){healthy=false;shared_healthy.store(0);fail_flight(*oldest);fail_pending();return false;}
       mark_ready(*oldest);return true;
     };
     auto runnable_now=[&](){
       for(int i=0;i<128;i++) if(!tile_fibers[i].done && !tile_fibers[i].waiting) return true;
       return false;
     };
     // Launch once 32 cells are waiting and take up to nb. The tail, smaller
     // than that, launches only when this thread has no runnable fiber and the
     // queue has stopped growing. Only worker 0 calls CUDA.
     auto pump=[&](bool allow_partial)->bool{
       if(worker!=0) return false;
       bool did=shared_poll();release_consumed();
       const int have=pending_have();
       if(have>=32){if(submit_shared(false)){did=true;stall=0;}}
       else if(allow_partial && have>0 && have==last_have){if(++stall>=8 && submit_shared(true)){did=true;stall=0;}}
       else if(have!=last_have) stall=0;
       last_have=pending_have();
       if(allow_partial && !did && shared_busy()){if(wait_oldest()) did=true;}
       return did;
     };
     for(;;){
       bool progressed=false,alive=false;
       if(worker==0 && !healthy){shared_poll();fail_pending();}
       for(int i=0;i<128;i++){
         auto &one=tile_fibers[i];
         if(reclaim_index(i)) progressed=true;
         if(one.done){
           if(one.global_cell>=0){++completed[worker];one.global_cell=-1;}
           if(exhausted && initial[i]<0)continue;
           const int cell=initial[i]>=0?initial[i]:next.fetch_add(1);
           initial[i]=-1;
           if(cell>=cells){exhausted=true;continue;}
           one.global_cell=cell;one.done=false;one.retry=false;one.backend=-2;
           one.waiting=nullptr;one.queued=false;
#ifdef SNRT_CHIMES_RHS_TESTING
           one.calls=0;
#endif
           fast_reset(one.ctx,one.stack.data(),one.stack.size(),worker_mxcsr,worker_fcw);
         }
         alive=true;
         if(one.waiting)continue;
         run.active=i;cell_backend=one.backend;retry_cell=one.retry;
#ifdef SNRT_CHIMES_RHS_TESTING
         cell_calls=one.calls;
#endif
         snrt_chimes_rate_cache_invalidate();
         snrt_chimes_fast_swap(&run.sched,&one.ctx);
         one.backend=cell_backend;one.retry=retry_cell;
#ifdef SNRT_CHIMES_RHS_TESTING
         one.calls=cell_calls;
#endif
         progressed=true;
         if(one.waiting && !one.queued){
           if(!shared_healthy.load(std::memory_order_acquire)){one.waiting->ok=0;one.waiting=nullptr;one.queued=false;}
           else publish_fiber(one);
         }
         if(worker==0 && pump(false)) progressed=true;
       }
       if(!alive){
         brokers_running.fetch_sub(1);
         if(worker==0){
           const auto drain_start=std::chrono::steady_clock::now();
           while(brokers_running.load()>0 || pending_have()>0 || shared_busy()){
             if(!pump(true)) std::this_thread::sleep_for(std::chrono::microseconds(50));
             if(std::chrono::steady_clock::now()-drain_start>std::chrono::seconds(10)){
               std::fprintf(stderr,"CHIMES shared drain timeout pending=%d brokers=%d\n",pending_have(),brokers_running.load());
               healthy=false;shared_healthy.store(0);shared_poll();fail_pending();
               break;
             }
           }
           const auto wait_consume=std::chrono::steady_clock::now();
           while(shared_busy() && std::chrono::steady_clock::now()-wait_consume<std::chrono::seconds(2)){
             release_consumed();std::this_thread::sleep_for(std::chrono::microseconds(50));
           }
         }
         break;
       }
       if(pump(!runnable_now())) progressed=true;
       if(!runnable_now() && !progressed) std::this_thread::sleep_for(std::chrono::microseconds(50));
     }
     tile_run=nullptr;
   }
   if(slot>=0){
     cudaStreamSynchronize(stream);
     for(auto event:events)if(event)cudaEventDestroy(event);
     cuda_release_stream(slot);
   }
   cpu_direct=saved_direct;cell_backend=saved_backend;retry_cell=saved_retry;
   snrt_chimes_rate_cache_invalidate();
 }
 unsigned long long total=0;for(auto count:completed)total+=count;
 std::fprintf(stderr,"CHIMES_LEVEL_QUEUE cells=%d completed=%llu workers=%zu gpu_brokers=%d active_limit_per_broker=128 batch_limit=%d buffers_per_broker=2\n",
     cells,total,completed.size(),active_brokers,nb);
 for(size_t i=0;i<completed.size();i++)std::fprintf(stderr,"CHIMES_LEVEL_WORKER worker=%zu cells=%llu\n",i,completed[i]);
 return total==(unsigned long long)cells?0:1;
}
#ifdef SNRT_CHIMES_RHS_TESTING
extern "C" void snrt_chimes_rhs_test_drop_after(int n){drop_after=n;}
extern "C" void snrt_chimes_rhs_test_submit_fail_after(int n){submit_fail_after=n;submit_attempts=0;}
extern "C" unsigned long long snrt_chimes_rhs_test_drops(void){return drops;}
extern "C" int snrt_chimes_rhs_test_retry_sticky(void){
 if(!enabled)return 0;
 const int saved_backend=cell_backend;const bool saved_retry=retry_cell;
 cell_backend=1;retry_cell=true;
 // An aborted trajectory is intentionally unavailable. No tables, stream
 // or UserData may be inspected once whole-cell CPU restart is requested.
 const int result=snrt_chimes_dark_rhs_gpu(nullptr);
 const bool sticky=result==-1 && retry_cell && cell_backend==1;
 cell_backend=saved_backend;retry_cell=saved_retry;
 return sticky;
}
#endif
extern "C" void snrt_chimes_rhs_counts(unsigned long long out[7]) {
 out[0]=used.load();out[1]=fallback.load();out[2]=batches.load();out[3]=errors.load();out[4]=uploads.load();
 out[5]=peak_active_batches.load();out[6]=(unsigned long long)__builtin_popcountll(used_stream_mask.load());
}
extern "C" void snrt_chimes_rhs_finalize(void){
 for(auto &b:buffers)b.clear();
 for(auto &b:broker_extra)b.clear();
}
extern "C" int snrt_chimes_gpu_net_cooling(double *net){
 if(!gpu_net_cooling_ready || !net || !std::isfinite(gpu_net_cooling)){gpu_net_cooling_ready=0;return 0;}
 *net=gpu_net_cooling;gpu_net_cooling_ready=0;return 1;
}
extern "C" int snrt_chimes_gpu_cooling_probe(UserData *d,double *net){
 if(!d || !net || !gpu_cooling || !build_tables() || !cool.enabled)return 0;
 cudaStream_t stream=nullptr;auto &b=buffers[0];
 if(!b.init(stream))return 0;
 Request q{};q.d=d;fill_input(q.input,d);
 std::vector<Request*> req{&q};
 b.hi[0]=q.input;
 bool ok=cudaMemcpyAsync(b.di,b.hi,sizeof(Input),cudaMemcpyHostToDevice,stream)==cudaSuccess;
 b.compact=select_compact(req);
 if(ok){reactions<<<1,256,0,stream>>>(host.t,b.v,b.r,b.entry,b.di,b.dout,b.compact?b.dp:nullptr);ok=cudaGetLastError()==cudaSuccess;}
 if(ok)ok=launch_dark_cooling(host.t,cool,b.v,b.entry,b.di,b.dout,b.compact?b.dp:nullptr,b.cr_heat,1,stream);
 if(ok)ok=copy_result_async(b,1,stream);
 if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
 consume_outputs(req,b,ok);
 if(!q.ok || !q.cooling_valid || !std::isfinite(q.net_cooling))return 0;
 *net=q.net_cooling;return 1;
}
