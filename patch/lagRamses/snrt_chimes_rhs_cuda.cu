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
static_assert(sizeof(ChimesFloat)==8 && sizeof(sunindextype)==8,
              "Device reaction slice requires FP64 / 64-bit SUNDIALS indices");
namespace {
constexpr int ns=157, nb=64, nrmax=1024, groups=8;
enum {TD,CON,REC,GRAIN,CR,H2D,H2C,CO};
enum {AXT,AXD,AXX,TDT,CONT,RECT,GRAINT,CRT,CRSEC,H2DT,KZERO,KLTE,CRITH,CRITH2,CRITHE,COT,NTAB};
struct Reaction {int kind,local,a[3];};
struct Tables {
 int nt,nd,nx,np,nr,count[groups][2],base[groups],offset[NTAB],secondary[2];
 int begin[2*ns+1]; // CSR: destruction first, creation second; duplicate nuclei retained
};
struct Input {double x[ns],t,td,nh,dust,boost,cr;int mol,ab[2];};
struct Output {double rate[nrmax],coefficient[nrmax],crit[3],destroy[ns],create[ns];};
struct HostTables {Tables t{};std::vector<double> v;std::vector<Reaction> r;std::vector<int> entries;bool valid=false;};
HostTables host;
bool enabled=false;
bool verify=false;
thread_local int cell_backend=-2; // -2 outside bridge, -1 undecided, 0 CPU, 1 GPU
thread_local bool retry_cell=false;
#ifdef SNRT_CHIMES_RHS_TESTING
int drop_after=0;thread_local int cell_calls=0;
std::atomic<unsigned long long> drops{0};
#endif
std::mutex verification_mutex;
double verify_rates=0,verify_species=0,verify_cooling=0;
std::atomic<unsigned long long> used{0},fallback{0},batches{0},errors{0},uploads{0};

bool build_tables() {
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
 h.valid=true;return true;
}
__device__ void index(const double *a,int n,double x,int &i,double &w){
 double den=(a[n-1]-a[0])/(n-1.0);
 if(x<=a[0]){i=0;w=0;}else if(x>=a[n-1]){i=n-2;w=1;}
 else{i=(int)floor((x-a[0])/den);w=(x-a[i])/den;}
}
__device__ double interp(const double *v,int i,double w){return (1.-w)*v[i]+w*v[i+1];}
__device__ double powten(double x){return exp(x*2.30258509299404568402);}
__global__ void reactions(Tables t,const double *v,const Reaction *r,const int *entry,const Input *in,Output *out){
 const Input &x=in[blockIdx.x];Output &y=out[blockIdx.x];
 __shared__ int ti,di,xi;__shared__ double tw,dw,xw,ncrit;
 if(threadIdx.x==0){
   index(v+t.offset[AXT],t.nt,log10(x.t),ti,tw);
   index(v+t.offset[AXD],t.nd,log10(x.td),di,dw);
   index(v+t.offset[AXX],t.nx,log10(fmax(x.x[sp_HII],DBL_MIN)),xi,xw);
   y.crit[0]=powten(interp(v+t.offset[CRITH],ti,tw));
   y.crit[1]=powten(interp(v+t.offset[CRITH2],ti,tw));
   y.crit[2]=powten(interp(v+t.offset[CRITHE],ti,tw));
   ncrit=x.x[sp_HI]/y.crit[0];ncrit+=2.*x.x[sp_H2]/y.crit[1];
   ncrit+=x.x[sp_HeI]/y.crit[2];ncrit*=x.nh;
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
 }
 __syncthreads();
 for(int s=threadIdx.x;s<ns;s+=blockDim.x){
   double d=0,c=0;
   for(int j=t.begin[s];j<t.begin[s+1];j++)d+=y.rate[entry[j]];
   for(int j=t.begin[ns+s];j<t.begin[ns+s+1];j++)c+=y.rate[entry[j]];
   y.destroy[s]=d;y.create[s]=c;
 }
}
struct Buffers {
 int device=-1;bool ready=false;
 double *v=nullptr;Reaction *r=nullptr;int *entry=nullptr;
 Input *hi=nullptr,*di=nullptr;Output *ho=nullptr,*dout=nullptr;
 void clear(){if(device>=0)cudaSetDevice(device);cudaFree(v);cudaFree(r);cudaFree(entry);cudaFree(di);cudaFree(dout);cudaFreeHost(hi);cudaFreeHost(ho);*this=Buffers{};}
 bool init(cudaStream_t stream){
   if(ready)return true;
   if(cudaGetDevice(&device)!=cudaSuccess)return false;
   auto &h=host;
   bool ok=cudaMalloc(&v,h.v.size()*sizeof(double))==cudaSuccess && cudaMalloc(&r,h.r.size()*sizeof(Reaction))==cudaSuccess &&
     cudaMalloc(&entry,h.entries.size()*sizeof(int))==cudaSuccess && cudaMallocHost(&hi,nb*sizeof(Input))==cudaSuccess &&
     cudaMallocHost(&ho,nb*sizeof(Output))==cudaSuccess && cudaMalloc(&di,nb*sizeof(Input))==cudaSuccess && cudaMalloc(&dout,nb*sizeof(Output))==cudaSuccess;
   if(ok)ok=cudaMemcpyAsync(v,h.v.data(),h.v.size()*sizeof(double),cudaMemcpyHostToDevice,stream)==cudaSuccess &&
     cudaMemcpyAsync(r,h.r.data(),h.r.size()*sizeof(Reaction),cudaMemcpyHostToDevice,stream)==cudaSuccess &&
     cudaMemcpyAsync(entry,h.entries.data(),h.entries.size()*sizeof(int),cudaMemcpyHostToDevice,stream)==cudaSuccess;
   if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
   if(!ok){clear();return false;}ready=true;uploads++;return true;
 }
};
Buffers buffers[MAX_CUDA_STREAMS];
struct Request {UserData *d;Input input;bool done=false;int ok=0;bool required=false;};
std::mutex mutex;std::condition_variable changed;std::vector<Request*> pending;bool working=false;
void execute(std::vector<Request*> &requests){
 // Keep GPU arithmetic through the tail of an already active GPU cell.
 // Switching within a finite-difference Jacobian amplifies roundoff.
 if(requests.size()<4 && std::none_of(requests.begin(),requests.end(),[](Request *r){return r->required;}))return;
 int slot=cuda_acquire_stream();if(slot<0)return;
 auto stream=cuda_get_stream_internal(slot);auto &b=buffers[slot];
 bool ok=build_tables() && b.init(stream);int n=(int)requests.size();
 if(ok){
   for(int i=0;i<n;i++)b.hi[i]=requests[i]->input;
   ok=cudaMemcpyAsync(b.di,b.hi,n*sizeof(Input),cudaMemcpyHostToDevice,stream)==cudaSuccess;
   if(ok){reactions<<<n,256,0,stream>>>(host.t,b.v,b.r,b.entry,b.di,b.dout);ok=cudaGetLastError()==cudaSuccess;}
   if(ok)ok=cudaMemcpyAsync(b.ho,b.dout,n*sizeof(Output),cudaMemcpyDeviceToHost,stream)==cudaSuccess;
 }
 if(cudaStreamSynchronize(stream)!=cudaSuccess)ok=false;
 if(ok){
   batches++;
   for(int i=0;i<n;i++){
     auto &o=b.ho[i];auto &d=*requests[i]->d;auto &t=host.t;auto *c=d.chimes_current_rates;
     bool valid=true;for(int k=0;k<t.nr;k++)if(!std::isfinite(o.rate[k])||!std::isfinite(o.coefficient[k]))valid=false;
     for(double v:o.crit)if(!std::isfinite(v))valid=false;
     for(int s=0;s<ns;s++)if(!std::isfinite(o.destroy[s])||!std::isfinite(o.create[s]))valid=false;
     if(!valid)continue; // original CPU RHS decides recoverable state status
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
     requests[i]->ok=1;used++;
   }
 }else{errors+=n;b.clear();}
 cuda_release_stream(slot);
}
}
extern "C" int snrt_chimes_rhs_configure(void){
 verify=std::getenv("SNRT_CHIMES_RHS_VERIFY")!=nullptr;
 const char *s=std::getenv("SNRT_CHIMES_RHS_BACKEND");
 if(!s||!*s||!std::strcmp(s,"cpu")){enabled=false;return 0;}
 if(std::strcmp(s,"cuda_dark_reactions")){std::fprintf(stderr,"Unknown SNRT_CHIMES_RHS_BACKEND=%s\n",s);return 1;}
 enabled=true;return 0;
}
extern "C" int snrt_chimes_dark_rhs_gpu(UserData *d){
 if(!enabled)return 0;
 if(cell_backend<=-2 || cell_backend==0){fallback++;return 0;}
#ifdef SNRT_CHIMES_RHS_TESTING
 if(cell_backend==1 && drop_after && ++cell_calls==drop_after){retry_cell=true;drops++;return -1;}
#endif
 auto decline=[&](){if(cell_backend==1){retry_cell=true;return -1;}cell_backend=0;fallback++;return 0;};
 auto *c=d->myGlobalVars;auto *g=d->myGasVars;
 int team=omp_in_parallel()?omp_get_num_threads():1;
 if(c->N_spectra!=0 || c->totalNumberOfSpecies!=ns || g->ThermEvolOn!=1 || team<4 || (d->mol_flag_index!=0 && d->mol_flag_index!=1))return decline();
 for(int s=0;s<ns;s++)if(c->speciesIndices[s]!=s)return decline();
 if(d->case_AB_index[0]<0||d->case_AB_index[0]>1||d->case_AB_index[1]<0||d->case_AB_index[1]>1)return decline();
 Request q{};q.d=d;q.required=cell_backend==1;auto &in=q.input;
 std::memcpy(in.x,g->abundances,sizeof(in.x));in.t=g->temperature;in.td=c->grain_temperature;
 in.nh=g->nH_tot;in.dust=g->dust_ratio;in.boost=g->dust_boost_factor;in.cr=g->cr_rate;in.mol=d->mol_flag_index;
 in.ab[0]=d->case_AB_index[0];in.ab[1]=d->case_AB_index[1];
 if(!(in.t>0) || !(in.td>0) || !(in.nh>0) || !std::isfinite(in.t))return decline();
 std::unique_lock<std::mutex> lock(mutex);pending.push_back(&q);changed.notify_all();
 while(!q.done){
   if(working){changed.wait(lock,[&]{return q.done||!working;});continue;}
   working=true;
   changed.wait_until(lock,std::chrono::steady_clock::now()+std::chrono::microseconds(50),[&]{return pending.size()>=(size_t)std::min(team,nb);});
   std::vector<Request*> batch;
   size_t n=std::min(pending.size(),(size_t)nb);batch.assign(pending.begin(),pending.begin()+n);pending.erase(pending.begin(),pending.begin()+n);
   lock.unlock();execute(batch);lock.lock();
   for(auto r:batch)r->done=true;
   working=false;changed.notify_all();
 }
 lock.unlock();
 if(!q.ok){
   if(cell_backend==1){retry_cell=true;return -1;}
   cell_backend=0;fallback++;
 }else cell_backend=1;
 if(q.ok && c->rt_update_flux){
   // Preserve host extension callbacks on their owning thread. In dark
   // mode the current receiver resets photo coefficients/adds zero sources.
   if(snrt_chimes_molecular_coefficients)snrt_chimes_molecular_coefficients(g,c,*d);
   if(snrt_chimes_secondary_rates)snrt_chimes_secondary_rates(d->species,g,c);
 }
 if(q.ok && verify){
   // One-off paired evaluation on the identical CVODE trial state. Return
   // the original CPU result to the integrator, so its trajectory is a
   // reference for examining device arithmetic independently of rejection.
   double rates[nrmax]={0},create[ns],destroy[ns];auto &t=host.t;
   auto *r=d->chimes_current_rates;
#define SAVE(field,k) std::memcpy(rates+t.base[k],r->field,t.count[k][d->mol_flag_index]*sizeof(double))
   SAVE(T_dependent_rate,TD);SAVE(constant_rate,CON);SAVE(recombination_AB_rate,REC);
   SAVE(grain_recombination_rate,GRAIN);SAVE(cosmic_ray_rate,CR);
   if(d->mol_flag_index){rates[t.base[H2D]]=r->H2_dust_formation_rate;SAVE(H2_collis_dissoc_rate,H2C);SAVE(CO_cosmic_ray_rate,CO);}
#undef SAVE
   for(int s=0;s<ns;s++){create[s]=d->species[s].creation_rate;destroy[s]=d->species[s].destruction_rate;}
   double cool=calculate_total_cooling_rate(g,c,*d,0);
   // Cache invalidation is a C function in the patched private CHIMES.
   snrt_chimes_rate_cache_invalidate();
   update_rate_coefficients(g,c,*d,g->ThermEvolOn);
   if(c->rt_update_flux)update_rt_photochemistry_coefficients(g,c,*d);
   update_rates(g,c,*d);update_rate_vector(d->species,g,c,*d);
   double dr=0,ds=0,dc=0;
   auto diff=[&](double a,double b){return std::abs(a-b)/std::max(std::abs(b),1e-60);};
#define CHECK(field,k) for(int j=0;j<t.count[k][d->mol_flag_index];j++)dr=std::max(dr,diff(rates[t.base[k]+j],r->field[j]))
   CHECK(T_dependent_rate,TD);CHECK(constant_rate,CON);CHECK(recombination_AB_rate,REC);
   CHECK(grain_recombination_rate,GRAIN);CHECK(cosmic_ray_rate,CR);
   if(d->mol_flag_index){dr=std::max(dr,diff(rates[t.base[H2D]],r->H2_dust_formation_rate));CHECK(H2_collis_dissoc_rate,H2C);CHECK(CO_cosmic_ray_rate,CO);}
#undef CHECK
   for(int s=0;s<ns;s++){
     ds=std::max(ds,diff(create[s],d->species[s].creation_rate));
     ds=std::max(ds,diff(destroy[s],d->species[s].destruction_rate));
   }
   dc=diff(cool,calculate_total_cooling_rate(g,c,*d,0));
   std::lock_guard<std::mutex> guard(verification_mutex);
   verify_rates=std::max(verify_rates,dr);verify_species=std::max(verify_species,ds);verify_cooling=std::max(verify_cooling,dc);
 }
 return q.ok;
}
extern "C" void snrt_chimes_rhs_verification(double out[3]){out[0]=verify_rates;out[1]=verify_species;out[2]=verify_cooling;}
extern "C" void snrt_chimes_rhs_cell_begin(void){cell_backend=-1;retry_cell=false;
#ifdef SNRT_CHIMES_RHS_TESTING
 cell_calls=0;
#endif
}
extern "C" int snrt_chimes_rhs_cell_retry(void){int retry=retry_cell;retry_cell=false;cell_backend=0;return retry;}
extern "C" void snrt_chimes_rhs_cell_end(void){cell_backend=-2;retry_cell=false;}
#ifdef SNRT_CHIMES_RHS_TESTING
extern "C" void snrt_chimes_rhs_test_drop_after(int n){drop_after=n;}
extern "C" unsigned long long snrt_chimes_rhs_test_drops(void){return drops;}
#endif
extern "C" void snrt_chimes_rhs_counts(unsigned long long out[5]) {out[0]=used;out[1]=fallback;out[2]=batches;out[3]=errors;out[4]=uploads;}
extern "C" void snrt_chimes_rhs_finalize(void){for(auto &b:buffers)b.clear();}
