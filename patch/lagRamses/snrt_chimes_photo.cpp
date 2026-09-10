#include "snrt_chimes_spectrum.h"
#include "snrt_chimes_spectrum_internal.h"
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <unordered_map>
#include <stdexcept>

extern "C" int snrt_chimes_fs_ready(void);
extern "C" int snrt_chimes_fs_grid(int,double *,double *,double *);
extern "C" int snrt_chimes_budget(const double *,double *,double *);

namespace {
constexpr int ns=157,ng=9,K=128,nledger=8+2*ng;
constexpr int ne=258;
using snrt_chimes_detail::Bank;
void require(bool ok){if(!ok)throw std::runtime_error("atomic photo integration failed");}
struct Term {int node,from,to,electrons,sample;double sigma,binding,excess;};
struct MolecularTerm {int node,from,to1,to2;double sigma,yield,heat;};
struct EnergyIndex {int index;double weight;};
struct Photo {
  double nh,cdt;
  double pumping=0;
  std::vector<int> node;
  std::vector<double> initial;
  std::vector<double> dust_alpha,node_energy;
  std::vector<double> survival;
  std::vector<Term> terms;
  std::vector<MolecularTerm> molecular;
  std::vector<EnergyIndex> samples;
  std::vector<std::array<double,5>> fractions;
  std::array<double,14> xi;
  std::array<double,ne> energy;
  std::array<std::array<double,84>,ne> grid;
};
struct Solver {
  N_Vector y=nullptr,tol=nullptr,constraints=nullptr;
  SUNLinearSolver linear=nullptr;
  void *cv=nullptr;
  ~Solver(){
    if(cv)CVodeFree(&cv);
    if(linear)SUNLinSolFree(linear);
    if(y)N_VDestroy(y);
    if(tol)N_VDestroy(tol);
    if(constraints)N_VDestroy(constraints);
  }
};
// Same raw FS2010 interpolation/target limiter as the existing grey adapter.
// Samples retain six raw channels; averaging already-normalized fractions
// at the xi endpoints would introduce a different interpolation law.
void partition(Photo &p,const double *state){
  const double ion=std::max(0.,std::min(1.,state[2]));
  const double xi=std::max(p.xi[0],std::min(p.xi[13],ion));
  int ix=0;while(ix<12 && xi>p.xi[ix+1])++ix;
  const double wx=(xi-p.xi[ix])/(p.xi[ix+1]-p.xi[ix]);
  const double xr=std::max(1e-4,std::min(.999,ion)),he=.24/(4*.76);
  const double available[3]={std::min(1.,std::max(0.,state[1])/(1-xr)),
    std::min(1.,std::max(0.,state[4])/(he*(1-xr))),
    std::min(1.,std::max(0.,state[5])/(he*xr))};
  double by_energy[ne][6];
  for(int e=0;e<ne;++e)for(int j=0;j<6;++j)
    by_energy[e][j]=(1-wx)*p.grid[e][ix*6+j]+wx*p.grid[e][(ix+1)*6+j];
  for(size_t k=0;k<p.samples.size();++k){
    const auto &s=p.samples[k];auto &f=p.fractions[k];double v[6],w[3];
    for(int j=0;j<6;++j)v[j]=std::max(0.,(1-s.weight)*by_energy[s.index][j]+s.weight*by_energy[s.index+1][j]);
    w[0]=v[3]*13.6;w[1]=v[4]*24.59;w[2]=v[5]*54.42;
    const double sum=w[0]+w[1]+w[2];f.fill(0);
    if(sum>0)for(int j=0;j<3;++j)f[j+1]=v[0]*w[j]/sum;
    const double total=v[1]+v[2]+f[1]+f[2]+f[3];require(total>0 && std::isfinite(total));
    for(int j=0;j<3;++j)f[j+1]=f[j+1]/total*available[j];
    f[4]=v[2]/total*available[0];f[0]=std::max(0.,1-f[1]-f[2]-f[3]-f[4]);
  }
}
int rhs(realtype,N_Vector y,N_Vector dy,void *context){
  try{
    auto &p=*static_cast<Photo*>(context);const int nn=p.node.size(),li=ns+nn;
    const double *v=N_VGetArrayPointer(y);double *out=N_VGetArrayPointer(dy);
    for(int j=0;j<li+nledger;++j)if(!std::isfinite(v[j]))return 1;
    std::fill(out,out+li+nledger,0.);partition(p,v);
    // Integrate accumulated optical depth, not the exponentially exhausted
    // fraction. d(tau)/ds = c*dt*opacity; survival=exp(-tau) cannot cross zero.
    // max is only for Newton trial tau; accepted negative states still reject.
    for(int k=0;k<nn;++k)p.survival[k]=std::exp(-std::max(0.,v[ns+k]));
    const double all_heat[5]={1,0,0,0,0};
    for(const auto &t:p.terms){
      // Positive trial rates for Newton iterations, not published clipping.
      const double opacity=p.cdt*t.sigma*std::max(0.,v[t.from]);
      const double factor=opacity*p.survival[t.node];
      const double rate=factor*p.initial[t.node]; // events per H over t/dt
      out[ns+t.node]+=opacity*p.nh;
      out[t.from]-=rate;out[t.to]+=rate;out[0]+=rate*t.electrons;
      const double *f=t.sample<0 ? all_heat : p.fractions[t.sample].data();
      const double power=rate*t.excess;
      for(int j=0;j<5;++j)out[li+j]+=power*f[j];
      const int from[3]={1,4,5},to[3]={2,5,6};
      const double threshold[3]={13.6,24.59,54.42};
      for(int j=0;j<3;++j){
        const double secondary=power*f[j+1]/threshold[j];
        out[from[j]]-=secondary;out[to[j]]+=secondary;out[0]+=secondary;
      }
      out[li+5]+=rate*t.binding;
      out[li+6]+=rate*(t.binding+t.excess);
      out[li+7]+=rate;
    }
    // Grain inventories are frozen within this photo step only and compete
    // with gas in the SAME accumulated optical depth and capture ledger.
    // alpha already includes grain mass density (cm^-1), not a cross section
    // per H. Normalize grain counters per H like the gas counters here.
    for(int k=0;k<nn;++k){
      const double loss=p.cdt*p.dust_alpha[k]*p.survival[k];
      const double rate=loss*p.initial[k]/p.nh;
      out[ns+k]+=p.cdt*p.dust_alpha[k];
      const int g=p.node[k]/K;
      out[li+8+g]+=rate*p.node_energy[k];out[li+8+ng+g]+=rate;
    }
    for(const auto &t:p.molecular){
      const double opacity=p.cdt*t.sigma*std::max(0.,v[t.from]);
      const double factor=opacity*p.survival[t.node];
      const double captures=factor*p.initial[t.node],events=captures*t.yield;
      out[ns+t.node]+=opacity*p.nh;
      out[t.from]-=events;out[t.to1]+=events;out[t.to2]+=events;
      const double power=captures*p.node_energy[t.node];
      // Native Draine-field pumping estimate, bounded by counted fluorescent
      // captures at 2 eV/pump. Density quenching is supplied by native CHIMES
      // critical-density tables, not invented in this spectral operator.
      const double pump=t.heat>0?std::min(events*2.7e-11/1.602176634e-12,
          captures*(1-t.yield)*2.)*p.pumping:0.;
      const double heat=events*t.heat+pump;
      out[li]+=heat;out[li+5]+=power-heat;
      out[li+6]+=power;out[li+7]+=captures;
    }
    for(int j=0;j<li+nledger;++j)if(!std::isfinite(out[j]))return 1;
    return 0;
  }catch(...){return -1;}
}
}

static int photo_step(void *handle,int nd,double nh,double dt,double chat,const double *dust_alpha,
    const snrt_chimes_detail::MoleculeBank *molecules,const double *shield,double pumping,
    const double *old,const double *number,const double *energy,double *next,
    double *next_number,double *next_energy,double *ledger,
    double *grain_number=nullptr,double *grain_energy=nullptr){
  if(!handle || nd<1 || nd>720 || !old || !number || !energy || !next || !next_number || !next_energy || !ledger)return 1;
  if(!std::isfinite(nh) || nh<=0 || !std::isfinite(dt) || dt<0 || !std::isfinite(chat) || chat<=0 || chat>2.99792458e10)return 2;
  const auto &b=*static_cast<const Bank*>(handle);
  double nuclei[11],charge;
  if(snrt_chimes_budget(old,nuclei,&charge) || std::abs(nuclei[0]-1)>1e-10)return 2;
  // A fixed external solid charge is allowed: preserve the incoming gas
  // charge, do not reinterpret it as a new electron reservoir.
  try{
    Photo p;p.nh=nh;p.cdt=chat*dt;p.pumping=pumping;require(std::isfinite(p.cdt));
    if(!std::isfinite(pumping) || pumping<0 || pumping>1)return 2;
    if(dust_alpha)for(int j=0;j<ng*K;++j)
      if(!std::isfinite(dust_alpha[j]) || dust_alpha[j]<0)return 2;
    if(molecules && (!shield || !std::isfinite(shield[0]) || !std::isfinite(shield[1]) ||
        shield[0]<0 || shield[0]>1 || shield[1]<0 || shield[1]>1))return 2;
    std::vector<double> rays(size_t(ng)*nd*K,0),initial(ng*K,0);
    double initial_energy=0;
    for(int g=0;g<ng;++g)for(int d=0;d<nd;++d){
      const size_t i=size_t(g)*nd+d;const double n=number[i],e=energy[i];
      if(!std::isfinite(n) || !std::isfinite(e) || n<0 || e<0 || (n==0 && e!=0))return 2;
      initial_energy+=e;
      if(n==0)continue;
      const double mean=e/n;std::array<double,K> weights;
      // Only FP64-sized endpoint roundoff, not the transport's FP32 margin.
      // N/E themselves are never projected; debit the actual absorbed nodes.
      if(mean<b.grids[g].e.front()*(1-2e-13) || mean>b.grids[g].e.back()*(1+2e-13) ||
         !b.grids[g].reconstruct(mean,weights))return 2;
      for(int j=0;j<K;++j){const double q=n*weights[j];rays[i*K+j]=q;initial[g*K+j]+=q;}
    }
    require(std::isfinite(initial_energy));
    auto identity=[&](){
      std::copy(old,old+ns,next);std::copy(number,number+nd*ng,next_number);
      std::copy(energy,energy+nd*ng,next_energy);std::fill(ledger,ledger+10,0.);
      if(grain_number)std::fill(grain_number,grain_number+ng,0.);
      if(grain_energy)std::fill(grain_energy,grain_energy+ng,0.);
    };
    if(dt==0 || initial_energy==0){identity();return 0;}
    require(snrt_chimes_fs_ready());
    std::unordered_map<double,int> sample_index;
    require(snrt_chimes_fs_grid(ne,p.energy.data(),p.xi.data(),p.grid[0].data())==0);
    for(int j=1;j<14;++j)require(p.xi[j]>p.xi[j-1]);
    for(int j=1;j<ne;++j)require(p.energy[j]>p.energy[j-1]);
    for(const auto &row:p.grid)for(double v:row)require(std::isfinite(v) && v>=0);
    for(int node=0;node<ng*K;++node){
      if(initial[node]==0)continue;
      const int active=p.node.size();bool present=false;
      for(size_t shell=0;shell<b.shell.size();++shell){
        const double sigma=b.sigma[shell*ng*K+node];if(sigma==0)continue;
        const int r=b.shell[shell];const int *map=b.reaction.data()+r*5;
        const double excess=std::max(0.,b.energy[node]-b.binding[shell]);
        int sample=-1;
        if(excess>=10){
          auto found=sample_index.find(excess);
          if(found!=sample_index.end())sample=found->second;
          else{
            sample=p.samples.size();
            const double ev=std::max(p.energy.front(),std::min(p.energy.back(),excess));
            const int index=std::min(ne-2,int(std::upper_bound(p.energy.begin(),p.energy.end(),ev)-p.energy.begin())-1);
            p.samples.push_back({index,(ev-p.energy[index])/(p.energy[index+1]-p.energy[index])});
            sample_index.emplace(excess,sample);
          }
        }
        p.terms.push_back({active,map[2],map[3],map[4],sample,sigma,b.binding[shell],excess});present=true;
      }
      const double alpha=dust_alpha?dust_alpha[node]:0.;
      if(molecules)for(int r=0;r<32;++r){
        const int offset=r*ng*K+node;
        const double ab=molecules->absorption[offset],di=molecules->dissociation[offset];
        const double attenuation=r<30?1.:shield[r-30];
        if(ab==0 || attenuation==0)continue;
        const int *map=molecules->reaction.data()+3*r;
        const double heat=r==30?6.4e-13/1.602176634e-12:0.;
        require(heat<=b.energy[node]);
        p.molecular.push_back({active,map[0],map[1],map[2],ab*attenuation,di/ab,heat});
        present=true;
      }
      if(present || alpha>0){
        p.node.push_back(node);p.initial.push_back(initial[node]);
        p.dust_alpha.push_back(alpha);p.node_energy.push_back(b.energy[node]);
      }
    }
    if(p.node.empty()){identity();return 0;}
    p.fractions.resize(p.samples.size());
    const int nn=p.node.size(),li=ns+nn,size=li+nledger;
    p.survival.resize(nn);
    Solver solver;solver.y=N_VNew_Serial(size);solver.tol=N_VNew_Serial(size);solver.constraints=N_VNew_Serial(size);
    require(solver.y && solver.tol && solver.constraints);
    double *y=N_VGetArrayPointer(solver.y),*tol=N_VGetArrayPointer(solver.tol);
    std::copy(old,old+ns,y);std::fill(y+ns,y+size,0.);
    std::fill(tol,tol+ns,1e-16);std::fill(tol+ns,tol+li,1e-14);
    const double scale=initial_energy/nh;
    std::fill(tol+li,tol+size,std::max(1e-30,scale*1e-14));
    tol[li+7]=std::max(1e-30,scale/b.grids.back().e.back()*1e-14);
    std::fill(tol+li+8+ng,tol+size,tol[li+7]);
    N_VConst(1.,solver.constraints); // nonnegative accepted species/tau/counters
    solver.cv=CVodeCreate(CV_BDF);require(solver.cv);
    require(CVodeInit(solver.cv,rhs,0.,solver.y)==0 && CVodeSetUserData(solver.cv,&p)==0);
    require(CVodeSVtolerances(solver.cv,1e-11,solver.tol)==0 && CVodeSetConstraints(solver.cv,solver.constraints)==0);
    require(CVodeSetMaxNumSteps(solver.cv,10000)==0);
    // Keep the requested endpoint; never interpolate backwards from a
    // later step with a different gas/grain capture budget.
    require(CVodeSetStopTime(solver.cv,1.)==0);
    solver.linear=SUNLinSol_SPGMR(solver.y,PREC_NONE,30);require(solver.linear);
    require(CVodeSetLinearSolver(solver.cv,solver.linear,nullptr)==0);
    double reached=0;
    const int result=CVode(solver.cv,1.,solver.y,&reached,CV_NORMAL);
    if((result!=CV_SUCCESS && result!=CV_TSTOP_RETURN) || reached!=1.)return 4;
    for(int j=0;j<size;++j)if(!std::isfinite(y[j]) || y[j]<0)return 51;
    double new_nuclei[11],new_charge;
    if(snrt_chimes_budget(y,new_nuclei,&new_charge))return 5;
    for(int e=0;e<11;++e)if(std::abs(new_nuclei[e]-nuclei[e])>1e-8*std::max(nuclei[e],1e-20))return 6;
    if(std::abs(new_charge-charge)>1e-8)return 6;
    std::vector<double> survival(ng*K,1),out_n(number,number+ng*nd),out_e(energy,energy+ng*nd);
    for(int k=0;k<nn;++k)survival[p.node[k]]=std::exp(-y[ns+k]);
    double absorbed_n=0,absorbed_e=0;
    for(int g=0;g<ng;++g)for(int d=0;d<nd;++d){
      const int i=g*nd+d;double total_n=0,total_e=0,survive_n=0,survive_e=0;
      bool changed=false;
      for(int j=0;j<K;++j){
        const double q=rays[i*K+j],s=survival[g*K+j],ev=b.energy[g*K+j];
        total_n+=q;total_e+=q*ev;survive_n+=q*s;survive_e+=q*s*ev;
        changed=changed || (q>0 && s<1);
      }
      // Recover the surviving population directly. Subtracting a rounded
      // almost-total capture from N can produce negative exhausted rays.
      // Ratios retain exact incoming N/E normalization of the reconstruction;
      // this is not clipping, and zero-opacity rays remain bitwise unchanged.
      if(changed){out_n[i]=number[i]*(survive_n/total_n);out_e[i]=energy[i]*(survive_e/total_e);}
      absorbed_n+=number[i]-out_n[i];absorbed_e+=energy[i]-out_e[i];
      if(!std::isfinite(out_n[i]) || !std::isfinite(out_e[i]) || out_n[i]<0 || out_e[i]<0)return 53;
    }
    std::array<double,10> budget{};
    std::array<double,ng> dust_n{},dust_e{};
    for(int j=0;j<8;++j)budget[j]=y[li+j]*nh;
    for(int g=0;g<ng;++g){
      dust_e[g]=y[li+8+g]*nh;dust_n[g]=y[li+8+ng+g]*nh;
      budget[8]+=dust_e[g];budget[9]+=dust_n[g];
      // Each grain receiver gets its actual band energy, never a fixed
      // representative energy times N. Reject counter drift before publish.
      const double lo=b.grids[g].e.front(),hi=b.grids[g].e.back();
      if(dust_e[g]<lo*dust_n[g]-1e-8*std::max(initial_energy,1e-30) ||
         dust_e[g]>hi*dust_n[g]+1e-8*std::max(initial_energy,1e-30))return 7;
    }
    const double e_scale=std::max(initial_energy,1e-30);
    double n_scale=0;for(int i=0;i<ng*nd;++i)n_scale+=number[i];
    if(std::abs(absorbed_n-budget[7]-budget[9])>1e-7*std::max(n_scale,1e-30) ||
       std::abs(absorbed_e-budget[6]-budget[8])>1e-7*e_scale)return 7;
    double accounted=0;for(int j=0;j<6;++j)accounted+=budget[j];
    if(std::abs(accounted-budget[6])>1e-8*e_scale)return 7;
    std::copy(y,y+ns,next);std::copy(out_n.begin(),out_n.end(),next_number);
    std::copy(out_e.begin(),out_e.end(),next_energy);std::copy(budget.begin(),budget.end(),ledger);
    if(grain_number)std::copy(dust_n.begin(),dust_n.end(),grain_number);
    if(grain_energy)std::copy(dust_e.begin(),dust_e.end(),grain_energy);
    return 0;
  }catch(...){return 3;}
}

extern "C" int snrt_chimes_band_photo_step(void *handle,int nd,double nh,double dt,double chat,
    const double *old,const double *number,const double *energy,double *next,
    double *next_number,double *next_energy,double *ledger){
  if(!ledger)return 1;
  double budget[10];
  const int status=photo_step(handle,nd,nh,dt,chat,nullptr,nullptr,nullptr,0,old,number,energy,next,next_number,next_energy,budget);
  if(status==0)std::copy(budget,budget+8,ledger);
  return status;
}

extern "C" int snrt_chimes_band_photo_dust_step(void *handle,int nd,double nh,double dt,double chat,
    const double *dust_alpha,const double *old,const double *number,const double *energy,double *next,
    double *next_number,double *next_energy,double *ledger){
  if(!dust_alpha)return 1;
  return photo_step(handle,nd,nh,dt,chat,dust_alpha,nullptr,nullptr,0,old,number,energy,next,next_number,next_energy,ledger);
}

extern "C" int snrt_chimes_band_photo_molecular_step(void *handle,void *molecules,int nd,
    double nh,double dt,double chat,const double *dust_alpha,const double *shield,double pumping,
    const double *old,const double *number,const double *energy,double *next,
    double *next_number,double *next_energy,double *ledger){
  if(!molecules || !dust_alpha)return 1;
  return photo_step(handle,nd,nh,dt,chat,dust_alpha,
      static_cast<const snrt_chimes_detail::MoleculeBank*>(molecules),shield,pumping,
      old,number,energy,next,next_number,next_energy,ledger);
}

extern "C" int snrt_chimes_band_photo_molecular_groups(void *handle,void *molecules,int nd,
    double nh,double dt,double chat,const double *dust_alpha,const double *shield,double pumping,
    const double *old,const double *number,const double *energy,double *next,
    double *next_number,double *next_energy,double *ledger,double *grain_number,double *grain_energy){
  if(!molecules || !dust_alpha || !grain_number || !grain_energy)return 1;
  return photo_step(handle,nd,nh,dt,chat,dust_alpha,
      static_cast<const snrt_chimes_detail::MoleculeBank*>(molecules),shield,pumping,
      old,number,energy,next,next_number,next_energy,ledger,grain_number,grain_energy);
}
