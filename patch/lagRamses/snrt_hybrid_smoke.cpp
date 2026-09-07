#include "snrt_hybrid.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
extern "C" int snrt_cuda_available_c();
extern "C" int snrt_dust_material_openmp_c(const double*,const double*,double*,int,int,int,int,double,double,double,double,int);
extern "C" int snrt_hybrid_dust_material_c(const double*,const double*,double*,int,int,int,int,double,double,double,double,int);

int main() {
  constexpr int n=1031,nw=n+3,nd=8,ng=9,g=n*ng,total=nw*nd*ng,batch=64;
  const bool gpu=snrt_cuda_available_c()>0;
  if(snrt_hybrid_configure_c(0,1,batch,4,1,gpu))return 1;
  double worst=0;
  for(int dusty=0;dusty<2;++dusty) {
    std::vector<float> input(total),dir(3*nd),tau(g),stau(3*g),dtau(g),budget(3*n);
    std::vector<int> neighbor(6*n);
    for(int i=0;i<total;++i)input[i]=0.001f*(1+(i*37)%51);
    for(int i=0;i<nd;++i)for(int d=0;d<3;++d)dir[3*i+d]=((i>>d)&1)?0.577350269f:-0.577350269f;
    for(int i=0;i<6*n;++i)neighbor[i]=i%11?1+(i*7)%nw:0;
    for(int i=0;i<3*n;++i)budget[i]=0.01f*(1+i%13);
    for(int j=0;j<ng;++j)for(int i=0;i<n;++i) {
      for(int s=0;s<3;++s)stau[s*g+j*n+i]=j>=4+s?0.03f*(1+(s+j+i)%17):0;
      dtau[j*n+i]=dusty?0.021f*(1+(j+i)%7):0;
      tau[j*n+i]=stau[j*n+i]+stau[g+j*n+i]+stau[2*g+j*n+i]+dtau[j*n+i];
    }
    std::vector<float> ref=input,atoms=budget,hhe(3*g),dust(g),ret(g),raw(g),abs_g(g),absorbed(n);
    if(snrt_openmp_species_dust_c(ref.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
        atoms.data(),hhe.data(),dust.data(),ret.data(),raw.data(),abs_g.data(),absorbed.data(),n,nw,nd,ng,.2f))return 2;
    for(int busy=0;busy<2;++busy) {
      const int held=(busy&&gpu)?cuda_acquire_stream():-1;
      if(busy&&gpu&&held<0)return 3;
      auto qs=input,at=budget;std::vector<float> hh(3*g),dd(g),rr(g),ww(g),aa(g),tt(n);
      int rc=snrt_hybrid_species_dust_c(qs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          at.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f);
      if(held>=0)cuda_release_stream(held);
      if(rc)return 4;
      int cpu_count,gpu_count;snrt_hybrid_counts_c(0,&cpu_count,&gpu_count);
      if(cpu_count+gpu_count!=1+(n-1)/batch)return 5;
      if((busy||!gpu)&&gpu_count)return 6;
      if(!busy&&gpu&&(cpu_count==0||gpu_count==0))return 7;
      auto compare=[&](const std::vector<float>& a,const std::vector<float>& b) {
        double scale=.5,error=0;
        for(size_t k=0;k<a.size();++k){scale=std::max(scale,double(std::abs(a[k])));error=std::max(error,double(std::abs(a[k]-b[k])));}
        worst=std::max(worst,error/scale);
        return error<3e-5*scale;
      };
      if(!compare(ref,qs)||!compare(atoms,at)||!compare(hhe,hh)||!compare(dust,dd)||!compare(ret,rr)||
          !compare(raw,ww)||!compare(abs_g,aa)||!compare(absorbed,tt))return 8;
      if((busy||!gpu)&&(ref!=qs||atoms!=at||hhe!=hh||dust!=dd||ret!=rr||raw!=ww||abs_g!=aa||absorbed!=tt))return 9;
      std::printf("PRIMARY hybrid dusty=%d held=%d CPU=%d GPU=%d PASS\n",dusty,busy,cpu_count,gpu_count);
      // Late-batch invalid opacity: successful sibling batches must not publish.
      const auto saved=qs,saved_at=at,saved_hh=hh,saved_dd=dd,saved_rr=rr,saved_ww=ww,saved_aa=aa,saved_tt=tt;
      const float old=tau[g-1];tau[g-1]=std::numeric_limits<float>::quiet_NaN();
      rc=snrt_hybrid_species_dust_c(qs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          at.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f);
      tau[g-1]=old;
      if(!rc||qs!=saved||at!=saved_at||hh!=saved_hh||dd!=saved_dd||rr!=saved_rr||ww!=saved_ww||aa!=saved_aa||tt!=saved_tt)return 10;
    }
  }
  constexpr int dc=1031,dg=2,nt=4;
  std::vector<double> table((dg+3)*nt),in(4*dc),ref((dg+2)*dc);
  const double temp[nt]={10,20,50,100},power[nt]={1,2,5,10};
  for(int k=0;k<nt;++k) {
    table[k]=std::log(temp[k]);table[nt+k]=power[k];table[2*nt+k]=temp[k]*temp[k];
    table[3*nt+k*dg]=.4*power[k];table[3*nt+k*dg+1]=.6*power[k];
  }
  for(int use_u=0;use_u<2;++use_u) {
    for(int i=0;i<dc;++i) {in[i]=.3;in[dc+i]=1;in[2*dc+i]=use_u?400:20;in[3*dc+i]=1;}
    if(snrt_dust_material_openmp_c(in.data(),table.data(),ref.data(),dc,dg,nt,use_u,.1,1,10,1e-9,4))return 11;
    for(int busy=0;busy<2;++busy) {
      const int held=(busy&&gpu)?cuda_acquire_stream():-1;
      if(busy&&gpu&&held<0)return 12;
      std::vector<double> trial(ref.size(),-1);
      int rc=snrt_hybrid_dust_material_c(in.data(),table.data(),trial.data(),dc,dg,nt,use_u,.1,1,10,1e-9,4);
      if(held>=0)cuda_release_stream(held);
      if(rc)return 13;
      int cpu_count,gpu_count;snrt_hybrid_counts_c(1,&cpu_count,&gpu_count);
      if(cpu_count+gpu_count!=1+(dc-1)/batch)return 14;
      if((busy||!gpu)&&gpu_count)return 15;
      if(!busy&&gpu&&(cpu_count==0||gpu_count==0))return 16;
      for(size_t i=0;i<ref.size();++i)if(std::abs(ref[i]-trial[i])>1e-12*std::max(std::abs(ref[i]),1.))return 17;
      std::printf("DUST hybrid material_u=%d held=%d CPU=%d GPU=%d PASS\n",use_u,busy,cpu_count,gpu_count);
      const auto saved=trial;
      in[dc-1]=1e99;
      rc=snrt_hybrid_dust_material_c(in.data(),table.data(),trial.data(),dc,dg,nt,use_u,.1,1,10,1e-9,4);
      in[dc-1]=.3;
      if(!rc||trial!=saved)return 18;
    }
  }
  std::printf("HYBRID_MIXED_BUSY_ROLLBACK_PASS gpu=%d primary_relative=%g\n",gpu,worst);
}
