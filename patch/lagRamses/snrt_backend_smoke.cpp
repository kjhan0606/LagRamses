#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
using step_fn=int(float*,const float*,const int*,const float*,const float*,const float*,float*,
    float*,float*,float*,float*,float*,float*,int,int,int,int,float);
extern "C" step_fn snrt_openmp_species_dust_c,snrt_cuda_multigroup_rt_step_species_dust_c;
extern "C" int snrt_cuda_available_c();
int main() {
  constexpr int n=17,nw=19,nd=8,ng=9,g=n*ng,total=nw*nd*ng;
  const bool gpu=snrt_cuda_available_c()>0;
  double worst=0;
  for(int dusty=0;dusty<2;++dusty) {
    std::vector<float> input(total),dir(3*nd),tau(g),stau(3*g),dtau(g),budget(3*n);
    std::vector<int> neighbor(6*n);
    for(int i=0;i<total;++i)input[i]=0.001f*(1+(i*37)%51);
    for(int i=0;i<nd;++i)for(int d=0;d<3;++d)dir[3*i+d]=((i>>d)&1)?0.577350269f:-0.577350269f;
    for(int i=0;i<6*n;++i)neighbor[i]=1+(i*7)%nw;
    for(int i=0;i<3*n;++i)budget[i]=0.01f*(1+i%13);
    for(int j=0;j<ng;++j)for(int i=0;i<n;++i) {
      for(int s=0;s<3;++s)stau[s*g+j*n+i]=j>=4+s?0.03f*(1+(s+j+i)%17):0;
      dtau[j*n+i]=dusty?0.021f*(1+(j+i)%7):0;
      tau[j*n+i]=stau[j*n+i]+stau[g+j*n+i]+stau[2*g+j*n+i]+dtau[j*n+i];
    }
    std::vector<float> state=input,atoms=budget,hhe(3*g),dust(g),returned(g),raw(g),abs_group(g),absorbed(n);
    int rc=snrt_openmp_species_dust_c(state.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
        atoms.data(),hhe.data(),dust.data(),returned.data(),raw.data(),abs_group.data(),absorbed.data(),n,nw,nd,ng,0.2f);
    if(rc) return 10+rc;
    for(int i=0;i<g;++i)if(fabsf(raw[i]-returned[i]-abs_group[i])>1e-6f)return 20;
    for(int s=0;s<3;++s)for(int i=0;i<n;++i) {
      float used=0;
      for(int j=0;j<ng;++j)used+=hhe[s*g+j*n+i];
      if(fabsf(budget[s*n+i]-atoms[s*n+i]-used)>1e-6f)return 21;
    }
    if(gpu) {
      auto compare=[&](const std::vector<float>& a,const std::vector<float>& b,bool state_array=false) {
        double scale=1e-30,error=0;
        for(size_t i=0;i<a.size();++i) {
          if(state_array && i%nw>=n)continue; // CUDA ghost workspace is not an output.
          scale=std::max(scale,double(fabsf(a[i])));
          error=std::max(error,double(fabsf(a[i]-b[i])));
        }
        // Nearly exhausted inventories are differences of large FP32 values.
        // Normalize their rounding floor to the initial photon/atom budget,
        // not the tiny remainder; retain a relative test for nonzero ledgers.
        const double inventory_scale=std::max(double(*std::max_element(budget.begin(),budget.end())),
            double(*std::max_element(input.begin(),input.end()))*nd);
        worst=std::max(worst,error/std::max(scale,inventory_scale));
        return error<=3e-5*scale+8*std::numeric_limits<float>::epsilon()*inventory_scale;
      };
      auto gs=input,ga=budget;std::vector<float> gh(3*g),gd(g),gr(g),gw(g),gg(g),gt(n);
      rc=snrt_cuda_multigroup_rt_step_species_dust_c(gs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          ga.data(),gh.data(),gd.data(),gr.data(),gw.data(),gg.data(),gt.data(),n,nw,nd,ng,0.2f);
      if(rc)return 30+rc;
      if(!compare(state,gs,true)||!compare(atoms,ga)||!compare(hhe,gh)||!compare(dust,gd)||
          !compare(returned,gr)||!compare(raw,gw)||!compare(abs_group,gg)||!compare(absorbed,gt))return 40;
    }
    const auto saved=state,saved_atoms=atoms;
    tau[0]=std::numeric_limits<float>::quiet_NaN();
    rc=snrt_openmp_species_dust_c(state.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),atoms.data(),
        hhe.data(),dust.data(),returned.data(),raw.data(),abs_group.data(),absorbed.data(),n,nw,nd,ng,0.2f);
    if(rc==0||state!=saved||atoms!=saved_atoms)return 50;
  }
  std::printf("PASS OpenMP closure/atomic reject; CUDA comparison=%s max-relative=%g\n",gpu?"PASS":"SKIP",worst);
}
