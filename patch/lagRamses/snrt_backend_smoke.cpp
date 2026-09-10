#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <omp.h>
#include "snrt_species_dust_cell.h"
#include "snrt_hybrid.h"
#include "snrt_band_spectrum.h"
extern "C" int snrt_openmp_band_energy_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*,double*,const double*,
    double*,double*,double*,const double*,const double*);
extern "C" int snrt_openmp_band_secondary_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*,double*,const double*,
    double*,double*,double*,const double*,const double*,snrt_band::Secondary,const double*,double*);
extern "C" int snrt_openmp_band_d03_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*,const double*,double*,double*,double*,
    const double*,const double*,snrt_band::Secondary,const double*,double*,const double*,const double*,const double*,
    const double*,const double*);
extern "C" int snrt_openmp_band_grains_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*,const double*,double*,double*,double*,
    const double*,const double*,snrt_band::Secondary,const double*,double*,const double*,const double*,const double*,
    const double*,const double*,int);

static int test_secondary(double e,double xi,double *f) {
  if(xi==.5)return 1; // deliberately fail AFTER absorption is staged
  f[0]=e/(e+100);f[1]=1-f[0];f[2]=f[3]=f[4]=0;
  return 0;
}

static int d03_tests() {
  constexpr int nd=2,ng=9,K=128;
  const double edges[10]={.01,1,5.6,11.2,13.6,24.59,54.42,500,2000,10000};
  const float dir[6]={1,0,0,-1,0,0};const int neighbors[6]={1,1,1,1,1,1};
  for(int mode=0;mode<5;++mode){
    float q[nd*ng]={},tau[ng]={},stau[3*ng]={},dtau[ng]={},budget[3]={},hh[3*ng]={};
    float dd[ng]={},ret[ng]={},raw[ng]={},ag[ng]={},at[1]={};
    double shift[nd*ng]={},ref[ng],he[3*ng]={},de[ng]={},em[3*ng]={},columns[3]={},xi[1]={.1},dep[8]={};
    double grains[4]={.7,0,0,0},a[4*ng*K]={},s[4*ng*K]={},ev[ng*K],weights[nd]={.25,.75};
    for(int g=0;g<ng;++g){
      snrt_band::Grid<K> grid(edges[g],edges[g+1],g==ng-1);
      std::copy(grid.e.begin(),grid.e.end(),ev+g*K);ref[g]=.5*(edges[g]+edges[g+1]);
      for(int j=0;j<K;++j){a[g*K+j]=1;s[g*K+j]=2;}
      q[g*nd]=4;
    }
    if(mode==1)a[50]=-1;
    if(mode==2)weights[0]=1;
    if(mode==3)ev[7]=ev[8];
    if(mode==4)grains[0]=0;
    const auto before=std::vector<float>(q,q+nd*ng);
    const int rc=snrt_openmp_band_d03_c(q,dir,neighbors,tau,stau,dtau,budget,hh,dd,ret,raw,ag,at,
        1,1,nd,ng,0,shift,ref,he,de,em,columns,edges,test_secondary,xi,dep,grains,a,s,ev,weights);
    if(mode>=1 && mode<=3){
      if(rc!=2 || !std::equal(q,q+nd*ng,before.begin()))return 240;
      for(double v:shift)if(v!=0)return 241;
      for(double v:de)if(v!=0)return 242;
      continue;
    }
    if(rc)return 243;
    // Four-bin ABI and six-bin zero-Fe path must be bitwise identical,
    // including all particle/energy/secondary ledgers, not merely totals.
    float q6[nd*ng],b6[3]={},hh6[3*ng]={},dd6[ng]={},r6[ng]={},raw6[ng]={},ag6[ng]={},at6[1]={};
    double sh6[nd*ng]={},he6[3*ng]={},de6[ng]={},em6[3*ng]={},dep6[8]={};
    std::copy(before.begin(),before.end(),q6);
    double col6[6]={grains[0],grains[1],grains[2],grains[3],0,0};
    std::vector<double> a6(6*ng*K,1),s6(6*ng*K,2);
    std::copy(a,a+4*ng*K,a6.begin());std::copy(s,s+4*ng*K,s6.begin());
    const int rc6=snrt_openmp_band_grains_c(q6,dir,neighbors,tau,stau,dtau,b6,hh6,dd6,r6,raw6,ag6,at6,
        1,1,nd,ng,0,sh6,ref,he6,de6,em6,columns,edges,test_secondary,xi,dep6,col6,a6.data(),s6.data(),ev,weights,6);
    if(rc6 || !std::equal(q,q+nd*ng,q6) || !std::equal(shift,shift+nd*ng,sh6) ||
        !std::equal(budget,budget+3,b6) || !std::equal(hh,hh+3*ng,hh6) || !std::equal(dd,dd+ng,dd6) ||
        !std::equal(ret,ret+ng,r6) || !std::equal(raw,raw+ng,raw6) || !std::equal(ag,ag+ng,ag6) ||
        at[0]!=at6[0] || !std::equal(he,he+3*ng,he6) || !std::equal(de,de+ng,de6) ||
        !std::equal(em,em+3*ng,em6) || !std::equal(dep,dep+8,dep6))return 247;
    for(int g=0;g<ng;++g){
      const double t=std::exp(-grains[0]),scat=std::exp(-2*grains[0]);
      if(std::abs(q[g*nd]-4*t*(scat+.25*(1-scat)))>3e-7 ||
          std::abs(q[g*nd+1]-4*t*.75*(1-scat))>3e-7)return 244;
      if(std::abs(de[g]-4*ref[g]*(1-t))>2e-12*ref[g])return 245;
      const double after=ref[g]*(q[g*nd]+double(q[g*nd+1]))+shift[g*nd]+shift[g*nd+1]+de[g];
      if(std::abs(after/(4*ref[g])-1)>2e-12)return 246;
    }
  }
  std::puts("PASS D03128 grey limit / elastic N-E conservation / invalid coefficient-weight-grid atomic reject");
  std::puts("PASS six-bin zero-Fe bitwise identity with original four-bin ABI");
  return 0;
}

static int secondary_tests() {
  constexpr int nc=4,nd=2,ng=9,gs=nc*ng,sz=gs*nd,g=6;
  const double edges[10]={.01,1,5.6,11.2,13.6,24.59,54.42,500,2000,10000};
  const float dir[6]={1,0,0,-1,0,0};
  const int saved_threads=omp_get_max_threads();
  std::vector<double> serial;
  for(int mode=0;mode<4;++mode) {
    omp_set_num_threads(mode==1?4:1);
    float q[sz]={},tau[gs]={},stau[3*gs]={},dtau[gs]={},budget[3*nc],hh[3*gs]={};
    float dd[gs]={},ret[gs]={},raw[gs]={},ag[gs]={},at[nc]={};
    int neighbor[6*nc];
    double shift[sz]={},ref[ng],he[3*gs]={},de[gs]={},em[3*gs]={},column[3*nc]={},xi[nc],out[8*nc];
    for(int c=0;c<nc;++c){
      xi[c]=mode==3?.5:.1;column[c]=1e18;
      for(int s=0;s<3;++s)budget[s*nc+c]=mode==2?.05f:100.f;
      for(int j=0;j<6;++j)neighbor[6*c+j]=c+1;
    }
    for(int b=0;b<ng;++b)ref[b]=.5*(edges[b]+edges[b+1]);
    for(int d=0;d<nd;++d)for(int c=0;c<nc;++c){
      const int k=(g*nd+d)*nc+c;
      q[k]=.4f+.1f*c+.2f*d;shift[k]=(100+300*d-ref[g])*q[k];
    }
    for(int c=0;c<nc;++c)tau[g*nc+c]=stau[g*nc+c]=float(column[c]*snrt_band::sigma(0,ref[g]));
    std::fill(out,out+8*nc,-17);
    const std::vector<float> before(q,q+sz),old_budget(budget,budget+3*nc);
    const std::vector<double> old_shift(shift,shift+sz);
    const int rc=snrt_openmp_band_secondary_c(q,dir,neighbor,tau,stau,dtau,budget,hh,dd,ret,raw,ag,at,
        nc,nc,nd,ng,0,nullptr,shift,ref,he,de,em,column,edges,test_secondary,xi,out);
    if(mode==3){
      if(rc!=9 || !std::equal(q,q+sz,before.begin()) || !std::equal(shift,shift+sz,old_shift.begin()) ||
          !std::equal(budget,budget+3*nc,old_budget.begin()))return 221;
      for(double v:out)if(v!=-17)return 222;
      for(double v:he)if(v!=0)return 223;
      for(float v:hh)if(v!=0)return 224;
      continue;
    }
    if(rc)return 225;
    if(mode==0)serial.assign(out,out+8*nc);
    if(mode==1 && !std::equal(out,out+8*nc,serial.begin()))return 226;
    snrt_band::Grid<> grid(edges[g],edges[g+1]);
    for(int c=0;c<nc;++c){
      double expected[8]={},energy=0;
      for(int d=0;d<nd;++d){
        std::array<double,snrt_band::nodes> p;
        if(!grid.reconstruct(100+300*d,p))return 227;
        for(int j=0;j<snrt_band::nodes;++j){
          const double capture=before[(g*nd+d)*nc+c]*p[j]*-std::expm1(-column[c]*grid.cross[0][j]);
          const double excess=grid.e[j]-13.6;double f[5];test_secondary(excess,.1,f);
          expected[0]+=capture;energy+=capture*grid.e[j];
          for(int k=0;k<5;++k)expected[3+k]+=capture*excess*f[k];
        }
      }
      const double cap=std::min(1.,old_budget[c]*(1-4*FLT_EPSILON)/expected[0]);
      for(int k=0;k<8;++k)if(std::abs(out[k*nc+c]-cap*expected[k])>2e-12*std::max(cap*expected[k],1e-100))return 228;
      const double deposited=13.6*out[c]+out[3*nc+c]+out[4*nc+c];
      if(std::abs(deposited-cap*energy)>2e-12*cap*energy || std::abs(deposited-he[g*nc+c])>2e-12*deposited)return 229;
    }
  }
  omp_set_num_threads(saved_threads);
  std::puts("PASS node deposition independent quadrature/finite cap/OpenMP exactness/callback rollback");
  return 0;
}

static int band_tests() {
  constexpr double edges[10]={.01,1,5.6,11.2,13.6,24.59,54.42,500,2000,10000};
  // Positive reconstruction and boundary convention, including FP32 edge roundoff.
  for(int g=0;g<9;++g) {
    snrt_band::Grid<> grid(edges[g],edges[g+1],g==8);
    for(double f:{0.,1e-8,.1,.5,.9,1-1e-8,1.}) {
      std::array<double,snrt_band::nodes> p;
      const double mean=edges[g]+f*(edges[g+1]-edges[g]);
      if(!grid.reconstruct(mean,p))return 201;
      double n=0,e=0;for(int j=0;j<snrt_band::nodes;++j){if(p[j]<0)return 202;n+=p[j];e+=p[j]*grid.e[j];}
      if(std::abs(n-1)>1e-13 || std::abs(e-mean)>1e-11*edges[g+1])return 203;
    }
    std::array<double,snrt_band::nodes> p;
    if(grid.reconstruct(edges[g]*.999,p) || grid.reconstruct(edges[g+1]*1.001,p))return 204;
    if(!grid.reconstruct(edges[g]*(1-FLT_EPSILON),p)||p[0]!=1)return 205;
  }
  if(snrt_band::sigma(0,13.59)!=0 || snrt_band::sigma(1,24.58)!=0 || snrt_band::sigma(2,54.41)!=0)return 206;
  // Independently evaluated Verner table/formula values, cm2.
  const double expected[3]={6.346296358990503e-18,7.434698694110650e-18,1.587280257538649e-18};
  for(int s=0;s<3;++s)if(std::abs(snrt_band::sigma(s,snrt_band::threshold[s])/expected[s]-1)>2e-12)return 207;
  snrt_band::Grid<> below(11.2,13.6);
  for(int j=0;j<snrt_band::nodes;++j)if(below.cross[0][j]!=0)return 208;
  double refinement_max=0;
  for(int g=4;g<9;++g)for(double column:{1e15,1e18,1e21}) {
    snrt_band::Grid<64> a(edges[g],edges[g+1],g==8);
    snrt_band::Grid<256> b(edges[g],edges[g+1],g==8);
    std::array<double,64> pa;std::array<double,256> pb;
    const double mean=edges[g]+.3*(edges[g+1]-edges[g]);
    if(!a.reconstruct(mean,pa)||!b.reconstruct(mean,pb))return 209;
    double na=0,nb=0,ea=0,eb=0;
    for(int j=0;j<64;++j){const double q=pa[j]*-std::expm1(-column*a.cross[0][j]);na+=q;ea+=q*a.e[j];}
    for(int j=0;j<256;++j){const double q=pb[j]*-std::expm1(-column*b.cross[0][j]);nb+=q;eb+=q*b.e[j];}
    refinement_max=std::max({refinement_max,std::abs(na/nb-1),std::abs(ea/eb-1)});
  }
  std::printf("band 64/256 quadrature max relative absorption moment error %.12e\n",refinement_max);
  if(refinement_max>.005)return 210;
  double max_n_error=0,max_e_error=0;
  for(int mode=0;mode<6;++mode) {
    constexpr int no=2,nd=2,ng=9,sz=no*nd*ng,gs=no*ng;
    float q[sz],dir[6]={1,0,0,-1,0,0},tau[gs]={},stau[3*gs]={},dtau[gs]={};
    float budget[6],hhe[3*gs]={},dust[gs]={},ret[gs]={},raw[gs]={},ag[gs]={},at[no]={};
    double shift[sz],ref[ng],he[3*gs]={},de[gs]={},em[3*gs]={},column[6];
    const int neighbor[12]={2,2,1,1,1,1,1,1,2,2,2,2};
    for(int s=0;s<3;++s)for(int c=0;c<no;++c){budget[s*no+c]=mode==2?.1f:1000.f;
      column[s*no+c]=mode==0?0:mode==4?1e25:1e18/(s+1);}
    for(int g=0;g<ng;++g) {
      ref[g]=.5*(edges[g]+edges[g+1]);
      for(int d=0;d<nd;++d)for(int c=0;c<no;++c) {
        const int k=(g*nd+d)*no+c;
        q[k]=.5f+.2f*c+.1f*d;
        // Normal FP32 packets can have SUBNORMAL upwind increments.
        // FTZ erases those increments in N but not FP64 E, producing an
        // out-of-band mean even though both input moments are admissible.
        if(mode==5)q[k]=g==7?q[k]*5e-38f:0;
        const double mean=mode==5?500.5:edges[g]+(.25+.25*d)*(edges[g+1]-edges[g]);
        shift[k]=(mean-ref[g])*q[k];
      }
      for(int s=0;s<3;++s)for(int c=0;c<no;++c){stau[s*gs+g*no+c]=float(column[s*no+c]*snrt_band::sigma(s,ref[g]));tau[g*no+c]+=stau[s*gs+g*no+c];}
    }
    const std::vector<float> initial(q,q+sz),initial_budget(budget,budget+6);
    const std::vector<double> initial_shift(shift,shift+sz);
    double before_n=0,before_e=0;
    for(int k=0;k<sz;++k){before_n+=q[k];before_e+=ref[k/(nd*no)]*q[k]+shift[k];}
    if(mode==3)shift[0]=-ref[0]*q[0]; // E=0 is not a valid positive-energy band moment.
    const int rc=snrt_openmp_band_energy_c(q,dir,neighbor,tau,stau,dtau,budget,hhe,dust,ret,raw,ag,at,
        no,no,nd,ng,.25f,nullptr,shift,ref,he,de,em,column,edges);
    if(mode==3) {
      if(rc!=2 || !std::equal(q,q+sz,initial.begin()) || !std::equal(budget,budget+6,initial_budget.begin()))return 211;
      for(float v:hhe)if(v!=0)return 212;
      continue;
    }
    if(rc){std::printf("FAIL band mode=%d rc=%d\n",mode,rc);return 213;}
    double after_n=0,after_e=0;
    for(int k=0;k<sz;++k){after_n+=q[k];after_e+=ref[k/(nd*no)]*q[k]+shift[k];}
    for(float v:hhe)after_n+=v;
    for(double v:he)after_e+=v;
    max_n_error=std::max(max_n_error,std::abs(after_n/before_n-1));
    max_e_error=std::max(max_e_error,std::abs(after_e/before_e-1));
    for(int s=0;s<3;++s)for(int c=0;c<no;++c){double used=0;for(int g=0;g<ng;++g)used+=hhe[s*gs+g*no+c];
      if(budget[s*no+c]<0 || used>initial_budget[s*no+c]*(1+2e-7))return 214;}
    if(mode==2 && (hhe[gs+5*no]<=0 || hhe[2*gs+6*no]<=0))return 218;
    if(mode==4)for(int k=4*nd*no;k<5*nd*no;++k)if(q[k]!=0 || shift[k]!=0)return 219;
    if(mode==1) { // Mean survivor energy rises under selective H absorption.
      const int k=4*nd*no;
      const double before=ref[4]+initial_shift[k]/initial[k],after=ref[4]+shift[k]/q[k];
      std::printf("band hardening mean before=%.9f after=%.9f eV\n",before,after);
      if(after<=before)return 215;
    }
    if(mode==0)for(int k=0;k<sz;++k) {
      const double e0=ref[k/(nd*no)]*initial[k]+initial_shift[k];
      const double e1=ref[k/(nd*no)]*initial[k^1]+initial_shift[k^1];
      if(std::abs(ref[k/(nd*no)]*q[k]+shift[k]-(.75*e0+.25*e1))>1e-10)return 216;
    }
  }
  std::printf("band conservative ledger max relative N=%.12e E=%.12e\n",max_n_error,max_e_error);
  if(max_n_error>2e-7 || max_e_error>2e-12)return 217;
  std::puts("PASS band reconstruction/selective absorption/finite inventory/rollback");
  return 0;
}
using step_fn=int(float*,const float*,const int*,const float*,const float*,const float*,float*,
    float*,float*,float*,float*,float*,float*,int,int,int,int,float);
extern "C" step_fn snrt_openmp_species_dust_c,snrt_cuda_multigroup_rt_step_species_dust_c;
using moment_fn=int(float*,const float*,const int*,const float*,const float*,const float*,float*,
    float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*);
extern "C" moment_fn snrt_openmp_species_dust_moment_c,snrt_cuda_species_dust_moment_c;
extern "C" int snrt_cuda_available_c();
int main() {
  const int band_status=band_tests();if(band_status){std::printf("band test failed %d\n",band_status);return band_status;}
  const int dust_status=d03_tests();if(dust_status){std::printf("D03 test failed %d\n",dust_status);return dust_status;}
  const int secondary_status=secondary_tests();
  if(secondary_status){std::printf("secondary test failed %d\n",secondary_status);return secondary_status;}
  // Two periodic cells, opposite rays, deliberately different actual photon
  // energies. Check the independent conservative stencil and accepted energy
  // ledger, rather than reconstructing absorbed energy from group means.
  for(int mode=0;mode<6;++mode) {
    float q[4]={4,2,1,3},dir[6]={1,0,0,-1,0,0};
    const int links[12]={2,2,1,1,1,1,1,1,2,2,2,2};
    const double reference[1]={20};
    double shift[4]={-40,40,10,-30},he[6],de[2],em[6],moment[6];
    if(mode==4)for(int k=0;k<4;++k)shift[k]=-20.*q[k]; // E=0 is admissible.
    const float initial[4]={q[0],q[1],q[2],q[3]};
    double initial_e[4];for(int k=0;k<4;++k)initial_e[k]=20.*q[k]+shift[k];
    float dtau[2]={mode==0?0.f:mode==5?100.f:.5f,mode==0?0.f:mode==5?100.f:.5f};
    float stau[6]={mode==2?1.f:0.f,mode==2?1.f:0.f,0,0,0,0};
    float tau[2]={dtau[0]+stau[0],dtau[1]+stau[1]};
    float atoms[6]={.1f,.2f,0,0,0,0},hh[6],dust[2],ret[2],raw[2],ag[2],at[2];
    const float cdt=mode==3?0.f:.25f;
    float legacy[4]={4,2,1,3},la[6]={.1f,.2f,0,0,0,0},lh[6],ld[2],lr[2],lw[2],lg[2],lt[2];
    double lm[6];
    if(snrt_openmp_species_dust_moment_c(legacy,dir,links,tau,stau,dtau,la,lh,ld,lr,lw,lg,lt,2,2,2,1,cdt,lm))return 100;
    int rc=snrt_openmp_species_dust_energy_c(q,dir,links,tau,stau,dtau,atoms,hh,dust,ret,raw,ag,at,
        2,2,2,1,cdt,moment,shift,reference,he,de,em);
    if(rc){std::printf("FAIL paired mode=%d rc=%d\n",mode,rc);return 101;}
    if(!std::equal(q,q+4,legacy)||!std::equal(atoms,atoms+6,la)||!std::equal(hh,hh+6,lh)||
        !std::equal(dust,dust+2,ld)||!std::equal(ret,ret+2,lr)||!std::equal(raw,raw+2,lw)||
        !std::equal(ag,ag+2,lg)||!std::equal(at,at+2,lt)||!std::equal(moment,moment+6,lm))return 102;
    double before=0,after=0;
    for(int k=0;k<4;++k) {
      const double e=20.*q[k]+shift[k];
      if(e<0 || !std::isfinite(e) || (q[k]==0 && shift[k]!=0))return 103;
      before+=initial_e[k];after+=e;
      if(mode==0 && std::abs(e-((1-cdt)*initial_e[k]+cdt*initial_e[k^1]))>1e-12)return 104;
    }
    for(double e:he)after+=e;
    for(double e:de)after+=e;
    if(std::abs(before-after)>2e-12*std::max(1.,before))return 105;
    for(int cell=0;cell<2;++cell) {
      if(std::abs(em[cell])>de[cell]*(1+1e-14) || em[2+cell]!=0 || em[4+cell]!=0)return 106;
      if(mode==2) {
        const double photons=(1-cdt)*(initial[cell]+initial[2+cell])+cdt*(initial[cell^1]+initial[(2+cell)^1]);
        const double energy=(1-cdt)*(initial_e[cell]+initial_e[2+cell])+cdt*(initial_e[cell^1]+initial_e[(2+cell)^1]);
        if(hh[cell]<=0 || std::abs(he[cell]-hh[cell]*energy/photons)>1e-6)return 113;
        if(cell==0 && std::abs(he[cell]-20.*hh[cell])<1e-3)return 114;
      }
      if(mode==3) {
        const double fraction=double(initial[cell]-q[cell])/initial[cell];
        const double expected=fraction*(initial_e[cell]-initial_e[2+cell]);
        if(std::abs(em[cell]-expected)>1e-6)return 107;
        if(std::abs(de[cell]-20.*dust[cell])<1e-3)return 108;
      }
    }
    // Device capability failure and invalid signed state leave every paired
    // output and inventory untouched.
    auto reject=[&](bool device) {
      const std::vector<float> saved_q(q,q+4),saved_a(atoms,atoms+6),saved_h(hh,hh+6),saved_d(dust,dust+2),
          saved_r(ret,ret+2),saved_w(raw,raw+2),saved_g(ag,ag+2),saved_t(at,at+2);
      const std::vector<double> ss(shift,shift+4),sh(he,he+6),sd(de,de+2),sm(em,em+6),sn(moment,moment+6);
      SnrtEnergyStep *fn=device?snrt_cuda_species_dust_energy_c:snrt_openmp_species_dust_energy_c;
      const int status=fn(q,dir,links,tau,stau,dtau,atoms,hh,dust,ret,raw,ag,at,2,2,2,1,cdt,moment,
          shift,reference,he,de,em);
      return status!=0 && (!device||status==8) && std::equal(q,q+4,saved_q.begin()) &&
          std::equal(atoms,atoms+6,saved_a.begin()) && std::equal(hh,hh+6,saved_h.begin()) &&
          std::equal(dust,dust+2,saved_d.begin()) && std::equal(ret,ret+2,saved_r.begin()) &&
          std::equal(raw,raw+2,saved_w.begin()) && std::equal(ag,ag+2,saved_g.begin()) &&
          std::equal(at,at+2,saved_t.begin()) && std::equal(shift,shift+4,ss.begin()) &&
          std::equal(he,he+6,sh.begin()) && std::equal(de,de+2,sd.begin()) &&
          std::equal(em,em+6,sm.begin()) && std::equal(moment,moment+6,sn.begin());
    };
    if(!reject(true))return 109;
    shift[3]=-20.*q[3]-1;
    if(!reject(false))return 110;
  }
  {
    // A read-only ghost brings a negative correction into a nominal-energy
    // owned cell. Both fields must use the ghost, without modifying it.
    float q[2]={2,4},dir[3]={1,0,0},tau=0,stau[3]={},dtau=0,atoms[3]={},hh[3],d,r,w,g,t;
    int links[6]={2,1,1,1,1,1};double shift[2]={0,-60},reference=20,he[3],de,em[3];
    if(snrt_openmp_species_dust_energy_c(q,dir,links,&tau,stau,&dtau,atoms,hh,&d,&r,&w,&g,&t,
        1,2,1,1,.25f,nullptr,shift,&reference,he,&de,em))return 111;
    if(q[0]!=2.5f || q[1]!=4 || shift[1]!=-60 || std::abs(20*q[0]+shift[0]-35)>1e-13)return 112;
    q[1]=0;shift[1]=1; // N=0 forbids a residual correction, even positive E.
    if(snrt_openmp_species_dust_energy_c(q,dir,links,&tau,stau,&dtau,atoms,hh,&d,&r,&w,&g,&t,
        1,2,1,1,.25f,nullptr,shift,&reference,he,&de,em)!=2 || q[0]!=2.5f || shift[1]!=1)return 115;
  }
  std::printf("PASS paired signed-energy stencil / actual absorption / energy moment / ghost / atomic reject / CUDA admission\n");
  // Fully absorbed cells must never return negative photons through a cap
  // rounded above one. Exercise zero and trace dust, with ample H/He inventory.
  uint32_t seed=173;
  auto sample=[&]() {seed=1664525u*seed+1013904223u;return float(1+(seed>>8)%10000)/10000;};
  for(int trial=0;trial<4096;++trial) {
    float state[8]={},removed[8],stau[3],dtau=trial%2?1e-12f:0.0f;
    for(float &v:removed)v=sample();
    for(float &v:stau)v=sample();
    float available[3]={10,10,10},hhe[3],dust,returned,raw,group,total;
    snrt_cap_species_dust_cell(state,removed,stau,&dtau,available,hhe,&dust,&returned,&raw,
        &group,&total,1,1,8,1,0);
    for(float v:state)if(v<0||!std::isfinite(v)) {
      std::printf("FAIL saturated photon cap trial=%d state=%g raw=%g assigned=%g\n",
          trial,v,raw,hhe[0]+hhe[1]+hhe[2]);return 60;
    }
    if(hhe[0]+hhe[1]+hhe[2]>raw||fabsf(raw-returned-group)>2e-6f||
        fabsf(group-hhe[0]-hhe[1]-hhe[2]-dust)>2e-6f) {
      std::printf("FAIL cap ledger trial=%d raw=%.9g HHe=%.9g group=%.9g dust=%.9g return=%.9g\n",
          trial,raw,hhe[0]+hhe[1]+hhe[2],group,dust,returned);return 61;
    }
  }
  constexpr int n=17,nw=19,nd=8,ng=9,g=n*ng,total=nw*nd*ng;
  const bool gpu=snrt_cuda_available_c()>0;
  double worst=0;
  // Accepted beam moments, isotropic cancellation and physical H/He excess
  // transfer. Photon bins already contain their angular weights.
  for(int mode=0;mode<4;++mode) {
    float q[2]={4.0f,mode==0?0.0f:mode==1?4.0f:1.0f};
    const float initial[2]={q[0],q[1]},direction[6]={1,0,0,-1,0,0};
    const int links[6]={1,1,1,1,1,1};
    float dust_tau=mode==3?0.0f:0.5f,stau[3]={mode>=2?1.0f:0.0f,0,0},tau=dust_tau+stau[0];
    float atoms[3]={0,0,0},hhe[3],dust,returned,raw,group,total;
    double moment[3]={17,18,19};
    int rc=snrt_openmp_species_dust_moment_c(q,direction,links,&tau,stau,&dust_tau,atoms,hhe,&dust,&returned,&raw,
        &group,&total,1,1,2,1,0,moment);
    const double expected=dust*double(initial[0]-initial[1])/(initial[0]+initial[1]);
    if(rc || std::abs(moment[0]-expected)>2e-7*std::max(1.0,double(dust)) || moment[1]!=0 || moment[2]!=0)return 62;
    if(mode>=2 && !(returned>0))return 63;
    if(mode==3 && (dust!=0 || moment[0]!=0))return 64;
  }
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
    auto ms=input,ma=budget;
    std::vector<float> mh(3*g),md(g),mr(g),mw(g),mg(g),mt(n);
    std::vector<double> moment(3*g,123.0);
    rc=snrt_openmp_species_dust_moment_c(ms.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
        ma.data(),mh.data(),md.data(),mr.data(),mw.data(),mg.data(),mt.data(),n,nw,nd,ng,0.2f,moment.data());
    if(rc||ms!=state||ma!=atoms||mh!=hhe||md!=dust||mr!=returned||mw!=raw||mg!=abs_group||mt!=absorbed)return 65;
    for(int k=0;k<g;++k) {
      double norm=0;for(int axis=0;axis<3;++axis)norm+=moment[axis*g+k]*moment[axis*g+k];
      if(!std::isfinite(norm)||std::sqrt(norm)>dust[k]*(1+1e-6))return 66;
    }
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
        const bool ok=error<=3e-5*scale+8*std::numeric_limits<float>::epsilon()*inventory_scale;
        if(!ok)std::printf("FAIL backend comparison dusty=%d count=%zu state=%d error=%g scale=%g inventory=%g\n",
            dusty,a.size(),state_array,error,scale,inventory_scale);
        return ok;
      };
      auto gs=input,ga=budget;std::vector<float> gh(3*g),gd(g),gr(g),gw(g),gg(g),gt(n);
      rc=snrt_cuda_multigroup_rt_step_species_dust_c(gs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          ga.data(),gh.data(),gd.data(),gr.data(),gw.data(),gg.data(),gt.data(),n,nw,nd,ng,0.2f);
      if(rc)return 30+rc;
      if(!compare(state,gs,true)||!compare(atoms,ga)||!compare(hhe,gh)||!compare(dust,gd)||
          !compare(returned,gr)||!compare(raw,gw)||!compare(abs_group,gg)||!compare(absorbed,gt))return 40;
      gs=input;ga=budget;std::vector<double> gm(3*g,124.0);
      rc=snrt_cuda_species_dust_moment_c(gs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          ga.data(),gh.data(),gd.data(),gr.data(),gw.data(),gg.data(),gt.data(),n,nw,nd,ng,0.2f,gm.data());
      if(rc||!compare(state,gs,true)||!compare(dust,gd))return 67;
      for(int k=0;k<3*g;++k)if(!std::isfinite(gm[k]) ||
          std::abs(gm[k]-moment[k])>3e-5*std::max(1e-20,double(raw[k%g])))return 68;
    }
    const auto saved=state,saved_atoms=atoms;
    tau[0]=std::numeric_limits<float>::quiet_NaN();
    rc=snrt_openmp_species_dust_c(state.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),atoms.data(),
        hhe.data(),dust.data(),returned.data(),raw.data(),abs_group.data(),absorbed.data(),n,nw,nd,ng,0.2f);
    if(rc==0||state!=saved||atoms!=saved_atoms)return 50;
    const auto saved_moment=moment;
    rc=snrt_openmp_species_dust_moment_c(state.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),atoms.data(),
        hhe.data(),dust.data(),returned.data(),raw.data(),abs_group.data(),absorbed.data(),n,nw,nd,ng,0.2f,moment.data());
    if(rc==0||state!=saved||atoms!=saved_atoms||moment!=saved_moment)return 69;
  }
  std::printf("PASS OpenMP closure/atomic reject; CUDA comparison=%s max-relative=%g\n",gpu?"PASS":"SKIP",worst);
  std::printf("PASS accepted dust angular moment / scalar ABI / beam / isotropic / returned photon exclusion\n");
}
