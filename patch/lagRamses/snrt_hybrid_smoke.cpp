#include "snrt_hybrid.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <limits>
using std::isfinite;
#include "snrt_dust_material_cell.h"
extern "C" int snrt_cuda_available_c();
extern "C" int snrt_dust_material_openmp_c(const double*,const double*,double*,int,int,int,int,double,double,double,double,int);
extern "C" int snrt_hybrid_dust_material_c(const double*,const double*,double*,int,int,int,int,double,double,double,double,int);

int main() {
  constexpr int n=1031,nw=n+3,nd=8,ng=9,g=n*ng,total=nw*nd*ng,batch=64;
  const bool gpu=snrt_cuda_available_c()>0;
  if(snrt_hybrid_configure_c(0,1,batch,4,1,gpu))return 1;
  {
    // Independent scalar BE residual and analytic fixed-Tdust ODE solution.
    // This exercises heating, cooling and weak/stiff limits, not CPU parity alone.
    for(double t0:{1.,20.,80.,1e4})for(double td:{10.,50.,100.})
      for(double kd:{0.,1e-12,.3,1e4,1e12}) {
        const double q=dust_gas_transfer(t0,1,kd,td),tg=t0-q;
        if(!std::isfinite(q)||tg<std::min(t0,td)-1e-11*t0||tg>std::max(t0,td)+1e-11*t0)return 91;
        const long double scaled_k=kd*sqrtl((long double)tg/t0);
        const long double r=scaled_k/(1+scaled_k);
        if(fabsl(q-r*(t0-td))>2e-12L*std::max(t0,td))return 92;
        if(kd==0&&q!=0)return 93;
      }
    for(double t0:{1.,1000.}) {
      constexpr double td=20.,alpha=.01,duration=10.;
      const double u=(sqrt(t0)-sqrt(td))/(sqrt(t0)+sqrt(td))*exp(-alpha*sqrt(td)*duration);
      const double exact=td*std::pow((1+u)/(1-u),2);
      double previous=0;
      for(int steps:{32,64,128,256}) {
        double tg=t0;
        for(int j=0;j<steps;++j)tg-=dust_gas_transfer(tg,1,alpha*sqrt(tg)*duration/steps,td);
        const double error=std::abs(tg-exact);
        if(!std::isfinite(error)||error<=0||(previous>0&&(previous/error<1.8||previous/error>2.2))) {
          std::printf("GAS_DUST_VARIABLE_SPEED refinement failure T0=%g steps=%d error=%g previous=%g\n",t0,steps,error,previous);
          return 94;
        }
        previous=error;
        std::printf("GAS_DUST_VARIABLE_SPEED ODE T0=%g steps=%d error=%.8g PASS\n",t0,steps,error);
      }
    }
    const double variable=dust_gas_transfer(1000,1,10,20),frozen=10./11*(1000-20);
    if(variable>=frozen||frozen-variable<50)return 95;
    std::printf("GAS_DUST_VARIABLE_SPEED residual/bounds/stiff/analytic_refinement PASS\n");
  }
  {
    constexpr int ec=1031,et=4;
    double table[2*et]={std::log(10.),std::log(20.),std::log(50.),std::log(100.),1e-23,2e-23,5e-23,1e-22};
    std::vector<double> in(5*ec),out(4*ec),reference(4*ec);
    for(int i=0;i<ec;++i){in[5*i]=(i%2?20:80)*1e-24;in[5*i+1]=1e-24;in[5*i+2]=(i%2?80:20)*1e-24;in[5*i+3]=1;in[5*i+4]=i%7?1e-24:0;}
    for(double dt: {0.,1e-4,10.,1e8}) {
      if(snrt_dust_exchange_c(in.data(),table,reference.data(),ec,et,dt,10,1))return 80;
      for(int mode=0;mode<=2;++mode)for(int busy=0;busy<2;++busy) {
        std::fill(out.begin(),out.end(),-1);
        int held=busy&&gpu?cuda_acquire_stream():-1;
        if(busy&&gpu&&held<0)return 81;
        int rc=snrt_dust_exchange_c(in.data(),table,out.data(),ec,et,dt,10,mode);
        if(held>=0)cuda_release_stream(held);
        if(mode==2&&(busy||!gpu)){if(!rc||std::any_of(out.begin(),out.end(),[](double v){return v!=-1;}))return 82;continue;}
        if(rc)return 83;
        for(int i=0;i<ec;++i) {
          const double eg=in[5*i],ed=in[5*i+2],cg=in[5*i+1],eg1=out[4*i],ed1=out[4*i+1],td1=out[4*i+2],q=out[4*i+3];
          if(std::abs((eg1+ed1)-(eg+ed))>1e-14*(eg+ed)||eg1<0||ed1<0||!std::isfinite(td1))return 84;
          if(dt==0||in[5*i+4]==0){if(eg1!=eg||ed1!=ed||q!=0)return 85;}
          else {
            // Independently check the BE equation with a scale safe at stiff dt.
            double r=dt*in[5*i+4]/(cg+dt*in[5*i+4]);
            if(std::abs(q-r*(eg-cg*td1))>1e-13*(eg+ed))return 86;
            if(i%2?q>0:q<0)return 87;
          }
          for(int k=0;k<4;++k)if(std::abs(out[4*i+k]-reference[4*i+k])>1e-12*std::max(std::abs(reference[4*i+k]),1e-24))return 88;
        }
        const auto saved=out;in[5*(ec-1)]=std::numeric_limits<double>::quiet_NaN();
        rc=snrt_dust_exchange_c(in.data(),table,out.data(),ec,et,dt,10,mode);in[5*(ec-1)]=80e-24;
        if(!rc||out!=saved)return 89;
      }
    }
    std::printf("GAS_DUST_EXCHANGE conservation/heating/cooling/zero/stiff/CPU/GPU/busy/rollback PASS\n");
  }
  // Native array order (direction,group,cell), unequal quadrature weights.
  // Independent analytic reference + half-step composition + atomic failure.
  {
    constexpr int sc=1031,sg=3,sd=4;
    double weights[sd]={1,2,3,4};
    std::vector<float> initial(sc*sg*sd),expected(initial.size());
    std::vector<double> tau(sc*sg),half(tau.size()),zero(tau.size(),0);
    const double optical_depth[4]={0,1e-8,.3,100};
    for(int i=0;i<sc;++i)for(int g=0;g<sg;++g) {
      tau[g*sc+i]=optical_depth[(i+g)%4];half[g*sc+i]=.5*tau[g*sc+i];
      double sum=0;
      for(int d=0;d<sd;++d){initial[(i*sg+g)*sd+d]=float((d==i%sd?8.:.1)*(g+1));sum+=initial[(i*sg+g)*sd+d];}
      const double remain=std::exp(-tau[g*sc+i]);
      for(int d=0;d<sd;++d)expected[(i*sg+g)*sd+d]=float(initial[(i*sg+g)*sd+d]*remain+sum*weights[d]/10*(1-remain));
    }
    auto identity=initial;
    if(snrt_isotropic_scatter_c(identity.data(),zero.data(),weights,sc,sg,sd,1)||identity!=initial)return 61;
    auto twice=initial;
    for(int repeat=0;repeat<2;++repeat)
      if(snrt_isotropic_scatter_c(twice.data(),half.data(),weights,sc,sg,sd,1))return 62;
    for(int mode=0;mode<=2;++mode)for(int busy=0;busy<2;++busy) {
      int held=(busy&&gpu)?cuda_acquire_stream():-1;
      if(busy&&gpu&&held<0)return 63;
      auto trial=initial;
      int rc=snrt_isotropic_scatter_c(trial.data(),tau.data(),weights,sc,sg,sd,mode);
      if(held>=0)cuda_release_stream(held);
      if(mode==2&&(busy||!gpu)) {if(!rc||trial!=initial)return 64;continue;}
      if(rc)return 65;
      int cpu,gpu_batches;snrt_hybrid_counts_c(4,&cpu,&gpu_batches);
      if(cpu+gpu_batches!=1+(sc-1)/batch)return 66;
      if((mode==1||busy||!gpu)&&gpu_batches)return 67;
      if(mode==2&&cpu)return 68;
      for(int i=0;i<sc;++i)for(int g=0;g<sg;++g) {
        double before=0,after=0;
        for(int d=0;d<sd;++d) {
          int k=(i*sg+g)*sd+d;before+=initial[k];after+=trial[k];
          if(!std::isfinite(trial[k])||trial[k]<0||std::abs(trial[k]-expected[k])>2e-6f||
              std::abs(trial[k]-twice[k])>2e-6f)return 69;
        }
        if(std::abs(after-before)>1e-7*before)return 70;
      }
      const auto saved=trial;const double old=tau.back();tau.back()=std::numeric_limits<double>::quiet_NaN();
      rc=snrt_isotropic_scatter_c(trial.data(),tau.data(),weights,sc,sg,sd,mode);tau.back()=old;
      if(!rc||trial!=saved)return 71;
      std::printf("SCATTER mode=%d held=%d CPU=%d GPU=%d conservation/analytic/rollback PASS\n",mode,busy,cpu,gpu_batches);
    }
  }
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
    std::vector<double> reference_moment(3*g);
    if(snrt_openmp_species_dust_moment_c(ref.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
        atoms.data(),hhe.data(),dust.data(),ret.data(),raw.data(),abs_g.data(),absorbed.data(),n,nw,nd,ng,.2f,
        reference_moment.data()))return 2;
    {
      auto eq=input,ea=budget;
      std::vector<double> reference(ng),shift(total),he(3*g),de(g),em(3*g),nm(3*g);
      std::vector<float> hh(3*g),dd(g),rr(g),ww(g),aa(g),tt(n);
      for(int j=0;j<ng;++j)reference[j]=20.+j;
      for(int k=0;k<total;++k)shift[k]=(k%2?-0.75:1.25)*reference[k/(nw*nd)]*input[k];
      const auto initial_shift=shift;
      if(snrt_openmp_species_dust_energy_c(eq.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          ea.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f,nm.data(),
          shift.data(),reference.data(),he.data(),de.data(),em.data()))return 120;
      if(eq!=ref||ea!=atoms||hh!=hhe||dd!=dust||rr!=ret||ww!=raw||aa!=abs_g||tt!=absorbed||nm!=reference_moment)return 121;
      for(int busy=0;busy<2;++busy) {
        const int held=(busy&&gpu)?cuda_acquire_stream():-1;
        auto q=input,a=budget;auto s=initial_shift;
        std::vector<double> h(3*g,-1),d(g,-1),m(3*g,-1),number_m(3*g,-1);
        const int rc=snrt_hybrid_species_dust_energy_c(q.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
            a.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f,number_m.data(),
            s.data(),reference.data(),h.data(),d.data(),m.data());
        if(held>=0)cuda_release_stream(held);
        if(rc||q!=eq||a!=ea||s!=shift||h!=he||d!=de||m!=em||number_m!=nm)return 122;
        int cpu_count,gpu_count;snrt_hybrid_counts_c(0,&cpu_count,&gpu_count);
        if(cpu_count!=1||gpu_count!=0)return 123;
        const auto saved_q=q,saved_a=a;
        const float old=tau[g-1];tau[g-1]=std::numeric_limits<float>::quiet_NaN();
        const int bad=snrt_hybrid_species_dust_energy_c(q.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
            a.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f,number_m.data(),
            s.data(),reference.data(),h.data(),d.data(),m.data());
        tau[g-1]=old;
        if(!bad||q!=saved_q||a!=saved_a||s!=shift||h!=he||d!=de||m!=em||number_m!=nm)return 124;
      }
      std::printf("PAIRED hybrid CPU admission / signed multigroup ghosts / scalar parity / failed trial rollback PASS\n");
    }
    for(int busy=0;busy<2;++busy) {
      const int held=(busy&&gpu)?cuda_acquire_stream():-1;
      if(busy&&gpu&&held<0)return 3;
      auto qs=input,at=budget;std::vector<float> hh(3*g),dd(g),rr(g),ww(g),aa(g),tt(n);
      std::vector<double> moment(3*g,123);
      int rc=snrt_hybrid_species_dust_moment_c(qs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          at.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f,moment.data());
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
      for(int k=0;k<3*g;++k)if(!std::isfinite(moment[k])||
          std::abs(moment[k]-reference_moment[k])>3e-5*std::max(1e-20,double(raw[k%g])))return 72;
      if((busy||!gpu)&&moment!=reference_moment)return 73;
      std::printf("PRIMARY hybrid dusty=%d held=%d CPU=%d GPU=%d PASS\n",dusty,busy,cpu_count,gpu_count);
      // Late-batch invalid opacity: successful sibling batches must not publish.
      const auto saved=qs,saved_at=at,saved_hh=hh,saved_dd=dd,saved_rr=rr,saved_ww=ww,saved_aa=aa,saved_tt=tt;
      const auto saved_moment=moment;
      const float old=tau[g-1];tau[g-1]=std::numeric_limits<float>::quiet_NaN();
      rc=snrt_hybrid_species_dust_moment_c(qs.data(),dir.data(),neighbor.data(),tau.data(),stau.data(),dtau.data(),
          at.data(),hh.data(),dd.data(),rr.data(),ww.data(),aa.data(),tt.data(),n,nw,nd,ng,.2f,moment.data());
      tau[g-1]=old;
      if(!rc||qs!=saved||at!=saved_at||hh!=saved_hh||dd!=saved_dd||rr!=saved_rr||ww!=saved_ww||aa!=saved_aa||tt!=saved_tt)return 10;
      if(moment!=saved_moment)return 74;
    }
  }
  constexpr int dc=1031,dg=2,nt=4;
  std::vector<double> table((dg+3)*nt),in(4*dc),ref((dg+2)*dc);
  const double temp[nt]={10,20,50,100},power[nt]={1,2,5,10};
  for(int k=0;k<nt;++k) {
    table[k]=std::log(temp[k]);table[nt+k]=power[k];table[2*nt+k]=temp[k]*temp[k];
    table[3*nt+k*dg]=.4*power[k];table[3*nt+k*dg+1]=.6*power[k];
  }
  for(int use_u=0;use_u<4;++use_u) {
    in.resize((use_u>=2?7:4)*dc);ref.resize((dg+2+(use_u>=2))*dc);
    for(int i=0;i<dc;++i) {in[i]=.3;in[dc+i]=1;in[2*dc+i]=use_u?400:20;in[3*dc+i]=1;}
    if(use_u>=2)for(int i=0;i<dc;++i){in[4*dc+i]=use_u==3&&i%2?100:400;in[5*dc+i]=10;in[6*dc+i]=.3;}
    if(use_u==3)for(int i=0;i<dc;++i)in[6*dc+i]=i%3==0?0:i%3==1?.3:1e8;
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
      if(use_u>=2)for(int i=0;i<dc;++i) {
        const double q=trial[(dg+2)*dc+i],ed=trial[(dg+1)*dc+i];
        double emitted=0;for(int g=0;g<dg;++g)emitted+=trial[i*dg+g]*.1;
        if(std::abs((ed-400)+emitted-.03-q)>1e-9)return 90;
        const double eg=in[4*dc+i],cg=in[5*dc+i],td=trial[dg*dc+i];
        const double kd=.1*in[6*dc+i]*(use_u==3?std::sqrt((eg-q)/eg):1);
        if(std::abs(q-kd/(cg+kd)*(eg-cg*td))>1e-10)return 96;
        if(use_u==2&&q<=0)return 97;
      }
      std::printf("DUST hybrid material_u=%d held=%d CPU=%d GPU=%d PASS\n",use_u,busy,cpu_count,gpu_count);
      const auto saved=trial;
      in[dc-1]=1e99;
      rc=snrt_hybrid_dust_material_c(in.data(),table.data(),trial.data(),dc,dg,nt,use_u,.1,1,10,1e-9,4);
      in[dc-1]=.3;
      if(!rc||trial!=saved)return 18;
    }
  }
  // IR uses independent frozen neighbor/ghost snapshots across batch edges.
  constexpr int ic=1031,ig=2,id=2,ir=ig*id,igroups=ic*ig;
  std::vector<double> iq(ic*ir),ghost(ir,.04),rho(ic),directions={1,0,0,-1,0,0};
  std::vector<double> sigma={1e-7,2.},weights={.5,.5},rate(igroups,.02);
  std::vector<int> links(6*ic),remote(6*ic),blocked(6*ic);
  for(int i=0;i<ic;++i) {
    rho[i]=.5+.01*(i%71);
    for(int k=0;k<ir;++k)iq[i*ir+k]=.01*(1+(i+k)%19);
    if(i)links[6*i]=i;
    if(i+1<ic)links[6*i+1]=i+2;
  }
  remote[0]=1;blocked[6*(ic-1)+1]=1;
  std::vector<double> tr(ic*ir),tx(igroups),lo(igroups),re(igroups),cand(ic*ir),abs(ic);
  auto transport=[&](auto& t,auto& x,auto& l,auto& r,int mode) {
    return snrt_ir_transport_c(iq.data(),ghost.data(),links.data(),remote.data(),blocked.data(),rho.data(),
        directions.data(),sigma.data(),t.data(),x.data(),l.data(),r.data(),ic,ig,id,1,.2,.2,mode);
  };
  auto absorb=[&](auto& c,auto& a,int mode) {
    return snrt_ir_absorb_c(tr.data(),tx.data(),lo.data(),re.data(),rate.data(),weights.data(),
        c.data(),a.data(),ic,ig,id,.1,1.,mode);
  };
  if(transport(tr,tx,lo,re,1)||absorb(cand,abs,1))return 20;
  for(int busy=0;busy<2;++busy) {
    const int held=busy&&gpu?cuda_acquire_stream():-1;
    if(busy&&gpu&&held<0)return 21;
    auto t=tr,x=tx,l=lo,r=re,c=cand,a=abs;
    int rc=transport(t,x,l,r,0);
    int cpu,gp;snrt_hybrid_counts_c(2,&cpu,&gp);
    if(rc||cpu+gp!=17||((busy||!gpu)&&gp)||(!busy&&gpu&&(!cpu||!gp)))return 22;
    auto equal=[](const auto& lhs,const auto& rhs) {
      for(size_t k=0;k<lhs.size();++k)if(!std::isfinite(rhs[k])||std::abs(lhs[k]-rhs[k])>1e-13*std::max(1.,std::abs(lhs[k])))return false;
      return true;
    };
    if(!equal(tr,t)||!equal(tx,x)||!equal(lo,l)||!equal(re,r))return 23;
    std::printf("IR transport held=%d CPU=%d GPU=%d PASS\n",busy,cpu,gp);
    rc=absorb(c,a,0);snrt_hybrid_counts_c(3,&cpu,&gp);
    if(rc||cpu+gp!=17||((busy||!gpu)&&gp)||(!busy&&gpu&&(!cpu||!gp)))return 24;
    if(!equal(cand,c)||!equal(abs,a))return 25;
    std::printf("IR absorption held=%d CPU=%d GPU=%d PASS\n",busy,cpu,gp);
    // Late-batch failure leaves every public output unchanged.
    const auto st=t,sx=x,sl=l,sr=r,sc=c,sa=a;
    const double saved_rho=rho.back();rho.back()=std::numeric_limits<double>::quiet_NaN();
    rc=transport(t,x,l,r,0);rho.back()=saved_rho;
    if(!rc||t!=st||x!=sx||l!=sl||r!=sr)return 26;
    rate.back()=std::numeric_limits<double>::quiet_NaN();
    rc=absorb(c,a,0);rate.back()=.02;
    if(!rc||c!=sc||a!=sa)return 27;
    // Explicit CUDA never silently becomes CPU when the only stream is held.
    if(busy||!gpu) {
      if(!transport(t,x,l,r,2)||!absorb(c,a,2)||t!=st||x!=sx||l!=sl||r!=sr||c!=sc||a!=sa)return 28;
    } else {
      if(transport(t,x,l,r,2)||absorb(c,a,2)||!equal(tr,t)||!equal(cand,c))return 29;
    }
    if(held>=0)cuda_release_stream(held);
  }
  std::printf("HYBRID_MIXED_BUSY_ROLLBACK_PASS gpu=%d primary_relative=%g\n",gpu,worst);
}
