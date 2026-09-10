// Host implementation of the production primary-photon operator. All writes
// are staged; invalid input/allocation failure leaves caller state untouched.
#include "snrt_species_dust_cell.h"
#include "snrt_band_spectrum.h"
#include <algorithm>
#include <limits>
#include <vector>
#include <new>
#include <cstdlib>
#include <omp.h>

namespace { int host_threads=0; }
extern "C" int snrt_openmp_configure_c(int local_size) {
  // Do not change the hydro/feedback team's global OpenMP settings. Respect
  // an explicit user team size; otherwise divide the available CPU budget
  // between local ranks (a scheduler per-task budget is already divided).
  const char *explicit_threads=std::getenv("OMP_NUM_THREADS");
  if(explicit_threads && *explicit_threads)host_threads=omp_get_max_threads();
  else {
    int available=omp_get_num_procs();
    const char *task_cpus=std::getenv("SLURM_CPUS_PER_TASK");
    char *end=nullptr;
    const long count=task_cpus?std::strtol(task_cpus,&end,10):0;
    if(count>0 && end && !*end)available=std::min(long(available),count);
    else available/=std::max(1,local_size);
    host_threads=std::max(1,available);
  }
  return host_threads;
}

static int snrt_host_species_dust_impl(
    float *state, const float *direction, const int *neighbor,
    const float *tau, const float *species_tau, const float *dust_tau,
    float *available, float *hhe, float *dust, float *returned, float *raw,
    float *group_absorbed, float *absorbed,
    int nowned, int nwork, int ndirection, int ngroup, float cdt, int threads,
    double *dust_moment=nullptr, double *shift=nullptr, const double *reference_ev=nullptr,
    double *hhe_energy=nullptr, double *dust_energy=nullptr, double *dust_energy_moment=nullptr,
    const double *columns=nullptr,const double *edges=nullptr,
    snrt_band::Secondary secondary=nullptr,const double *xi=nullptr,double *deposition=nullptr,
    const snrt_band::Dust *grains=nullptr) {
  if (!state || !direction || !neighbor || !tau || !species_tau || !dust_tau ||
      !available || !hhe || !dust || !returned || !raw || !group_absorbed || !absorbed ||
      nowned<=0 || nwork<nowned || ndirection<=0 || ngroup<=0 || !std::isfinite(cdt) || cdt<0)
    return 1;
  const size_t limit=std::numeric_limits<size_t>::max()/sizeof(float);
  if (size_t(nwork)>limit/size_t(ndirection)/size_t(ngroup) ||
      size_t(nowned)>limit/size_t(ngroup)/3) return 1;
  if(grains){
    if(grains->bins!=4 && grains->bins!=6)return 1;
    if(!columns || !secondary || dust_moment || !grains->columns || !grains->abs ||
        !grains->transport || !grains->energy || !grains->weights || ngroup!=9 || !edges)return 1;
    for(int b=0;b<grains->bins;++b)for(int c=0;c<nowned;++c)
      if(!std::isfinite(grains->columns[b*nowned+c]) || grains->columns[b*nowned+c]<0)return 2;
    double norm=0;
    for(int d=0;d<ndirection;++d){
      if(!std::isfinite(grains->weights[d]) || grains->weights[d]<=0)return 2;
      norm+=grains->weights[d];
    }
    if(std::abs(norm-1)>2e-14)return 2;
    for(int g=0;g<=ngroup;++g)
      if(!std::isfinite(edges[g]) || edges[g]<=0 || (g>0 && edges[g]<=edges[g-1]) || edges[g]>1e4)return 2;
    constexpr int K=snrt_band::dust_nodes;
    for(int g=0;g<ngroup;++g){
      const snrt_band::Grid<K> reference_grid(edges[g],edges[g+1],g==ngroup-1);
      for(int j=0;j<K;++j){
        const double e=grains->energy[g*K+j];
        if(!std::isfinite(e) || std::abs(e/reference_grid.e[j]-1)>2e-13)return 2;
        if(j>0 && e<=grains->energy[g*K+j-1])return 2;
        for(int b=0;b<grains->bins;++b){
          const size_t k=(b*ngroup+g)*K+j;
          if(!std::isfinite(grains->abs[k]) || grains->abs[k]<0 ||
              !std::isfinite(grains->transport[k]) || grains->transport[k]<0)return 2;
        }
      }
      if(grains->energy[g*K]!=edges[g] || grains->energy[(g+1)*K-1]!=edges[g+1])return 2;
    }
  }
  const size_t total=size_t(nwork)*ndirection*ngroup, groups=size_t(nowned)*ngroup;
  if(secondary) {
    if(!columns || !xi || !deposition)return 1;
    for(int c=0;c<nowned;++c)if(!std::isfinite(xi[c]) || xi[c]<0 || xi[c]>1)return 2;
  } else if(xi || deposition)return 1;
  if(columns) {
    if(!shift || !edges || ngroup!=9)return 1;
    for(int g=0;g<=ngroup;++g)
      if(!std::isfinite(edges[g]) || edges[g]<=0 || (g>0 && edges[g]<=edges[g-1]) || edges[g]>1e4)return 2;
    for(int s=0;s<3;++s)for(int c=0;c<nowned;++c)
      if(!std::isfinite(columns[s*nowned+c]) || columns[s*nowned+c]<0)return 2;
    for(size_t i=0;i<groups;++i)if(dust_tau[i]!=0)return 8;
  } else if(edges)return 1;
  if(shift) {
    if(!reference_ev || !hhe_energy || !dust_energy || !dust_energy_moment ||
        total>std::numeric_limits<size_t>::max()/sizeof(double) ||
        groups>std::numeric_limits<size_t>::max()/sizeof(double)/3)return 1;
    for(int g=0;g<ngroup;++g)if(!std::isfinite(reference_ev[g]) || reference_ev[g]<=0)return 2;
    for(size_t i=0;i<total;++i) {
      const double e=reference_ev[i/(size_t(nwork)*ndirection)]*state[i]+shift[i];
      if(!std::isfinite(shift[i]) || !std::isfinite(e) || e<0 || (state[i]==0 && shift[i]!=0))return 2;
      if(columns && state[i]>0) {
        const int g=i/(size_t(nwork)*ndirection);
        if(e<edges[g]*state[i]*(1-snrt_band::moment_tolerance) ||
            e>edges[g+1]*state[i]*(1+snrt_band::moment_tolerance))return 2;
      }
    }
    for(int d=0;d<ndirection;++d)
      if(double(cdt)*(fabsf(direction[3*d])+double(fabsf(direction[3*d+1]))+fabsf(direction[3*d+2]))>1)return 2;
  } else if(reference_ev || hhe_energy || dust_energy || dust_energy_moment)return 1;
  if(dust_moment && groups>std::numeric_limits<size_t>::max()/sizeof(double)/3)return 1;
  for(size_t i=0;i<total;++i) if(!std::isfinite(state[i]) || state[i]<0) return 2;
  for(size_t i=0;i<size_t(3)*ndirection;++i) if(!std::isfinite(direction[i])) return 2;
  for(size_t i=0;i<size_t(6)*nowned;++i) if(neighbor[i]<0 || neighbor[i]>nwork) return 2;
  for(size_t i=0;i<3*groups;++i) if(!std::isfinite(species_tau[i]) || species_tau[i]<0) return 2;
  for(size_t i=0;i<size_t(3)*nowned;++i) if(!std::isfinite(available[i]) || available[i]<0) return 2;
  for(size_t i=0;i<groups;++i) {
    const float sum=species_tau[i]+species_tau[groups+i]+species_tau[2*groups+i]+dust_tau[i];
    const float scale=fmaxf(fmaxf(fabsf(tau[i]),fabsf(sum)),FLT_MIN);
    if(!std::isfinite(tau[i]) || tau[i]<0 || !std::isfinite(dust_tau[i]) || dust_tau[i]<0 ||
       !std::isfinite(sum) || fabsf(tau[i]-sum)>8*FLT_EPSILON*scale) return 2;
  }
  try {
    std::vector<snrt_band::Grid<>> grids;
    std::vector<snrt_band::Grid<snrt_band::dust_nodes>> grain_grids;
    if(grains)for(int g=0;g<ngroup;++g)
      grain_grids.emplace_back(edges[g],edges[g+1],g==ngroup-1,grains->energy+g*snrt_band::dust_nodes);
    else if(columns)for(int g=0;g<ngroup;++g)grids.emplace_back(edges[g],edges[g+1],g==ngroup-1);
    const int team=std::min(nowned,threads);
    std::vector<float> next(state,state+total), removed(total,0), budget(available,available+size_t(3)*nowned);
    std::vector<float> hh(3*groups),dd(groups),rr(groups),raw_stage(groups),aa(groups),sum(nowned);
    std::vector<double> moment(dust_moment?3*groups:0);
    std::vector<double> secondary_stage(secondary?size_t(8)*nowned:0);
    std::vector<double> shift_next,energy,he(shift?3*groups:0),de(shift?groups:0),em(shift?3*groups:0);
    std::vector<float> transported(shift?total:0);
    if(shift) {shift_next.assign(shift,shift+total);energy.resize(total);}
    int invalid=0;
    #pragma omp parallel for num_threads(team) reduction(max:invalid) schedule(static)
    for(size_t linear=0;linear<total;++linear) {
      const int cell=linear%nwork;
      if(cell>=nowned) continue;
      const int idir=(linear/nwork)%ndirection, group=linear/(size_t(nwork)*ndirection);
      const size_t base=linear-cell;
      const float q=state[linear],mx=direction[3*idir],my=direction[3*idir+1],mz=direction[3*idir+2];
      const int x=neighbor[6*cell+(mx>=0?0:1)],y=neighbor[6*cell+(my>=0?2:3)],z=neighbor[6*cell+(mz>=0?4:5)];
      const float qx=x>0?state[base+x-1]:q,qy=y>0?state[base+y-1]:q,qz=z>0?state[base+z-1]:q;
      const float before=q-cdt*(fabsf(mx)*(q-qx)+fabsf(my)*(q-qy)+fabsf(mz)*(q-qz));
      if(!std::isfinite(before) || before<0) invalid=1;
      next[linear]=columns?before:before*expf(-fmaxf(0,tau[size_t(group)*nowned+cell]));
      removed[linear]=before-next[linear];
      if(shift) {
        // Advect the signed correction with the photon stencil. Compensate
        // FP32 photon rounding when forming actual transported energy.
        const double s=shift[linear],sx=x>0?shift[base+x-1]:s,
            sy=y>0?shift[base+y-1]:s,sz=z>0?shift[base+z-1]:s;
        const double nb=double(q)-double(cdt)*(fabsf(mx)*(double(q)-qx)+
            fabsf(my)*(double(q)-qy)+fabsf(mz)*(double(q)-qz));
        const double sb=s-double(cdt)*(fabsf(mx)*(s-sx)+fabsf(my)*(s-sy)+fabsf(mz)*(s-sz));
        double eb=reference_ev[group]*nb+sb;
        // A convex energy stencil is the stable evaluation near E=0,
        // where nominal energy and signed correction nearly cancel.
        const double e=reference_ev[group]*q+s;
        const double ex=reference_ev[group]*qx+sx,ey=reference_ev[group]*qy+sy,ez=reference_ev[group]*qz+sz;
        const double wx=double(cdt)*fabsf(mx),wy=double(cdt)*fabsf(my),wz=double(cdt)*fabsf(mz);
        if(eb<1e-8*(reference_ev[group]*fabs(nb)+fabs(sb)))
          eb=(1-wx-wy-wz)*e+wx*ex+wy*ey+wz*ez;
        if(!std::isfinite(eb) || eb<0 || (before==0 && eb!=0)){
          invalid=1;
          if(columns)std::fprintf(stderr,"SNRT band invalid advection cell=%d group=%d N=%.17g E=%.17g\n",cell,group,double(before),eb);
        }
        energy[linear]=eb;transported[linear]=before;
        shift_next[linear]=eb-reference_ev[group]*next[linear];
      }
    }
    if(invalid) return 3;
    #pragma omp parallel for num_threads(team) reduction(max:invalid) schedule(static)
    for(int cell=0;cell<nowned;++cell) {
      if(columns) {
        try {
          if(grains){
            int status=snrt_band::absorb_cell(next.data(),shift_next.data(),transported.data(),
                energy.data(),reference_ev,columns,grain_grids,budget.data(),hh.data(),he.data(),rr.data(),
                raw_stage.data(),aa.data(),sum.data(),nowned,nwork,ndirection,ngroup,cell,
                secondary,xi,secondary_stage.data(),grains,dd.data(),de.data());
            if(status==0)status=snrt_band::scatter_cell(next.data(),shift_next.data(),reference_ev,
                grain_grids,*grains,nowned,nwork,ndirection,ngroup,cell);
            invalid=std::max(invalid,status);
          }else
          invalid=std::max(invalid,snrt_band::absorb_cell(next.data(),shift_next.data(),transported.data(),
              energy.data(),reference_ev,columns,grids,budget.data(),hh.data(),he.data(),rr.data(),
              raw_stage.data(),aa.data(),sum.data(),nowned,nwork,ndirection,ngroup,cell,
              secondary,xi,secondary_stage.data()));
        } catch(const std::bad_alloc&) {invalid=4;}
      } else snrt_cap_species_dust_cell(next.data(),removed.data(),species_tau,dust_tau,budget.data(),
          hh.data(),dd.data(),rr.data(),raw_stage.data(),aa.data(),sum.data(),nowned,nwork,ndirection,ngroup,cell,
          direction,dust_moment?moment.data():nullptr,shift?shift_next.data():nullptr,
          energy.data(),transported.data(),reference_ev,he.data(),de.data(),em.data());
    }
    if(invalid)return invalid;
    // Detect arithmetic overflow before publishing either photons or atoms.
    for(float v:next) if(!std::isfinite(v) || v<0) return 3;
    for(float v:sum) if(!std::isfinite(v) || v<0) return 3;
    for(float v:raw_stage) if(!std::isfinite(v) || v<0) return 3;
    for(double v:moment) if(!std::isfinite(v))return 3;
    for(double v:secondary_stage)if(!std::isfinite(v) || v<0)return 3;
    if(shift) {
      for(size_t i=0;i<total;++i) {
        const double e=reference_ev[i/(size_t(nwork)*ndirection)]*next[i]+shift_next[i];
        if(!std::isfinite(shift_next[i]) || !std::isfinite(e) || e<0 || (next[i]==0 && shift_next[i]!=0))return 3;
      }
      for(double v:he)if(!std::isfinite(v) || v<0)return 3;
      for(double v:de)if(!std::isfinite(v) || v<0)return 3;
      for(double v:em)if(!std::isfinite(v))return 3;
    }
    std::copy(next.begin(),next.end(),state); std::copy(budget.begin(),budget.end(),available);
    std::copy(hh.begin(),hh.end(),hhe); std::copy(dd.begin(),dd.end(),dust);
    std::copy(rr.begin(),rr.end(),returned); std::copy(raw_stage.begin(),raw_stage.end(),raw);
    std::copy(aa.begin(),aa.end(),group_absorbed); std::copy(sum.begin(),sum.end(),absorbed);
    if(dust_moment)std::copy(moment.begin(),moment.end(),dust_moment);
    if(secondary)std::copy(secondary_stage.begin(),secondary_stage.end(),deposition);
    if(shift) {
      std::copy(shift_next.begin(),shift_next.end(),shift);
      std::copy(he.begin(),he.end(),hhe_energy);std::copy(de.begin(),de.end(),dust_energy);
      std::copy(em.begin(),em.end(),dust_energy_moment);
    }
    return 0;
  } catch(const std::bad_alloc&) { return 4; }
}

extern "C" int snrt_openmp_band_d03_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,
    double *shift,const double *reference_ev,double *hhe_energy,double *dust_energy,double *energy_moment,
    const double *columns,const double *edges,snrt_band::Secondary secondary,const double *xi,double *deposition,
    const double *grain_columns,const double *kabs,const double *ksca,const double *node_energy,const double *weights) {
  if(!columns || !secondary)return 1;
  const snrt_band::Dust grains{grain_columns,kabs,ksca,node_energy,weights};
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads(),nullptr,
      shift,reference_ev,hhe_energy,dust_energy,energy_moment,columns,edges,secondary,xi,deposition,&grains);
}

// Separate entry keeps the original four-bin ABI and its callers intact.
extern "C" int snrt_openmp_band_grains_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *stau,const float *dtau,float *budget,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int no,int nw,int nd,int ng,float cdt,
    double *shift,const double *reference,double *he,double *de,double *em,const double *columns,const double *edges,
    snrt_band::Secondary secondary,const double *xi,double *deposition,const double *grain_columns,
    const double *kabs,const double *ksca,const double *node_energy,const double *weights,int bins) {
  const snrt_band::Dust grains{grain_columns,kabs,ksca,node_energy,weights,bins};
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,stau,dtau,budget,hhe,dust,returned,raw,
      group_absorbed,absorbed,no,nw,nd,ng,cdt,host_threads>0?host_threads:omp_get_max_threads(),nullptr,
      shift,reference,he,de,em,columns,edges,secondary,xi,deposition,&grains);
}

extern "C" int snrt_openmp_band_energy_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,double *moment,
    double *shift,const double *reference_ev,double *hhe_energy,double *dust_energy,double *energy_moment,
    const double *columns,const double *edges) {
  if(!columns)return 1;
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads(),moment,
      shift,reference_ev,hhe_energy,dust_energy,energy_moment,columns,edges);
}

extern "C" int snrt_openmp_band_secondary_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,double *moment,
    double *shift,const double *reference_ev,double *hhe_energy,double *dust_energy,double *energy_moment,
    const double *columns,const double *edges,snrt_band::Secondary secondary,const double *xi,double *deposition) {
  if(!columns || !secondary)return 1;
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads(),moment,
      shift,reference_ev,hhe_energy,dust_energy,energy_moment,columns,edges,secondary,xi,deposition);
}

extern "C" int snrt_openmp_species_dust_c(float *state,const float *direction,const int *neighbor,const float *tau,const float *species_tau,
    const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt) {
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads());
}
extern "C" int snrt_serial_species_dust_c(float *state,const float *direction,const int *neighbor,const float *tau,const float *species_tau,
    const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt) {
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,1);
}

extern "C" int snrt_openmp_species_dust_moment_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,double *moment) {
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads(),moment);
}
extern "C" int snrt_serial_species_dust_moment_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,double *moment) {
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,1,moment);
}

extern "C" int snrt_openmp_species_dust_energy_c(float *state,const float *direction,const int *neighbor,const float *tau,
    const float *species_tau,const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt,double *moment,
    double *shift,const double *reference_ev,double *hhe_energy,double *dust_energy,double *energy_moment) {
  if(!shift)return 1;
  return snrt_host_species_dust_impl(state,direction,neighbor,tau,species_tau,dust_tau,available,hhe,dust,returned,raw,
      group_absorbed,absorbed,nowned,nwork,ndirection,ngroup,cdt,host_threads>0?host_threads:omp_get_max_threads(),moment,
      shift,reference_ev,hhe_energy,dust_energy,energy_moment);
}
