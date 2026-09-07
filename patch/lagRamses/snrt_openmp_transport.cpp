// Host implementation of the production primary-photon operator. All writes
// are staged; invalid input/allocation failure leaves caller state untouched.
#include "snrt_species_dust_cell.h"
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

extern "C" int snrt_openmp_species_dust_c(
    float *state, const float *direction, const int *neighbor,
    const float *tau, const float *species_tau, const float *dust_tau,
    float *available, float *hhe, float *dust, float *returned, float *raw,
    float *group_absorbed, float *absorbed,
    int nowned, int nwork, int ndirection, int ngroup, float cdt) {
  if (!state || !direction || !neighbor || !tau || !species_tau || !dust_tau ||
      !available || !hhe || !dust || !returned || !raw || !group_absorbed || !absorbed ||
      nowned<=0 || nwork<nowned || ndirection<=0 || ngroup<=0 || !std::isfinite(cdt) || cdt<0)
    return 1;
  const size_t limit=std::numeric_limits<size_t>::max()/sizeof(float);
  if (size_t(nwork)>limit/size_t(ndirection)/size_t(ngroup) ||
      size_t(nowned)>limit/size_t(ngroup)/3) return 1;
  const size_t total=size_t(nwork)*ndirection*ngroup, groups=size_t(nowned)*ngroup;
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
    const int team=std::min(nowned,host_threads>0?host_threads:omp_get_max_threads());
    std::vector<float> next(state,state+total), removed(total,0), budget(available,available+size_t(3)*nowned);
    std::vector<float> hh(3*groups),dd(groups),rr(groups),raw_stage(groups),aa(groups),sum(nowned);
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
      next[linear]=before*expf(-fmaxf(0,tau[size_t(group)*nowned+cell]));
      removed[linear]=before-next[linear];
    }
    if(invalid) return 3;
    #pragma omp parallel for num_threads(team) schedule(static)
    for(int cell=0;cell<nowned;++cell)
      snrt_cap_species_dust_cell(next.data(),removed.data(),species_tau,dust_tau,budget.data(),
          hh.data(),dd.data(),rr.data(),raw_stage.data(),aa.data(),sum.data(),nowned,nwork,ndirection,ngroup,cell);
    // Detect arithmetic overflow before publishing either photons or atoms.
    for(float v:next) if(!std::isfinite(v) || v<0) return 3;
    for(float v:sum) if(!std::isfinite(v) || v<0) return 3;
    for(float v:raw_stage) if(!std::isfinite(v) || v<0) return 3;
    std::copy(next.begin(),next.end(),state); std::copy(budget.begin(),budget.end(),available);
    std::copy(hh.begin(),hh.end(),hhe); std::copy(dd.begin(),dd.end(),dust);
    std::copy(rr.begin(),rr.end(),returned); std::copy(raw_stage.begin(),raw_stage.end(),raw);
    std::copy(aa.begin(),aa.end(),group_absorbed); std::copy(sum.begin(),sum.end(),absorbed);
    return 0;
  } catch(const std::bad_alloc&) { return 4; }
}
