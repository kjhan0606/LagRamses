// Nonblocking shared-stream leases: each OMP worker tries a cell batch once,
// otherwise computes it itself. Inputs stay immutable until all batches finish.
#include "snrt_hybrid.h"
#include <cmath>
using std::isfinite;
#include "snrt_dust_material_cell.h"
#include <algorithm>
#include <vector>
#include <new>
#include <limits>
#include <cstdio>
#include <omp.h>

namespace {
int batch_cells=256,team_size=1,device_sharers=1;
bool gpu_enabled=false,reported[2]={false,false};
int cpu_count[2]={0,0},gpu_count[2]={0,0};
struct Lease {
  int slot;
  explicit Lease(long long bytes):slot(gpu_enabled?snrt_hybrid_try_acquire_c(bytes,device_sharers):-1) {}
  ~Lease(){if(slot>=0)cuda_release_stream(slot);}
};
void counts(int op,int cpu,int gpu) {
  cpu_count[op]=cpu;gpu_count[op]=gpu;
  if(!reported[op]) {
    std::printf(" SNRT hybrid %s batches CPU=%d GPU=%d batch_cells=%d\n",
        op==0?"primary":"dust",cpu,gpu,batch_cells);
    reported[op]=true;
  }
}
}
extern "C" int snrt_hybrid_configure_c(int rank,int streams,int cells,int threads,int sharers,int gpu) {
  if(streams<1||cells<1||cells>1048576||threads<1||sharers<1)return 1;
  batch_cells=cells;team_size=threads;device_sharers=sharers;
  if(gpu)cuda_pool_init(rank,streams);
  gpu_enabled=gpu && cuda_pool_is_initialized();
  return gpu && !gpu_enabled?2:0;
}
extern "C" void snrt_hybrid_counts_c(int op,int *cpu,int *gpu) {
  if(op<0||op>1){*cpu=*gpu=-1;return;}
  *cpu=cpu_count[op];*gpu=gpu_count[op];
}

extern "C" int snrt_hybrid_species_dust_c(
    float *state,const float *direction,const int *neighbor,const float *tau,const float *species_tau,
    const float *dust_tau,float *available,float *hhe,float *dust,float *returned,float *raw,
    float *group_absorbed,float *absorbed,int nowned,int nwork,int ndirection,int ngroup,float cdt) {
  if(!state||!direction||!neighbor||!tau||!species_tau||!dust_tau||!available||!hhe||!dust||
      !returned||!raw||!group_absorbed||!absorbed||nowned<1||nwork<nowned||ndirection<1||ngroup<1||
      !std::isfinite(cdt)||cdt<0)return 1;
  const size_t limit=std::numeric_limits<size_t>::max()/sizeof(float);
  if(size_t(nwork)>limit/ndirection/ngroup||size_t(nowned)>limit/ngroup/3)return 1;
  const size_t total=size_t(nwork)*ndirection*ngroup,groups=size_t(nowned)*ngroup;
  // Validate indices before gathering. Unused ghosts must also satisfy the
  // same input contract as the original full-array operator.
  for(size_t k=0;k<6*size_t(nowned);++k)if(neighbor[k]<0||neighbor[k]>nwork)return 2;
  for(size_t k=0;k<total;++k)if(!std::isfinite(state[k])||state[k]<0)return 2;
  try {
    std::vector<float> next(state,state+total),atoms(available,available+3*size_t(nowned));
    std::vector<float> hh(3*groups),dd(groups),rr(groups),raw_out(groups),aa(groups),sum(nowned);
    const int nbatch=1+(nowned-1)/batch_cells;
    int error=0,cpu=0,gpu=0;
    #pragma omp parallel for num_threads(std::min(team_size,nbatch)) schedule(dynamic,1) reduction(max:error) reduction(+:cpu,gpu)
    for(int ib=0;ib<nbatch;++ib) {
      try {
        const int first=ib*batch_cells,n=std::min(batch_cells,nowned-first),nw=7*n;
        const size_t g=size_t(n)*ngroup,t=size_t(nw)*ndirection*ngroup;
        // Complete six-neighbor snapshot per owned cell (including remote
        // ghosts). No batch ever reads another batch's updated photon state.
        std::vector<float> q(t),ta(g),st(3*g),dt(g),at(3*size_t(n));
        std::vector<float> bh(3*g),bd(g),br(g),bw(g),ba(g),bt(n);
        std::vector<int> links(6*size_t(n));
        for(int i=0;i<n;++i) {
          const int cell=first+i;
          for(int face=0;face<6;++face)links[6*i+face]=neighbor[6*cell+face]?n+face*n+i+1:0;
          for(int s=0;s<3;++s)at[s*n+i]=available[size_t(s)*nowned+cell];
          for(int j=0;j<ngroup;++j) {
            ta[size_t(j)*n+i]=tau[size_t(j)*nowned+cell];
            dt[size_t(j)*n+i]=dust_tau[size_t(j)*nowned+cell];
            for(int s=0;s<3;++s)st[size_t(s)*g+size_t(j)*n+i]=species_tau[size_t(s)*groups+size_t(j)*nowned+cell];
            for(int d=0;d<ndirection;++d) {
              const size_t src=(size_t(j)*ndirection+d)*nwork,dst=(size_t(j)*ndirection+d)*nw;
              q[dst+i]=state[src+cell];
              for(int face=0;face<6;++face) {
                const int other=neighbor[6*cell+face];
                q[dst+n+face*n+i]=state[src+(other?other-1:cell)];
              }
            }
          }
        }
        // Mirrors wrapper arrays, plus headroom. This is admission, not a
        // reservation against other processes; errors after launch reject.
        Lease lease(4LL*(2*t+3LL*ndirection+10LL*n+12*g+1)+16777216LL);
        int rc;
        if(lease.slot>=0) {
          ++gpu;
          rc=snrt_cuda_species_dust_batch_c(q.data(),direction,links.data(),ta.data(),st.data(),dt.data(),
              at.data(),bh.data(),bd.data(),br.data(),bw.data(),ba.data(),bt.data(),n,nw,ndirection,ngroup,cdt,lease.slot);
        } else {
          ++cpu;
          rc=snrt_serial_species_dust_c(q.data(),direction,links.data(),ta.data(),st.data(),dt.data(),
              at.data(),bh.data(),bd.data(),br.data(),bw.data(),ba.data(),bt.data(),n,nw,ndirection,ngroup,cdt);
        }
        if(rc){error=std::max(error,rc);continue;}
        for(int i=0;i<n;++i) {
          const int cell=first+i;
          sum[cell]=bt[i];
          for(int s=0;s<3;++s)atoms[size_t(s)*nowned+cell]=at[s*n+i];
          for(int j=0;j<ngroup;++j) {
            const size_t dst=size_t(j)*nowned+cell,src=size_t(j)*n+i;
            dd[dst]=bd[src];rr[dst]=br[src];raw_out[dst]=bw[src];aa[dst]=ba[src];
            for(int s=0;s<3;++s)hh[size_t(s)*groups+dst]=bh[size_t(s)*g+src];
            for(int d=0;d<ndirection;++d)
              next[(size_t(j)*ndirection+d)*nwork+cell]=q[(size_t(j)*ndirection+d)*nw+i];
          }
        }
      } catch(const std::bad_alloc&) {error=4;}
    }
    counts(0,cpu,gpu);
    if(error)return error;
    for(float v:next)if(!std::isfinite(v)||v<0)return 3;
    for(float v:sum)if(!std::isfinite(v)||v<0)return 3;
    for(float v:raw_out)if(!std::isfinite(v)||v<0)return 3;
    std::copy(next.begin(),next.end(),state);std::copy(atoms.begin(),atoms.end(),available);
    std::copy(hh.begin(),hh.end(),hhe);std::copy(dd.begin(),dd.end(),dust);
    std::copy(rr.begin(),rr.end(),returned);std::copy(raw_out.begin(),raw_out.end(),raw);
    std::copy(aa.begin(),aa.end(),group_absorbed);std::copy(sum.begin(),sum.end(),absorbed);
    return 0;
  } catch(const std::bad_alloc&) {return 4;}
}

extern "C" int snrt_hybrid_dust_material_c(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int) {
  if(!input||!table||!output||nc<1||ng<1||nt<2||!std::isfinite(dt)||dt<=0)return 7;
  try {
    std::vector<double> trial(size_t(ng+2)*nc);
    const int nbatch=1+(nc-1)/batch_cells;
    int error=0,cpu=0,gpu=0;
    #pragma omp parallel for num_threads(std::min(team_size,nbatch)) schedule(dynamic,1) reduction(max:error) reduction(+:cpu,gpu)
    for(int ib=0;ib<nbatch;++ib) {
      try {
        const int first=ib*batch_cells,n=std::min(batch_cells,nc-first);
        Lease lease(8LL*((ng+6LL)*n+(ng+3LL)*nt)+16777216LL);
        if(lease.slot<0) {
          ++cpu;
          for(int i=first;i<first+n;++i)error=std::max(error,dust_material_cell(input,table,trial.data(),
              nc,ng,nt,use_u,dt,background,bath,tolerance,i));
        } else {
          ++gpu;
          std::vector<double> in(4*size_t(n)),out(size_t(ng+2)*n);
          for(int f=0;f<4;++f)std::copy_n(input+size_t(f)*nc+first,n,in.data()+size_t(f)*n);
          const int rc=snrt_dust_material_batch_c(in.data(),table,out.data(),n,ng,nt,use_u,
              dt,background,bath,tolerance,lease.slot);
          if(rc){error=std::max(error,rc);continue;}
          std::copy_n(out.data(),size_t(ng)*n,trial.data()+size_t(ng)*first);
          for(int f=0;f<2;++f)std::copy_n(out.data()+size_t(ng+f)*n,n,trial.data()+size_t(ng+f)*nc+first);
        }
      } catch(const std::bad_alloc&) {error=7;}
    }
    counts(1,cpu,gpu);
    if(error)return error;
    std::copy(trial.begin(),trial.end(),output);
    return 0;
  } catch(const std::bad_alloc&) {return 7;}
}
