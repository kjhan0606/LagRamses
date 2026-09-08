#include <cstddef>
#include <cmath>
using std::isfinite;
#include "snrt_dust_material_cell.h"
#include <vector>
#include <algorithm>
#include <new>
#include <omp.h>

extern "C" int snrt_dust_material_openmp_c(const double *input,const double *table,double *output,
    int nc,int ng,int nt,int use_u,double dt,double background,double bath,double tolerance,int threads) {
  if(nc<1||ng<1||nt<2||use_u<0||use_u>15||use_u==4||
      ((use_u&8)&&(!(use_u&1)||nt>256))||!std::isfinite(dt)||dt<=0||threads<1)return 7;
  try {
    std::vector<double> trial(size_t(ng+2+(use_u%4>=2))*nc);
    int error=0;
    #pragma omp parallel for schedule(static) num_threads(std::min(nc,threads)) reduction(max:error)
    for(int i=0;i<nc;++i)
      error=std::max(error,dust_material_cell(input,table,trial.data(),nc,ng,nt,use_u,dt,background,bath,tolerance,i));
    if(error)return error;
    std::copy(trial.begin(),trial.end(),output);
    return 0;
  } catch(const std::bad_alloc&) {return 7;}
}
