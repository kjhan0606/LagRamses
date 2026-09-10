// Toolkit-free SNRT link. These are unavailable-device ABI endpoints, NOT
// surrogate physics kernels. Auto/OpenMP use the existing host operators;
// every explicit GPU compute entry rejects without changing caller arrays.
#include "snrt_hybrid.h"
#include <algorithm>
namespace { constexpr int unavailable=7; }
extern "C" {
int snrt_cuda_available_c() { return 0; }
int snrt_cuda_species_dust_energy_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*,double*,const double*,
    double*,double*,double*) { return 8; }
int snrt_cuda_prepare_c(int, char *uuid) {
  if(uuid)std::fill_n(uuid,33,' ');
  return unavailable;
}
long long snrt_cuda_free_bytes_c() { return 0; }
int snrt_scatter_batch_c(const double*,const double*,double*,int,int,int,double,int) { return unavailable; }
int snrt_exchange_batch_c(const double*,const double*,double*,int,int,double,double,int) { return unavailable; }
void cuda_pool_init(int,int) {}
int cuda_pool_is_initialized() { return 0; }
int cuda_acquire_stream() { return -1; }
void cuda_release_stream(int) {}
int snrt_hybrid_try_acquire_c(long long,int) { return -1; }
int snrt_cuda_angular_reduce_tf32_c(const float*,const float*,float*,int,int,int) { return unavailable; }
int snrt_cuda_weighted_sum_fp32_c(const float*,const float*,float*,int,int) { return unavailable; }
int snrt_cuda_upwind_periodic_c(float*,const float*,int,int,int,int,float) { return unavailable; }
int snrt_cuda_upwind_sparse_c(float*,const float*,const int*,int,int,float) { return unavailable; }
int snrt_cuda_absorb_c(float*,const float*,float*,int,int) { return unavailable; }
int snrt_cuda_transport_absorb_c(float*,const float*,const int*,const float*,float*,int,int,float) {
  return unavailable;
}
int snrt_cuda_transport_absorb_limited_c(float*,const float*,const int*,const float*,const float*,
    float*,int,int,float) { return unavailable; }
int snrt_cuda_multigroup_rt_step_c(float*,const float*,const int*,const float*,const float*,
    float*,float*,int,int,int,float) { return unavailable; }
int snrt_cuda_multigroup_rt_step_owned_c(float*,const float*,const int*,const float*,const float*,
    float*,float*,int,int,int,int,float) { return unavailable; }
int snrt_cuda_multigroup_rt_step_species_c(float*,const float*,const int*,const float*,const float*,
    float*,float*,float*,int,int,int,int,float) { return unavailable; }
int snrt_cuda_multigroup_rt_step_species_dust_c(float*,const float*,const int*,const float*,const float*,
    const float*,float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float) { return unavailable; }
int snrt_cuda_species_dust_batch_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,int) { return unavailable; }
int snrt_cuda_species_dust_moment_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,double*) { return unavailable; }
int snrt_cuda_species_dust_moment_batch_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,int,double*) { return unavailable; }
int snrt_dust_material_cuda_c(const double*,const double*,double*,int,int,int,int,
    double,double,double,double,int) { return unavailable; }
int snrt_dust_material_batch_c(const double*,const double*,double*,int,int,int,int,
    double,double,double,double,int) { return unavailable; }
int snrt_ir_batch_c(const double*,const int*,const double*,double*,int,int,int,double,double,int,int) {
  return unavailable;
}
}
