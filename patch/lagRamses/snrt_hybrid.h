#pragma once
#include <cstddef>
using SnrtStep = int(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float);
extern "C" {
int snrt_scatter_batch_c(const double*,const double*,double*,int,int,int,double,int);
int snrt_isotropic_scatter_c(float*,const double*,const double*,int,int,int,int);
int snrt_exchange_batch_c(const double*,const double*,double*,int,int,double,double,int);
int snrt_dust_exchange_c(const double*,const double*,double*,int,int,double,double,int);
int snrt_ir_batch_c(const double*,const int*,const double*,double*,int,int,int,double,double,int,int);
int snrt_ir_transport_c(const double*,const double*,const int*,const int*,const int*,const double*,
    const double*,const double*,double*,double*,double*,double*,int,int,int,int,double,double,int);
int snrt_ir_absorb_c(const double*,const double*,const double*,const double*,const double*,const double*,
    double*,double*,int,int,int,double,double,int);
SnrtStep snrt_openmp_species_dust_c, snrt_serial_species_dust_c, snrt_hybrid_species_dust_c;
int snrt_cuda_species_dust_batch_c(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float,int);
int snrt_dust_material_batch_c(const double*,const double*,double*,int,int,int,int,
    double,double,double,double,int);
int snrt_hybrid_try_acquire_c(long long bytes,int sharers);
void cuda_pool_init(int local_rank,int streams);
int cuda_pool_is_initialized();
int cuda_acquire_stream();
void cuda_release_stream(int slot);
int snrt_hybrid_configure_c(int rank,int streams,int cells,int threads,int sharers,int gpu);
void snrt_hybrid_counts_c(int op,int *cpu,int *gpu);
}
