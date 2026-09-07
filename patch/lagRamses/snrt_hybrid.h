#pragma once
#include <cstddef>
using SnrtStep = int(float*,const float*,const int*,const float*,const float*,const float*,
    float*,float*,float*,float*,float*,float*,float*,int,int,int,int,float);
extern "C" {
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
