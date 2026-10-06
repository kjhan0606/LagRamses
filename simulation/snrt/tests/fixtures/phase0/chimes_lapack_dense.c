/* Optimized double-precision direct LU, experimental native comparison.
 * Standard partial-pivot LU; no iterative tolerance or physics changes.
 * Unlike exact-zero reuse this can change rounding, so endpoint comparisons
 * must apply the existing species/energy criteria, not a speed-only check.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <mkl.h>
#include <sundials/sundials_dense.h>
static _Thread_local long factors,solves;
sunindextype denseGETRF(realtype **a,sunindextype m,sunindextype n,sunindextype *p)
{
    if(m>256 || n>256 || m<n || n<1)abort();
    for(int j=0;j<n;j++)if(a[j]!=a[0]+j*m)abort();
    MKL_INT pivots[256];
    int previous=mkl_set_num_threads_local(1);
    MKL_INT info=LAPACKE_dgetrf_work(LAPACK_COL_MAJOR,m,n,a[0],m,pivots);
    mkl_set_num_threads_local(previous);
    if(info<0)abort();
    for(int i=0;i<n;i++)p[i]=pivots[i]-1;
    factors++;return info;
}
void denseGETRS(realtype **a,sunindextype n,sunindextype *p,realtype *b)
{
    if(n>256 || n<1)abort();
    for(int j=0;j<n;j++)if(a[j]!=a[0]+j*n)abort();
    MKL_INT pivots[256];for(int i=0;i<n;i++)pivots[i]=p[i]+1;
    int previous=mkl_set_num_threads_local(1);
    MKL_INT info=LAPACKE_dgetrs_work(LAPACK_COL_MAJOR,'N',n,1,a[0],n,pivots,b,n);
    mkl_set_num_threads_local(previous);
    if(info)abort();solves++;
}
__attribute__((destructor)) static void report(void)
{
    if(factors || solves)fprintf(stderr,"LAPACK_DENSE factors=%ld solves=%ld inner_threads=1\n",factors,solves);
}
