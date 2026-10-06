/* SPDX-License-Identifier: BSD-3-Clause
 * Portions derived from SUNDIALS:
 * Copyright (c) 2002-2021, Lawrence Livermore National Security and Southern
 * Methodist University. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/* Exact-zero sparse triangular solves for the existing SUNDIALS dense LU.
 * Retains pivot selection and factorization; no numerical drop tolerance.
 * Experimental LD_PRELOAD adapter. Single-thread per cell, TLS across cells.
 * Bounded to the CHIMES matrix size. Nonfinite RHS/factors use upstream.
 */
#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>
#include <sundials/sundials_dense.h>
#define LIMIT 256
typedef sunindextype (*factor_fn)(realtype **,sunindextype,sunindextype,sunindextype *);
typedef void (*solve_fn)(realtype **,sunindextype,sunindextype *,realtype *);
static _Thread_local factor_fn original_factor;
static _Thread_local solve_fn original_solve;
static _Thread_local realtype **matrix;
static _Thread_local sunindextype *pivots,dimension;
static _Thread_local int lower[LIMIT+1],upper[LIMIT+1],indices[LIMIT*LIMIT];
static _Thread_local realtype values[LIMIT*LIMIT];
static _Thread_local long factors,solves,verified,fallback,entries,dense_entries;
static _Thread_local long factor_checked,factor_sparse_updates,factor_dense_updates;
static void *symbol(const char *name)
{
    void *p=dlsym(RTLD_NEXT,name);if(!p){fprintf(stderr,"missing %s\n",name);abort();}return p;
}
/* Same right-looking partial-pivot LU and arithmetic order as SUNDIALS 5.8
 * sundials_dense.c, denseGETRF (LLNL/SMU, BSD-3-Clause). Only exactly-zero
 * multipliers are skipped, and only when the operands are finite. */
static sunindextype zero_lu(realtype **a,int m,int n,sunindextype *p)
{
    int rows[LIMIT];
    for(int k=0;k<n;k++) {
        realtype *col=a[k];int pivot=k;
        for(int i=k+1;i<m;i++)if(fabs(col[i])>fabs(col[pivot]))pivot=i;
        p[k]=pivot;
        if(col[pivot]==0)return k+1;
        if(pivot!=k)for(int j=0;j<n;j++){realtype tmp=a[j][pivot];a[j][pivot]=a[j][k];a[j][k]=tmp;}
        realtype mult=1.0/col[k];int count=0,finite=1;
        for(int i=k+1;i<m;i++) {
            col[i]*=mult;
            if(!isfinite(col[i]))finite=0;
            if(col[i]!=0)rows[count++]=i;
        }
        for(int j=k+1;j<n;j++) {
            realtype *dest=a[j],v=dest[k];
            if(v==0)continue;
            /* A cost switch only, never a physical/numerical drop tolerance. */
            if(finite && isfinite(v) && count*2<m-k-1) {
                for(int r=0;r<count;r++){int i=rows[r];dest[i]-=v*col[i];}
                factor_sparse_updates++;
            } else {
                for(int i=k+1;i<m;i++)dest[i]-=v*col[i];
                factor_dense_updates++;
            }
        }
    }
    return 0;
}
sunindextype denseGETRF(realtype **a,sunindextype m,sunindextype n,sunindextype *p)
{
    if(!original_factor)original_factor=(factor_fn)symbol("denseGETRF");
    matrix=NULL;
    int sparse=getenv("SNRT_CHIMES_SPARSE_FACTOR") && m>=n && m<=LIMIT && n>=1;
    int check=sparse && getenv("SNRT_CHIMES_JAC_VERIFY");
    realtype *ref=NULL,*cols[LIMIT];sunindextype refp[LIMIT];
    if(check) {
        ref=malloc(m*n*sizeof(*ref));if(!ref)abort();
        for(int j=0;j<n;j++){cols[j]=ref+j*m;memcpy(cols[j],a[j],m*sizeof(*ref));}
    }
    sunindextype rc=sparse?zero_lu(a,m,n,p):original_factor(a,m,n,p);
    if(check) {
        sunindextype expected=original_factor(cols,m,n,refp);
        if(rc!=expected)abort();
        for(int j=0;j<(rc?rc:n);j++)if(p[j]!=refp[j])abort();
        for(int j=0;j<n;j++)for(int i=0;i<m;i++)
            if(a[j][i]!=cols[j][i] && !(isnan(a[j][i]) && isnan(cols[j][i]))) {
                fprintf(stderr,"SPARSE_LU_MISMATCH row=%d col=%d\n",i,j);abort();
            }
        factor_checked++;free(ref);
    }
    if(rc || m!=n || n<1 || n>LIMIT)return rc;
    int count=0;
    for(int i=0;i<n;i++) {
        lower[i]=count;
        for(int k=0;k<i;k++) {
            if(!isfinite(a[k][i]))return rc;
            if(a[k][i]!=0){indices[count]=k;values[count++]=a[k][i];}
        }
    }
    lower[n]=count;
    for(int i=0;i<n;i++) {
        upper[i]=count;
        if(!isfinite(a[i][i]))return rc;
        for(int k=n-1;k>i;k--) {
            if(!isfinite(a[k][i]))return rc;
            if(a[k][i]!=0){indices[count]=k;values[count++]=a[k][i];}
        }
    }
    upper[n]=count;matrix=a;pivots=p;dimension=n;factors++;
    return rc;
}
void denseGETRS(realtype **a,sunindextype n,sunindextype *p,realtype *b)
{
    if(!original_solve)original_solve=(solve_fn)symbol("denseGETRS");
    int eligible=matrix==a && pivots==p && n==dimension;
    if(eligible)for(int i=0;i<n;i++)if(!isfinite(b[i])){eligible=0;break;}
    if(!eligible){fallback++;original_solve(a,n,p,b);return;}
    realtype input[LIMIT];memcpy(input,b,n*sizeof(*b));
    int check=getenv("SNRT_CHIMES_JAC_VERIFY")!=NULL;
    for(int k=0;k<n;k++)if(p[k]!=k){realtype tmp=b[k];b[k]=b[p[k]];b[p[k]]=tmp;}
    /* Gather each row in exactly the subtraction order of upstream's
     * column scatter. The accumulator stays in a register. */
    for(int i=1;i<n;i++) {
        realtype sum=b[i];
        for(int j=lower[i];j<lower[i+1];j++)sum-=values[j]*b[indices[j]];
        b[i]=sum;
    }
    for(int i=n-1;i>=0;i--) {
        realtype sum=b[i];
        for(int j=upper[i];j<upper[i+1];j++)sum-=values[j]*b[indices[j]];
        b[i]=sum/a[i][i];
    }
    /* Overflow would make omitted zero*Inf terms observable: re-evaluate the
     * original path rather than silently changing nonfinite propagation. */
    int finite=1;for(int i=0;i<n;i++)if(!isfinite(b[i]))finite=0;
    if(!finite){memcpy(b,input,n*sizeof(*b));fallback++;original_solve(a,n,p,b);return;}
    if(check) {
        original_solve(a,n,p,input);
        for(int i=0;i<n;i++)if(b[i]!=input[i]) {
            fprintf(stderr,"SPARSE_TRIANGULAR_MISMATCH i=%d trial=%.17g ref=%.17g\n",i,b[i],input[i]);abort();
        }
        verified++;
    }
    solves++;entries+=upper[n];dense_entries+=n*(n-1);
}
__attribute__((destructor)) static void report(void)
{
    if(factors || solves)fprintf(stderr,"SPARSE_TRIANGULAR factors=%ld solves=%ld verified=%ld fallback=%ld entries=%ld dense_entries=%ld\n",
                                factors,solves,verified,fallback,entries,dense_entries);
    if(factor_checked || factor_sparse_updates)fprintf(stderr,"SPARSE_LU verified=%ld sparse_updates=%ld dense_updates=%ld\n",
        factor_checked,factor_sparse_updates,factor_dense_updates);
}
