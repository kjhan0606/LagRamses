#ifdef AMR_INDEX_CHECK
use amr_index, only: icell_of, igrid_of, ichild_of
#else
use amr_parameters, only: twotondim, amr_block_size
#endif
#ifndef AMR_INDEX_MACROS_DEFINED
#define AMR_INDEX_MACROS_DEFINED
#ifdef AMR_INDEX_CHECK
#  define ICELL_OF(g,c)   icell_of((g),(c))
#  define IGRID_OF(k)     igrid_of((k))
#  define ICHILD_OF(k)    ichild_of((k))
#else
#  define ICELL_OF(g,c)   (ncoarse+(((g)-1)/amr_block_size)*(twotondim*amr_block_size)+((c)-1)*amr_block_size+mod((g)-1,amr_block_size)+1)
#  define IGRID_OF(k)     (((((k)-ncoarse-1)/(twotondim*amr_block_size))*amr_block_size)+mod(mod((k)-ncoarse-1,twotondim*amr_block_size),amr_block_size)+1)
#  define ICHILD_OF(k)    ((mod((k)-ncoarse-1,(twotondim*amr_block_size))/amr_block_size)+1)
#endif
#endif
