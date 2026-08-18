#ifdef AMR_INDEX_CHECK
use amr_index, only: icell_of, igrid_of, ichild_of
#endif
#ifndef AMR_INDEX_MACROS_DEFINED
#define AMR_INDEX_MACROS_DEFINED
#ifdef AMR_INDEX_CHECK
#  define ICELL_OF(g,c)   icell_of((g),(c))
#  define IGRID_OF(k)     igrid_of((k))
#  define ICHILD_OF(k)    ichild_of((k))
#else
#  define ICELL_OF(g,c)   (ncoarse+((c)-1)*ngridmax+(g))
#  define IGRID_OF(k)     (mod((k)-ncoarse-1,ngridmax)+1)
#  define ICHILD_OF(k)    (((k)-ncoarse-1)/ngridmax+1)
#endif
#endif
