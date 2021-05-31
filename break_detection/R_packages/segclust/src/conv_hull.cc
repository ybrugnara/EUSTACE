#include "conv_hull.h"

#include "helperfuncs.h"

namespace segclust {

gsl_vector *
conv_hull(gsl_vector *J, gsl_vector *pen)
{


  //K = length(J)
  long K=J->size;

  //k = 1
  long k=1;

  //kv = c()
  gsl_vector *kv=gsl_vector_calloc(K);
  long kv_next=0;

  //pv = c()

  //  gsl_vector *pv=gsl_vector_calloc(K);
  //long pv_next=0;

  //while (k<K){
  while (k<K) {
    //pk = (J[(k+1):K]-J[k]) / (pen[k]-pen[(k+1):K])

    gsl_vector *Jrange=copy_range_gsl_vector(J,k+1-1,K-1);
    gsl_vector_add_constant(Jrange,-1.0*gsl_vector_get(J,k-1));

    gsl_vector *penrange=copy_range_gsl_vector(pen,k+1-1,K-1);
    gsl_vector_scale(penrange,-1.0);
    gsl_vector_add_constant(penrange,gsl_vector_get(pen,k-1));

    gsl_vector *pk=gsl_vector_alloc(Jrange->size);
    gsl_vector_memcpy(pk,Jrange);
    gsl_vector_free(Jrange);
    gsl_vector_div(pk,penrange);
    gsl_vector_free(penrange);

    //dm = which.max(pk)
    long dm=checked_max_index_gsl_vector(pk)+1;
    gsl_vector_free(pk);

    //kv = c(kv,k)
    gsl_vector_set(kv,kv_next,k);
    kv_next++;

    //k  = k + dm 
    k=k+dm;

    //}
  }
 
  //kv = c(kv,K)
  gsl_vector_set(kv,kv_next,K);
  kv_next++;

  //invisible(kv=kv)
  gsl_vector *kv_ret=gsl_vector_alloc(kv_next);
  for (int i=0;i<kv_next;i++)
    gsl_vector_set(kv_ret,i,gsl_vector_get(kv,i));

  gsl_vector_free(kv);

  return kv_ret;

}

}
