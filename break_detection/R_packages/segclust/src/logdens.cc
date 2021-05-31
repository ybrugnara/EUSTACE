#include <cmath>
using namespace std;

#include "logdens.h"

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375
#endif 

namespace segclust {

gsl_matrix *
logdens(gsl_vector *xk, int P, gsl_vector *phi) {

  static double sqrt_2_pi=sqrt(2*M_PI);

  int i,p;
  double ld,sp,mp,sum,sumterm;

  int m_index=0;
  int s_index=P;
  int nk=xk->size;

  gsl_matrix *tmp=gsl_matrix_alloc(1,P);
  
  for (p=0;p<P;p++) {
    sp=gsl_vector_get(phi,s_index+p);
    mp=gsl_vector_get(phi,m_index+p);

    sum=0.0;
    for (i=0;i<xk->size;i++) {
      sumterm=gsl_vector_get(xk,i)-mp;
      sum+=sumterm*sumterm;
    }
    ld=-nk*log(sqrt_2_pi*sp)-0.5*sum/(sp*sp);

    gsl_matrix_set(tmp,0,p,ld);
  }
  
  return tmp;
}


}
