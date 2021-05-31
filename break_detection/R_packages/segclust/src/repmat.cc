#include "repmat.h"

namespace segclust {

gsl_matrix *
repmat(gsl_matrix *a, int n, int m)
{
  gsl_matrix *res=gsl_matrix_alloc(n*a->size1,m*a->size2);
  int i,j,r,c;

  for (i=0;i<n;i++) 
    for (j=0;j<m;j++)
      for (r=0;r<a->size1;r++)
	for (c=0;c<a->size2;c++)
	  gsl_matrix_set(res,i*a->size1+r,j*a->size2+c,gsl_matrix_get(a,r,c));

  return res;
}

}
