#include "Estep.h"

#include <cmath>
#include <iostream>
using namespace std;

#include "helperfuncs.h"
#include "repmat.h"

namespace segclust {

void
Estep(gsl_matrix *logdensity, gsl_vector *phi,
      gsl_matrix **tau_res, double *lvinc_res)
{

  //K = nrow(logdensity)
  int K=logdensity->size1;

  //P = ncol(logdensity)
  int P=logdensity->size2;

  //tau     = repmat( t (log( phi[(2*P+1):(3*P)] )),K,1)+logdensity
  int cmin=2*P+1;
  int cmax=3*P;
  gsl_vector *logphi=copy_range_gsl_vector(phi,cmin-1,cmax-1);
  gsl_vector *logphi_tmp=apply_basicfunc_gsl_vector(logphi,log);
  gsl_vector_free(logphi);
  logphi=logphi_tmp;

  gsl_matrix *logphi_mat=gsl_matrix_alloc(1,logphi->size);
  gsl_matrix_set_row(logphi_mat,0,logphi);
  gsl_vector_free(logphi);

  gsl_matrix *logphi_repmat=repmat(logphi_mat,K,1);
  gsl_matrix_free(logphi_mat);

  gsl_matrix *tau=gsl_matrix_calloc(logphi_repmat->size1,logphi_repmat->size2);
  gsl_matrix_memcpy(tau,logphi_repmat);
  gsl_matrix_free(logphi_repmat);
  gsl_matrix_add(tau,logdensity);


  //tau_max = apply(tau,1,max)
  gsl_vector *tau_max=rowmax_gsl_matrix(tau);

  //tau     = exp(tau-repmat(tau_max,1,P))
  gsl_matrix *tau_max_mat=gsl_matrix_calloc(tau_max->size,1);
  gsl_matrix_set_col(tau_max_mat,0,tau_max);
  gsl_matrix *tau_max_repmat=repmat(tau_max_mat,1,P);
  gsl_matrix_free(tau_max_mat);
  gsl_matrix_sub(tau,tau_max_repmat);
  gsl_matrix_free(tau_max_repmat);
  gsl_matrix *tau_tmp=apply_basicfunc_gsl_matrix(tau,exp);
  gsl_matrix_free(tau);
  tau=tau_tmp;

  //lvinc  = sum(log( apply(tau,1,sum)) + tau_max)
  gsl_vector *tau_rowsum=rowsum_gsl_matrix(tau);
  gsl_vector *tau_rowsumtmp=apply_basicfunc_gsl_vector(tau_rowsum,log);
  gsl_vector_free(tau_rowsum);
  tau_rowsum=tau_rowsumtmp;
  gsl_vector_add(tau_rowsum,tau_max);
  double lvinc=sum_gsl_vector(tau_rowsum);
  gsl_vector_free(tau_max);
  gsl_vector_free(tau_rowsum);

  //tau     = tau / repmat( apply(tau,1,sum) ,1,P)
  tau_rowsum=rowsum_gsl_matrix(tau);
  gsl_matrix *tau_rowsum_mat=gsl_matrix_calloc(tau_rowsum->size,1);
  gsl_matrix_set_col(tau_rowsum_mat,0,tau_rowsum);
  gsl_vector_free(tau_rowsum);
  gsl_matrix *tau_rowsum_repmat=repmat(tau_rowsum_mat,1,P);
  gsl_matrix_free(tau_rowsum_mat);
  gsl_matrix_div_elements(tau,tau_rowsum_repmat);
  gsl_matrix_free(tau_rowsum_repmat);
  
  //invisible(list(tau,lvinc))
  if (*tau_res)
    gsl_matrix_free(*tau_res);
  *tau_res=tau;

  *lvinc_res=lvinc;
}

}
