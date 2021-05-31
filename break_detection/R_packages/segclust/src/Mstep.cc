#include "Mstep.h"

#include <cmath>
#include <iostream>
using namespace std;

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

#include "helperfuncs.h"
#include "repmat.h"

namespace segclust {

gsl_vector *
Mstep(gsl_vector *x, gsl_matrix *rupt, gsl_matrix *tau, gsl_vector *phi,
      bool vh)
{
  //K = nrow(tau)
  int K=tau->size1;

  //P = ncol(tau)
  int P=tau->size2;

  // m = matrix(nrow=P,ncol=1)
  gsl_vector *m=gsl_vector_calloc(P);
  
  //s = matrix(nrow=P,ncol=1)
  gsl_vector *s=gsl_vector_calloc(P);

  //prop = matrix(nrow=P,ncol=1)
  gsl_vector *prop=gsl_vector_calloc(P);

  //Yk = apply(rupt,1,FUN=function(y) sum(x[y[1]:y[2]]))
  gsl_vector *Yk=gsl_vector_calloc(rupt->size1);
  for (int r=0;r<rupt->size1;r++) {
    int y1=(int)gsl_matrix_get(rupt,r,1-1);
    int y2=(int)gsl_matrix_get(rupt,r,2-1);
    if (y1>y2) {
      int tmp_y=y1;
      y1=y2;
      y2=tmp_y;
    }
    double sum=0.0;
    for (int index=y1; index<=y2;index++)
      sum+=gsl_vector_get(x,index-1);
    gsl_vector_set(Yk,r,sum);
  }

  

  //nk = rupt[,2]-rupt[,1]+1
  gsl_vector *nk=gsl_vector_calloc(rupt->size1);

  gsl_vector *rupt2=gsl_vector_calloc(rupt->size1);
  gsl_matrix_get_col(rupt2,rupt,2-1);

  gsl_vector *rupt1=gsl_vector_calloc(rupt->size1);
  gsl_matrix_get_col(rupt1,rupt,1-1);

  gsl_vector_memcpy(nk,rupt2);
  gsl_vector_sub(nk,rupt1);
  gsl_vector_add_constant(nk,+1.0);

  gsl_vector_free(rupt1);
  gsl_vector_free(rupt2);

  //n  = sum(nk)
  double n=sum_gsl_vector(nk);

  //#homogeneous variances

  if (vh==true) {

    //for (p in (1:P)){
    for (int p=1;p<=P;p++) {
      
      // np    = nk %*% tau[,p]
      gsl_vector *taup=gsl_vector_calloc(tau->size1);
      gsl_matrix_get_col(taup,tau,p-1);
      double np=0.0;
      gsl_blas_ddot(nk,taup,&np);
      
      //m[p]  = Yk %*% tau[,p]/np
      gsl_vector *taup_div_np=gsl_vector_calloc(taup->size);
      gsl_vector_memcpy(taup_div_np,taup);
      gsl_vector_scale(taup_div_np,1.0/np);
      double prod=0.0;
      gsl_blas_ddot(Yk,taup_div_np,&prod);
      gsl_vector_free(taup_div_np);
      
      gsl_vector_set(m,p-1,prod);
      
      //tmp   = tau[,p] * ( apply(rupt,1,FUN=function(y) sum(    (x[y[1]:y[2]]-m[p])^2   )))  
      gsl_vector *res_apply=gsl_vector_calloc(taup->size);
      for (int r=0;r<rupt->size1;r++) {
	int y1=(int)gsl_matrix_get(rupt,r,1-1);
	int y2=(int)gsl_matrix_get(rupt,r,2-1);
	if (y1>y2) {
	  int tmp_y=y1;
	  y1=y2;
	  y2=tmp_y;
	}
	gsl_vector *diff_vec=copy_range_gsl_vector(x,y1-1,y2-1);
	gsl_vector_add_constant(diff_vec,-1.0*gsl_vector_get(m,p-1));
	gsl_vector_mul(diff_vec,diff_vec);
	double sum=sum_gsl_vector(diff_vec);
	gsl_vector_free(diff_vec);
	gsl_vector_set(res_apply,r,sum);
      }
      gsl_vector *tmp=gsl_vector_calloc(taup->size);
      gsl_vector_memcpy(tmp,taup);
      gsl_vector_mul(tmp,res_apply);
      gsl_vector_free(res_apply);
      
      //s[p]=sum(tmp)
      double sum=sum_gsl_vector(tmp);
      gsl_vector_free(tmp);
      gsl_vector_set(s,p-1,sum);
      
      gsl_vector_free(taup);
      //} for p
    }
    
    gsl_vector_free(Yk);
    gsl_vector_free(nk);
    
    
    //s = repmat(sum(s)/n,P,1)
    double sum=sum_gsl_vector(s)/(1.0*n);
    for (int index=0;index<s->size;index++)
      gsl_vector_set(s,index,sum);

  } else {
    
    gsl_vector *taup=0;
    gsl_vector *taup_div_np=0;
    gsl_vector *tmp=0;
    gsl_vector *ruptk1_k2=0;
    //for (p in (1:P)){
    for (int p=1;p<=P;p++) {
      
      // np    = nk %*% tau[,p]
      if (taup)
	gsl_vector_free(taup);
      taup=gsl_vector_calloc(tau->size1);
      gsl_matrix_get_col(taup,tau,p-1);
      double np=0.0;
      gsl_blas_ddot(nk,taup,&np);
      
      //m[p]  = Yk %*% tau[,p]/np
      if (taup_div_np)
	gsl_vector_free(taup_div_np);
      taup_div_np=gsl_vector_calloc(taup->size);
      gsl_vector_memcpy(taup_div_np,taup);
      gsl_vector_scale(taup_div_np,1.0/np);
      double prod=0.0;
      gsl_blas_ddot(Yk,taup_div_np,&prod);
      gsl_vector_free(taup_div_np);
      taup_div_np=0;

      gsl_vector_set(m,p-1,prod);

      //tmp = rep(Inf,K)
      if (tmp)
	gsl_vector_free(tmp);
      tmp=gsl_vector_calloc(K);
      gsl_vector_set_all(tmp,GSL_POSINF);

      //for (k in (1:K)) {
      for (int k=1;k<=K;k++) {

	// tmp[k] = tau[k,p] *  sum ( (x[(rupt[k,2]:rupt[k,1])]-m[p])^2)
	int ruptk2=(int)gsl_matrix_get(rupt,k-1,2-1);
	int ruptk1=(int)gsl_matrix_get(rupt,k-1,1-1);
	if (ruptk1_k2)
	  gsl_vector_free(ruptk1_k2);
	ruptk1_k2=copy_range_gsl_vector(x,ruptk2-1,ruptk1-1);
	gsl_vector_add_constant(ruptk1_k2,-1.0*gsl_vector_get(m,p-1));
	gsl_vector_mul(ruptk1_k2,ruptk1_k2);
	double sum=sum_gsl_vector(ruptk1_k2);
	gsl_vector_free(ruptk1_k2);
	ruptk1_k2=0;
	gsl_vector_set(tmp,k-1,gsl_matrix_get(tau,k-1,p-1)*sum);

      //}
      }

      // s[p] = sum(tmp)/np
      gsl_vector_set(s,p-1,sum_gsl_vector(tmp)/np);
      gsl_vector_free(tmp);
      tmp=0;
      //}
      
    }
      
  }





  //s = sqrt(s)
  gsl_vector *tmp_s=apply_basicfunc_gsl_vector(s,sqrt);
  gsl_vector_free(s);
  s=tmp_s;

  //prop = apply(tau,2,sum)/K
  gsl_vector *colsum=colsum_gsl_matrix(tau);
  gsl_vector_scale(colsum,1.0/K);
  gsl_vector_memcpy(prop,colsum);
  gsl_vector_free(colsum);

  /**
   * !!! WARNING !!!
   * Due to numerical precision limits, the sort order of m
   * in the C++ code can be different form the order in the R code.
  **/

  //b    = order(m)
  gsl_permutation *b=gsl_permutation_calloc(m->size);
  gsl_sort_vector_index(b,m);

  //m    = sort(m)
  gsl_sort_vector(m);

  //s    = s[b]
  tmp_s=gsl_vector_calloc(s->size);
  for (int index=0; index<s->size;index++)
    gsl_vector_set(tmp_s,index,gsl_vector_get(s,gsl_permutation_get(b,index)));
  gsl_vector_free(s);
  s=tmp_s;

  //prop = prop[b]
  gsl_vector *tmp_prop=gsl_vector_alloc(prop->size);
  for (int index=0;index<prop->size;index++)
    gsl_vector_set(tmp_prop,index,gsl_vector_get(prop,gsl_permutation_get(b,index)));
  gsl_vector_free(prop);
  prop=tmp_prop;

  gsl_permutation_free(b);

  // phi  = c(m,s,prop)
  gsl_vector *phi_res=gsl_vector_calloc(m->size+s->size+prop->size);
  for (int index=0;index<m->size;index++)
    gsl_vector_set(phi_res,index,gsl_vector_get(m,index));
  for (int index=0;index<s->size;index++)
    gsl_vector_set(phi_res,m->size+index,gsl_vector_get(s,index));
  for (int index=0;index<prop->size;index++)
    gsl_vector_set(phi_res,m->size+s->size+index,gsl_vector_get(prop,index));

  gsl_vector_free(prop);
  gsl_vector_free(s);
  gsl_vector_free(m);

  return phi_res;
}

}
