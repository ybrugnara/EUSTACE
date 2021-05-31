#include "EM_algo.h"

#include <cmath>
#include <iostream>
using namespace std;

#include <gsl/gsl_math.h>

#include "Estep.h"
#include "Mstep.h"
#include "helperfuncs.h"
#include "logdens.h"

namespace segclust {

void
EM_algo(gsl_vector *x, gsl_matrix *rupt, int P, gsl_vector *phi0,
	bool vh,
	gsl_vector **phi_res, gsl_matrix **tau_res,
	double *lvinc_res, int *empty_res, int *dv_res)
{
  //K     = nrow(rupt)
  int K = rupt->size1;

  //delta = 1
  double delta=1.0;

  //empty = 0
  int empty=0;

  //dv    = 0
  int dv=0;

  //tau   = matrix(1,nrow = K,ncol = P)
  gsl_matrix *tau=gsl_matrix_calloc(K,P);
  gsl_matrix_set_all(tau,1.0);

  //phi   = phi0
  gsl_vector *phi=gsl_vector_calloc(phi0->size);
  gsl_vector_memcpy(phi,phi0);

  //iter  = 0
  int iter=0;

  //np    = apply(tau,2,sum)
  gsl_vector *np=colsum_gsl_matrix(tau);


  //eps   = 10e-10
  double eps=10e-4;
  


  double lvinc=0.0;


  gsl_vector *phi_temp=0;
  gsl_matrix *logdensity=0;
  gsl_vector *x_tmp=0;
  gsl_matrix *logdens_mat=0;
  gsl_vector *logdens_mat_row=0;
  gsl_vector *phi_tmp=0;
  gsl_vector *absdiff=0;
  gsl_vector *absdiff_tmp=0;
  //while ( (delta>=1e-4) & (min(np)>eps) & (iter<=5000) ){
  while ((delta>=1e-4) && (checked_min_gsl_vector(np)>eps) && (iter<=5000)) {


    //iter       = iter+1
    iter++;

    //phi_temp   = phi
    if (phi_temp)
      gsl_vector_free(phi_temp);
    phi_temp=gsl_vector_alloc(phi->size);
    gsl_vector_memcpy(phi_temp,phi);


    //logdensity = t( apply(rupt,1,FUN=function(y) logdens(   x[ y[1]:y[2] ] ,P,phi)))
    if (logdensity)
      gsl_matrix_free(logdensity);
    logdensity=gsl_matrix_calloc(rupt->size1,P);
    for (int r=0;r<rupt->size1;r++) {
      int y1=(int)gsl_matrix_get(rupt,r,1-1);
      int y2=(int)gsl_matrix_get(rupt,r,2-1);
      if (x_tmp)
	gsl_vector_free(x_tmp);
      x_tmp=copy_range_gsl_vector(x,y1-1,y2-1);
      if (logdens_mat)
	gsl_matrix_free(logdens_mat);
      logdens_mat=logdens(x_tmp,P,phi);
      gsl_vector_free(x_tmp);
      x_tmp=0;
      if (logdens_mat_row)
	gsl_vector_free(logdens_mat_row);
      logdens_mat_row=gsl_vector_calloc(P);
      gsl_matrix_get_row(logdens_mat_row,logdens_mat,0);
      gsl_matrix_free(logdens_mat);
      logdens_mat=0;
      gsl_matrix_set_row(logdensity,r,logdens_mat_row);
      gsl_vector_free(logdens_mat_row);
      logdens_mat_row=0;
      
    }

    //Estepout   = Estep(logdensity,phi)
    //tau        = Estepout[[1]]
    //lvinc      = Estepout[[2]]
    Estep(logdensity,phi,&tau,&lvinc);
    gsl_matrix_free(logdensity);
    logdensity=0;
    
    //phi        = Mstep(x,rupt,tau,phi)
    if (phi_tmp)
      gsl_vector_free(phi_tmp);
    phi_tmp=Mstep(x,rupt,tau,phi,vh);
    gsl_vector_free(phi);
    phi=phi_tmp;
    phi_tmp=0;

    //np         = apply(tau,2,sum)
    gsl_vector_free(np);
    np=colsum_gsl_matrix(tau);
    
    // delta      = max(abs(phi_temp-phi)/phi)
    if (absdiff)
      gsl_vector_free(absdiff);
    absdiff=gsl_vector_alloc(phi->size);
    gsl_vector_memcpy(absdiff,phi_temp);
    gsl_vector_free(phi_temp);
    phi_temp=0;
    gsl_vector_sub(absdiff,phi);
    if (absdiff_tmp)
      gsl_vector_free(absdiff_tmp);
    absdiff_tmp=apply_basicfunc_gsl_vector(absdiff,fabs);
    gsl_vector_free(absdiff);
    absdiff=absdiff_tmp;
    absdiff_tmp=0;
    gsl_vector_div(absdiff,phi);
    delta=checked_max_gsl_vector(absdiff);
    gsl_vector_free(absdiff);
    absdiff=0;

    //}
  }


  //if (min(np)<eps){
  if (checked_min_gsl_vector(np)<eps) {
    //empty = 1
    empty=1;

    //lvinc = -Inf
    lvinc=GSL_NEGINF;
    
  //}
  }

  //if (iter>5000){
  if (iter>5000) {

    //dv    = 2
    dv=2;
    //lvinc = -Inf
    lvinc=GSL_NEGINF;

  //}
  }

  gsl_vector_free(np);
  np=0;

  //rm(delta,logdensity)
  
  //invisible(list(phi = phi,tau = tau,lvinc = lvinc,empty = empty,dv = dv))
  if (*phi_res)
    gsl_vector_free(*phi_res);
  *phi_res=phi;
  
  if (*tau_res)
    gsl_matrix_free(*tau_res);
  *tau_res=tau;

  *lvinc_res=lvinc;
  *empty_res=empty;
  *dv_res=dv;

}

}
