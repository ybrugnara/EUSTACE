#include "neighbors.h"

#include <cstdio>
#include <iostream>
using namespace std;

#include <gsl/gsl_math.h>



#include "EM_algo.h"
#include "Segmentation_mixt.h" 
#include "helperfuncs.h"

namespace segclust {

void
neighbors(gsl_vector *x, long Kmax,
	  gsl_vector *L, long k, 
	  param_struct param[], long param_size,
	  long P, long lmin, long lmax, bool vh,
	  gsl_vector **L_res, 
	  param_struct *param_res[], long *param_res_size)
{

  if (*param_res) {
    for (int index=0;index<*param_res_size;index++) {
      if ((*param_res)[index].phi)
	gsl_vector_free((*param_res)[index].phi);
      if ((*param_res)[index].rupt) {
	gsl_matrix_free((*param_res)[index].rupt);
      }
    }
    delete[] *param_res;
  }

  *param_res=new param_struct[param_size];
  *param_res_size=param_size;
  for (int index=0;index<param_size;index++) {
    gsl_vector *phi=0;
    if (param[index].phi) {
      phi=gsl_vector_calloc(param[index].phi->size);
      gsl_vector_memcpy(phi,param[index].phi);
    }
    (*param_res)[index].phi=phi;

    gsl_matrix *rupt=0;
    if (param[index].rupt) {
      rupt=gsl_matrix_calloc(param[index].rupt->size1,
			     param[index].rupt->size2);
      gsl_matrix_memcpy(rupt,param[index].rupt);
    }
    (*param_res)[index].rupt=rupt;
  }

  if (*L_res)
    gsl_vector_free(*L_res);
  *L_res=gsl_vector_calloc(L->size);
  gsl_vector_memcpy(*L_res,L);



  //# left neighbor
  //a  = max(L[1:(k-1)])
  //K1 = which.max(L[1:(k-1)])
  double a=GSL_NAN;
  int K1=-1;
  //  fprintf(stdout,"k=%d\n",k);
  if (1<=k-1) {
    gsl_vector *L_range=copy_range_gsl_vector(L,1-1,k-1-1);
    a=checked_max_gsl_vector(L_range);
    //    fprintf(stdout,"L_range\n");
    //    gsl_vector_fprintf(stdout,L_range,"%lf");
    K1=checked_max_index_gsl_vector(L_range)+1;
    gsl_vector_free(L_range);
  }

  //  fprintf(stdout,"a=%lf\nK1=%d\n",a,K1);
  gsl_vector *phi1=0;

  gsl_matrix *rupt1=0;
  gsl_vector *out_EM1_phi=0;
  gsl_matrix *out_EM1_tau=0;
  double out_EM1_lvinc=0.0;
  int out_EM1_empty=0;
  int out_EM1_dv=0;

  //if ( (is.null(K1)) || (a==-Inf) || (is.na(a) ) ){
  if ((K1<=0) || (gsl_isinf(a)==-1) || gsl_isnan(a)) {

    //  phi1          = rep(-Inf,3*P)
    phi1=gsl_vector_calloc(3*P);
    gsl_vector_set_all(phi1,GSL_NEGINF);

    //  out.EM1$lvinc = list(lvinc=- Inf)
    out_EM1_lvinc=GSL_NEGINF;

  //}    else { 
  } else {

    //phi1                     = param[[K1]]$phi     
    if (phi1)
      gsl_vector_free(phi1);
    phi1=gsl_vector_calloc(param[K1-1].phi->size);
    gsl_vector_memcpy(phi1,param[K1-1].phi);
    
    //G                        = Gmixt(x,lmin=lmin,phi1,P)
    gsl_matrix *t_est=0;

    compute_segmentation_mixt(x,k,lmin,lmax,phi1,P,NULL,&t_est);

    //rupt1      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    if (rupt1) {
      gsl_matrix_free(rupt1);
    }
    rupt1=gsl_matrix_calloc(t_est->size2,2);
    gsl_matrix_set(rupt1,0,0,1.0);
    for (int index=1;index<=k-1;index++)
      gsl_matrix_set(rupt1,index,0,gsl_matrix_get(t_est,k-1,index-1)+1.0);
    gsl_vector *rowvec=gsl_vector_calloc(t_est->size2);
    gsl_matrix_get_row(rowvec,t_est,k-1);
    gsl_matrix_free(t_est);
    gsl_matrix_set_col(rupt1,1,rowvec);
    gsl_vector_free(rowvec);

    //out.EM1    = EM.algo(x,rupt1,P,phi1,vh)
    EM_algo(x,rupt1,P,phi1,vh,&out_EM1_phi,&out_EM1_tau,&out_EM1_lvinc,
	    &out_EM1_empty,&out_EM1_dv);
    
    gsl_matrix_free(out_EM1_tau);

    // } #end else
  }
  gsl_vector_free(phi1);
  phi1=0;

  //# right neighbor 
  //a  = max(L[(k+1):Kmax])  
  //K2 = which.max(L[(k+1):Kmax])
  //K2 = K2 + k
  a=GSL_NAN;
  int K2=-1;
  if (k+1<=Kmax) {
    gsl_vector *L_range=copy_range_gsl_vector(L,k+1-1,Kmax-1);
    a=checked_max_gsl_vector(L_range);

    K2=checked_max_index_gsl_vector(L_range)+1;
    gsl_vector_free(L_range);

    K2=K2+k;
  }

  gsl_vector *phi2=0;
  gsl_matrix *rupt2=0;
  gsl_vector *out_EM2_phi=0;
  gsl_matrix *out_EM2_tau=0;
  double out_EM2_lvinc=0.0;
  int out_EM2_empty=0;
  int out_EM2_dv=0;


  //if ( (K2==0) || (a==-Inf) || (is.na(a)) ){
  if ((K2<=0) || (gsl_isinf(a)==-1) || (gsl_isnan(a))) {




    //  phi2          = rep(-Inf,3*P)
    phi2=gsl_vector_calloc(3*P);
    gsl_vector_set_all(phi2,GSL_NEGINF);

    //  out.EM2       = list(lvinc=-Inf)
    out_EM2_lvinc=GSL_NEGINF;

    //}   else {
  } else {

    //phi2                     = param[[K2]]$phi
    phi2=gsl_vector_calloc(param[K2-1].phi->size);
    gsl_vector_memcpy(phi2,param[K2-1].phi);

    //G                        = Gmixt(x,lmin=lmin,phi2,P)
    gsl_matrix *t_est=0;

    compute_segmentation_mixt(x,k,lmin,lmax,phi2,P,NULL,&t_est);


    //rupt2      = matrix(ncol=2,c(c(1,t.est[k,1:(k-1)]+1),t.est[k,]))
    if (rupt2)
      gsl_matrix_free(rupt2);
    rupt2=gsl_matrix_calloc(t_est->size2,2);
    gsl_matrix_set(rupt2,0,0,1.0);
    for (int index=1;index<=k-1;index++)
      gsl_matrix_set(rupt2,index,0,gsl_matrix_get(t_est,k-1,index-1)+1.0);
    gsl_vector *rowvec=gsl_vector_calloc(t_est->size2);
    gsl_matrix_get_row(rowvec,t_est,k-1);
    gsl_matrix_free(t_est);
    gsl_matrix_set_col(rupt2,1,rowvec);
    gsl_vector_free(rowvec);
    
    //out.EM2    = EM.algo(x,rupt2,P,phi2,vh)
    EM_algo(x,rupt2,P,phi2,vh,&out_EM2_phi,&out_EM2_tau,&out_EM2_lvinc,
	    &out_EM2_empty,&out_EM2_dv);
    
    gsl_matrix_free(out_EM2_tau);
    out_EM2_tau=0;
  
    //} #end else
  }

  gsl_vector_free(phi2);
  phi2=0;
  

    
  //# choice of the best likelihood
  //a = which.max( c(out.EM1$lvinc, out.EM2$lvinc,  L[k]) ) 
  gsl_vector *maxvec=gsl_vector_calloc(3);
  gsl_vector_set(maxvec,0,out_EM1_lvinc);
  gsl_vector_set(maxvec,1,out_EM2_lvinc);
  gsl_vector_set(maxvec,2,gsl_vector_get(L,k-1));
  int aindex=checked_max_index_gsl_vector(maxvec)+1;
  //  fprintf(stdout,"maxvec=\n");
  //  gsl_vector_fprintf(stdout,maxvec,"%lf");
  //  fprintf(stdout,"max=%d\n",aindex);
    gsl_vector_free(maxvec);
  

  //# parameter update
  //if (length(a)==0) {
  if (aindex==-1) {

    //      L[k] = L
    gsl_vector_set(*L_res,k-1,gsl_vector_get(*L_res,0));

    //     param[[k]] = param
    if ((*param_res)[k-1].phi)
      gsl_vector_free((*param_res)[k-1].phi);
    (*param_res)[k-1].phi=0;
    if ((*param_res)[0].phi) {
      (*param_res)[k-1].phi=gsl_vector_calloc((*param_res)[0].phi->size);
      gsl_vector_memcpy((*param_res)[k-1].phi,(*param_res)[0].phi);
    }

    if ((*param_res)[k-1].rupt)
      gsl_matrix_free((*param_res)[k-1].rupt);
    (*param_res)[k-1].rupt=0;
    if ((*param_res)[0].rupt) {
      (*param_res)[k-1].rupt=gsl_matrix_calloc((*param_res)[0].rupt->size1,
					       (*param_res)[0].rupt->size2);
      gsl_matrix_memcpy((*param_res)[k-1].rupt,(*param_res)[0].rupt);
    }

    if (out_EM1_phi)
      gsl_vector_free(out_EM1_phi);
    if (out_EM2_phi)
      gsl_vector_free(out_EM2_phi);
    if (rupt1)
      gsl_matrix_free(rupt1);
    rupt1=0;
    if (rupt2)
      gsl_matrix_free(rupt2);
    rupt2=0;

  //} else {
  } else {
  //if (a==1){
    if (aindex==1) {

      //  param[[k]] = list(phi = out.EM1$phi, rupt = rupt1)
      if ((*param_res)[k-1].phi)
	gsl_vector_free((*param_res)[k-1].phi);
      (*param_res)[k-1].phi=out_EM1_phi;
      if (out_EM2_phi)
	gsl_vector_free(out_EM2_phi);
      out_EM2_phi=0;

      if ((*param_res)[k-1].rupt)
	gsl_matrix_free((*param_res)[k-1].rupt);
      (*param_res)[k-1].rupt=rupt1;
      if (rupt2)
	gsl_matrix_free(rupt2);
      rupt2=0;
      
      //  L[k] = out.EM1$lvinc}
      gsl_vector_set(*L_res,k-1,out_EM1_lvinc);
      
      
    }
    //if (a==2) {
    if (aindex==2) {
      //  param[[k]] = list(phi = out.EM2$phi,rupt = rupt2)
      if ((*param_res)[k-1].phi)
	gsl_vector_free((*param_res)[k-1].phi);
      (*param_res)[k-1].phi=out_EM2_phi;
      if (out_EM1_phi)
	gsl_vector_free(out_EM1_phi);
      out_EM1_phi=0;

      if ((*param_res)[k-1].rupt)
	gsl_matrix_free((*param_res)[k-1].rupt);
      (*param_res)[k-1].rupt=rupt2;
      if(rupt1)
	gsl_matrix_free(rupt1);
      rupt1=0;
      //  L[k] = out.EM2$lvinc}
      gsl_vector_set(*L_res,k-1,out_EM2_lvinc);
    }

    if (aindex==3) {
      if (out_EM1_phi)
	gsl_vector_free(out_EM1_phi);
      if (out_EM2_phi)
	gsl_vector_free(out_EM2_phi);
      if (rupt1)
	gsl_matrix_free(rupt1);
      rupt1=0;
      if (rupt2)
	gsl_matrix_free(rupt2);
      rupt2=0;
    }

  }
// invisible(list(L=L,param=param))  
    
}

}
