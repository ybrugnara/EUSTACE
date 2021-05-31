#include "hybrid.h"

#include <cassert>
#include <cstdio>
#include <iostream>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_vector.h>

#include "EM_algo.h"
#include "EM_init.h"
#include "Segmentation_mean.h" 
#include "Segmentation_mixt.h" 
#include "conv_hull.h"
#include "helperfuncs.h"
#include "logdens.h"
#include "neighbors.h"

namespace segclust {

void
hybrid(gsl_vector *x, long P, long Kmax, long lmin, long lmax,
       bool vh,  gsl_matrix **Linc_res, param_struct *param_res[],
       long *param_res_size)
{

  //Linc  = matrix(-Inf,nrow=Kmax,ncol=1)
  gsl_matrix *Linc=gsl_matrix_calloc(Kmax,1);
  gsl_matrix_set_all(Linc,GSL_NEGINF);


  //n     = length(x)
  long n=x->size;

  // param = list()
  long param_size=Kmax;
  param_struct *param=new param_struct[param_size];
  for (int index=0;index<param_size;index++) {
    param[index].phi=0;
    param[index].rupt=0;
  }
    
  long Kmin=P;

  // if (P==1){
  if (P == 1) {
    
    
    //   rupt               = c(1, length(x))
    gsl_matrix *rupt=gsl_matrix_calloc(1,2);
    gsl_matrix_set(rupt,0,0,1.0);
    gsl_matrix_set(rupt,0,1,x->size);

    //   phi                = c(mean(x),sd(x),1)
    gsl_vector *phi=gsl_vector_calloc(3);
    gsl_vector_set(phi,0,mean_gsl_vector(x));
    gsl_vector_set(phi,1,stddev_gsl_vector(x));
    gsl_vector_set(phi,2,1.0);

    //   Linc[Kmin:Kmax]    = logdens(x,P,phi)
    // P==1 -> logdensity returns a 1x1 matrix.
    gsl_matrix *tmp=logdens(x,P,phi);
    for (long Kindex=Kmin;Kindex<=Kmax;Kindex++) {
      gsl_matrix_set(Linc,Kindex-1,0,gsl_matrix_get(tmp,0,0));
    }
    gsl_matrix_free(tmp);
    tmp=0;
    
    //   param[[1]]         = list(phi=phi,rupt=rupt)

    //    delete[] param;
    //    param=new param_struct[1];
    param[0].phi=phi;
    param[0].rupt=rupt;

    // } else {
  } else {

    cout << "estimation for P = " << P << endl;  
    gsl_matrix *t_est=0;

    gsl_vector *phi_temp=0;
    gsl_vector *phi_temp2=0;
    gsl_vector *tmp_th=0;
    gsl_vector *th=0;

    gsl_matrix *rupt=0;
    gsl_vector *tmpvec=0;

    gsl_vector *phi=0;

    gsl_vector *phi_res=0;
    gsl_matrix *tau_res=0;
    //   for (K in (Kmin:Kmax)){
    for (long K=Kmin;K<=Kmax;K++) {
    
      //     cat("P",P,"K", K ,"\n")
      //cout << "P" << P << "K" << K << endl;

      //     j      = 0
      long j=0;

      //     delta  = Inf
      double delta=GSL_POSINF;

      //     empty  = 0
      int empty=0;


      //     dv     = 0
      int dv=0;

	
      if (t_est)
	gsl_matrix_free(t_est);
      t_est=0;
    
      compute_segmentation_mean(x,K,lmin,lmax,vh,NULL,&t_est);
  
      //     th     = t.est[K,1:K]
      if (tmp_th)
	gsl_vector_free(tmp_th);
      tmp_th=gsl_vector_calloc(t_est->size2);
      gsl_matrix_get_row(tmp_th,t_est,K-1);
      if (th)
	gsl_vector_free(th);
      th=copy_range_gsl_vector(tmp_th,1-1,K-1);
      gsl_vector_free(tmp_th);
      tmp_th=0;
      if (t_est)
	gsl_matrix_free(t_est);
      t_est=0;

      //     rupt   = matrix(ncol=2,c(c(1,th[1:(K-1)]+1),th)) 
      if (rupt)
	gsl_matrix_free(rupt);
      rupt=gsl_matrix_calloc(th->size,2);
      if (tmpvec)
	gsl_vector_free(tmpvec);
      tmpvec=gsl_vector_calloc(th->size);
      gsl_vector_set(tmpvec,0,1.0);
      for (int Kindex=1;Kindex<=K-1;Kindex++)
	gsl_vector_set(tmpvec,Kindex,gsl_vector_get(th,Kindex-1)+1.0);
      gsl_matrix_set_col(rupt,0,tmpvec);
      gsl_vector_free(tmpvec);
      tmpvec=0;
      gsl_matrix_set_col(rupt,1,th);
      gsl_vector_free(th);
      th=0;

      //     phi    = EM.init(x,rupt,K,P,vh)
      if (phi)
	gsl_vector_free(phi);
      phi=EM_init(x,rupt,K,P,vh);

      //     out.EM = EM.algo(x,rupt,P,phi,vh)
      if (phi_res)
	gsl_vector_free(phi_res);
      phi_res=0;
      if (tau_res)
	gsl_matrix_free(tau_res);
      tau_res=0;
      double lvinc=0.0;
      EM_algo(x,rupt,P,phi,vh,&phi_res,&tau_res,&lvinc,&empty,&dv);
      if (tau_res)
	gsl_matrix_free(tau_res);
      tau_res=0;
      
      //     phi    = out.EM$phi
      gsl_vector_free(phi);
      phi=phi_res;
      phi_res=0;

      double lvinc_mixt=0.0;


      //     while ( (delta>1e-4) & (empty==0) & (dv==0) & (j<=100)){           
      while ( (delta>1e-4) && (empty==0) && (dv==0) & (j<=100)) {
      
	//       j          = j+1 
	j++;
	//       phi.temp   = phi
	if (phi_temp)
	  gsl_vector_free(phi_temp);
	phi_temp=gsl_vector_calloc(phi->size);
	gsl_vector_memcpy(phi_temp,phi);

	//       G          = Gmixt(x,lmin=lmin,phi,P)
	compute_segmentation_mixt(x,K,lmin,lmax,phi,P,NULL,&t_est);

	
	//       rupt       = matrix(ncol=2,c(c(1,t.est[K,1:(K-1)]+1),t.est[K,]))
	if (rupt)
	  gsl_matrix_free(rupt);
	rupt=gsl_matrix_calloc(t_est->size2,2);
	gsl_matrix_set(rupt,0,0,1.0);
	for (int Kindex=1;Kindex<=K-1;Kindex++)
	  gsl_matrix_set(rupt,Kindex,0,gsl_matrix_get(t_est,K-1,Kindex-1)+1.0);
	for (int r=0;r<rupt->size1;r++)
	  gsl_matrix_set(rupt,r,1,gsl_matrix_get(t_est,K-1,r));
	gsl_matrix_free(t_est);
	  t_est=0;
	

	//       out.EM     = EM.algo(x,rupt,P,phi.temp,vh)
	//       phi        = out.EM$phi
	//       empty      = out.EM$empty
	//       dv         = out.EM$dv
	//       lvinc.mixt = out.EM$lvinc
	EM_algo(x,rupt,P,phi_temp,vh,&phi,&tau_res,&lvinc_mixt,&empty,&dv);
	if (tau_res)
	  gsl_matrix_free(tau_res);
	tau_res=0;
	//       delta      = max(abs(phi.temp-phi)/phi)
	gsl_vector_sub(phi_temp,phi);
	if (phi_temp2)
	  gsl_vector_free(phi_temp2);
	phi_temp2=apply_basicfunc_gsl_vector(phi_temp,fabs);
	gsl_vector_free(phi_temp);
	phi_temp=phi_temp2;
	phi_temp2=0;
	gsl_vector_div(phi_temp,phi);
	
	delta=checked_max_gsl_vector(phi_temp);

	gsl_vector_free(phi_temp);
	phi_temp=0;

	//     }#end while
      }
      
      
      //     Linc[K]=lvinc.mixt
      gsl_matrix_set(Linc,K-1,0,lvinc_mixt);

      //     param[[K]] = list(phi=phi,rupt=rupt)
      param[K-1].phi=phi;
      phi=0;
      param[K-1].rupt=rupt;
      rupt=0;
      
      //   } #end K
    }
 
    //   Ltmp= rep(-Inf,Kmax)
    gsl_matrix *Ltmp=gsl_matrix_calloc(Kmax,1);
    gsl_matrix_set_all(Ltmp,GSL_NEGINF);

    //   cat("tracking local maxima for P =",P,"\n")
    //cout << "tracking local maxima for P = " << P << endl;

    //    fprintf(stdout,"tracking local maxima for %d\n",P);


    gsl_matrix *L_inc=0;
    gsl_vector *Linc_col=0;
    gsl_vector *Linc_range=0;

    gsl_vector *L_res=0;
    gsl_vector *Linc_vec=0;

    gsl_vector *kvfinite=0;
    gsl_vector *Lfinite=0;
    gsl_vector *a=0;
    gsl_vector *a_tmp=0;
    gsl_vector *Kconc=0;
    gsl_vector *Kconc_keep=0;
    gsl_vector *Kconc_tmp=0;




    //   while (sum(Ltmp!=Linc)>=1) {
    while (sum_diffs_gsl_matrices(Ltmp,Linc) >=1) {

      //     # find the convex hull of the likelihood
      //     Ltmp     = Linc  
      gsl_matrix_free(Ltmp);
      Ltmp=gsl_matrix_calloc(Linc->size1,Linc->size2);
      gsl_matrix_memcpy(Ltmp,Linc);

      //     kvfinite = which(is.finite(Linc[P:Kmax]))+P-1
      if (Linc_col)
	gsl_vector_free(Linc_col);
      Linc_col=gsl_vector_calloc(Linc->size1);
      gsl_matrix_get_col(Linc_col,Linc,0);
      if (Linc_range)
	gsl_vector_free(Linc_range);
      Linc_range=copy_range_gsl_vector(Linc_col,P-1,Kmax-1);
      gsl_vector_free(Linc_col);
      Linc_col=0;
      int nfinite=0;
      for (int index=0;index<Linc_range->size;index++) {
	if (gsl_finite(gsl_vector_get(Linc_range,index)))
	  nfinite++;
      }

      if (kvfinite)
	gsl_vector_free(kvfinite);
      kvfinite=gsl_vector_calloc(nfinite);
      
      for (int index=0,offset=0;index<Linc_range->size;index++) 
	if (gsl_finite(gsl_vector_get(Linc_range,index))) {
	  // kvfinite is R-indexed.
	  gsl_vector_set(kvfinite,offset,index+1);
	  offset++;
	}
      gsl_vector_add_constant(kvfinite,P-1.0);
      gsl_vector_free(Linc_range);
      Linc_range=0;
      

      //     Lfinite  = Linc[kvfinite]
      if (Lfinite)
	gsl_vector_free(Lfinite);
      Lfinite=gsl_vector_calloc(kvfinite->size);
      for (int index=0;index<kvfinite->size;index++)
	gsl_vector_set(Lfinite,index,gsl_matrix_get(Linc,(int)gsl_vector_get(kvfinite,index)-1,0));
      

      //     a        = conv.hull(-Lfinite,kvfinite)
      gsl_vector_scale(Lfinite,-1.0);
      if (a)
	gsl_vector_free(a);
      a=conv_hull(Lfinite,kvfinite);
      gsl_vector_free(Lfinite);
      Lfinite=0;

      
      //     a        = kvfinite[a]
      if (a_tmp)
	gsl_vector_free(a_tmp);
      a_tmp=gsl_vector_calloc(a->size);
      for (int index=0;index<a->size;index++) {
	int a_index=(int)gsl_vector_get(a,index);
	gsl_vector_set(a_tmp,index,gsl_vector_get(kvfinite,a_index-1));
      }
      if (a)
	gsl_vector_free(a);
      a=a_tmp;
      a_tmp=0;

      //Kminfin = min(a)
      double Kminfin=checked_min_gsl_vector(a);

      //Kmaxfin = max(a)
      double Kmaxfin=checked_max_gsl_vector(a);

      //     oumin    = which(a==(Kminfin))
      int oumin=0;
      for (int index=0;index<a->size && !oumin;index++)
	if (gsl_vector_get(a,index)==Kminfin) {
	  oumin=index+1;
	}
       
      assert(oumin>0);

      //     oumax    = which(a==(Kmaxfin))
      int oumax=0;
      for (int index=0;index<a->size && !oumax;index++)
	if (gsl_vector_get(a,index)==Kmaxfin) {
	  oumax=index+1;
	}

      assert(oumax>0);

      //     a        = a[oumin:oumax]
      if (a_tmp)
	gsl_vector_free(a_tmp);
      a_tmp=copy_range_gsl_vector(a,oumin-1,oumax-1);
      gsl_vector_free(a);
      a=a_tmp;
      a_tmp=0;

      
      //     kvfinite = sort(a)
      if (kvfinite)
	gsl_vector_free(kvfinite);
      kvfinite=gsl_vector_alloc(a->size);
      gsl_vector_memcpy(kvfinite,a);
      gsl_sort_vector(kvfinite);
      gsl_vector_free(a);
      a=0;
      
      //     # find the coordinates of points out of the convex hull
      //     Kconc    = c(1:Kmax)
      if (Kconc)
	gsl_vector_free(Kconc);
      Kconc=build_range_gsl_vector(1,Kmax);
      //      printf("Kconc 1 : \n");
      //      gsl_vector_fprintf(stdout,Kconc,"%lf");

      //     Kconc    = Kconc[-which(Kconc %in% c(kvfinite))]
      int n_keep=Kconc->size;
      if (Kconc_keep)
	gsl_vector_free(Kconc_keep);
      Kconc_keep=gsl_vector_alloc(Kconc->size);
      gsl_vector_set_all(Kconc_keep,1.0);
      for (int index=0;index<Kconc->size;index++) {
	bool remove=false;
	for (int index2=0;index2<kvfinite->size && !remove;index2++)
	  if (gsl_vector_get(Kconc,index)==gsl_vector_get(kvfinite,index2)) {
	    remove=true;
	    gsl_vector_set(Kconc_keep,index,0.0);
	    n_keep--;
	  }
      }

      if (kvfinite)
	gsl_vector_free(kvfinite);
      kvfinite=0;

      if (Kconc_tmp)
	gsl_vector_free(Kconc_tmp);
      Kconc_tmp=gsl_vector_calloc(n_keep);
      for (int index=0,offset=0;index<Kconc_keep->size && offset<Kconc_tmp->size;index++)
	if (gsl_vector_get(Kconc_keep,index) != 0.0) {
	  gsl_vector_set(Kconc_tmp,offset,gsl_vector_get(Kconc,index));
	  offset++;
	}
      gsl_vector_free(Kconc_keep);
      Kconc_keep=0;
      gsl_vector_free(Kconc);
      Kconc=Kconc_tmp;
      Kconc_tmp=0;
      //printf("Kconc 2 : \n");
      //gsl_vector_fprintf(stdout,Kconc,"%lf");

      //     Kconc    = Kconc[Kconc>=Kmin]
      n_keep=Kconc->size;

      Kconc_keep=gsl_vector_calloc(Kconc->size);
      gsl_vector_set_all(Kconc_keep,1.0);

      for (int index=0;index<Kconc->size;index++) {
	if (gsl_vector_get(Kconc,index)<Kmin) {
	  gsl_vector_set(Kconc_keep,index,0.0);
	  n_keep--;
	}
      }
      
      
      if (n_keep<=0) {
	gsl_vector_free(Kconc_keep);
      } else {
	Kconc_tmp=gsl_vector_alloc(n_keep);
	for (int index=0,offset=0;index<Kconc->size && offset<Kconc_keep->size;index++) {
	  if (gsl_vector_get(Kconc_keep,index) != 0.0) {
	    gsl_vector_set(Kconc_tmp,offset,gsl_vector_get(Kconc,index));
	    offset++;
	  }	  
	}
	
	gsl_vector_free(Kconc_keep);
	Kconc_keep=0;
	gsl_vector_free(Kconc);
	Kconc=Kconc_tmp;
	Kconc_tmp=0;
		
	//     for (k in Kconc){
	for (int Kconc_index=0;Kconc_index<Kconc->size;Kconc_index++) {
	  int k=(int)gsl_vector_get(Kconc,Kconc_index);
	  //	fprintf(stdout,"k=\n%d\n",k);
	  
	  //	fprintf(stdout,"Linc before neighbors=\n");
	  //	gsl_matrix_fprintf(stdout,Linc,"%lf");
	  
	  //       out.neighbors  = neighbors(L=Linc,k=k,param=param,P=P,lmin=lmin,vh=vh)
	  if (L_res)
	  gsl_vector_free(L_res);
	  L_res=0;
	  if (Linc_vec)
	    gsl_vector_free(Linc_vec);
	  Linc_vec=gsl_vector_calloc(Linc->size1);
	  gsl_matrix_get_col(Linc_vec,Linc,0);
	  //	fprintf(stdout,"call 1\n");
	  
	  neighbors(x,k,Linc_vec,k,param,param_size,P,lmin,lmax,vh,&L_res,param_res,param_res_size);
	  
	  //       param          = out.neighbors$param
	  for (int index=0;index<param_size;index++) {
	    if (param[index].phi)
	      gsl_vector_free(param[index].phi);
	    if (param[index].rupt)
	      gsl_matrix_free(param[index].rupt);
	  }
	  delete[] param;
	  param=0;
	  
	  param=new param_struct[*param_res_size];
	  for (int index=0;index<*param_res_size;index++) {
	    param[index].phi=(*param_res)[index].phi;
	    (*param_res)[index].phi=0;
	    param[index].rupt=(*param_res)[index].rupt;
	    (*param_res)[index].rupt=0;
	  }
	  param_size=*param_res_size;
	  delete[] (*param_res);
	  *param_res=0;
	  *param_res_size=0;
	  
	  //       Linc           = out.neighbors$L
	  
	  gsl_matrix_free(Linc);
	  Linc=gsl_matrix_calloc(L_res->size,1);
	  gsl_matrix_set_col(Linc,0,L_res);
	  gsl_vector_free(L_res);
	  L_res=0;
	  //fprintf(stdout,"Linc after neighbors=\n");
	  //gsl_matrix_fprintf(stdout,Linc,"%lf");
	  
	  //       #plot(1:length(Linc),Linc,col=1)     lines(1:length(Ltmp),Ltmp,col=2)     lines(k,Linc[k],col=3)
	  
	  //     } # end k
	}
	//      fclose(fres);

	//////////// Commentaire pour Mark ////////////////	
	//////////// le else de la condition n_keep ne se finit pas ici car on ne peut pas faire tourner neighbors quand n_keep=0 ////////////////	
	if (Kconc) 
	  gsl_vector_free(Kconc);
	Kconc=0;
	
	//     out.neighbors  = neighbors(L=Linc,k=Kmax,param=param,P=P,lmin=lmin,vh)
	if (L_res)
	  gsl_vector_free(L_res);
	L_res=0;
	if (Linc_vec)
	  gsl_vector_free(Linc_vec);
	Linc_vec=gsl_vector_calloc(Linc->size1);
	gsl_matrix_get_col(Linc_vec,Linc,0);
	//      fprintf(stdout,"call 2\n");
	     	
	neighbors(x,Kmax,Linc_vec,Kmax,param,param_size,P,lmin,lmax,vh,&L_res,param_res,param_res_size);
	gsl_vector_free(Linc_vec);
	Linc_vec=0;
	
	//     param          = out.neighbors$param
	for (int index=0;index<param_size;index++) {
	  if (param[index].phi)
	    gsl_vector_free(param[index].phi);
	  if (param[index].rupt)
	    gsl_matrix_free(param[index].rupt);
	}
	delete[] param;
	
	param_size=*param_res_size;
	param=new param_struct[param_size];
	for (int index=0;index<param_size;index++) {
	  param[index].phi=(*param_res)[index].phi;
	  (*param_res)[index].phi=0;
	  param[index].rupt=(*param_res)[index].rupt;
	  (*param_res)[index].rupt=0;
	}
	delete[] (*param_res);
	*param_res=0;
	*param_res_size=0;
	
	//     Linc           = out.neighbors$L
	gsl_matrix_free(Linc);
	Linc=gsl_matrix_calloc(L_res->size,1);
	gsl_matrix_set_col(Linc,0,L_res);
	gsl_vector_free(L_res);
	L_res=0;
	//      fprintf(stdout,"Linc res =\n");
	//gsl_matrix_fprintf(stdout,Linc,"%lf");
       
      } // end else n_keep  
//////////////////////////////// le end du else n_keep est ici !! //////////////////////////////////
      //   } # end while
    }
    
    if (Ltmp)
      gsl_matrix_free(Ltmp);
    Ltmp=0;
    // } # end else Pmin ==1
  }
  
  
  // invisible(list(Linc=Linc,param=param))
  *Linc_res=Linc;
  *param_res=param;
  *param_res_size=param_size;
  
// } #end function


}

}
