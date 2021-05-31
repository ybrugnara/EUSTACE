#include "EM_init.h"

#include <iostream>
using namespace std;

#include <gsl/gsl_math.h>

#include "helperfuncs.h"

namespace segclust {

double 
compute_dist_var_homo(double nk, double vk, double nr, double vr, double varpool)
{
  // Dist[k,r] = -n[k]*v[k] - n[r]*v[r] + (n[k]+n[r])*varpool
  return (-1.0*nk*vk  - nr*vr + (nk+nr)*varpool);
}

double 
compute_dist_var_inhomo(double nk, double vk, double nr, double vr, double varpool)
{
  // Dist[k,r] = -n[k]*log(v[k])-n[r]*log(v[r])+(n[k]+n[r])*log(varpool)
  return ( -1.0*nk*log(vk) - nr*log(vr) + (nk+nr)*log(varpool));
}



gsl_vector *
EM_init(gsl_vector *x, gsl_matrix *rupt, long K, long P, bool vh)
{
  double (*compute_dist)(double, double, double, double, double);

  //  FILE *f=fdopen(1,"w");

  // m = apply(rupt,1,FUN=function(y) mean(x[y[1]:y[2]]))
  gsl_vector *m=gsl_vector_calloc(rupt->size1);
  for (int r=0;r<rupt->size1;r++) {
    long y1=(long)gsl_matrix_get(rupt,r,1-1);
    long y2=(long)gsl_matrix_get(rupt,r,2-1);
    gsl_vector *meanvec=copy_range_gsl_vector(x,y1-1,y2-1);
    double mean=mean_gsl_vector(meanvec);
    gsl_vector_free(meanvec);
    gsl_vector_set(m,r,mean);
  }

  //cout << "m" << endl;
  //gsl_vector_fprintf(stdout,m,"%f");
  
  // v =  apply(rupt,1,FUN=function(y) var(x[y[1]:y[2]]))
  gsl_vector *v=gsl_vector_calloc(rupt->size1);
  for (int r=0;r<rupt->size1;r++) {
    long y1=(long)gsl_matrix_get(rupt,r,1-1);
    long y2=(long)gsl_matrix_get(rupt,r,2-1);
    gsl_vector *varvec=copy_range_gsl_vector(x,y1-1,y2-1);
    double var=var_gsl_vector(varvec);
    gsl_vector_free(varvec);
    gsl_vector_set(v,r,var);
  }
  
  //  cout << "v" << endl;
  //gsl_vector_fprintf(stdout,v,"%f");

  // v[is.na(v)]=0
  for (int index=0;index<v->size;index++)
    if (gsl_isnan(gsl_vector_get(v,index)))
      gsl_vector_set(v,index,0);

  // n = apply(rupt,1,diff)+1
  /** rupt has only 2 columns => the diff matrix, n, is a vector ! **/
  gsl_vector *n=gsl_vector_calloc(rupt->size1);
  for (int r=0;r<rupt->size1;r++) {
      gsl_vector_set(n,r,
		     gsl_matrix_get(rupt,r,1)-gsl_matrix_get(rupt,r,0)+1.0);
  }
  

  // v=(n-1)/n * v
  for (int index=0;index<v->size;index++) {
    double vi=gsl_vector_get(v,index);
    double ni=gsl_vector_get(n,index);
    gsl_vector_set(v,index,(ni-1)/ni*vi);
  }

  // cout << "n" << endl;
  //gsl_vector_fprintf(stdout,n,"%f");

  // Dist = matrix(Inf,K,K)
  gsl_matrix *Dist=gsl_matrix_calloc(K,K);
  gsl_matrix_set_all(Dist,GSL_POSINF);


  // if (vh==TRUE) {
  if (vh==true) {
    compute_dist=compute_dist_var_homo;
  //} else {
  } else {
    compute_dist=compute_dist_var_inhomo;
  //}
  }

  // for (k in (1: (K-1))){
  for (long k=1; k<=K-1;k++) {
  
    //    for (r in ( (k+1) :K)) {
    for (long r=k+1; r<=K;r++) {

      double mk=GSL_NAN;
      if (k-1<m->size)
	mk=gsl_vector_get(m,k-1);

      double mr=GSL_NAN;
      if (r-1<m->size)
	mr=gsl_vector_get(m,r-1);

      double vk=GSL_NAN;
      if (k-1<v->size)
	vk=gsl_vector_get(v,k-1);

      double vr=GSL_NAN;
      if (r-1<v->size)
	vr=gsl_vector_get(v,r-1);

      double nk=GSL_NAN;
      if (k-1<n->size)
	nk=gsl_vector_get(n,k-1);

      double nr=GSL_NAN;
      if (r-1<n->size)
	nr=gsl_vector_get(n,r-1);

      //       ybar      = (n[k]*m[k]+n[r]*m[r])/(n[k]+n[r])
      double ybar=(nk*mk+nr*mr)/(nk+nr);

      //       varpool   = (  n[r]*v[r] + n[k]*v[k] + n[r]*(m[r]-ybar)^2 + n[k]*(m[k]-ybar)^2 ) / (n[r]+n[k])
      double varpool_1=nr*vr+nk*vk+nr*(mr-ybar)*(mr-ybar);
      double varpool_2=mk-ybar;

      double varpool=(varpool_1+nk*varpool_2*varpool_2)/(nr+nk);


      // eval(parse(text=stringD))
      double Distkr=compute_dist(nk,vk,nr,vr,varpool);
      gsl_matrix_set(Dist,k-1,r-1,Distkr);


      //       }
    }

  // }
  }
  
  // LR = NULL
  
  // if (nrow(Dist)>P) {
  if (Dist->size1 > P) {

    // for (indice in (1:K)){
    for (long indice=1;indice<=K;indice++) {

      // #while (length(Dist) > P){
      
      //    if (nrow(as.matrix(Dist)) <= P) {
      if (Dist->size1 <=P) {

      //    break 
	break;

      //    }
      }
      
      //    out = which(Dist==min(Dist),arr.ind=TRUE)
      //    imin = out[1]
      //    jmin = out[2]
      //    Dmin = min(Dist)
      long imin=0;
      long jmin=0;
      double Dmin=gsl_matrix_get(Dist,0,0);
      for (int r=0;r<Dist->size1;r++)
	for (int c=0;c<Dist->size2;c++) {
	  double Distval=gsl_matrix_get(Dist,r,c);
	  if (Distval<Dmin) {
	    Dmin=Distval;
	    imin=r;
	    jmin=c;
	  }
	}
      
      //    LR[nrow(Dist)] = Dmin

      double nimin=gsl_vector_get(n,imin);
      double njmin=gsl_vector_get(n,jmin);

      double mimin=gsl_vector_get(m,imin);
      double mjmin=gsl_vector_get(m,jmin);

      double vimin=gsl_vector_get(v,imin);
      double vjmin=gsl_vector_get(v,jmin);

      //    ntmp = n[imin]+n[jmin]
      double ntmp=nimin+njmin;

      //    mtmp = (n[imin]*m[imin]+n[jmin]*m[jmin])/ntmp
      double mtmp=(nimin*mimin+njmin*mjmin)/ntmp;

      //    vtmp = (  n[imin]*v[imin] + n[jmin]*v[jmin] + n[imin]*(m[imin]-mtmp)^2 + n[jmin]*(m[jmin]-mtmp)^2 ) / ntmp
      double vtmp=(nimin*vimin+njmin*vjmin+nimin*(mimin-mtmp)*(mimin-mtmp)+
		   njmin*(mjmin-mtmp)*(mjmin-mtmp))/ntmp;

      
      //    Dist <- Dist[-c(imin,jmin), -c(imin,jmin)]
      int n_remove=2;
      if (imin==jmin)
	n_remove=1;
      gsl_matrix *Dist_tmp=gsl_matrix_calloc(Dist->size1-n_remove,
					   Dist->size2-n_remove);
      for (int r1=0,r2=0;r1<Dist->size1;r1++) {
	if (r1!=imin && r1!=jmin) {
	  for (int c1=0,c2=0;c1<Dist->size2;c1++) {
	    if (c1!=imin && c1!=jmin) {
	      gsl_matrix_set(Dist_tmp,r2,c2,gsl_matrix_get(Dist,r1,c1));
	      c2++;
	    }
	  }
	  r2++;
	}
      }
      gsl_matrix_free(Dist);
      Dist=Dist_tmp;

      //fprintf(stdout,"Dist");
      //gsl_matrix_fprintf(stdout,Dist,"%lf");
	

      
      //    m = m[-c(imin,jmin)]
      gsl_vector *m_tmp=gsl_vector_calloc(m->size-n_remove);
      for (int i=0,index=0;i<m->size;i++) {
	if (i!=imin && i!=jmin) {
	  gsl_vector_set(m_tmp,index,gsl_vector_get(m,i));
	  index++;
	}
      }
      gsl_vector_free(m);
      m=m_tmp;


      //    v = v[-c(imin,jmin)]
      gsl_vector *v_tmp=gsl_vector_calloc(v->size-n_remove);
      for (int i=0,index=0;i<v->size;i++) {
	if (i!=imin && i!=jmin) {
	  gsl_vector_set(v_tmp,index,gsl_vector_get(v,i));
	  index++;
	}
      }
      gsl_vector_free(v);
      v=v_tmp;


      //    n = n[-c(imin,jmin)]
      gsl_vector *n_tmp=gsl_vector_calloc(n->size-n_remove);
      for (int i=0,index=0;i<n->size;i++) {
	if (i!=imin && i!=jmin) {
	  gsl_vector_set(n_tmp,index,gsl_vector_get(n,i));
	  index++;
	}
      }
      gsl_vector_free(n);
      n=n_tmp;

      //    # update Dist N Nplus LogL
      //    m = c(m , mtmp)
      m_tmp=gsl_vector_calloc(m->size+1);
      for (int index=0;index<m->size;index++)
	gsl_vector_set(m_tmp,index,gsl_vector_get(m,index));
      gsl_vector_free(m);
      m=m_tmp;
      gsl_vector_set(m,m->size-1,mtmp);


      //    n = c(n, ntmp)
      n_tmp=gsl_vector_calloc(n->size+1);
      for (int index=0;index<n->size;index++)
	gsl_vector_set(n_tmp,index,gsl_vector_get(n,index));
      gsl_vector_free(n);
      n=n_tmp;
      gsl_vector_set(n,n->size-1,ntmp);


      //    v = c(v, vtmp)
      v_tmp=gsl_vector_calloc(v->size+1);
      for (int index=0;index<v->size;index++)
	gsl_vector_set(v_tmp,index,gsl_vector_get(v,index));
      gsl_vector_free(v);
      v=v_tmp;
      gsl_vector_set(v,v->size-1,vtmp);

      //    Dtmp = rep(Inf,times=nrow(as.matrix(Dist)))
      gsl_vector *Dtmp=gsl_vector_calloc(Dist->size1);
      gsl_vector_set_all(Dtmp,GSL_POSINF);

      //    for (k  in (1:nrow(as.matrix(Dist)))){
      for (int k=1;k<=Dist->size1;k++) {

	double nk=GSL_NAN;
	if (k-1<n->size)
	  nk=gsl_vector_get(n,k-1);

	double mk=GSL_NAN;
	if (k-1<m->size)
	  mk=gsl_vector_get(m,k-1);

	double vk=GSL_NAN;
	if (k-1<v->size)
	  vk=gsl_vector_get(v,k-1);

	//       ybar      = (n[k]*m[k]+ntmp*mtmp)/(n[k]+ntmp)
	double ybar=(nk*mk+ntmp*mtmp)/(nk+ntmp);


	//       varpool   = (  n[k]*v[k] + ntmp*vtmp + n[k]*(m[k]-ybar)^2 + ntmp*(mtmp-ybar)^2  ) / (n[k]+ntmp)
	double varpool=(nk*vk+ntmp*vtmp+nk*(mk-ybar)*(mk-ybar)+
			ntmp*(mtmp-ybar)*(mtmp-ybar))/(nk+ntmp);

	// eval(parse(text=stringDtmp))
	double DtmpVal=compute_dist(nk,vk,ntmp,vtmp,varpool);
	gsl_vector_set(Dtmp,k-1,DtmpVal);

	//       }
      }

      //       Dist = rbind(cbind(Dist,Dtmp),rep(Inf,times=(nrow(as.matrix(Dist))+1)))
      gsl_vector *repInf=gsl_vector_calloc(Dist->size1+1);
      gsl_vector_set_all(repInf,GSL_POSINF);
      
      gsl_matrix *cbind=gsl_matrix_calloc(Dist->size1,Dist->size2+1);
      for (int c=0;c<Dist->size2;c++) {
	gsl_vector *colvec=gsl_vector_calloc(Dist->size2);
	gsl_matrix_get_col(colvec,Dist,c);
	gsl_matrix_set_col(cbind,c,colvec);
	gsl_vector_free(colvec);
      }
      gsl_matrix_set_col(cbind,cbind->size2-1,Dtmp);
      gsl_vector_free(Dtmp);

      gsl_matrix *rbind=gsl_matrix_calloc(cbind->size1+1,cbind->size2);
      for (int r=0;r<cbind->size1;r++) {
	gsl_vector *rowvec=gsl_vector_calloc(cbind->size2);
	gsl_matrix_get_row(rowvec,cbind,r);
	gsl_matrix_set_row(rbind,r,rowvec);
	gsl_vector_free(rowvec);
      }
      gsl_matrix_free(cbind);

      gsl_matrix_set_row(rbind,rbind->size1-1,repInf);
      gsl_vector_free(repInf);
      gsl_matrix_free(Dist);
      Dist=rbind;
      
	// }
    }
      
      
      // }
  }

  gsl_matrix_free(Dist);

  gsl_vector *phi0=gsl_vector_calloc(m->size+v->size+n->size);
  // if (vh==TRUE) {
  if (vh==true) {

    //tau = matrix(0,ncol=P,nrow=K)
    gsl_matrix *tau=gsl_matrix_calloc(K,P);

    //s = rep(0,P);
    gsl_vector *s=gsl_vector_calloc(P);


    gsl_vector *tmp=gsl_vector_calloc(P);
    //for (k in (1:K)) {
    for (int k=1; k<=K;k++) {

      // tmp = rep(Inf,P);
      gsl_vector_set_all(tmp,GSL_POSINF);
      
      // for (p in (1:P)) {
      for (int p=1; p<=P;p++) {

	// tmp[p] = sum ((x[(rupt[k,2]:rupt[k,1])]-m[p])^2)
	int k1=(int)gsl_matrix_get(rupt,k-1,1-1);
	int k2=(int)gsl_matrix_get(rupt,k-1,2-1);
	gsl_vector *ruptk1k2=copy_range_gsl_vector(x,k1-1,k2-1);
	gsl_vector_add_constant(ruptk1k2,-1.0*gsl_vector_get(m,p-1));
	gsl_vector_mul(ruptk1k2,ruptk1k2);

	double sum=sum_gsl_vector(ruptk1k2);
	
	gsl_vector_free(ruptk1k2);
	
	gsl_vector_set(tmp,p-1,sum);

	//}
      }

      //tau[k,which.min(tmp)]=1
      int minindex=checked_min_index_gsl_vector(tmp);
      gsl_matrix_set(tau,k-1,minindex,1);

      //}
    }

    //    fprintf(stdout,"tau=");
    //gsl_matrix_fprintf(stdout,tau,"%lf");

    gsl_vector_free(tmp);

    tmp=gsl_vector_alloc(tau->size1);
    //for (p in (1:P)) {
    for (int p=1;p<=P;p++) {
      gsl_vector_set_all(tmp,0.0);
      // tmp = tau[,p] * (apply(rupt,1,FUN=function(y) sum( (x[y[1]:y[2]]-m[p])^2 )))
      gsl_matrix_get_col(tmp,tau,p-1);
      for (int l=0;l<rupt->size1;l++) {
	int y1=(int)gsl_matrix_get(rupt,l,1-1);
	int y2=(int)gsl_matrix_get(rupt,l,2-1);
	gsl_vector *xtmp=copy_range_gsl_vector(x,y1-1,y2-1);
	gsl_vector_add_constant(xtmp,-1.0*gsl_vector_get(m,p-1));
	gsl_vector_mul(xtmp,xtmp);
	double sum=sum_gsl_vector(xtmp);
	gsl_vector_set(tmp,l,gsl_matrix_get(tau,l,p-1)*sum);
	gsl_vector_free(xtmp);
      }
      //fprintf(stdout,"tmp=");
     // gsl_vector_fprintf(stdout,tmp,"%lf");
      
      // s[p] = sum(tmp)
      gsl_vector_set(s,p-1,sum_gsl_vector(tmp));

      //}
    }

    //    fprintf(stdout,"s=");
    //gsl_vector_fprintf(stdout,s,"%lf");

    // s = repmat(sum(s)/sum(n),P,1);
    double sums=sum_gsl_vector(s);
    double sumn=sum_gsl_vector(n);
    gsl_vector_set_all(s,sums/sumn);
    
    
    //phi0 = c(m,sqrt(s),n/sum(n));
    for (int i=0;i<m->size;i++)
      gsl_vector_set(phi0,i,gsl_vector_get(m,i));
    
    gsl_vector *tmp_s=apply_basicfunc_gsl_vector(s,sqrt);
    for (int i=0;i<s->size;i++)
      gsl_vector_set(phi0,m->size+i,gsl_vector_get(tmp_s,i));
    gsl_vector_free(tmp_s);
    
    for (int i=0;i<n->size;i++)
      gsl_vector_set(phi0,m->size+v->size+i,gsl_vector_get(n,i)/sumn);

    gsl_vector_free(tmp);
    gsl_vector_free(s);

   // } else {
    gsl_matrix_free(tau);
  } else {

    // phi0= c(m,sqrt(v),n/sum(n))
    for (int i=0;i<m->size;i++)
      gsl_vector_set(phi0,i,gsl_vector_get(m,i));
    
    gsl_vector *tmp_v=apply_basicfunc_gsl_vector(v,sqrt);
    for (int i=0;i<v->size;i++)
      gsl_vector_set(phi0,m->size+i,gsl_vector_get(tmp_v,i));
    gsl_vector_free(tmp_v);
    
    double sum=sum_gsl_vector(n);
    gsl_vector_scale(n,1.0/sum);
    
    for (int i=0;i<n->size;i++)
      gsl_vector_set(phi0,m->size+v->size+i,gsl_vector_get(n,i));
    
    //}
  }
  gsl_vector_free(m);
  gsl_vector_free(v);
  gsl_vector_free(n);
  // invisible(phi0)
  return phi0;
}

}
