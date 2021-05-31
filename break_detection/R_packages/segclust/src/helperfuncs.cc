#include "helperfuncs.h"

#include <cmath>
#include <cassert>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

namespace segclust {

gsl_matrix * 
apply_basicfunc_gsl_matrix(gsl_matrix *mat, basicfunc bf) {

  gsl_matrix *res=gsl_matrix_calloc(mat->size1,mat->size2);
  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      gsl_matrix_set(res,r,c,bf(gsl_matrix_get(mat,r,c)));

  return res;
  
}

gsl_vector * 
apply_basicfunc_gsl_vector(gsl_vector *v, basicfunc bf) {

  gsl_vector *res=gsl_vector_calloc(v->size);
  for (int i=0;i<v->size;i++)
      gsl_vector_set(res,i,bf(gsl_vector_get(v,i)));

  return res;
  
}


gsl_vector *
copy_range_gsl_vector(gsl_vector *v, int min, int max)
{
  if (min>max) {
    int tmp=min;
    min=max;
    max=min;
  }

  gsl_vector *res=gsl_vector_calloc((max-min)+1);

  for (int i=min;i<=max;i++)
    gsl_vector_set(res,i-min,gsl_vector_get(v,i));

  return res;
}


gsl_vector *
build_range_gsl_vector(int min, int max)
{
  int inc=1;
  int size=(max-min)+1;

  if (min>max) {
    size=(min-max)+1;
    inc=-1;
  }

  gsl_vector *res=gsl_vector_calloc(max-min+1);
  for (int index=0,val=min;index<size;index++,val+=inc) {
    gsl_vector_set(res,index,val);
  }

  return res;
}

gsl_vector *
rowsum_gsl_matrix(gsl_matrix *mat)
{
  gsl_vector *vec=gsl_vector_calloc(mat->size1);

  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      gsl_vector_set(vec,r,gsl_vector_get(vec,r)+gsl_matrix_get(mat,r,c));

  return vec;
}

gsl_vector *
colsum_gsl_matrix(gsl_matrix *mat)
{
  gsl_vector *vec=gsl_vector_calloc(mat->size2);

    for (int c=0;c<mat->size2;c++)
      for (int r=0;r<mat->size1;r++)
	gsl_vector_set(vec,c,gsl_vector_get(vec,c)+gsl_matrix_get(mat,r,c));

  return vec;
}

double
sum_gsl_matrix(gsl_matrix *mat)
{
  double sum=0.0;

  for (int r=0;r<mat->size1;r++)
    for (int c=0;c<mat->size2;c++)
      sum+=gsl_matrix_get(mat,r,c);

  return sum;
}

double
sum_gsl_vector(gsl_vector *v)
{
  double sum=0.0;

  for (int i=0;i<v->size;i++)
    sum+=gsl_vector_get(v,i);

  return sum;
}

gsl_vector *
cumsum_gsl_matrix(gsl_matrix *m)
{
  gsl_vector *v=gsl_vector_calloc(m->size1*m->size2);

  int index=0;
  for (int r=0;r<m->size1;r++)
    for (int c=0;c<m->size2;c++) {
      if (r == 0 && c == 0) {
	gsl_vector_set(v,index,gsl_matrix_get(m,r,c));
      } else {
	gsl_vector_set(v,index,gsl_vector_get(v,index-1)+gsl_matrix_get(m,r,c));
      }
      index++;
    }
  return v;
}

gsl_vector *
cumsum_gsl_vector(gsl_vector *v)
{
  gsl_vector *res=gsl_vector_calloc(v->size);
  gsl_vector_set(res,0,gsl_vector_get(v,0));
  for (int index=1;index<v->size;index++)
    gsl_vector_set(res,index,gsl_vector_get(res,index-1)+gsl_vector_get(v,index));
  return res;
  
}

gsl_vector *
colmax_gsl_matrix(gsl_matrix *m)
{
  gsl_vector *res=gsl_vector_calloc(m->size2);
  for (int c=0;c<m->size2;c++) {
    double colmax=gsl_matrix_get(m,0,c);
    for (int r=1;r<m->size1;r++) {
      colmax=GSL_MAX(colmax,gsl_matrix_get(m,r,c));
    }
    gsl_vector_set(res,c,colmax);
  }

  return res;
}

gsl_vector *
rowmax_gsl_matrix(gsl_matrix *m)
{
  gsl_vector *res=gsl_vector_calloc(m->size1);
  for (int r=0;r<m->size1;r++) {
    double rowmax=gsl_matrix_get(m,r,0);
    for (int c=1;c<m->size2;c++) {
      rowmax=GSL_MAX(rowmax,gsl_matrix_get(m,r,c));
    }
    gsl_vector_set(res,r,rowmax);
  }

  return res;
}

gsl_vector *
colmin_gsl_matrix(gsl_matrix *m)
{
  gsl_vector *res=gsl_vector_calloc(m->size2);
  for (int c=0;c<m->size2;c++) {
    double colmin=gsl_matrix_get(m,0,c);
    for (int r=1;r<m->size1;r++) {
      colmin=GSL_MIN(colmin,gsl_matrix_get(m,r,c));
    }
    gsl_vector_set(res,c,colmin);
  }

  return res;
}

gsl_vector *
rowmin_gsl_matrix(gsl_matrix *m)
{
  gsl_vector *res=gsl_vector_calloc(m->size1);
  for (int r=0;r<m->size1;r++) {
    double rowmin=gsl_matrix_get(m,r,0);
    for (int c=1;c<m->size2;c++) {
      rowmin=GSL_MIN(rowmin,gsl_matrix_get(m,r,c));
    }
    gsl_vector_set(res,r,rowmin);
  }

  return res;
}

double
mean_gsl_vector(gsl_vector *v) {
  double mean=GSL_NAN;
  if (v->size) {
    mean=0.0;
    for (int i=0;i<v->size;i++)
      mean+=gsl_vector_get(v,i);
    mean/=(1.0*v->size);
  }
  return mean;
  
}

double
var_gsl_vector(gsl_vector *v) {
  double var=GSL_NAN;
  if (v->size) {
    var=0.0;
    double mean=mean_gsl_vector(v);
    for (int i=0;i<v->size;i++) {
      var+=(gsl_vector_get(v,i)*gsl_vector_get(v,i))-mean*mean;
    }
    if (v->size>1)
      var/=(v->size-1.0);
  }

  return var;

}

double
stddev_gsl_vector(gsl_vector *v) {

  double stddev=GSL_NAN;
  if (v->size)
    stddev=sqrt(var_gsl_vector(v));
  
  return stddev;

}

gsl_vector *
diff_gsl_vector(gsl_vector *v) {

  gsl_vector *res=gsl_vector_calloc(v->size-1);
  for (int i=1; i<v->size;i++)
    gsl_vector_set(res,i-1,gsl_vector_get(v,i)-gsl_vector_get(v,i-1));

  return res;
}

int count_diffs_gsl_vectors(gsl_vector *v1, gsl_vector *v2)
{
  gsl_vector *shortest=v1;
  gsl_vector *longest=v2;
  
  if (shortest->size > longest->size) {
    gsl_vector *tmp=shortest;
    shortest=longest;
    longest=tmp;
  }

  int n=longest->size-shortest->size;

  for (int i=0;i<shortest->size;i++)
    if (gsl_vector_get(shortest,i) != gsl_vector_get(longest,i))
      n++;

  return n;
}

int
sum_diffs_gsl_matrices(gsl_matrix *m1, gsl_matrix *m2)
{
  assert(m1->size1 == m2->size1);
  assert(m1->size2 == m2->size2);

  int n=0;

  for (int r=0;r<m1->size1;r++)
    for (int c=0;c<m1->size2;c++)
      if (gsl_matrix_get(m1,r,c) != gsl_matrix_get(m2,r,c))
	n++;

  return n;
  
}

double checked_max_gsl_vector(gsl_vector *v)
{
  double maxval=GSL_NEGINF;
  for (int index=0;index<v->size;index++)
    if (!isnan(gsl_vector_get(v,index)) && gsl_vector_get(v,index)>maxval)
      maxval=gsl_vector_get(v,index);


  return maxval;
}

int checked_max_index_gsl_vector(gsl_vector *v)
{
  double maxval=GSL_NEGINF;
  int index=-1;

  for (int i=0;i<v->size;i++) {
    if ((!isnan(gsl_vector_get(v,i))) && gsl_vector_get(v,i)>maxval) {
      maxval=gsl_vector_get(v,i);
      index=i;
    }
  }

  return index;
}

double checked_min_gsl_vector(gsl_vector *v)
{
  double minval=GSL_POSINF;
  for (int index=0;index<v->size;index++)
    if (!isnan(gsl_vector_get(v,index)) && gsl_vector_get(v,index)<minval)
      minval=gsl_vector_get(v,index);

  return minval;
}

int checked_min_index_gsl_vector(gsl_vector *v)
{
  double minval=GSL_POSINF;
  int index=-1;

  for (int i=0;i<v->size;i++)
    if ((!isnan(gsl_vector_get(v,i))) && gsl_vector_get(v,i)<minval) {
      minval=gsl_vector_get(v,i);
      index=i;
    }

  return index;
}

}
