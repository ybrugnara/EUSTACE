#ifndef SEGCLUST_HELPERFUNCS_H
#define SEGCLUST_HELPERFUNCS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

extern "C" {

  typedef double (*basicfunc)(const double);

}

namespace segclust {



gsl_matrix * apply_basicfunc_gsl_matrix(gsl_matrix *, basicfunc);

gsl_vector * apply_basicfunc_gsl_vector(gsl_vector *, basicfunc);

gsl_vector *copy_range_gsl_vector(gsl_vector *, int, int);

gsl_vector *build_range_gsl_vector(int, int);


gsl_vector * rowsum_gsl_matrix(gsl_matrix *);

gsl_vector * colsum_gsl_matrix(gsl_matrix *);

double sum_gsl_matrix(gsl_matrix *);

double sum_gsl_vector(gsl_vector *);

gsl_vector * cumsum_gsl_matrix(gsl_matrix *);

gsl_vector * cumsum_gsl_vector(gsl_vector *);

gsl_vector *colmax_gsl_matrix(gsl_matrix *);

gsl_vector *rowmax_gsl_matrix(gsl_matrix *);

gsl_vector *colmin_gsl_matrix(gsl_matrix *);

gsl_vector *rowmin_gsl_matrix(gsl_matrix *);


double mean_gsl_vector(gsl_vector *);

double var_gsl_vector(gsl_vector *);

double stddev_gsl_vector(gsl_vector *);

gsl_vector *diff_gsl_vector(gsl_vector *);

int sum_diffs_gsl_vectors(gsl_vector *, gsl_vector *);

int sum_diffs_gsl_matrices(gsl_matrix *, gsl_matrix *);

double checked_max_gsl_vector(gsl_vector *);

int checked_max_index_gsl_vector(gsl_vector *);


double checked_min_gsl_vector(gsl_vector *);

int checked_min_index_gsl_vector(gsl_vector *);

}
#endif
