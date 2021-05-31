#ifndef SEGCLUST_EM_ALGO_H
#define SEGCLUST_EM_ALGO_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace segclust {

void EM_algo(gsl_vector *x, gsl_matrix *rupt, int P, gsl_vector *phi0,
	     bool vh,
	     gsl_vector **phi_res, gsl_matrix **tau_res,
	     double *lvinc, int *empty, int *dv);

}
#endif
