#ifndef SEGCLUST_NEIGHBORS_H
#define SEGCLUST_NEIGHBORS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "segclust_types.h"

namespace segclust {

void neighbors(gsl_vector *x, long Kmax,
	       gsl_vector *L, long k, 
	       param_struct param[], long param_size,
	       long P, long lmin, long lmax, bool vh,
	       gsl_vector **L_res, 
	       param_struct *param_res[],long *param_res_size);

}
#endif
