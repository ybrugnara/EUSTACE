#ifndef SEGCLUST_HYBRID_H
#define SEGCLUST_HYBRID_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "segclust_types.h"

namespace segclust {

  void hybrid(gsl_vector *x, long P, long Kmax, long lmin, long lmax,
	    bool vh, gsl_matrix **Linc_res, param_struct *param_res[],
	    long *param_res_size);
}
#endif
