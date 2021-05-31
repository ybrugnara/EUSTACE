#ifndef SEGCLUST_MSTEP_H
#define SEGCLUST_MSTEP_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace segclust {

gsl_vector * Mstep(gsl_vector *x, gsl_matrix *rupt, gsl_matrix *tau,
		   gsl_vector *phi, bool vh);

}

#endif
