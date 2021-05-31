#ifndef SEGCLUST_ESTEP_H
#define SEGCLUST_ESTEP_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace segclust {

void Estep(gsl_matrix *logdensity, gsl_vector *phi,
	   gsl_matrix **tau, double *lvinc);

}
#endif

