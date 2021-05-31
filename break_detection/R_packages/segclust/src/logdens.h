#ifndef SEGCLUST_LOGDENS_H
#define SEGCLUST_LOGDENS_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace segclust {

gsl_matrix *
logdens(gsl_vector *xk, int P, gsl_vector *phi);

}
#endif
